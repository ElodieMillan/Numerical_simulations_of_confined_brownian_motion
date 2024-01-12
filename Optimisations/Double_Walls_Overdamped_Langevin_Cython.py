"""
Élodie Millan
Janvier 2022
Particule confinée entre deux mur à z=+/- Hp.
Le centre de la particule ne peux pas aller plus loin que +H = Hp-a et -H = -(Hp-a).
-> comme le set up expérimentale de Maxime/Yacine/Nicolas.
"""

import numpy as np
from tqdm import tqdm
import sys
from scipy import interpolate
from Functions_cythonised import trajectory_cython

sys.path.append(r"../Overdamped_Langevin_Python")
from Overdamped_Langevin import Langevin

sys.path.append(r"../Double_Walls_Overdamped_Langevin_Python")
from Trajectory_Calculation_Functions import gamma_x_tot, gamma_z_tot


class DoubleWallsLangevin_Cy( Langevin ):
    def __init__(self, dt, Nt, a, H, lD, lB,
                 B=4.8, Nt_sub=1, eta0=0.001, T=300, R0=None):
        """
        :param dt: Numerical time step (s).
        :param Nt: Number of points.
        :param a: Particle radius (m).
        :param H: 2Hp = 2H + 2, is the gap between the two walls.
        :param lD: Boltzman length (m).
        :param lB: Debye length (m).
        :param B: Dimensionless constant characteristic of surface charge interactions.
        :param Nt_sub: Modulation of the number of points recorded.
                       Exemple: if Nt_sub=10, then points are calculated every 10 steps.
        :param eta0: Fluid viscodity (default = 0.001 (Pa.s)).
        :param T: Temperature of the system (default = 300 (K)).
        :param R0: Initial position vector |X0, Z0|,
                   (default = |0,Zeq| with Zeq shoot on equilibrium distribution Peq(z) (m,m)).
        """
        super().__init__(dt, Nt, a, eta0, T)
        self.a = a
        self.H = H
        self.Nt_sub = Nt_sub
        self.eta0 = eta0
        self.T = T
        self.kBT = 1.38e-23*self.T
        self.Hp = self.H + self.a
        self.D0 = self.kBT/(6*np.pi*self.eta0*self.a)
        self.lD = lD
        self.B = B
        self.lB = lB
        self.sample_f = self.sample()
        self.R0 = R0
        if R0 == None:
            # Random selection of altitude Z0
            # on the equilibrium distribution P_eq(z)
            self.R0 = (0.0, self.return_samples(1))

        del self.t #to save some space memory
        # because it could be compute with sel.dt and self.Nt

    def trajectory(self, output=False):
        """
        :param output: Booléen - Si VRAI, return (Par défaut : FAUX).
        :return: Trajectoire (X_n, Z_n).
        """
        Rn = np.zeros((2, self.Nt))
        # Initial condition.
        Rn[0,0] = self.R0[0]
        Rn[1,0] = self.R0[1]

        Rn = np.asarray(trajectory_cython(self.Nt, self.Nt_sub,
                                   Rn,
                                   self.dt,
                                   self.a, self.eta0,
                                   self.kBT, self.H, self.lB, self.lD, self.B))

        self.Xn = Rn[0,:]
        self.Zn = Rn[1,:]
        if output:
            return Rn

    ## SOME ANALYSIS FUNCTIONS
    def MSD(self, axis):
        """
        :param axis: Choose between "x", "y" and "z".

        :return (time, MSD(axis))
        """
        if axis == "x":
            position = self.Xn
        elif axis == "z":
            position = self.Zn
        else:
            raise ValueError('WRONG AXIS : choose between "x" and "z" !')

        list_Ntau_MSD = np.array([], dtype=int)
        for i in range(len(str(self.Nt)) - 1):
            # Take just 10 points by decade.
            list_Ntau_MSD = np.concatenate(
                (
                    list_Ntau_MSD,
                    np.arange(10 ** i, 10 ** (i + 1), 10 ** i, dtype=int),
                )
            )
        # -----------------------
        NumberOfMSDPoint = len(list_Ntau_MSD)
        msd = np.zeros(NumberOfMSDPoint)
        for k, i in enumerate(tqdm(list_Ntau_MSD)):
            if i == 0:
                msd[k] = 0
                continue
            msd[k] = np.mean((position[i:]-position[:-i])**2)

        return self.dt*self.Nt_sub*list_Ntau_MSD, msd

    def logarithmic_hist(self, position: np.ndarray, begin: float, stop: float, num: int = 50, base: int = 10) -> (np.ndarray, np.ndarray, np.ndarray):
        """Function to compute a pdf using  logspaced bins
        :param axis: "x" or "z".
        :param begin: Start value along axis.
        :param stop: End value along axis.
        :param num: Number of points (Default = 50).
        :param base: Logarithmic base value (Default = 10).
        :return: bins_center, widths, hist
        """

        if begin == 0:
            beg = stop / num
            bins = np.logspace(
                np.log(beg) / np.log(base),
                np.log(stop) / np.log(base),
                num - 1,
                base=base,
            )
            widths = bins[1:] - bins[:-1]
            bins = np.cumsum(widths[::-1])
            bins = np.concatenate(([0], bins))
            widths = bins[1:] - bins[:-1]
        else:
            bins = np.logspace(
                np.log(begin) / np.log(base),
                np.log(stop) / np.log(base),
                num,
                base=base,
            )
            widths = bins[1:] - bins[:-1]

        hist, bins = np.histogram(position, bins=bins, density=True)
        # normalize by bin width
        bins_center = (bins[1:] + bins[:-1]) / 2

        return bins_center, widths, hist

    def Cumulant4(self, axis):
        """

        :param axis: Choose between "x" and "z".
        :return (time, C4)
        """
        if axis == "x":
            position = self.Xn
        elif axis == "z":
            position = self.Zn
        else:
            raise ValueError('WRONG AXIS : choose between "x" and "z" !')

        list_Nt_c4 = np.array([], dtype=int)
        for i in range(len(str(self.Nt)) - 1):
            # Take just 10 points by decade.
            list_Nt_c4 = np.concatenate(
                (
                    list_Nt_c4,
                    np.arange(10 ** i, 10 ** (i + 1), 10 ** i, dtype=int),
                )
            )
        c4 = np.zeros(len(list_Nt_c4))
        # Compute fourth cumulant
        for k, i in enumerate(list_Nt_c4):
            if i == 0:
                c4[k] = 0
                continue
            deltaX = position[i:] - position[:-i]
            c4[k] = (np.mean(deltaX**4) - 3 * (np.mean(deltaX**2))**2)

        return self.dt*self.Nt_sub * list_Nt_c4, c4

    '''
    Some fonctions of the problem usefull.
    '''
    def D_x(self, z):
        """
        :param z: Altitude of the particle.
        :return: Effective diffusion coeficient parallele to walls.
        """
        D = [self.kBT / gamma_x_tot(i, self.a, self.eta0, self.H) for i in z]
        return np.asarray(D)

    def D_z(self, z):
        """
        :param z: Altitude of the particle.
        :return: Effective diffusion coeficient perpendiculare to walls.
        """
        D = [self.kBT / gamma_z_tot(i, self.a, self.eta0, self.H) for i in z]
        return np.asarray(D)

    def P_eq_Z(self, z):
        """
        :param z: Altitude of the particle.
        :return: Equilibrium distribution P_eq(z).
        """
        if type(z) != np.ndarray:  # for z : float
            if (z > self.H) or (z < -self.H):
                return 0
            return np.exp(-self.B * np.exp(-self.H / self.lD) * (np.exp(-z / self.lD) + np.exp(z / self.lD)) - (
                        self.H + z) / self.lB)
        # for z : np.ndarray
        Pz = lambda z: np.exp(
            -self.B * np.exp(-self.H / self.lD) * (np.exp(-z / self.lD) + np.exp(z / self.lD)) - (self.H + z) / self.lB)
        P = np.array([Pz(zz) for zz in z])
        P[z < -self.H] = 0
        P[z > self.H] = 0
        return P

    """
    Random draw on the equilibrium distribution P_eq(z) to start z_0 at equilibrium.
    """
    def f(self, z):
        # :param z: Array of altitude of the particle.
        # :return: Equilibrium distribution P_eq(z). No need to be normalized.
        zz = np.zeros_like(z)
        for n, i in enumerate(z):
            if (i < -self.H) or (i > self.H):
                zz[n] = 0
            else:
                zz[n] = np.exp(
                    -self.B * np.exp(-self.H / self.lD) * (np.exp(-i / self.lD) + np.exp(i / self.lD)) - i / self.lB)
        return zz

    def sample(self):
        # :return: inverse cumulative distribution function of P_eq(z).
        z = np.linspace(-self.H, self.H, 1000)
        y = self.f(z)  # Probability density function, pdf
        cdf_y = np.cumsum(y)  # Cumulative distribution function, cdf
        cdf_y = cdf_y / cdf_y.max()  # Takes care of normalizing cdf to 1.0
        inverse_cdf = interpolate.interp1d(cdf_y, z)  # This is a function
        return inverse_cdf

    def return_samples(self, N=1):
        # :param N: Number of random number wanted (Default = 1).
        # :return: N random altitude z generated on distribution f(z).
        try:
            uniform_samples = np.random.random(int(N))
            required_samples = self.sample_f(uniform_samples)
            return required_samples
        except ValueError:
            self.return_samples(N)


"""
FIN CLASSE
"""

def test():
    import matplotlib.pyplot as plt
    import seaborn as sns
    custom_params = {
        "xtick.direction": "in",
        "ytick.direction": "in",
        "lines.markeredgecolor": "k",
        "lines.markeredgewidth": 0.3,
        "figure.dpi": 200,
        "text.usetex": True,
        "font.family": "serif",
    }
    sns.set_theme(context="paper", style="ticks", rc=custom_params)
    a = 1.519e-6
    H = 40e-6
    lD = 88.0e-9
    T = 300
    lB = 526e-9
    B = 5.0
    eta0 = 0.001
    simu = DoubleWallsLangevin_Cy(dt=0.01, Nt=100_000,
                               a=a, H=H, lD=lD, lB=lB, B=B,
                               Nt_sub=1, eta0=eta0,
                               T=T, R0=None)
    simu.trajectory()
    ###Plot MSD
    MSDz, tau_z = simu.MSD(axis="z")
    dt_theo = np.linspace(simu.dt, simu.dt * simu.Nt, 100)
    plt.figure(figsize=(1.5 * 3.375, 1.5 * 3.375 / 1.68), tight_layout=True)
    plt.loglog(tau_z, MSDz, "o")
    plt.plot(tau_z, 2 * simu.D0 * tau_z, "k-")
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\langle X_t^2 \rangle$")
    plt.show()

if __name__ == '__main__':
    test()