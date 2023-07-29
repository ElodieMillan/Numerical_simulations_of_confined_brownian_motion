#cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp, sqrt, log, cosh, sinh
from libc.stdlib cimport rand, RAND_MAX, srand
import time

cdef double pi = np.pi

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double gamma_x_tot(double zn, double a, double eta0, double H):
    """
    [Libshaber]
    @return: Parallel mobility of the particle depending of Zn due to walls.
    """
    # Wall of the top at z = +Hp
    cdef double zn_T = a / ((H-zn) + a)
    cdef double gam_x_T = (
        6.
        * pi
        * a
        * eta0
        *
        (
            1.
            - 9./16. * zn_T
            + 1./8. * zn_T**3.
            - 45./256. * zn_T**4.
            - 1./16. * zn_T** 5.
        )
        ** (-1)
    )
    # Wall of the bottom at z = -Hp
    cdef double zn_B = a / ((H+zn) + a)
    cdef double gam_x_B = (
        6
        * pi
        * a
        * eta0
        * (
            1.
            - 9./16. * zn_B
            + 1./8. * zn_B**3
            - 45./256. * zn_B**4
            - 1./16. * zn_B** 5
        )
        ** (-1)
    )
    cdef double gam_0 = 6 * pi * a * eta0

    return (gam_x_T + gam_x_B - gam_0)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double gamma_z_tot(double zn, double a, double eta0, double H):
    """
    [Padé approximation]
    @return: Perpendicular mobility of the particle depending of Zn due to walls.
    """
    # Wall of the top at z = +Hp
    cdef double gam_z_T = (
        6
        * pi
        * a
        * eta0
        * (
            (
                (6 * (H-zn)**2 + 9*a*(H-zn) + 2*a**2)
                / (6 * (H-zn)**2 + 2*a*(H-zn))
            )
        )
    )
    # Wall of the bottom at z = -Hp
    cdef double gam_z_B = (
        6
        * pi
        * a
        * eta0
        * ((6*(H+zn)**2 + 9*a*(H+zn) + 2*a**2)/ (6 * (H+zn)**2 + 2*a*(H+zn)))
    )
    cdef double gam_0 = 6 * pi * a * eta0

    return (gam_z_T + gam_z_B - gam_0)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double Noise(double gamma, double kBT):
    """
    :return: Noise amplitude of brownian motion.
    """
    cdef double noise = sqrt(2 * kBT / gamma)
    return noise

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double next_Xn(double xn, double zn, double Wn, double dt, double a,
                               double eta0, double kBT, double H):
    """
    :return: Parallel position at time tn+dt
    """
    cdef double gamma = gamma_x_tot(zn, a, eta0, H)  #gamma effectif avec 2 murs
    cdef double nextXn = xn + Noise(gamma, kBT) * Wn * sqrt(dt)
    return nextXn

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double Dprime_z(double zn, double kBT, double eta0, double a, double H):
    """
    :return: Spurious force to correct overdamping. (Author Dr. Maxime Lavaud).
    """
    cdef double eta_B_primes = -(a * eta0 * (2 * a ** 2 + 12 * a * (H + zn) + 21 * (H + zn) ** 2)) / (
        2 * (H + zn) ** 2 * (a + 3 * (H + zn)) ** 2
    )

    cdef double eta_T_primes = (
        a
        * eta0
        * (2 * a ** 2 + 12 * a * (H-zn) + 21 * (H-zn) ** 2)
        / (2 * (a + 3*H - 3*zn) ** 2*(H-zn) ** 2)
    )

    cdef double eta_eff_primes = eta_B_primes + eta_T_primes

    cdef double eta_B = eta0 * (6*(H+zn)**2 + 9*a*(H+zn) + 2*a**2) / (6*(H+zn)**2 + 2*a*(H+zn))
    cdef double eta_T = eta0 * (6*(H-zn)**2 + 9*a*(H-zn) + 2*a**2) / (6*(H-zn)**2 + 2*a*(H-zn))

    cdef double eta_eff = eta_B + eta_T - eta0

    return  - kBT / (6*pi*a) * eta_eff_primes / eta_eff**2


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double Forces(double zn, double H, double kBT, double B, double lD, double lB):
    """
    @return: Total extern force on particle (without friction).
    """
    cdef double Felec = B * kBT/lD * exp(-H/lD) * (exp(-zn/lD) - exp(zn/lD))
    cdef double Fgrav = -kBT/lB
    return Felec + Fgrav


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double next_Zn(double zn, double Wn, double dt, double a,
                    double eta0, double kBT, double H, double lB, double lD, double B):
    """
    @return: Perpendicular position Zn+1 at time tn+dt
    """
    cdef double gamma = gamma_z_tot(zn, a, eta0, H) #gamma effectif avec 2 murs

    cdef double nextZn = zn  + Dprime_z(zn, kBT, eta0, a, H )*dt \
                     + Forces(zn, H, kBT, B, lD, lB)*dt /gamma \
                     + Noise(gamma, kBT)*Wn*sqrt(dt)
    if nextZn < -(H):
        nextZn = -2*H - nextZn
    if nextZn > H:
        nextZn =  2*H - nextZn

    return nextZn

## Random uniform classic
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double random_uniform():
    cdef double r = rand()
    return r / RAND_MAX

'''
## MERSENNE-TWISTER ALGORITHM IN CYTHON
## Copyright (c) 2017 Anandh Swaminathan
'''
cdef unsigned NN = 312
cdef unsigned MM = 156
cdef unsigned long long MATRIX_A = 0xB5026F5AA96619E9ULL
cdef unsigned long long UM = 0xFFFFFFFF80000000ULL
cdef unsigned long long LM = 0x7FFFFFFFULL
cdef unsigned long long mt[312]
cdef unsigned mti = NN + 1
cdef unsigned long long mag01[2]

@cython.nonecheck(False)
cdef mt_seed(unsigned long long seed):
    global mt
    global mti
    global mag01
    global NN
    global MATRIX_A
    mt[0] = seed
    for mti in range(1,NN):
        mt[mti] = (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti)
    mag01[0] = 0ULL
    mag01[1] = MATRIX_A
    mti = NN

@cython.nonecheck(False)
cdef unsigned long long genrand64():
    cdef int i
    cdef unsigned long long x
    global mag01
    global mti
    global mt
    global NN
    global MM
    global UM
    global LM
    if mti >= NN:
        for i in range(NN-MM):
            x = (mt[i]&UM) | (mt[i+1]&LM)
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[int(x&1ULL)]
        for i in range(NN-MM, NN-1):
            x = (mt[i]&UM)|(mt[i+1]&LM)
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[int(x&1ULL)]
        x = (mt[NN-1]&UM)|(mt[0]&LM)
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[int(x&1ULL)]
        mti = 0
    x = mt[mti]
    mti += 1
    x ^= (x >> 29) & 0x5555555555555555ULL
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL
    x ^= (x << 37) & 0xFFF7EEE000000000ULL
    x ^= (x >> 43);
    return x

# Seed the random number generator
@cython.nonecheck(False)
cdef seed_random(unsigned long long seed):
    """
    Seed the C random number generator with the current system time.
    :return: none
    """
    if seed == 0:
        mt_seed(time.time())
    else:
        mt_seed(seed)

@cython.nonecheck(False)
@cython.cdivision(True)
cdef double uniform_rv():
    """
    :return: (double) Generate a uniform random variable in [0,1].
    """
    return (genrand64() >> 11) * (1.0/9007199254740991.0)

'''
## BOX-MULLER ALGORITHM 
'''
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double random_gaussian():
    cdef double x1, x2, w

    w = 2.0
    while (w >= 1.0):
        x1 = 2.0 * uniform_rv() - 1.0
        x2 = 2.0 * uniform_rv() - 1.0
        w = x1 * x1 + x2 * x2
    w = ((-2.0 * log(w)) / w) ** 0.5
    return x1 * w

'''
Trajectory compute (Xn, Zn).
'''
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef double[:,:] trajectory_cython(unsigned long int Nt,
                                   unsigned long int Nt_sub,
                                   double[:,:] Rn,
                                   double dt,
                                   double a, double eta0,
                                   double kBT, double H, double lB, double lD, double B):
    """
    @param Nt: Number of points.
    @param Nt_sub: Modulation of the number of points recorded.
                Exemple: if Nt_sub=10, then points are calculated every 10 steps.
    @param Rn: Total trajectory vector of size [2, Nt] (m, m).
    @param dt: Numerical time step (s).
    @param a: Particle radius (m).
    @param eta0: Fluid viscodity (Pa.s).
    @param kBT: Thermal energie with kB: Boltzman constante et T: Temperature (K).
    @param H: 2Hp = 2H + 2, is the gap between the two walls.
    @param lB: Boltzman length.
    @param lD: Debye length.
    @param B: Dimensionless constant characteristic of surface charge interactions.
    @return: Trajectory (Xn, Zn) for all time [0,Nt]*dt
    """
    cdef unsigned long int i
    cdef unsigned long int j
    cdef double Xn = Rn[0,0]
    cdef double Zn = Rn[1,0]
    seed_random(0)

    for i in range(1, Nt):
        for j in range(0, Nt_sub):
            Xn = next_Xn(Xn, Zn, random_gaussian(), dt, a, eta0, kBT, H)
            Zn = next_Zn(Zn, random_gaussian(), dt, a, eta0, kBT, H, lB, lD, B)

        Rn[0,i] = Xn
        Rn[1,i] = Zn

    return Rn