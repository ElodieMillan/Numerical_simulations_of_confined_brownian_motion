import numpy as np
import time
import sys
from tqdm import tqdm
import multiprocessing
sys.path.append('../Optimisations')
from Double_Walls_Overdamped_Langevin_Cython import DoubleWallsLangevin_Cy


t1 = time.time()
input1 = sys.argv[1] # Name of file to save data (automatized in batch_multiprocessing.py)

'''Fixed parameter'''
bins = 100  # Choose the number of bins to compute PDF
Nsigma = 6  # Choose the width of PDF by unit of STD
Ns = 100 # Choose the number of trajectory wanted

a = 1.519e-6
H = 40e-6
lD = 88.0e-9
lB = 526e-9
B = 5.0
dt = 1e-2
Nt = 1_000_000

list_Ntau = np.logspace(np.log(1)/np.log(10), np.log(1000001)/np.log(10),120).astype(int)
list_Ntau = np.unique(list_Ntau)
List_tau = np.array(list_Ntau) * dt #list of lag time to cumpue C4

def calcul(_):
    '''Initialization'''
    PDFS = np.zeros((len(List_tau), bins))
    simu = DoubleWallsLangevin_Cy(dt=dt, Nt=Nt, a=a, H=H, lD=lD, lB=lB, B=B, Nt_sub=1)
    Sigma0 = np.sqrt(2 * simu.D0 * List_tau) #List of all gaussian STD in function of lag time tau (order of magnitude)
    '''Compute trajectory'''
    simu.trajectory()
    for n, i in enumerate(list_Ntau):
        dX  = simu.Xn[i:] - simu.Xn[:-i]
        STD = Sigma0[n]
        X, _ = np.histogram(dX, bins = bins, range=[-STD*Nsigma, STD*Nsigma]) #Compute histogramme
        PDFS[n,:] = X
    del simu    
    return PDFS

def calcul_z(_):
    '''Initialization'''
    PDFS = np.zeros((len(List_tau), bins))
    simu = DoubleWallsLangevin_Cy(dt=dt, Nt=Nt, a=a, H=H, lD=lD, lB=lB, B=B, Nt_sub=1)
    Sigma0 = np.sqrt(2 * simu.D0 * List_tau) #List of all gaussian STD in function of lag time tau (order of magnitude)
    '''Compute trajectory'''
    simu.trajectory()
    for n, i in enumerate(list_Ntau):
        dZ  = simu.Zn[i:] - simu.Zn[:-i]
        STD = Sigma0[n]
        tau = simu.dt*i
        if tau<50:
            Z, _ = np.histogram(dZ, bins=bins, range=[-STD * Nsigma, STD * Nsigma])
        elif tau >=50: # Adapt bins for equilibrium
            Z, _ = np.histogram(dZ, bins=bins, range=[-15e-6, 15e-6])
        PDFS[n,:] = Z
    del simu
    return PDFS

pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()) #count number of CPU
# run the "calcul" function in parallel on the CPUs Ns times
result = list(tqdm(pool.imap(calcul_z, range(Ns)), total=Ns))
#result = pool.map(calcul_z, range(Ns))
PDFS = np.zeros((len(List_tau), bins)) #initialisation
for i in result:
    PDFS += i[1] # sum all the PDFs for the final result

np.save(input1, PDFS)
print(time.time() - t1)

