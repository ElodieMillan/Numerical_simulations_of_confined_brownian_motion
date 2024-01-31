#!/usr/bin/env python
import os

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

'''Parameter: should be identical to Multiprocessing.py file (defined here to save in proper name)'''
bins = 200  # Choose the number of bins to compute PDF
Nsigma = 6  # Choose the width of PDF by unit of STD
Ns = 10_000 # Choose the number of trajectory wanted
Njob = 200 # Choose the number total of job wanted (Nb total of simulation = Ns*Njob)

'''Defined jobs and save directories/files'''
# With "os.getcwd()", the current working directory of a process
job_directory = "%s/.job" %os.getcwd()
data_dir = "{}/DatasPDFs_simuMultipross_dXbins{}_Nsigma{}_Nsimu{}_Njob{}".format(os.getcwd(), bins, Sigma, Ns, Njob)
mkdir_p(data_dir)
jobs = [str(i) for i in range(0, Njob)]

for j in jobs:
    job_file = os.path.join(job_directory,"%s.job" %j)
    data = os.path.join(data_dir, j)
    # Create lizard directories
    f = data + ".npy" #save data as .npy

    with open(job_file, "w") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=%s.job\n" % j)
        fh.writelines("#SBATCH --output=.out/%s.out\n" % j)
        fh.writelines("#SBATCH --error=.out/%s.err\n" % j)
        fh.writelines("#SBATCH --time=4:30:00\n")
        #fh.writelines("#SBATCH -n 64\n")
        fh.writelines("#SBATCH -N 1\n")
        fh.writelines("#SBATCH --ntasks-per-node=32\n")
        fh.writelines("#SBATCH --mem 30G\n")
        fh.writelines("python Multiprocessing.py %s\n" % f )

    os.system("sbatch %s" %job_file)
