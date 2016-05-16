# choose templates for the summing based on S/N and flux density

import numpy as np, sys, glob, os, subprocess
import matplotlib.pyplot as plt

allPSRs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
           'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y',
           'Z', 'aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'ai']
skip = set(['A', 'P', 'ad', 'W', 'ae', 'O'])
PSRs = sorted(set(allPSRs) - skip)

def get_snr(f):
    call = "pdv -fFTp %s" %f
    # print(call)
    out = subprocess.check_output(call, shell=True)
    snr = float(out.split()[-1])
    return snr    

def get_length(f):
    call = "psrstat -c length %s" %f
    # print(call)
    out = subprocess.check_output(call, shell=True)
    length = float(out.split("=")[-1])    
    return length

def calc(psr):
    files = glob.glob("Lband/Ter5%s/*.bscr" %psr)
    snr = np.array([get_snr(f) for f in files])
    length = np.array([get_length(f) for f in files])
    expected_snr = np.sqrt((sum(length)/length))*snr
    print(psr + " " +str(np.mean(expected_snr)))


if __name__ == '__main__':
    calc('C')
    # for psr in PSRs:
    #     calc(psr)
