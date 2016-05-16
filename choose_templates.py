# choose templates for the summing based on S/N and flux density

import numpy as np, sys, glob, os, subprocess
import matplotlib.pyplot as plt

allPSRs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
           'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y',
           'Z', 'aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'ai']
skip = set(['A', 'P', 'ad', 'W', 'ae', 'O'])
PSRs = sorted(set(allPSRs) - skip)

def get_vals(f):
    call = "pdv -fFTp %s" %f
    print(call)
    out = subprocess.check_output(call, shell=True)
    snr = float(out.split()[-1])
    flux = float(out.split()[16])
    return snr, flux    

def plot(psr):
    files = glob.glob("Sband/Ter5%s/*.bscr" %psr)
    a = np.array([get_vals(f) for f in files])
    S = a[:,1] # flux density
    snr = a[:,0]

    fig, ax = plt.subplots()
    ax.scatter(S, snr, c='k')
    plt.title("Ter5%s Observations" %psr)
    plt.xlabel("Flux Density (mJy)")
    plt.ylabel("SNR")
    for i, txt in enumerate(files):
        ax.annotate(txt, (S[i], snr[i]))
    plt.show()

if __name__ == '__main__':
    plot("X")
    # for psr in ['P', 'ad']:
    #      plot(psr)    
