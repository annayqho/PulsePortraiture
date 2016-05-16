import os, glob, sys, shlex
import numpy as np
from ppalign import align_archive
import subprocess as sub

def get_data(psr):
    if np.logical_or(psr == "P", psr == "ad"):
        script = "get_data_Sband_P_ad.sh"
    else:
        script = "get_data_Sband.sh"
    call = "sh %s %s zap" %(script, psr)
    print(call)
    os.system(call)
    call = "sh %s %s calib" %(script, psr)
    print(call)
    os.system(call)
    # when .zap exists, delete .calib
    files = glob.glob("Sband/Ter5%s/*.zap" %psr)
    for f in files:
        call = "rm " + f.split('.')[0] + ".calib"
        print(call)
        os.system(call)
    

def prep_data(psr, nbin, extension):
    call = "sh prep_data_Sband.sh %s %s %s" %(psr, nbin, extension)
    print(call)
    os.system(call)


def update_ephemeris(psr):
    if psr == 'ad':
        print("updating eph for Ter5ad")
        eph_file = "/nimrod1/GBT/Ter5/GUPPI/parfiles/1748-2446%s.par" %psr
    elif psr == 'P':
        print("updating eph for Ter5P")
        files = glob.glob("Sband/Ter5P/*.zap")
        direc = "/nimrod2/bprager/ECLIPSES/TER5/TER5P/"
        dates = np.array([f.split('_')[2] for f in files])
        eph_files = np.array([glob.glob(direc+d+"/FLUXES/*/%s.par" %d)[0] for d in dates])
        for i,f in enumerate(files):
            call = "pam -e zap_ne -E %s %s" %(eph_files[i],f)
            print(call)
            os.system(call)
    else:
        eph_file = "/nimrod1/GBT/Ter5/GUPPI/2015_timing_update/1748-2446%s.par" %psr
    if psr != 'P':
        call = "pam -m -E %s Sband/Ter5%s/*.zap" %(eph_file,psr)
        print(call)
        os.system(call)


def align_data(psr, tscr = 1):
    toalign = glob.glob("Sband/Ter5%s/*.zap" %psr)
    for archive in toalign:
        vap_cmd = "vap -c nbin %s" %archive 
        info = sub.Popen(
            shlex.split(vap_cmd), stdout=sub.PIPE).stdout.readlines()[-1]
        nbin = info.split(" ")[-1].replace("\n","")
        #if psr == 'ad':
        if psr == 'ah': 
           template = "/nimrod1/GBT/Ter5/GUPPI/PSRs/templates/Ter5%s_%s_gaussians.template" %('G', nbin)
        elif psr == 'P':
            template = "Sband/Ter5P/Ter5P.template"
            archive = archive.split(".")[0] + ".eclipse_zap"
        else:
            template = "/nimrod1/GBT/Ter5/GUPPI/PSRs/templates/Ter5%s_%s_gaussians.template" %(psr, nbin)
        print("aligning %s") %archive
        align_archive(archive, template, archive+"_aligned", tfac= tscr)


def write_files_to_align(psr):
    files = glob.glob("Sband/Ter5%s/*.bscr" %psr)
    inputf = open("files_to_align", "w")
    for f in files:
        inputf.write(f+'\n')
    inputf.close()


def sum_archives(psr, template):
    write_files_to_align(psr)
    tf = "Sband/Ter5%s/GUPPI_Ter5%s_%s_0001.bscr" %(psr, psr, template)
    print("one iteration")
    call="python ppalign.py -I %s -M files_to_align -o GUPPI_Ter5%s_summed_numit1" %(tf, psr)
    print(call)
    os.system(call)
    call="python ppalign.py -I %s -M files_to_align -o GUPPI_Ter5%s_summed_numit2 --niter 2" %(tf, psr)
    print(call)
    os.system(call)
    call="python ppalign.py -I %s -M files_to_align -o GUPPI_Ter5%s_summed_numit3 --niter 3" %(tf, psr)
    print(call)
    os.system(call)
    call="python ppalign.py -I %s -M files_to_align -o GUPPI_Ter5%s_summed_numit4 --niter 4" %(tf, psr)
    print(call)
    os.system(call)


if __name__ == "__main__":
    sys.path.append("/nimrod1/aho/alignment/PulsePortraiture")    
    allPSRs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
               'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y',
               'Z', 'aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'ai']
    bins = [256, 256, 256, 256, 256, 256, 64, 256, 256, 256, 128, 256, 
            256, 256, 256, 256, 128, 128, 128, 128, 64, 256, 128, 128,
            128, 128, 64, 64, 256, 128, 128, 64, 64, 256]
    a = np.array(open("templates_S.txt", "r").read().split('\n'))
    a = np.delete(a, -1) # last element is empty
    p = [b.split(' ')[0] for b in a]
    temp = [b.split(' ')[-1] for b in a]
    templates = dict(zip(p, temp))
    nbins = dict(zip(allPSRs, bins))
    skip = set(['A', 'P', 'ad', 'W', 'ae', 'O'])
    PSRs = sorted(set(allPSRs)-skip)

    # sum the isolated pulsars
    # for psr in PSRs:
        # get_data(psr)
        # update_ephemeris(psr)
        # prep_data(psr, nbins[psr], ".zap")
    #     sum_archives(psr, templates[psr]) 
    # prep_data('ah', nbins['ah'], ".zap")
    # sum_archives('ah', templates['ah'])

    # sum the fast binaries
    psr = 'W'
    # get_data(psr) WARNING don't do this again b/c you will overwrite
    # all my beautifully-zapped files
    # update_ephemeris(psr)
    # align_data(psr, tscr=2)

    psr = 'ae'
    # get_data(psr) WARNING don't do this again b/c you will overwrite
    # all my beautifully-zapped files
    # update_ephemeris(psr)
    # align_data(psr, tscr=4)

    psr = 'O'
    # get_data(psr) WARNING don't do this again b/c you will overwrite
    # all my beautifully-zapped files
    # update_ephemeris(psr)
    # align_data(psr, tscr=1)
    
    # for psr in ['W', 'ae', 'O']:
        # prep_data(psr, nbins[psr], ".zap")
        # make sure you've run choose_template
    #     sum_archives(psr, templates[psr])

    # update_ephemeris('P')    

    # align_data('ad', tscr=1)
    # align_data('P', tscr=1)

    # prep_data('ad', nbins['ad'], ".zap_aligned")
    # prep_data('P', nbins['P'], ".eclipse_zap")
    # prep_data('ah', nbins['ah'], ".zap_aligned")
    # sum_archives('ah', templates['ah'])

    # for psr in ['P']:
    #sum_archives('P', templates['P'])
    #align_data('X', tscr=1)
    #prep_data('X', nbins['X'], ".zap_aligned")
    sum_archives('X', templates['X'])
