#!/usr/bin/env python

#########
#ppalign#
#########

#ppalign is a command-line program used to align homogeneous data (i.e. from
#    the receiver, with the same center frequency, bandwidth, and number of
#    channels).  This is useful for making averaged portraits to either pass to
#    ppgauss.py with -M to make a Gaussian model, or to smooth and use as a
#    model with pptoas.py.

#Written by Timothy T. Pennucci (TTP; pennucci@virginia.edu).

#Need option for constant Gaussian initial guess.

import os, shlex
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sub
import psrchive
from pptoas import *

def psradd_archives(metafile, outfile, palign=False):
    """Add together archives using psradd.

    This function will call psradd with an option to pass -P and can be used to
    make an initial guess for align_archives.

    metafile is a file containing PSRFITS archive names to be averaged.
    outfile is the name of the output archive.
    palign=True passes -P to psradd, which phase-aligns the archives, intead of
        using the ephemeris (maybe?).
    """
    psradd_cmd = "psradd "
    if palign:
        psradd_cmd += "-P "
    psradd_cmd += "-T -o %s -M %s"%(outfile, metafile)
    psradd_call = sub.Popen(shlex.split(psradd_cmd))
    psradd_call.wait()


def psrsmooth_archive(archive, options="-W"):
    """Smooth an archive using psrsmooth.

    This function will call psrsmooth with options to smooth an output archive
    from align_archives.

    archive is the PSRFITS archive to be smoothed.
    options are the options passed to psrsmooth.
    """
    psrsmooth_cmd = "psrsmooth " + options + " %s"%archive
    psrsmooth_call = sub.Popen(shlex.split(psrsmooth_cmd))
    psrsmooth_call.wait()


def align_archives(metafile, initial_guess, outfile=None, rot_phase=0.0,
        place=None, niter=1, quiet=False):
    """Iteratively align and average archives.

    Each archive is fitted for a phase, a DM, and channel amplitudes against
    initial_guess.  The average is weighted by the fitted channel amplitudes
    and channel S/N.  The average becomes the new initial alignment template
    for additional iterations.  The output archive will have a 0 DM value and
    dmc=0.

    metafile is a file containing PSRFITS archive names to be averaged.
    initial_guess is the PSRFITS archive providing the initial alignment guess.
    outfile is the name of the output archive; defaults to
        <metafile>.algnd.fits.
    rot_phase is an overall rotation to be applied to the final output archive.
    niter is the number of iterations to complete.  1-5 seems to work ok.
    quiet=True suppresses output.

    """
    datafiles = [datafile[:-1] for datafile in open(metafile, "r").readlines()]
    if outfile is None:
        outfile = metafile + ".algnd.fits"

    # all of the files should be identical
    datasize = load_data(datafiles[0]).subints.shape
    nfiles = len(datafiles)
    npol = datasize[1]
    nchan = datasize[2]
    nbin = datasize[3]
    nsub = 1

    model_data = load_data(initial_guess, dedisperse=True, dededisperse=False,
            tscrunch=True, pscrunch=True, fscrunch=False, rm_baseline=True,
            flux_prof=False, refresh_arch=True, return_arch=True, quiet=quiet)
    model_port = (model_data.masks * model_data.subints)[0,0]
    count = 1
    while(niter):
        load_quiet = quiet
        # 1 is here b/c there should only be 1 subint
        aligned_subint = np.zeros((nsub, npol, nchan, nbin))
        total_weights = np.zeros(np.shape(aligned_subint))
        for ifile in xrange(len(datafiles)):
            data_tot = load_data(datafiles[ifile], dedisperse=False,
                        tscrunch=True, pscrunch=True, fscrunch=False,
                        rm_baseline=True, quiet=load_quiet)
            data = load_data(datafiles[ifile], dedisperse=False,
                    tscrunch=True, pscrunch=False, fscrunch=False,
                    rm_baseline=True, quiet=load_quiet)
            DM_guess = data_tot.DM
            for isub in data_tot.ok_isubs:
                print("subint %s" %isub)
                # 0 is here because we have pscrunched! 
                port = data_tot.subints[isub,0,data_tot.ok_ichans[isub]]
                freqs = data_tot.freqs[isub,data_tot.ok_ichans[isub]]
                model = model_port[data_tot.ok_ichans[isub]]
                P = data_tot.Ps[isub]
                SNRs = data_tot.SNRs[isub,0,data_tot.ok_ichans[isub]]
                errs = data_tot.noise_stds[isub,0,data_tot.ok_ichans[isub]]
                nu_fit = guess_fit_freq(freqs, SNRs)
                rot_port = rotate_data(port, 0.0, DM_guess, P, freqs, nu_fit)
                phase_guess = fit_phase_shift(rot_port.mean(axis=0),
                        model.mean(axis=0)).phase
                if len(freqs) > 1:
                    print("starting fit")
                    results = fit_portrait(port, model,
                            np.array([phase_guess, DM_guess]), P, freqs,
                            nu_fit, None, errs, quiet=False)
                else:  #1-channel hack
                    results = fit_phase_shift(port[0], model[0], errs[0])
                    results.DM = data.DM
                    results.DM_err = 0.0
                    results.nu_ref = freqs[0]
                    results.nfeval = 0
                    results.return_code = -2
                    results.scales = np.array([results.scale])
                    results.scale_errs = np.array([results.scale_error])
                    results.covariance = 0.0
                weights = np.outer(results.scales / errs**2, np.ones(nbin))
                for i in range(0, npol):
                    choose = data.subints[isub,i,data.ok_ichans[isub]]
                    aligned_subint[isub,i,data.ok_ichans[isub]] += weights * \
                        rotate_data(choose,results.phase,results.DM,P,freqs,results.nu_ref)
                    total_weights[isub, i, data.ok_ichans[isub]] +=  weights
            load_quiet = True
        for ipol in range(0, npol):
            aligned_subint[0,ipol,np.where(total_weights[0,ipol] > 0)[0]] /= \
                total_weights[0,ipol,np.where(total_weights[0,ipol] > 0)[0]]
        aligned_port = np.sum(np.sum(aligned_subint, axis=0), axis=0)
        model_port = aligned_port
        niter -= 1
        count += 1
    if rot_phase:
        aligned_subint = rotate_data(aligned_subint, rot_phase)
    if place is not None:
        prof = aligned_port.mean(axis=0)
        delta = prof.max() * gaussian_profile(len(prof), place, 0.0001)
        phase = fit_phase_shift(prof, delta).phase
        aligned_subint = rotate_data(aligned_subint, rot_phase)
    # create template for the output
    tmp = "psradd.tmp.fits"
    psradd_archives(metafile, outfile=tmp, palign=True)
    template_out = load_data(tmp, dedisperse=True, dededisperse=False,
            tscrunch=True, pscrunch=True, fscrunch=False, rm_baseline=True,
            flux_prof=False, refresh_arch=True, return_arch=True, quiet=quiet)
    arch = template_out.arch
    print("tscrunching")
    arch.tscrunch()
    arch.set_dispersion_measure(0.0)
    for isub,subint in enumerate(arch):
        for ipol in range(npol):
            for ichan in range(nchan):
                prof = subint.get_Profile(ipol, ichan)
                prof.get_amps()[:] = aligned_subint[isub, ipol, ichan]
                if total_weights[isub, ipol, ichan].sum() == 0.0:
                    subint.set_weight(ichan, 0.0)
    arch.unload(outfile)
    if not quiet: print "\nUnloaded %s.\n"%outfile


def align_archive(filename, template, outfile, tfac=1):
    """ Align an archive across subintegrations
    
    Dedisperses and fscrunches the archive, then fits for the shift between 
    the profile at each subint and the input template. 
    Applies this shift to each subint of the non-fscrunched, non-dedispersed archive.
    Saves this shifted archive.

    Parameters
    ----------
    filename: name of the input archive to be aligned
    template: file to align each subint against
    outfile: name of output file to save
    tfac: factor to tscrunch by
    """
    arch = psrchive.Archive_load(filename)        
    arch.convert_state('Stokes')
    if arch.get_dedispersed() is True:
        print("Already dedispersed")
    else:
        if arch.get_dispersion_measure() == 0:
            print("Bad dispersion measure")
        else: arch.dedisperse()
    if arch.get_nchan() == 1:
        print("already fscrunched")
    else:
        arch.fscrunch()
    tfac = int(tfac)
    arch.tscrunch(tfac)

    arch_template = psrchive.Archive_load(template)
    if arch_template.get_npol() > 1:
        print("Warning: template has > 1 pols")
    if arch_template.get_nchan() > 1:
        print("Warning: template has > 1 channels")
    if arch_template.get_nsubint() > 1:
        print("Warning: template has > 1 subints")
    # careful: this assumes that the template has data shape (1,1,1,nbins)
    tmpl_prof = arch_template.get_Profile(0,0,0)

    shifts = np.zeros(arch.get_nsubint())
    var_shifts = np.zeros(shifts.shape)
    psf = psrchive.ProfileShiftFit()
    psf.set_standard(tmpl_prof)
    for subint in range(arch.get_nsubint()):
        data_prof = arch.get_Profile(subint,0,0)
        psf.set_Profile(data_prof)
        (shift, var_shift) = psf.get_shift()
        shifts[subint] = shift
        var_shifts[subint] = var_shift
    shifts[np.isnan(var_shifts)] = 0.
    arch_to_change = psrchive.Archive_load(filename)
    arch_to_change.tscrunch(tfac)
    for isub,subint in enumerate(arch_to_change):
        subint.rotate_phase(shifts[isub])
    arch_to_change.unload(outfile)


if __name__ == "__main__":

    from optparse import OptionParser

    usage = "Usage: %prog -M <metafile> [options]"
    parser = OptionParser(usage)
    #parser.add_option("-h", "--help",
    #                  action="store_true", dest="help", default=False,
    #                  help="Show this help message and exit.")
    parser.add_option("-M", "--metafile",
                      default=None,
                      action="store", metavar="metafile", dest="metafile",
                      help="Metafile of archives to average together.")
    parser.add_option("-I", "--init",
                      default=None,
                      action="store", metavar="initial_guess",
                      dest="initial_guess",
                      help="Archive containing initial alignment guess.  psradd is used if -I is not used. [default=None]")
    parser.add_option("-o", "--outfile",
                      default=None,
                      action="store", metavar="outfile", dest="outfile",
                      help="Name of averaged output archive. [default=metafile.algnd.fits]")
    parser.add_option("-P", "--palign",
                      default=False,
                      action="store_true", dest="palign",
                      help="Passes -P to psradd if -I is not used. [default=False]")
    parser.add_option("-s", "--smooth",
                      default=False,
                      action="store_true", dest="smooth",
                      help="Smooth the output average (second output archive) with psrsmooth -W. [default=False]")
    parser.add_option("-r", "--rot",
                      default=0.0,
                      action="store", metavar="phase", dest="rot_phase",
                      help="Additional rotation to add to averaged archive. [default=0.0]")
    parser.add_option("-t", "--templatefile",
                      default=None,
                      action="store", dest="template",
                      help="Template file used to compare and align the profile in each subintegration. [default=None]")
    parser.add_option("--tfac",
                      action="store", metavar="int", dest="tfac", default=1,
                      help="Before aligning subint, tscrunch by this factor. [default=1]")
    parser.add_option("--place",
                      default=None,
                      action="store", metavar="place", dest="place",
                      help="Roughly place pulse to be at the phase of the provided argument.  Overrides --rot. [default=None]")
    parser.add_option("--niter",
                      action="store", metavar="int", dest="niter", default=1,
                      help="Number of iterations to complete. [default=1]")
    parser.add_option("--verbose",
                      action="store_false", dest="quiet", default=True,
                      help="More to stdout.")

    (options, args) = parser.parse_args()

    if options.metafile is None or not options.niter:
        print "\nppalign.py - Aligns and averages homogeneous archives by fitting DMs and phases\n"
        parser.print_help()
        print ""
        parser.exit()

    metafile = options.metafile
    initial_guess = options.initial_guess
    outfile = options.outfile
    palign = options.palign
    smooth = options.smooth
    rot_phase = np.float64(options.rot_phase)
    if options.template is not None:
        template = options.template
        tfac = options.tfac
    else:
        template = None
        tfac = None
    if options.place is not None:
        place = np.float64(options.place) 
        rot_phase = 0.0
    else:
        place = None
    niter = int(options.niter)
    quiet = options.quiet

    rm = False
    if initial_guess is None:
        tmp_file = "ppalign.tmp.fits"
        print(metafile)
        psradd_archives(metafile, outfile=tmp_file, palign=palign)
        initial_guess = tmp_file
        #rm = True
    
    if template is not None:
        print("template file provided, aligning subints")
        new_metafile = "%s_aligned" %metafile
        datafiles = [datafile[:-1] for datafile in open(metafile, "r").readlines()]
        datafiles_new = ["%s_aligned" %datafile for datafile in datafiles]
        np.savetxt(new_metafile, datafiles_new, fmt='%s', delimiter='\n')
        for ifile in xrange(len(datafiles)):
            align_archive(datafiles[ifile], template, datafiles_new[ifile], tfac=tfac)
        metafile = new_metafile
        
    print("Using metafile %s" %metafile)
    align_archives(metafile, initial_guess=initial_guess, outfile=outfile,
            rot_phase=rot_phase, place=place, niter=niter, quiet=quiet)
    if smooth:
        if outfile is None:
            outfile = metafile + ".algnd.fits"
        psrsmooth_archive(outfile, options="-W")
    if rm:
        rm_cmd = "rm -f %s"%tmp_file
        rm_call = sub.Popen(shlex.split(rm_cmd))
        rm_call.wait()
