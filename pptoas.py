#!/usr/bin/env python

#Calls PulsePortraiture to generate TOAs and DM corrections

from PulsePortraiture import *
from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-d", "--datafile",
                  action="store", metavar="ARCHIVE", dest="datafile",
                  help="PSRCHIVE archive from which to generate TOAs.")
parser.add_option("-m", "--modelfile",
                  action="store", metavar="MODEL", dest="modelfile",
                  help=".model file created by ppgauss. psrsmooth models soon to be accepted.")
parser.add_option("-o", "--outfile",
                  action="store", metavar="TIMFILE", dest="outfile", default=None,
                  help="Name of output .tim file name. Will append. [default=ARCHIVE.tim]")
parser.add_option("--showplot",
                  action="store_true", dest="showplot", default=False,
                  help="Plot fit results. Only useful if nsubint > 1. [default=False]")
parser.add_option("--quiet",
                  action="store_true", dest="quiet", default=False,
                  help="Minimal to stdout.")

(options, args) = parser.parse_args()

if options.datafile is None or options.modelfile is None:
    parser.print_help()
    sys.exit()

#Make template
datafile = str(options.datafile)
modelfile = str(options.modelfile)
outfile = str(options.outfile)
showplot = options.showplot
quiet = options.quiet
#quiet = options.quiet

gt = GetTOAs(datafile,modelfile,outfile,quiet=quiet)
if showplot: gt.show_results()
