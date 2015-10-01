import glob
from ppalign import align_archive

filenames = glob.glob("/nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*Ter5O*.zap")
template = "/nimrod1/GBT/Ter5/GUPPI/PSRs/templates/Ter5O_256_gaussians.template"
outfiles = [f+"_aligned" for f in filenames]

for f,o in zip(filenames, outfiles):
    align_archive(f, template, o, tfac=1)

