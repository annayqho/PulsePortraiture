import glob
from ppalign import align_archive


# template = "/nimrod1/GBT/Ter5/GUPPI/PSRs/templates/Ter5P_256_gaussians.template"

template = "/nimrod1/GBT/Ter5/GUPPI/Sband_tscr/GUPPI_Ter5ad_140501_0001.template"

# for one file
filename = "/nimrod1/GBT/Ter5/GUPPI/Sband_tscr/GUPPI_Ter5ad_150407_0001.zap"
outfile = filename+"_aligned_tscr1"
align_archive(filename, template, outfile, tfac=1)

# for multiple files

# filenames = glob.glob("/nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*Ter5O*.zap")
# outfiles = [f+"_aligned_noStokes" for f in filenames]

# for f,o in zip(filenames, outfiles):
#     align_archive(f, template, o, tfac=1)

