#!/bin/bash

# for isolated pulsars, just pull in the .tscr
# sh /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/get_data.sh $1 tscr
# sh /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/prep_data.sh

# pick the brightest profile to serve as the template
python ppalign.py -I /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*$2*.fscr -M files_to_align -o GUPPI_Ter5$1_summed
python ppalign.py -I /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*$2*.fscr -M files_to_al
ign -o GUPPI_Ter5$1_summed_numit2 --niter 2
python ppalign.py -I /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*$2*.fscr -M files_to_al
ign -o GUPPI_Ter5$1_summed_numit3 --niter 3
python ppalign.py -I /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*$2*.fscr -M files_to_al
ign -o GUPPI_Ter5$1_summed_numit4 --niter 4

python ppalign.py -I GUPPI_Ter5$1_summed -M files_to_align -o GUPPI_Ter5$1_summed_2
python ppalign.py -I GUPPI_Ter5$1_summed_numit2 -M files_to_align -o GUPPI_Ter5$1_summed_2_numit2
python ppalign.py -I GUPPI_Ter5$1_summed_numit3 -M files_to_align -o GUPPI_Ter5$1_summed_2_numit3
python ppalign.py -I GUPPI_Ter5$1_summed_numit4 -M files_to_align -o GUPPI_Ter5$1_summed_2_numit4
