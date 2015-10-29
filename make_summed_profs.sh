#!/bin/bash
python ppalign.py -I /nimrod1/GBT/Ter5/GUPPI/Sband_tscr/*$2*.fscr -M files_to_align -o GUPPI_Ter5$1_summed
python ppalign.py -I GUPPI_Ter5$1_summed -M files_to_align -o GUPPI_Ter5$1_summed_2 
