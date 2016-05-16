#!/bin/bash

# pam --mult 0.05 -m Sband/Ter5$1/*100815*$3
pam -T -e tscr Sband/Ter5$1/*$3
pam -R 180 -e rm_tscr Sband/Ter5$1/*.tscr
pam --setnchn 8 -e fscr Sband/Ter5$1/*.rm_tscr
pam --setnbin $2 -e bscr Sband/Ter5$1/*.fscr
