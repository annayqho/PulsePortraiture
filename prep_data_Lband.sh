#!/bin/bash

pam --mult 0.05 -m Lband/Ter5$1/*100501*$3
pam -T -e tscr Lband/Ter5$1/*$3
pam -R 180 -e rm_tscr Lband/Ter5$1/*.tscr
pam --setnchn 8 -e fscr Lband/Ter5$1/*.rm_tscr
pam --setnbin $2 -e bscr Lband/Ter5$1/*.fscr
