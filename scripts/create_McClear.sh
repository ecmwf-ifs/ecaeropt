#!/bin/bash

#SBATCH --job-name=ecaeropt-ifs-test
#SBATCH --output=logs/calc.ifs_CY49R1_McClear.%j.out

module load python3/3.8.8-01
module load nco

source scripts/create_ifs_prepare.sh

./ecaeropt -s settings/IFS_CY49R1_McClear.toml

