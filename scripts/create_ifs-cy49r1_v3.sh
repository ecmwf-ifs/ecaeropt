#!/bin/bash

#SBATCH --job-name=ecaeropt-ifs-test
#SBATCH --output=logs/calc.ifs_CY49R1v3.%j.out

module load python3
module load nco

source scripts/create_ifs_prepare.sh

./ecaeropt -s settings/IFS_CY49R1_v3.toml


