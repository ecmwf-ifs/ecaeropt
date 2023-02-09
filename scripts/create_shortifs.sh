#!/bin/bash

#SBATCH --job-name=ecaeropt-ifs-test
#SBATCH --output=logs/calc.short_ifs_example.%j.out

module load python3

source scripts/create_ifs_prepare.sh

./ecaeropt -s tests/settings/short_ifs_example.toml


