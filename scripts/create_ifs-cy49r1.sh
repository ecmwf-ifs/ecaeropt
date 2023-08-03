#!/bin/bash

################################################################################
# scripts/create_ifs_cy49r1.sh
#
#   (C) Copyright 2022- ECMWF.
#  
#   This software is licensed under the terms of the Apache Licence Version 2.0
#   which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
#   In applying this licence, ECMWF does not waive the privileges and immunities
#   granted to it by virtue of its status as an intergovernmental organisation
#   nor does it submit to any jurisdiction.
#
#  Author:
#     Ramiro Checa-Garcia. ECMWF
# 
#################################################################################


#SBATCH --job-name=ecaeropt-ifs-test
#SBATCH --output=logs/calc.ifs_CY49R1.%j.out

module load python3/old
module load nco

source scripts/create_ifs_prepare.sh

./ecaeropt -s settings/IFS_CY49R1_v4.toml


