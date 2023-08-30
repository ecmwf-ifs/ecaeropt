#  +----------------------------------------------------------------------------------------+
#  | Makefike                                                                               |
#  |                                                                                        |
#  |   (C) Copyright 2022- ECMWF.                                                           |
#  |                                                                                        |
#  |   This software is licensed under the terms of the Apache Licence Version 2.0          |
#  |   which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                 |
#  |                                                                                        |
#  |   In applying this licence, ECMWF does not waive the privileges and immunities         |
#  |   granted to it by virtue of its status as an intergovernmental organisation           |
#  |   nor does it submit to any jurisdiction.                                              |
#  |                                                                                        |
#  |  Author:                                                                               |
#  |     Ramiro Checa-Garcia. ECMWF                                                         |
#  |                                                                                        |
#  |  Modifications:                                                                        |
#  |     10-Dec-2022   Ramiro Checa-Garcia    1st. version                                  |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+

.PHONY : clean

check :
	scripts/check_dependencies.sh	

docum :
	scripts/build_docs.sh 

build : check docum
	scripts/build.sh
clean : 
	scripts/clean.sh
test :
	scripts/tests.sh
	
IFS-example:
	@echo ""
	@echo " Example of calculation of IFS aerosol optical properties (few wavelengths)"
	@echo " -> submitted with sbatch: "
	@echo "    1. Check status with squeue -u username"
	@echo "    2. Output stored in a logs/calc.short_ifs_example.JOBID.out file."
	@echo ""
	sbatch scripts/create_shortifs.sh

IFS-sbatch:
	@echo ""
	@echo " Calculation of IFS aerosol optical properties at high-resolution spectral"
	@echo " -> submitted with sbatch: "
	@echo "    1. Check status with squeue -u username"
	@echo "    2. Output stored in a logs/calc.IFS-version.JOBID.out file."
	@echo ""
	
IFS-CY46R1: IFS-sbatch
	@echo " === IFS-CY46R1 ==="
	sbatch scripts/create_ifs-cy46r1.sh

IFS-CY48R1: IFS-sbatch
	@echo " === IFS-CY48R1 ==="
	sbatch scripts/create_ifs-cy48r1.sh

IFS-CY49R1: IFS-sbatch
	@echo " === IFS-CY49R1 ==="
	sbatch scripts/create_ifs-cy49r1.sh 

IFS-CY49R1-McClear: IFS-sbatch
	@echo " === IFS-CY49R1-McClear ==="
	sbatch scripts/create_ifs_McClear.sh
	


