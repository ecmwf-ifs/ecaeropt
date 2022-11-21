
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
	sbatch scripts/create_shortifs.sh

IFS-CY46R1:
	@echo ""
	@echo " Calculation of aerosol optical properties for CY48R1 high-resolution spectral"
	@echo " -> submitted with sbatch: "
	@echo "    1. Check status with squeue -u username"
	@echo "    2. Output stored in a slurm-JOBID.log file."
	@echo ""
	sbatch scripts/create_ifs-cy46r1.sh

IFS-CY48R1:
	@echo ""
	@echo " Calculation of aerosol optical properties for CY48R1 high-resolution spectral"
	@echo " -> submitted with sbatch: "
	@echo "    1. Check status with squeue -u username"
	@echo "    2. Output stored in a slurm-JOBID.log file."
	@echo ""
	sbatch scripts/create_ifs-cy48r1.sh

ifs48R2:

ifs49R1:
