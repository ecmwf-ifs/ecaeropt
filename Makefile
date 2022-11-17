
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
	
IFS-CY46R1:
	./ecaeropt -s settings/IFS_CY46R1.toml

IFS-CY46R1S-CY48R1:
	./ecaeropt -s settings/IFS_CY48R1.toml

ifs48R2:

ifs49R1:
