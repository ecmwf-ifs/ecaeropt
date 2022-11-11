
.PHONY : clean

check :
	scripts/check_dependencies.sh	

build : check docum
	scripts/build.sh
clean : 
	scripts/clean.sh
test :
	scripts/tests.sh
	
docum :
	scripts/build_docs.sh 

ifs48R1:

ifs48R2:

ifs49R1:
