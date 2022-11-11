
.PHONY : clean
check :
	scripts/check_dependencies.sh	
build : check
	scripts/build.sh
clean : 
	scripts/clean.sh
test :
	scripts/tests.sh
	
ifs48R1:

ifs48R2:

ifs49R1:
