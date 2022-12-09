
f77:
	@bash compilef77.sh > buildrun_f77.log 2>&1

build:
	@echo ""
	@echo "Compilation of source code in build.log"
	@bash compile.sh > build.log 2>&1

run: build
	@echo ""
	@echo "Running the compiled code... results in test file "
	@./bin/tmd.out   > run.log

show: f77 run
	@echo ""
	@echo "Results of refactored code"
	head --l 10 test 
	@echo ""
	@echo "Results of original code"
	@head --l 10 ../src_f77/test
	@echo ""
	@echo ""
	@diff -y --suppress-common-lines -W 200 test ../src_f77/test

