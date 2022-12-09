

 gfortran -g -c src_tmat/tmd_params.lp.f90      -Jmod -o objs/test_params_del.o
 gfortran -g -c src_tmat/tmd_aux.lp.f90         -Jmod -o objs/test_aux_del.o
 gfortran -g -c src_tmat/tmd_CGcoeff.lp.f90     -Jmod -o objs/test_CGcoeff_del.o
 gfortran -g -c src_tmat/tmd_bessel.lp.f90      -Jmod -o objs/test_bessel_del.o
 gfortran -g -c src_tmat/tmd_gsp.lp.f90         -Jmod -o objs/test_gsp_del.o
 gfortran -g -c src_tmat/tmd_tmat_proced.lp.f90 -Jmod -o objs/test_proc_del.o
 gfortran -g -c src_tmat/tmd_tmatrix.lp.f90     -Jmod -o objs/test_tmatrix_del.o

 gfortran -g tmd.lp.f90 src_lapack/lpd.f -Jmod -o bin/tmd.out objs/test_proc_del.o objs/test_aux_del.o objs/test_CGcoeff_del.o objs/test_params_del.o objs/test_bessel_del.o objs/test_gsp_del.o objs/test_tmatrix_del.o
