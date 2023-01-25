#/bin/bash

# ================================================================
# Building and installation script for engine libraries
#
# Date: 2022-Nov-09
# Author: Ramiro Checa-Garcia
# Contact: ramiro.checa-garcia at ecmwf.int
#
# INFO: Code is organized with functions, each function is doing
#       the build. TODO change this for  gnu-make
#
# ================================================================


# (0) Variables needed for this script.
#     We assume that user did not change directories structure.

path_ecaeropt=$(pwd)
date_now=$(date "+%Y%m%d_%H%M")

path_mie_BB=$path_ecaeropt"/engines/mie_Boucher_Bozzo"
path_Tmatrx=$path_ecaeropt"/engines/tmatrix_mischenko"


echo ""
echo "..... Checking directories logs and libs"

     if [ ! -d "logs" ]; then
         echo "       -> directory logs not present."
	 echo "       -> creating it."
         mkdir logs
     fi

     if [ ! -d $path_mie_BB"/libs" ]; then
         echo "       -> directory libs for mie_BB not present."
         echo "       -> creating it."
         mkdir $path_mie_BB"/libs/"
     fi


     if [ ! -d "outputnc" ]; then
         echo "       -> directory outputnc not present."
	 echo "       -> creating it."
         mkdir outputnc
     fi


     if [ ! -d "outputplt" ]; then
     echo "       -> directory outputplt not present."
	 echo "       -> creating it."
         mkdir outputplt
     fi




# (1) Mie-BB code

function mieBB {
    #
    # Needs three arguments
    #       $1=path of library
    #       $2=path_ecaeropt
    #       $3=date string for log file
    
    logfile=$2"/logs/build_mieBB"$3".log"
    
    cd $1"/src"


    # new_mie_aerosols.f90            => legacy version
    # new_mie_aerosols_struc.f90      => new version structured
    # new_mie_aerosols_openmp_vec.f90 => dev version with OpemMP

    f2py3 -c -m mie_BB parkind1.F90 new_mie_aerosols.f90 interp_ri.f90 &> $logfile
    #f2py3 -c -m mie_BB parkind1.F90 new_mie_aerosols_struc.f90 interp_ri.f90 &> $logfile
    #f2py3 --f90flags=-fopenmp -lgomp -c -m mie_BB parkind1.F90 new_mie_aerosols_openmp_dev.f90 interp_ri.f90 &> $logfile

    status=$?

    case $status in
         0)
         echo "      -> Successfully built mie_BB in current computer."
         mv mie_BB*.so $path_mie_BB/libs/
         ;;
         127)
         echo "      -> There was a problem."
         echo "         Please check that command f2py3 is installed this computer."
         ;;
         *)
         echo "      -> There was a problem."
         echo "         Please check the output log file :"
	 echo "         "$logfile
         ;;
    esac
    cd $2

}
echo ""
echo "..... Building the library mie_BB [using f2py3]"
mieBB $path_mie_BB $path_ecaeropt $date_now
echo  ""

