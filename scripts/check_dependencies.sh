#/bin/bash

# ================================================================
# Dependencies checking script for engine libraries
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

path_mie_BB=$path_ecaeropt"/engines/mie_Boucher_Bozzo"
path_Tmatrx=$path_ecaeropt"/engines/tmatrix_mischenko"


echo ""
echo "..... Checking directories logs, tmp, outputnc, libs"

     if [ ! -d "logs" ]; then
         echo "       -> directory logs not present."
	 echo "       -> creating it."
         mkdir logs
     fi


     if [ ! -d "tmp" ]; then
         echo "       -> directory tmp not present."
	 echo "       -> creating it."
         mkdir tmp
     fi

     if [ ! -d $path_mie_BB"/libs" ]; then
         echo "       -> directory libs for mie_BB not present."
         echo "       -> creating it."
         mkdir $path_mie_BB"/libs/"
     fi


     if [ ! -d $path_Tmatrx"/libs" ]; then
         echo "       -> directory libs for T-matrix Mischenko not present."
         echo "       -> creating it."
         mkdir $path_Tmatrx"/libs/"
     fi

     if [ ! -d "outputnc" ]; then
         echo "       -> directory outputnc not present."
	 echo "       -> creating it."
         mkdir outputnc
     fi

function cmd_present {
 
    eval $1 &> tmp/checking_deps.log

    status=$?
    case $status in
         0)
         echo "                 * "$1" it's detected!"
	 ;;
         127)
         echo "                 * Command "$1" might be not present."
	 echo "                 * Check something like "$2
         ;;
         *)
         echo "                 * There was a problem running "$1
         echo "                 * Please check configuration of "$1
         ;;
    esac

}

echo ""
echo "..... Checking dependencies on external tools & libs"
echo ""
echo "      -> f2py3   ... "; cmd_present "f2py3"       "module load python3"
echo "      -> nco     ... "; cmd_present "which ncks"  "module load nco"
echo "      -> gfortran... "; cmd_present "gfortran -v" "module load gfotran?"
echo ""
rm tmp/checking_deps.log



