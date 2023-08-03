#/bin/bash

################################################################################
# scripts/check_dependencies.sh
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
#  Modifications:
#     10-Dec-2022   Ramiro Checa-Garcia    1st. version
#
#  Info: 
#       * Dependencies checking script for engine libraries
#       * Code is organized with functions,
#################################################################################

#
# (0) Variables needed for this script.
#     We assume that user did not change directories structure.

path_ecaeropt=$(pwd)

path_mie=$path_ecaeropt"/engines/mie_code"
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

     if [ ! -d $path_mie"/libs" ]; then
         echo "       -> directory libs for Mie code not present."
         echo "       -> creating it."
         mkdir $path_mie"/libs/"
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



