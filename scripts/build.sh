#/bin/bash

################################################################################
# scripts/build_docs.sh
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
#       * Code is organized with functions, each function is doing one build 
#         for external optical libraries.
#
#################################################################################

# (0) Variables needed for this script.
#     We assume that user did not change directories structure.

path_ecaeropt=$(pwd)
date_now=$(date "+%Y%m%d_%H%M")

path_mie=$path_ecaeropt"/engines/mie_code/"
path_Tmatrx=$path_ecaeropt"/engines/tmatrix_mischenko"


echo ""
echo "..... Checking directories logs and libs"

     if [ ! -d "logs" ]; then
         echo "       -> directory logs not present."
	 echo "       -> creating it."
         mkdir logs
     fi

     if [ ! -d $path_mie"/libs" ]; then
         echo "       -> directory libs for Mie code not present."
         echo "       -> creating it."
         mkdir $path_mie"/libs/"
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




# (1) Mie-code

function mielib {
    #
    # Needs three arguments
    #       $1=path of library
    #       $2=path_ecaeropt
    #       $3=date string for log file
    
    logfile=$2"/logs/build_mie_code"$3".log"
    
    cd $1"/src"

    # mie_aerosols_struc.f90      => new version structured

    f2py3 -c -m mie parkind1.F90 mie_aerosols_struc.f90 interp_ri.f90 &> $logfile

    status=$?

    case $status in
         0)
         echo "      -> Successfully built Mie library in current computer."
         mv mie*.so $path_mie/libs/
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
echo "..... Building the Mie code library [using f2py3]"
mielib $path_mie $path_ecaeropt $date_now
echo  ""
echo "..... Giving exec. permission to ecaeropt and comparing scripts"
chmod +x ecaeropt
chmod +x aeropt/compare_secure.sh
chmod +x aeropt/compare_secure_noangle.sh
