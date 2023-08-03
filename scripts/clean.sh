##################################################################################
# scripts/clean.sh
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
#################################################################################


echo ""
echo "**************************************************************"
echo "**** Current path:"$(pwd)
echo ""
echo "     Cleaning tmp files"
rm tmp/*
echo "     Cleaning log files"
rm logs/*
echo "     Cleaning previous compiled engine libs"
rm engines/*/libs/*
echo "     For security reasons outputnc should be cleaned manually."
echo ""
echo "**************************************************************"
echo ""

