#!/bin/bash

################################################################################
# scripts/create_ifs_prepare.sh
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
#       * Code is organized with functions
#
#################################################################################


function cleanDir {
  if [ -d "$1" ]; then
	  if [ "$(ls -A $1)" ]; then
	    rm  $1/*.nc
      rm  $1/*.tmp
      rm  $1/*.log
	  fi
  fi
}

function createDir {
  if [ -d "$1" ]; then
    echo $1" present"
  else
    mkdir $1
    echo $1" created"
  fi
}

cleanDir  tmp
cleanDir  outputnc
createDir outputnc/store_ifsnc



