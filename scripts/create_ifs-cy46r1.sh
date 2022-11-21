#!/bin/bash

module load python3

function cleanDir {
  if [ -d "$1" ]; then
	  if [ "$(ls -A $1)" ]; then
	    rm $1"/*"
	  fi
  fi
}

cleanDir tmp
cleanDir outputnc/*

./ecaeropt -s settings/IFS_CY46R1.toml


