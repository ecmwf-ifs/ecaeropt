#!/bin/bash

module load python3

function cleanDir {
  pwd
  echo $1
  if [ -d "$1" ]; then
	  if [ "$(ls -A $1)" ]; then
	    rm -v $1/*
	  fi
  fi
}

cleanDir tmp
cleanDir outputnc

./ecaeropt -s settings/IFS_CY48R1.toml


