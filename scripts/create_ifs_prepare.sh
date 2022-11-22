#!/bin/bash

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



