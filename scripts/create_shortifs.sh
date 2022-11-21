#!/bin/bash

module load python3

rm tmp/*
rm outputnc/*
rm test_ifs.nc

./ecaeropt -s settings/short_ifs_example.toml


