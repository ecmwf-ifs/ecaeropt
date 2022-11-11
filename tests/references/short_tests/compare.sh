#!/usr/bin/env sh

tmp="../../../tmp/"
testname=$(basename "$1" .nc)
refsname=$(basename "$2" .nc)
test_sort=$tmp$testname"_reorder.nc"
refs_sort=$tmp$refsname"_reorder.nc"
ncpdq -O -a wavelength,rel_hum,size_bin,angle  $1 $test_sort
ncpdq -O -a wavelength,rel_hum,size_bin,angle  $2 $refs_sort
ncdiff $test_sort $refs_sort $tmp$testname"_minus_"$refsname".nc"

