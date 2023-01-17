#!/usr/bin/env sh

tmp="tmp/"

testname=$(basename "$1" .nc)
refsname=$(basename "$2" .nc)

testname_noangle=$tmp$testname"_noangle.nc"
refsname_noangle=$tmp$refsname"_noangle.nc"

ncks -O -x -v angle,phase_function $1 $testname_noangle
ncks -O -x -v angle,phase_function $2 $refsname_noangle

test_sort=$tmp$testname"_reorder.nc"
refs_sort=$tmp$refsname"_reorder.nc"

ncpdq -O -a wavelength,rel_hum,size_bin,component  $testname_noangle $test_sort
ncpdq -O -a wavelength,rel_hum,size_bin,component  $refsname_noangle $refs_sort

ncdiff $test_sort $refs_sort $tmp$testname"_minus_"$refsname".nc"

ncdump $tmp$testname"_minus_"$refsname".nc" > $3
