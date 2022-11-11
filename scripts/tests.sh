
echo "**************************************************************"
echo "**   Running several tests to compare with stored netcdf files"
echo "**    ... tests will take about 30 min.                       "
echo "**    ... please see info provided at the end                 "
echo "**    ... tests include one hydrophobic/insolube aerosol      "
echo "**    ...               one hydrophilic/solube   aerosol      "
echo "**    ..                one externally mixed     aerosol      "
./ecaeropt -t tests/test_single_mixed_original.toml
echo "**** DONE"

