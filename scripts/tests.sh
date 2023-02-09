
echo ""
echo "    ******************************************************************************"
echo "    **  TESTs"
echo "    **        Calculating several aerosols to compare with stored netcdf files"
echo "    **    ... tests will take about 10 min.                       "
echo "    **    ... please see info provided at the end                 "
echo "    **    ... tests include one hydrophobic/insolube aerosol      "
echo "    **    ...               one hydrophilic/solube   aerosol      "
echo "    **    ...               one externally mixed     aerosol      "
echo "    **"
echo "    ******************************************************************************"

rm tmp/*
rm outputnc/*

./ecaeropt -t tests/settings/test_single_mixed_original.toml
#echo "**** DONE"

