
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

