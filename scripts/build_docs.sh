

echo ""
echo "... Building the documentation ..."
cd docs
echo "    * documentation in html "
make html     &>   ../logs/build_docum.log 
echo "    * documentation in pdf  "
make latexpdf &>   ../logs/build_docum_pdf.log
cd ..
rm docs.html
ln -s docs/build/html/index.html docs.html
echo "    * created a link named docs.html"
rm manual_ecaeropt.pdf
ln -s docs/build/latex/main.pdf manual_ecaeropt.pdf
echo "    * created a link named manual_ecaeropt.pdf"
echo ""
