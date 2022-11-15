

echo ""
echo "... Building the documentation ..."
cd docs
echo "    * documentation in html "
make html &> ../logs/build_docum.log 
cd ..
rm docs.html
ln -s docs/build/html/index.html docs.html
echo "    * created a link named docs.html"
echo ""
