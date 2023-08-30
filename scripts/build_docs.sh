#  +----------------------------------------------------------------------------------------+
#  | scripts/build_docs.sh                                                                  |
#  |                                                                                        |
#  |   (C) Copyright 2022- ECMWF.                                                           |
#  |                                                                                        |
#  |   This software is licensed under the terms of the Apache Licence Version 2.0          |
#  |   which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                 |
#  |                                                                                        |
#  |   In applying this licence, ECMWF does not waive the privileges and immunities         |
#  |   granted to it by virtue of its status as an intergovernmental organisation           |
#  |   nor does it submit to any jurisdiction.                                              |
#  |                                                                                        |
#  |  Author:                                                                               |
#  |     Ramiro Checa-Garcia. ECMWF                                                         |
#  |                                                                                        |
#  |  Modifications:                                                                        |
#  |     10-Dec-2022   Ramiro Checa-Garcia    1st. version                                  |
#  |                                                                                        |
#  +----------------------------------------------------------------------------------------+


echo ""
echo ".... Building the documentation ...."
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
