

#  +------------------------------------------------------------------------------------+
#  | docs/source/conf.py                                                                |
#  |                                                                                    |
#  | (C) Copyright 2022- ECMWF.                                                         |
#  |                                                                                    |
#  | This software is licensed under the terms of the Apache Licence Version 2.0        |
#  | which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.               |
#  |                                                                                    |
#  | In applying this licence, ECMWF does not waive the privileges and immunities       |
#  | granted to it by virtue of its status as an intergovernmental organisation         |
#  | nor does it submit to any jurisdiction.                                            |
#  |                                                                                    |
#  | Author:                                                                            |
#  |    Ramiro Checa-Garcia. ECMWF                                                      |
#  |                                                                                    |
#  | Modifications:                                                                     |
#  |    10-Dec-2022   Ramiro Checa-Garcia    1st version                                |
#  |                                                                                    |
#  | Info:                                                                              |
#  |                                                                                    |
#  +------------------------------------------------------------------------------------+


# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ecaeropt'
copyright = '2022- ECMWF. (Under Apache License Version 2.0)'
author = 'Ramiro Checa-Garcia'
release = '1.16'
master_doc = 'index'
source_suffix = '.rst'
sys.path.append("/home/parc/_ecmwf/codes/ecaeropt/")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax'
]


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'nature'
#html_static_path = ['_static']


# -- Options for LaTeX output -------------------------------------------------

latex_engine = 'pdflatex'
latex_theme  = 'manual'
latex_table_style = ['booktabs', 'colorrows']

latex_documents = [
    (master_doc, 'main.tex', 'EC-AER-OPT offline tool for aerosols optical properties',
     'Ramiro Checa-Garcia', 'report') ]

latex_logo = "ECMWF_logo.png"
latex_elements = {
        'atendofbody': 'ECWMF',
        'papersize': 'a4paper',
        'releasename':" v1.16 ECMWF",
        'pointsize': '10pt'    }

