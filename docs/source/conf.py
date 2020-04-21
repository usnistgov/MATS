# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0,"C:\\Users\\ema3\\Documents\\MATS\\MATS")

# NIST required stuff
def setup(app):
    #NIST header and footer
    app.add_stylesheet('https://pages.nist.gov/nist-header-footer/css/nist-combined.css')
    app.add_javascript('https://pages.nist.gov/nist-header-footer/js/jquery-1.9.0.min.js')
    app.add_javascript('https://pages.nist.gov/nist-header-footer/js/nist-header-footer.js')
    
    #Google analytics
    app.add_javascript('https://dap.digitalgov.gov/Universal-Federated-Analytics-Min.js?agency=NIST&subagency=github&pua=UA-42404149-54&yt=true&exts=ppsx,pps,f90,sch,rtf,wrl,txz,m1v,xlsm,msi,xsd,f,tif,eps,mpg,xml,pl,xlt,c',id='_fed_an_ua_tag')
    
    #LeaveNotice
    app.add_javascript('https://code.jquery.com/jquery.1.12.4.min.js')
    app.add_javascript('https://pages.nist.gov/leaveNotice/js/jquery.leaveNotice-nist.min.js')
    app.add_stylesheet('https://pages.nist.gov/leaveNotice/css/jquery.leaveNotice.css')
    app.add_javascript('leave_notice.js')
    return

# -- Project information -----------------------------------------------------

project = 'MATS'
copyright = ''
author = 'Erin M. Adkins'
html_show_copyright = False

# The full version, including alpha/beta/rc tags
release = '1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.todo', 'sphinx.ext.viewcode', 'sphinx.ext.autodoc','sphinx.ext.napoleon'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']