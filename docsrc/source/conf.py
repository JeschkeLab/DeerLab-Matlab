# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------

import sys
import os
import re

sys.path.append(os.path.abspath('.'))

# Import Read-the-Docs (RTD) theme
from sphinx.locale import _
from sphinx_rtd_theme import __version__

# Project details
project = 'DeerLab'
copyright = '2019-2020, Luis Fábregas Ibáñez, Stefan Stoll, and others'
author = 'Fabregas Ibanez'
language = 'en'

# Add sphinx extensions
extensions = [
    'sphinx.ext.intersphinx',
    'sphinxcontrib.matlab',
    'sphinx.ext.viewcode',
    'sphinxcontrib.httpdomain',
    'sphinx.ext.imgmath',
]

# Render Latex math equations as svg instead of rendering with JavaScript
imgmath_image_format = 'svg'
imgmath_dvisvgm = 'dvisvgm'
imgmath_latex_preamble = r'''
\newcommand{\mr}[1]{\mathrm{#1}}
\newcommand{\mx}[1]{\boldsymbol{#1}}
\newcommand{\vc}[1]{\boldsymbol{#1}}
\DeclareMathOperator*{\argmin}{\arg\!\min}
'''

# Setup template stuff
templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'
intersphinx_mapping = {
    'rtd': ('https://docs.readthedocs.io/en/latest/', None),
    'sphinx': ('http://www.sphinx-doc.org/en/stable/', None),
}
html_theme = 'sphinx_rtd_theme'

# Integrate version control system
# -------------------------------------------------------------
html_context = {
    "display_github": False, # Integrate GitHub
    "github_user": "JeschkeLab", # Username
    "github_repo": "DeerLab", # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/source/", # Path in the checkout to the docs root
}


# Read-the-Docs options configuration
# -------------------------------------------------------------
html_theme_options = {
    'sticky_navigation': False,
    'titles_only': False,
    'collapse_navigation': False,
    'logo_only': True,
    'navigation_depth':1,
    #'style_nav_header_background': '#2a7bf8'
}
html_copy_source = False
html_theme_path = ["../.."]
html_logo = "demo/static/logo-wordmark-light.svg"
html_show_sourcelink = True
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']

# Extensions to theme docs
def setup(app):
    from sphinx.domains.python import PyField
    from sphinx.util.docfields import Field
    app.add_css_file('/source/_static/custom.css')
    app.add_stylesheet('/source/_static/theme_override.css')
    app.add_object_type(
        'confval',
        'confval',
        objname='configuration value',
        indextemplate='pair: %s; configuration value',
        doc_field_types=[
            PyField(
                'type',
                label=_('Type'),
                has_arg=False,
                names=('type',),
                bodyrolename='class'
            ),
            Field(
                'default',
                label=_('Default'),
                has_arg=False,
                names=('default',),
            ),
        ]
    )
    # Patch the MATLAB lexer to correct wrong highlighting
    from sphinx.highlighting import lexers
    from pygments.lexers import MatlabLexer
    from pygments.token import Name
    from pygments.filters import NameHighlightFilter
    matlab_lexer = MatlabLexer()
    matlab_lexer.add_filter(NameHighlightFilter(
            names=['function', 'end', 'if', 'else', 'elseif', 'switch', 'case', 'return', 'otherwise'],
            tokentype=Name.Function,
            ))
    app.add_lexer('matlab', matlab_lexer)


# These folders are copied to the documentation's HTML output
html_static_path = ['_static']

# Add path to custom CSS file to overwrite some of the default CSS settings
html_css_files = [
    'custom.css',
    'table_util.css',
    'table_main.css',
    'theme_override.css'
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Default role
default_role = 'math'  # with this, :math:`\psi` can be written simply as `\psi`


# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes", ]
html_static_path = ['_static']
highlight_language = 'matlab'
primary_domain = 'mat'
html_logo = '_static/logo.png'


# Design pygments patch for MATLAB true code highlighting
# --------------------------------------------------------

# Import the pygments python library
from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, Number, Operator, Generic, Text, Other, Comment, Whitespace

# Define custom style for the MATLAB highlighting

class MyFancyStyle(Style):
    background_color = "#f7f7f7"
    default_style = "default"
    styles = {
        Comment:                'italic #22924a',
        Keyword:                '#0000FF',
        Name:                   '#000',
        Name.Function:          'noinherit #0000FF',
        Name.Class:             'noinherit #0000FF',
        String:                 '#a022f0',
    }

# Create patch for applying the style 	
def pygments_monkeypatch_style(mod_name, cls):
    import sys
    import pygments.styles
    cls_name = cls.__name__
    mod = type(__import__("os"))(mod_name)
    setattr(mod, cls_name, cls)
    setattr(pygments.styles, mod_name, mod)
    sys.modules["pygments.styles." + mod_name] = mod
    from pygments.styles import STYLE_MAP
    STYLE_MAP[mod_name] = mod_name + "::" + cls_name

pygments_monkeypatch_style("my_fancy_style", MyFancyStyle)
pygments_style = "my_fancy_style"

