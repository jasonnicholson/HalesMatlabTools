import os
import sys
from datetime import datetime

# Add project source to sys.path for autodoc/matlabdomain
sys.path.insert(0, os.path.abspath(os.path.join('..', 'source')))

project = 'HalesMatlabTools'
author = 'Jason Nicholson'
current_year = datetime.now().year
copyright = f'{current_year}, {author}'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages',
    'sphinxcontrib.matlab',
]

master_doc = 'index'
source_suffix = {'.rst': 'restructuredtext'}

html_theme = 'alabaster'
html_static_path = ['_static']

autodoc_member_order = 'bysource'

# matlab domain configuration
matlab_src_dir = os.path.abspath(os.path.join('..', 'source'))
matlab_short_links = True
