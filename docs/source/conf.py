# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'KRED'
copyright = '2025, Knox Long & Sean Points'
author = 'Knox Long & Sean Points'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'autoapi.extension',
         'sphinx.ext.viewcode',  # Add [source] links
        ]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

# Path to your Python source code (relative to conf.py)
autoapi_dirs = ['../../py_progs']  # Adjust this path to your source directory

# Type of documentation to generate
autoapi_type = 'python'

# Where to put the generated documentation (relative to your source dir)
autoapi_root = 'api'

# Optional: Add a template directory if you want custom templates
# autoapi_template_dir = '_templates/autoapi'


# AutoAPI configuration
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
    'special-members',
    'imported-members',
]

# Sort members alphabetically within each module
autoapi_member_order = 'alphabetical'  # or 'groupwise' or 'bysource'

# Skip certain files or patterns (optional)
autoapi_ignore = [
    '*/tests/*',
    '*/test_*.py',
    '*/__pycache__/*',
]

# Keep the generated files for inspection (optional, useful for debugging)
autoapi_keep_files = True

# Add any modules that should be mocked (like your missing modules)
autodoc_mock_imports = [
    'lvmdrp',
    'sdss_access',
    # Add other missing dependencies here
]

import os

def sort_autoapi_toctree(app, exception):
    """Sort the autoapi index.rst toctree alphabetically after build."""
    if exception is not None:
        return  # Build failed, don't process
    
    index_path = os.path.join(app.srcdir, 'api', 'index.rst')
    if not os.path.exists(index_path):
        return
    
    with open(index_path, 'r') as f:
        lines = f.readlines()
    
    # Find and sort the toctree entries
    new_lines = []
    in_toctree = False
    toctree_entries = []
    indent = ''
    
    for line in lines:
        if '.. toctree::' in line:
            in_toctree = True
            new_lines.append(line)
        elif in_toctree:
            stripped = line.lstrip()
            if stripped.startswith(':') or not stripped:
                # toctree option or blank line
                new_lines.append(line)
            elif line[0] not in (' ', '\t'):
                # End of toctree
                in_toctree = False
                # Sort and add collected entries
                toctree_entries.sort(key=lambda x: x.strip().lower())
                new_lines.extend(toctree_entries)
                toctree_entries = []
                new_lines.append(line)
            else:
                # This is a toctree entry
                if not indent and line != line.lstrip():
                    indent = line[:len(line) - len(line.lstrip())]
                toctree_entries.append(line)
        else:
            new_lines.append(line)
    
    # Don't forget remaining entries
    if toctree_entries:
        toctree_entries.sort(key=lambda x: x.strip().lower())
        new_lines.extend(toctree_entries)
    
    # Write back
    with open(index_path, 'w') as f:
        f.writelines(new_lines)

def setup(sphinx):
    """Sphinx setup hook to sort autoapi index."""
    sphinx.connect('build-finished', sort_autoapi_toctree)

