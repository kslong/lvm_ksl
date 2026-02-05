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

def sort_autoapi_toctree(index_path):
    """Sort the autoapi index.rst toctree alphabetically."""
    if not os.path.exists(index_path):
        return

    with open(index_path, 'r') as f:
        content = f.read()

    # Split into lines for processing
    lines = content.split('\n')

    # Find and sort the toctree entries
    new_lines = []
    in_toctree = False
    toctree_entries = []
    past_options = False

    for line in lines:
        if '.. toctree::' in line:
            in_toctree = True
            past_options = False
            new_lines.append(line)
        elif in_toctree:
            stripped = line.strip()
            # Check for toctree options (start with :)
            if stripped.startswith(':'):
                new_lines.append(line)
                continue
            # Empty line after options signals start of entries
            if not stripped and not past_options:
                new_lines.append(line)
                past_options = True
                continue
            # If line starts with whitespace and we're past options, it's an entry
            if line and line[0] in (' ', '\t') and stripped:
                toctree_entries.append(line)
            elif not line or not line[0] in (' ', '\t'):
                # End of toctree - sort and add entries
                in_toctree = False
                toctree_entries.sort(key=lambda x: x.strip().lower())
                new_lines.extend(toctree_entries)
                toctree_entries = []
                new_lines.append(line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Don't forget remaining entries at end of file
    if toctree_entries:
        toctree_entries.sort(key=lambda x: x.strip().lower())
        new_lines.extend(toctree_entries)

    # Write back
    with open(index_path, 'w') as f:
        f.write('\n'.join(new_lines))

def sort_before_read(app, docname, source):
    """Sort the autoapi index.rst before Sphinx reads it."""
    if docname == 'api/index':
        index_path = os.path.join(app.srcdir, 'api', 'index.rst')
        sort_autoapi_toctree(index_path)
        # Re-read the sorted content
        with open(index_path, 'r') as f:
            source[0] = f.read()

def setup(sphinx):
    """Sphinx setup hook to sort autoapi index."""
    sphinx.connect('source-read', sort_before_read)

