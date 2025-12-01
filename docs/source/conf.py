# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'lvm_ksl'
copyright = '2025, Knox Long'
author = 'Knox Long'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'autoapi.extension',

        ]

templates_path = ['_templates']
exclude_patterns = []

language = 'python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# Path to your Python source code (relative to conf.py)
autoapi_dirs = ['../../py_progs']  # Adjust this path to your source directory

# Type of documentation to generate
autoapi_type = 'python'

# Where to put the generated documentation (relative to your source dir)
autoapi_root = 'api'

# Optional: Add a template directory if you want custom templates
# autoapi_template_dir = '_templates/autoapi'

# Optional but recommended settings:
autoapi_options = [
    'members',           # Document all members
    'undoc-members',     # Include members without docstrings
    'show-inheritance',  # Show inheritance diagrams
    'show-module-summary',  # Show module summary
    'imported-members',  # Document imported members
]

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
