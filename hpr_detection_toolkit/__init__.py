# __init__.py
"""
This is the initialization file for your Python package.  It is
executed when the package is imported.
"""

# Optional:  Import commonly used functions or classes for easier access.
# Example:
# from .module1 import my_function
# from .module2 import MyClass

# Optional: Define package-level variables.  Good for versioning or
# metadata.
__author__ = "Bob Oeyen"
__email__ = "bob.oeyen@outlook.com"

__version__ = "0.1.0"  Important for package management

'''try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

# obtain version number from pyproject.toml (developer version)
parent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
toml_file = os.path.join(parent_dir, 'pyproject.toml')
if os.path.isfile(toml_file):
    import toml
    toml_dict = toml.load(toml_file)
    try:
        if toml_dict['tool']['poetry']['name'] == "ToolkitDetectionHPR": # check this is the right pyproject.toml
            __version__ = toml_dict['tool']['poetry']['version']
    except KeyError:
        pass'''

# Optional:  Perform any necessary package-level initialization.
# Be careful with this; keep it lean.  Avoid lengthy operations here.
# Example:

import logging
logging.basicConfig(level=logging.INFO)
#logging.setLoggerClass(ToolkitHPRLogger)
#_setup_logger(name="ToolkitDetectionHPR")

# Optional:  If your package has subpackages, you might want to
# import them here to make them available at the top level.  This
# is less common in modern Python, but can be useful in some cases.
# Example:
# from . import subpackage1
# from . import subpackage2

# Best Practice:  Explicitly control what gets imported with
# from my_package import *
#__all__ = [
#    'my_function',  # If you imported it above
#    'MyClass',      # If you imported it above
    # 'subpackage1',  # If you have subpackages and want them in *
    #  ... any other names you want exposed at the package level
#]
