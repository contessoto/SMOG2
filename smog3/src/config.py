import math

VERSION = "2.7beta"
PI = 3.14159265358979

# OpenSMOG restricted names
OS_RESTRICT = {
    name: 0 for name in [
        "q1", "q2", "theta", "r", "r_c", "i", "j", "k", "l", "m", "n",
        "type1", "type2", "sqrt", "exp", "log", "sin", "cos", "sec", "csc",
        "tan", "cot", "asin", "acos", "atan", "sinh", "cosh", "tanh", "erf",
        "erfc", "min", "max", "abs", "floor", "ceil", "step", "delta", "select", "null"
    ]
}

# Global variables to track state
NB_TYPES_PRESENT = {}

import os
import sys

# Determine the base directory where templates are located.
# If running from source (pip install -e . or direct execution), it might be relative to this file.
# If installed, it should be in the package data or a shared location.

def get_template_path():
    # Try to find the 'share' directory relative to the package installation
    # This assumes structure:
    #   <prefix>/lib/pythonX.Y/site-packages/smog3/src/config.py
    #   <prefix>/share/smog3/templates
    # OR
    #   <source>/smog3/src/config.py
    #   <source>/share/templates

    # Check for local source structure first
    this_dir = os.path.dirname(os.path.abspath(__file__))
    source_share = os.path.abspath(os.path.join(this_dir, '..', '..', 'share', 'templates'))
    if os.path.isdir(source_share):
        return source_share

    # Check for installed package structure (sys.prefix)
    # Typically <sys.prefix>/share/smog3/templates
    installed_share = os.path.join(sys.prefix, 'share', 'smog3', 'templates')
    if os.path.isdir(installed_share):
        return installed_share

    # Fallback: check CWD/share/templates (useful if running from repo root)
    cwd_share = os.path.join(os.getcwd(), 'share', 'templates')
    if os.path.isdir(cwd_share):
        return cwd_share

    # Fallback to current directory (legacy behavior)
    return os.getcwd()

DEFAULT_TEMPLATE_PATH = get_template_path()
