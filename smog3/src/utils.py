import sys
import re
import os
import math
import numpy as np
from .config import VERSION

class SmogError(Exception):
    pass

def smog_quit(message, warn=False):
    """
    Exits the program with an error message or prints a warning.
    """
    if warn:
        print(f"\nWARNING: {message}\n")
    else:
        print(f"\n\nFATAL ERROR: {message}\n\nFor more information about specific errors, you can check the FAQ page on smog-server.org,\nthe SMOG 2 manual, or you can email us at info@smog-server.org.\n")
        sys.exit(1)

def smog_note(message):
    print(f"\nNOTE: {message}\n")

def check_comment(line):
    """
    Splits a line into content and comment.
    """
    line = line.strip()
    if ';' in line:
        content, comment = line.split(';', 1)
        return content.strip(), ';' + comment
    return line, ""

def has_content(line):
    """
    Checks if a line has non-comment content.
    """
    line, _ = check_comment(line)
    return len(line) > 0

def trim(s):
    return s.strip()

def what_am_i(s):
    """
    Returns 1 for integer, 2 for float, 3 for other.
    """
    if re.match(r'^[+-]?[0-9]+[0-9,eE+-]*$', s):
        return 1
    if re.match(r'^[+-]?[0-9]*\.[0-9]*[eE]?[+-]?[0-9]*$', s): # Simplified regex
        return 2
    return 3

def check_suffix(name, suf):
    if not name.endswith(suf):
        return name + suf
    return name

def check_already_exists(filen):
    maxbu = 10
    if os.path.exists(filen):
        for bu in range(1, maxbu + 1):
            buname = f"{filen}.bu{bu}"
            if not os.path.exists(buname):
                print(f"{filen} already exists. Backing up to {buname}")
                os.rename(filen, buname)
                return
        smog_quit(f"Already backed up {maxbu} copies of {filen}.")

def eval_sub(expression, substitutions):
    """
    Evaluates an expression with variable substitutions.
    substitutions: dict of var -> value
    """
    if not isinstance(substitutions, dict):
        # Backward compat if value passed directly (assumed for '?')
        # But wait, logic changed in bonded.py to pass dict.
        pass

    # Create safe context
    context = {"__builtins__": None, "sqrt": math.sqrt, "exp": math.exp, "log": math.log, "sin": math.sin, "cos": math.cos, "tan": math.tan, "abs": abs, "min": min, "max": max, "pi": math.pi, "PI": math.pi}

    # Add substitutions to context
    context.update(substitutions)

    try:
        return eval(expression, context)
    except Exception as e:
        # print(f"Eval error: {e}")
        return None

# Vector utilities
def norm(vec):
    return np.linalg.norm(vec)

def dot(vec1, vec2):
    return np.dot(vec1, vec2)

def cross(vec1, vec2):
    return np.cross(vec1, vec2)

def sub(vec1, vec2):
    return vec1 - vec2

def add(vec1, vec2):
    return vec1 + vec2

def mult(vec, scalar):
    return vec * scalar
