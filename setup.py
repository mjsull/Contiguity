# cx_Freeze setup file

import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["khmer"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "Contiguity",
        version = "1.0.4",
        description = "Assembly graph construction and visualisation.",
        options = {"build_exe": build_exe_options},
        executables = [Executable("Contiguity.py", base=base)])
