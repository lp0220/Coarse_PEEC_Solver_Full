import math
import os

# Physical constants
E0 = 8.85418782e-12       # Vacuum permittivity (F/m)
U0 = 4 * math.pi * 1e-7   # Vacuum permeability (H/m)

# Unit scale: input coordinates in mils, convert to meters
# 1 meter = 39370 mils
UNIT_SCALE = 39370

# L/P computation parameters
DEFAULT_GAUSS_ORDER = 3
DEFAULT_DIV_X = 1
DEFAULT_DIV_Y = 1

# Directory paths
INPUT_DIR = 'input'
WORKSPACE_DIR = 'workspace'

# Improved packaged solver executable (relative to this project directory)
CPP_SOLVER_EXE = os.path.join('..', 'Solver_exe - Improve', 'Quasi_Static_Solver.exe')
