# Coarse PEEC Solver

Coarse PEEC Solver is an integrated PEEC workflow for extracting partial
inductance and coefficient-of-potential matrices from a coarse node/branch
model, then running a packaged C++ frequency-domain solver to produce S, Z,
and Y parameter results.

The Python side reads geometry and port definitions, computes the PEEC `L` and
`P` matrices, writes the input files expected by the C++ solver, stages the
packaged executable into a workspace-local runtime directory, and copies the
solver outputs back into the project workspace.

## Project Layout

```text
.
|-- Integrated_PEEC_Solver/
|   |-- run.py                 # Main Python entry point
|   |-- config.py              # Constants, paths, integration settings
|   |-- core/                  # PEEC kernels and self-inductance integration
|   |-- utils/                 # Geometry and file conversion helpers
|   |-- input/                 # Input model, material, and solver settings
|   `-- workspace/             # Generated intermediate and output files
|-- Solver_exe - Improve/
|   |-- Quasi_Static_Solver.exe
|   `-- Data/                  # Packaged solver data directory
|-- coarse solver documentation (.docx)
`-- simulation comparison images
```

## Requirements

- Windows, because the bundled frequency-domain solver is
  `Quasi_Static_Solver.exe`.
- Python 3.10 or newer.
- Python packages:
  - `numpy`
  - `matplotlib` for the optional plotting script
- No separate Intel MKL runtime configuration is required for the bundled
  `Quasi_Static_Solver.exe` in this repository.

## Quick Start

From the Python solver directory, run:

```powershell
cd Integrated_PEEC_Solver
python run.py
```

The current entry point is configured to read this sample model:

```text
input/68/Node_mod_sqrt1.txt
input/68/Branch.txt
input/68/Port.txt
```

To run another model, update the three input paths in `run.py`.

## Input Files

### Node File

The first line is the node count. Each following row must use one of these
formats:

```text
x y z b
x y z b c
```

- `x y z`: node coordinates.
- `b`: node characteristic size.
- `c`: optional self-potential distance modifier.

Coordinates and `b` are divided by `config.UNIT_SCALE`. The default assumes
input values are in mils:

```python
UNIT_SCALE = 39370
```

A single node file must not mix four-column and five-column rows. The fifth
column currently affects only the self coefficient of potential.

### Branch File

The first line is the branch count. Each following row is:

```text
node_start node_end
```

Branch nodes are 1-based in the branch input. The generated `B2N.txt` file is
converted to 0-based indexing for the C++ solver. Each branch width is computed
as the average characteristic size of its two endpoint nodes.

### Port File

The first line is the port count. Each following row is:

```text
node_pos node_neg Zin
```

`Zin` is the reference impedance. The current workflow preserves the port node
indices as provided in the input file; the bundled sample uses solver-side
0-based node indices.

### Additional Solver Inputs

These files are copied from `Integrated_PEEC_Solver/input/` into the generated
solver runtime directory when present:

```text
set.txt
DIELECTRIC.txt
Source.txt
SourceConfig.txt
```

Common fields in `set.txt` include:

- `FS`: start frequency.
- `FE`: stop frequency.
- `N_FP`: number of frequency points.
- `DIM`: compatibility parameter used by the C++ solver.

## Computation Flow

1. Load node, branch, port, dielectric, and frequency settings.
2. Compute the coefficient-of-potential matrix `P` in Python.
3. Compute the partial inductance matrix `L` in Python.
4. Write `L`, `P`, branch-to-node connectivity, and port files into
   `Integrated_PEEC_Solver/workspace/`.
5. Create `workspace/_solver_runtime/`.
6. Copy `Solver_exe - Improve/Data` and `Quasi_Static_Solver.exe` into the
   runtime directory.
7. Copy matrix and configuration files into the fixed paths expected by the
   packaged solver.
8. Execute the C++ frequency-domain solver.
9. Copy result files back into `Integrated_PEEC_Solver/workspace/`.

## PEEC Extraction Notes

### Coefficient-of-Potential Matrix

The Python code computes the coefficient-of-potential matrix `P`. The C++
solver reads `P` and inverts it internally to obtain the equivalent capacitance
matrix.

- Self terms use the node characteristic size `b`.
- If the node file has a fifth column, the self term uses both `b` and `c`.
- Mutual terms are based on node-to-node distance in the current code.

### Partial Inductance Matrix

Each branch is represented as a rectangular region derived from its endpoint
nodes.

- If the branch is mainly aligned with `x`, length is assigned to `x` and
  width to `y`.
- If the branch is mainly aligned with `y`, length is assigned to `y` and
  width to `x`.
- Self inductance uses the `A1` term from the implemented integration routine.
  `A2` and `A3` are computed by the helper but are not added in the current
  main workflow.
- Mutual inductance uses a center-line approximation and is capped against the
  geometric mean of the corresponding self-inductance terms.

## Generated Files

The Python workflow writes these intermediate solver inputs:

```text
Integrated_PEEC_Solver/workspace/L.txt
Integrated_PEEC_Solver/workspace/P.txt
Integrated_PEEC_Solver/workspace/LL_OO.txt
Integrated_PEEC_Solver/workspace/PP_OO.txt
Integrated_PEEC_Solver/workspace/B2N.txt
Integrated_PEEC_Solver/workspace/PORT.txt
```

The main result files are:

```text
Integrated_PEEC_Solver/workspace/map_Final.txt   # S parameters
Integrated_PEEC_Solver/workspace/Z_in_Final.txt  # Z parameters
Integrated_PEEC_Solver/workspace/Y_in_Final.txt  # Y parameters
```

The packaged solver may also produce:

```text
map.txt
Z_in.txt
Y_in.txt
output_voltage.txt
1.png
```

## Plotting

After a run, the optional plotting script in `Integrated_PEEC_Solver/workspace`
can plot selected S-parameter columns from `map.txt`:

```powershell
cd Integrated_PEEC_Solver\workspace
python Plot.py
```

## Current Limitations

- The sample model paths are hard-coded in `Integrated_PEEC_Solver/run.py`.
- The fifth node column currently affects only self potential terms.
- Self inductance currently writes only the `A1` term into `L[i, i]`.
- The packaged solver expects fixed relative paths under `./Data/Output`, so
  `run.py` creates a workspace-local runtime tree before execution.

## Reference Material

The original Chinese usage document is included in the repository as a Word
document, together with simulation comparison images for the coarse PEEC
examples.
