import numpy as np
import os


def load_node_data(filepath, scale):
    """
    Load node data.
    Format: first line = node count, then either:
      'x y z size'
      'x y z size distance_modifier'
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError("Node file not found: " + filepath)

    with open(filepath, 'r') as f:
        lines = f.readlines()
        n_n = int(lines[0].split()[-1])
        points = np.zeros((n_n, 3))
        node_sizes = np.zeros(n_n)
        d_mod = None
        has_modifier = None

        for i, line in enumerate(lines[1:n_n + 1]):
            parts = line.split()
            if len(parts) not in (4, 5):
                raise ValueError("Node row {} must have 4 or 5 columns: {}".format(i + 1, line.strip()))
            row_has_modifier = len(parts) == 5
            if has_modifier is None:
                has_modifier = row_has_modifier
            elif row_has_modifier != has_modifier:
                raise ValueError("Node data cannot mix 4-column and 5-column rows")

            points[i] = np.array(parts[:3], dtype='float64') / scale
            node_sizes[i] = float(parts[3]) / scale

            if row_has_modifier:
                if d_mod is None:
                    d_mod = np.zeros(n_n)
                d_mod[i] = float(parts[4])

    return points, node_sizes, d_mod


def load_branch_data(filepath, node_sizes):
    """
    Load branch data.
    Format: first line = branch count, then 'node1 node2' (1-based) per line.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError("Branch file not found: " + filepath)

    with open(filepath, 'r') as f:
        lines = f.readlines()
        n_b = int(lines[0].split()[-1])
        connects = np.zeros((n_b, 2), dtype="int")
        branch_sizes = np.zeros(n_b)

        for i, line in enumerate(lines[1:n_b + 1]):
            parts = line.split()
            n1_idx, n2_idx = int(parts[0]), int(parts[1])
            connects[i] = [n1_idx, n2_idx]
            branch_sizes[i] = (node_sizes[n1_idx - 1] + node_sizes[n2_idx - 1]) / 2.0

    return connects, branch_sizes


def save_matrix_cpp_format(matrix, filepath):
    """Save matrix in C++ solver format: first line 'rows cols', then data rows."""
    rows, cols = matrix.shape
    with open(filepath, 'w') as f:
        f.write("{} {}\n".format(rows, cols))
        for i in range(rows):
            line = " ".join("{:.18e}".format(matrix[i, j]) for j in range(cols))
            f.write(line + "\n")


def save_b2n(connects, filepath):
    """Save branch-to-node connectivity in 0-based format for C++ solver."""
    n_b = len(connects)
    with open(filepath, 'w') as f:
        f.write("{} 2\n".format(n_b))
        for c in connects:
            f.write("{} {}\n".format(c[0] - 1, c[1] - 1))


def convert_port_file(port_file, output_file):
    """
    Convert Port.txt from user format (1-based node indices)
    to C++ format (0-based node indices).
    """
    with open(port_file, 'r') as f:
        lines = f.readlines()

    n_port = int(lines[0].strip())

    with open(output_file, 'w') as f:
        f.write("{}\n".format(n_port))
        for line in lines[1:n_port + 1]:
            parts = line.split()
            n1 = int(parts[0])
            n2 = int(parts[1])
            zin = parts[2]
            f.write("{} {} {}\n".format(n1, n2, zin))
