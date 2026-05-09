"""
Integrated PEEC Solver
======================
Combines Python L/P matrix extraction with C++ PEEC frequency domain solver.

Input:  input/Node.txt, input/Branch.txt, input/Port.txt
Output: workspace/ directory (S/Z/Y parameter result files)

Usage:
    python run.py
"""

import os
import sys
import shutil
import subprocess
import math
import numpy as np

import config
from core.kernels import cal_p_self, cal_p_oth, cal_l_oth_approx
from core.solver import PEECSolver
from utils.geometry import get_branch_info
from utils.io_handler import (load_node_data, load_branch_data,
                               save_matrix_cpp_format, save_b2n,
                               convert_port_file)


def cap_mutual_by_self_terms(mutual, self_l1, self_l2):
    self_product = self_l1 * self_l2
    if self_product > 0.0 and mutual * mutual > self_product:
        return math.copysign(0.9 * math.sqrt(self_product), mutual)
    return mutual


def prepare_solver_runtime(workspace, solver_exe):
    """
    Stage the improved packaged solver into a workspace-local runtime tree.
    The improved executable reads fixed paths under ./Data/Output.
    """
    runtime_dir = os.path.join(workspace, '_solver_runtime')
    data_dir = os.path.join(runtime_dir, 'Data')
    packaged_data_dir = os.path.join(os.path.dirname(solver_exe), 'Data')

    if os.path.isdir(runtime_dir):
        shutil.rmtree(runtime_dir)

    runtime_exe = os.path.join(runtime_dir, os.path.basename(solver_exe))
    shutil.copytree(packaged_data_dir, data_dir)
    shutil.copy2(solver_exe, runtime_exe)

    output_dir = os.path.join(data_dir, 'Output')

    # Matrix and topology inputs expected by the improved executable.
    output_inputs = {
        'LL_OO.txt': 'LL_OO.txt',
        'PP_OO.txt': 'PP_OO.txt',
        'B2N.txt': 'B2N.txt',
        'PORT.txt': 'PORT.txt',
    }
    for src_name, dst_name in output_inputs.items():
        shutil.copy2(os.path.join(workspace, src_name), os.path.join(output_dir, dst_name))

    # Solver configuration lives one level above Output.
    for fname in ['set.txt', 'DIELECTRIC.txt', 'Source.txt', 'SourceConfig.txt']:
        src = os.path.join(workspace, fname)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(data_dir, fname))

    return runtime_dir, runtime_exe


def collect_solver_outputs(runtime_dir, workspace):
    """
    Copy improved solver outputs back into the workspace, preserving both
    the packaged names and the legacy *_Final names used by current scripts.
    """
    output_dir = os.path.join(runtime_dir, 'Data', 'Output')
    output_aliases = {
        'map.txt': ['map.txt', 'map_Final.txt'],
        'Y_in.txt': ['Y_in.txt', 'Y_in_Final.txt'],
        'Z_in.txt': ['Z_in.txt', 'Z_in_Final.txt'],
        'output_voltage.txt': ['output_voltage.txt'],
        '1.png': ['1.png'],
    }

    for src_name, dst_names in output_aliases.items():
        src = os.path.join(output_dir, src_name)
        if not os.path.exists(src):
            continue
        for dst_name in dst_names:
            shutil.copy2(src, os.path.join(workspace, dst_name))


def main():
    print("=" * 60)
    print("  Integrated PEEC Solver")
    print("  Node/Branch -> L/P Extraction -> Circuit Simulation")
    print("=" * 60)

    # Create workspace directory
    os.makedirs(config.WORKSPACE_DIR, exist_ok=True)

    input_dir = config.INPUT_DIR
    workspace = config.WORKSPACE_DIR

    node_file = os.path.join(input_dir, '68/Node_mod_sqrt1.txt')
    branch_file = os.path.join(input_dir, '68/Branch.txt')
    port_file = os.path.join(input_dir, '68/Port.txt')

    # Verify input files exist
    for fpath, name in [(node_file, '68/Node_mod_sqrt1.txt'), (branch_file, '68/Branch.txt'),
                         (port_file, '68/Port.txt')]:
        if not os.path.exists(fpath):
            print("  Error: {} not found at {}".format(name, os.path.abspath(fpath)))
            return

    # ================================================================
    # Step 1: Compute L and P Matrices
    # ================================================================
    print("\n[Step 1/3] Computing L and P matrices...")
    scale = config.UNIT_SCALE
    points, node_sizes, d_mod = load_node_data(node_file, scale)
    connects, branch_sizes = load_branch_data(branch_file, node_sizes)

    num_nodes = len(points)
    num_branches = len(connects)
    print("  Mesh: {} nodes, {} branches".format(num_nodes, num_branches))

    # P matrix (coefficient of potential)
    print("  Computing P matrix...")
    P = np.zeros((num_nodes, num_nodes))
    for i in range(num_nodes):
        p_self_modifier = None if d_mod is None else d_mod[i]
        P[i, i] = cal_p_self(node_sizes[i], p_self_modifier)
        for j in range(i + 1, num_nodes):
            p_mutual_modifier_i = None if d_mod is None else d_mod[i]
            p_mutual_modifier_j = None if d_mod is None else d_mod[j]
            P[i, j] = P[j, i] = cal_p_oth(
                points[i], points[j], p_mutual_modifier_i, p_mutual_modifier_j)

    # L matrix (partial inductance)
    print("  Computing L matrix...")
    L = np.zeros((num_branches, num_branches))
    peec_solver = PEECSolver(gauss_order=config.DEFAULT_GAUSS_ORDER)

    branch_regions = []
    branch_axes = []
    branch_dims = []
    for i in range(num_branches):
        region, axis, dims = get_branch_info(i, connects, branch_sizes, points)
        branch_regions.append(region)
        branch_axes.append(axis)
        branch_dims.append(dims)

    print("  Computing L self terms...")
    for i in range(num_branches):
        _, t_a1, t_a2, t_a3 = peec_solver.compute_pair_integral(
            branch_regions[i], branch_regions[i],
            branch_axes[i], branch_axes[i],
            branch_dims[i], branch_dims[i])
        L[i, i] = t_a1 #+ t_a2 + t_a3

        sys.stdout.write("\r  Self progress: {}/{}".format(i + 1, num_branches))
        sys.stdout.flush()
    print()

    print("  Computing L mutual terms...")
    capped_mutual_count = 0
    for i in range(num_branches):
        for j in range(i + 1, num_branches):
            mutual = cal_l_oth_approx(points, connects[i], connects[j])
            capped_mutual = cap_mutual_by_self_terms(mutual, L[i, i], L[j, j])
            if capped_mutual != mutual:
                capped_mutual_count += 1
            L[i, j] = L[j, i] = capped_mutual

        sys.stdout.write("\r  Mutual progress: {}/{}".format(i + 1, num_branches))
        sys.stdout.flush()
    print()
    print("  Mutual terms capped by self inductance: {}".format(capped_mutual_count))

    # ================================================================
    # Step 2: Format Output for C++ Solver
    # ================================================================
    print("\n[Step 2/3] Writing solver input files...")
    save_matrix_cpp_format(L, os.path.join(workspace, 'L.txt'))
    save_matrix_cpp_format(P, os.path.join(workspace, 'P.txt'))
    save_matrix_cpp_format(L, os.path.join(workspace, 'LL_OO.txt'))
    save_matrix_cpp_format(P, os.path.join(workspace, 'PP_OO.txt'))
    save_b2n(connects, os.path.join(workspace, 'B2N.txt'))
    convert_port_file(port_file, os.path.join(workspace, 'PORT.txt'))

    # Copy solver config files to workspace
    for fname in ['set.txt', 'DIELECTRIC.txt', 'Source.txt', 'SourceConfig.txt']:
        src = os.path.join(input_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(workspace, fname))

    print("  Files written to {}".format(os.path.abspath(workspace)))

    # ================================================================
    # Step 3: Run C++ Frequency Solver
    # ================================================================
    print("\n[Step 3/3] Running C++ frequency domain solver...")
    solver_exe = os.path.abspath(config.CPP_SOLVER_EXE)
    workspace_abs = os.path.abspath(workspace)

    if not os.path.exists(solver_exe):
        print("  Warning: packaged solver executable not found at:")
        print("    " + solver_exe)
        print("  Please make sure 'Solver_exe - Improve' is present.")
        return

    runtime_dir, runtime_exe = prepare_solver_runtime(workspace_abs, solver_exe)

    env = os.environ.copy()
    extra_paths = []
    if hasattr(config, 'MKL_DLL_DIR') and os.path.isdir(config.MKL_DLL_DIR):
        extra_paths.append(config.MKL_DLL_DIR)
    extra_paths.append(runtime_dir)
    extra_paths.append(os.path.dirname(solver_exe))
    env["PATH"] = os.pathsep.join(extra_paths) + os.pathsep + env.get("PATH", "")

    cmd = [runtime_exe]
    print("  Executing: {} (cwd={})".format(os.path.basename(runtime_exe), runtime_dir))
    print("-" * 60)

    result = subprocess.run(cmd, env=env, cwd=runtime_dir)

    print("-" * 60)
    if result.returncode == 0:
        collect_solver_outputs(runtime_dir, workspace_abs)
        print("\n" + "=" * 60)
        print("  Simulation completed successfully!")
        print("  Results directory: " + workspace_abs)
        print("  Output files:")
        print("    - map_Final.txt   (S parameters)")
        print("    - Z_in_Final.txt  (Z parameters)")
        print("    - Y_in_Final.txt  (Y parameters)")
        print("=" * 60)
    else:
        print("\n  Solver exited with return code {}".format(result.returncode))


if __name__ == "__main__":
    main()
