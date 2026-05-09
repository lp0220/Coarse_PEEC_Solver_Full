import numpy as np


def get_branch_info(branch_idx, connect, branch_size, points):
    n1_raw, n2_raw = connect[branch_idx]
    p1, p2 = points[n1_raw - 1], points[n2_raw - 1]
    width = branch_size[branch_idx]

    center = (p1 + p2) / 2.0
    vec = p2 - p1
    length = np.linalg.norm(vec)

    long_axis = 0 if abs(vec[0]) >= abs(vec[1]) else 1
    dx = length if long_axis == 0 else width
    dy = width if long_axis == 0 else length

    region = [
        center[0] - dx / 2, center[0] + dx / 2,
        center[1] - dy / 2, center[1] + dy / 2,
        center[2], center[2]
    ]
    return region, long_axis, [dx, dy, 0]


def split_grid(region, x_div, y_div):
    xmin, xmax, ymin, ymax, zmin, zmax = region
    dx_step = (xmax - xmin) / max(1, x_div)
    dy_step = (ymax - ymin) / max(1, y_div)

    subgrids = []
    for i in range(x_div):
        for j in range(y_div):
            subgrids.append([
                xmin + i * dx_step, xmin + (i + 1) * dx_step,
                ymin + j * dy_step, ymin + (j + 1) * dy_step,
                zmin, zmax
            ])
    return subgrids
