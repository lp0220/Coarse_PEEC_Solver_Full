import math

import numpy as np

from config import E0, U0


def cal_p_self(b, c=None):
    coeff = 2 * math.log(1 + math.sqrt(2)) / (math.pi * E0)
    if c is None:
        return coeff / b
    return coeff / (math.sqrt(b) * c)


def cal_p_oth(p1, p2, c1=None, c2=None):
    """Calculate mutual coefficient of potential between two nodes."""
    dist = np.linalg.norm(p1 - p2)
    modifier = 1.0
    if c1 is not None and c2 is not None:
        modifier = math.sqrt(c1 * c2)
    return 1 / (dist * 4 * math.pi * E0)


def cal_l_oth_approx(points, c1, c2):
    pa1, pa2 = points[c1[0] - 1], points[c1[1] - 1]
    pb1, pb2 = points[c2[0] - 1], points[c2[1] - 1]

    mid1 = (pa1 + pa2) / 2
    mid2 = (pb1 + pb2) / 2
    dist = np.linalg.norm(mid1 - mid2)

    v1 = pa2 - pa1
    v2 = pb2 - pb1
    return U0 * np.dot(v1, v2) / (dist * 4 * math.pi)
