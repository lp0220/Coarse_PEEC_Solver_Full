import numpy as np
from numpy.polynomial.legendre import leggauss
from config import U0, E0


class PEECSolver:
    def __init__(self, gauss_order=3):
        self.gauss_order = gauss_order
        self.math_eps = 1e-20

    def get_gauss_weights(self, a, b):
        pts, wts = leggauss(self.gauss_order)
        mid = (a + b) / 2.0
        hw = (b - a) / 2.0
        return mid + hw * pts, hw * wts

    def rect_pulse(self, r, region, dims, long_axis):
        x, y, _ = r
        xmin, xmax, ymin, ymax, _, _ = region
        if (xmin <= x <= xmax) and (ymin <= y <= ymax):
            width = dims[1] if long_axis == 0 else dims[0]
            val = 1.0 / width if width > self.math_eps else 0.0
            vec = np.zeros(3)
            vec[long_axis] = val
            return vec
        return np.zeros(3)

    def calculate_mt_value(self, grid_i, grid_m, axis_i, axis_m):
        """Calculate analytical Mt integral value."""
        xi_min, xi_max, yi_min, yi_max, zi_min, zi_max = grid_i
        xm_min, xm_max, ym_min, ym_max, zm_min, zm_max = grid_m

        a, d = abs(yi_max - yi_min), abs(ym_max - ym_min)
        l1, l2 = abs(xi_max - xi_min), abs(xm_max - xm_min)
        l3, P, E = xm_min - xi_min, zm_min - zi_min, ym_min - yi_min

        def mt_kernel(x, y, p_val):
            sq = np.sqrt(y ** 2 + p_val ** 2 + x ** 2)
            t1 = ((y ** 2 - p_val ** 2) / 2.0) * x * np.log(x + sq + self.math_eps)
            t2 = ((x ** 2 - p_val ** 2) / 2.0) * y * np.log(y + sq + self.math_eps)
            t3 = -(1.0 / 6.0) * (y ** 2 - 2 * p_val ** 2 + x ** 2) * sq
            t4 = -y * p_val * x * np.arctan2(y * x, p_val * sq + self.math_eps)
            return t1 + t2 + t3 + t4

        y_intervals = [E - a, E + d - a, E + d, E]
        x_intervals = [l3 - l1, l3 + l2 - l1, l3 + l2, l3]

        total = 0
        for i in range(4):
            for k in range(4):
                val = mt_kernel(x_intervals[k], y_intervals[i], P)
                total += ((-1) ** (i + k)) * val

        # axis=0 means the branch runs along x, so its strip width is in y.
        # axis=1 means the branch runs along y, so its strip width is in x.
        width_i = a if axis_i == 0 else l1
        width_m = d if axis_m == 0 else l2
        denom = width_i * width_m * 4 * np.pi
        if abs(denom) < self.math_eps:
            return 0.0

        return total / denom

    def compute_pair_integral(self, g_i, g_m, axis_i, axis_m, dims_i, dims_m, eps_r=1.0):

        # 1. A1 Part 1 (Mt)
        Mt = self.calculate_mt_value(g_i, g_m, axis_i, axis_m)
        if np.isnan(Mt): Mt = 0.0

        zi = (g_i[4] + g_i[5]) / 2.0
        zm = (g_m[4] + g_m[5]) / 2.0
        ci = [(g_i[0] + g_i[1]) / 2, (g_i[2] + g_i[3]) / 2, zi]
        cm = [(g_m[0] + g_m[1]) / 2, (g_m[2] + g_m[3]) / 2, zm]

        b_vec_i = self.rect_pulse(ci, g_i, dims_i, axis_i)
        b_vec_m = self.rect_pulse(cm, g_m, dims_m, axis_m)

        A1_p1 = U0 * Mt

        # 2. A1 Part 2, A2, A3 (Gauss quadrature)
        xi_p, xi_w = self.get_gauss_weights(g_i[0], g_i[1])
        yi_p, yi_w = self.get_gauss_weights(g_i[2], g_i[3])
        xm_p, xm_w = self.get_gauss_weights(g_m[0], g_m[1])
        ym_p, ym_w = self.get_gauss_weights(g_m[2], g_m[3])

        sum_A1_p2 = 0.0
        sum_A2 = 0.0
        sum_A3 = 0.0
        eps = eps_r * E0

        for xi, wxi in zip(xi_p, xi_w):
            for yi, wyi in zip(yi_p, yi_w):
                ri = np.array([xi, yi, zi])
                bi = self.rect_pulse(ri, g_i, dims_i, axis_i)

                for xm, wxm in zip(xm_p, xm_w):
                    for ym, wym in zip(ym_p, ym_w):
                        rm = np.array([xm, ym, zm])
                        bm = self.rect_pulse(rm, g_m, dims_m, axis_m)

                        weight = wxi * wyi * wxm * wym

                        R_vec = rm - ri
                        R_val = np.linalg.norm(R_vec)

                        if R_val < self.math_eps:
                            R_hat = np.zeros(3)
                            R_eff = self.math_eps
                        else:
                            R_eff = R_val
                            R_hat = R_vec / R_eff

                        # A1 Part 2
                        dot_r_bi = np.dot(R_hat, bi)
                        dot_r_bm = np.dot(R_hat, bm)
                        term_a1 = (dot_r_bi * dot_r_bm) / (8 * np.pi * R_eff)
                        sum_A1_p2 += U0 * term_a1 * weight

                        # A2
                        term_a2 = np.dot(bi, bm) / (6 * np.pi)
                        sum_A2 += term_a2 * weight

                        # A3
                        dot_bb = np.dot(bi, bm)
                        term_a3 = (3 * dot_bb - dot_r_bi * dot_r_bm) / (32 * np.pi)
                        sum_A3 += term_a3 * R_val * weight

        res_A1 = A1_p1 + sum_A1_p2
        res_A2 = sum_A2 * np.sqrt(U0 ** 3 * eps)
        res_A3 = -sum_A3 * (U0 ** 2 * eps)

        return Mt, res_A1, res_A2, res_A3
