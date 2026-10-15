#!/usr/bin/env python3
"""Small finite contour consistency check; it is not the uniform proof."""

import json
import platform
import mpmath as mp


mp.mp.dps = 65
rows = []
I = [0, 1]
J = [2, 3]

for xv, yv in ((mp.mpf("2.3"), mp.mpf("3.7")),
                (mp.mpf("1.4"), mp.mpf("4.2"))):
    det_i = (1 / (xv + I[0]) / (yv + I[1]) -
             1 / (yv + I[0]) / (xv + I[1]))
    det_j = (1 / (xv + J[0]) / (yv + J[1]) -
             1 / (yv + J[0]) / (xv + J[1]))
    left = (xv - 1) ** 3 * (yv - 1) ** 3 * det_i * det_j
    right = ((xv - 1) ** 3 * (yv - 1) ** 3 * (yv - xv) ** 2 /
             mp.fprod([xv + k for k in range(4)]) /
             mp.fprod([yv + k for k in range(4)]))
    assert abs(left - right) < mp.mpf("1e-55")

q = mp.mpf(31)
kappa = mp.mpf(4)
c = mp.mpf("1.5")
series_values = []
contour_values = []
for i in I:
    srow = []
    crow = []
    for j in J:
        function = lambda zz: (zz - 1) ** 3 / ((zz + i) * (zz + j))
        assert abs(mp.diff(function, 1, 2)) < mp.mpf("1e-55")
        direct = mp.fsum([(-q) ** (-r) * function(mp.mpf(r))
                          for r in range(2, 101)])

        def integrand(t):
            zz = c + kappa * t + 1j * t
            value = mp.exp(-mp.log(q) * zz) * function(zz) / mp.sin(mp.pi * zz)
            return -(value.real + kappa * value.imag)

        contour = mp.quad(integrand,
                          [0, mp.mpf(1) / 16, mp.mpf(1) / 4, 1, 4, mp.inf])
        assert abs(direct - contour) < mp.mpf("1e-55")
        srow.append(direct)
        crow.append(contour)
    series_values.append(srow)
    contour_values.append(crow)

bound = max(mp.sqrt(2) * mp.exp(mp.pi / (4 * kappa)),
            2 / (1 - mp.exp(-mp.pi / (2 * kappa))))
for t in (0, mp.mpf(".001"), mp.mpf(1) / 16, mp.mpf(".125"), 1, 10):
    zz = c + kappa * t + 1j * t
    assert abs(abs(mp.sin(mp.pi * zz)) ** 2 -
               (mp.cos(mp.pi * kappa * t) ** 2 + mp.sinh(mp.pi * t) ** 2)) < mp.mpf("1e-30")
    assert 1 / abs(mp.sin(mp.pi * zz)) <= bound * mp.exp(-mp.pi * t)

print(json.dumps({
    "exact_cauchy_quotient": True,
    "triple_zero_second_derivative": True,
    "sine_identity_and_bound_samples": True,
    "determinant_series": mp.nstr(mp.det(mp.matrix(series_values)), 25),
    "determinant_contour": mp.nstr(mp.det(mp.matrix(contour_values)), 25),
    "warning": "finite consistency check only; the uniform estimate is proved in the paper",
    "python": platform.python_version(),
    "mpmath": mp.__version__,
}, indent=2))
