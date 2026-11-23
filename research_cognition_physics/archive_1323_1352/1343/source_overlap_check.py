"""Finite overlap and clock tolerance for the color-source comparison."""
import math
import numpy as np


def excite(arm, eta, theta):
    if arm == "unitary":
        xi = np.array([-.5, .25])
        w = np.array([3, 6])/9
    else:
        xi = np.array([-1., -.5, .5])
        w = np.array([1, 3, 5])/9
    return w @ np.sin(theta*np.sqrt(2+2*eta*xi))**2


for group, dim in [("unitary", 8), ("rotation", 3)]:
    assert abs(excite(group, 0, 1.) - math.sin(math.sqrt(2))**2) < 1e-14
    for eta in [-1., 0., .5, 1.]:
        theta = .02
        q = (4+4*eta**2/dim)/3
        # Sine series, using normalized squared transfer x <= 4.
        assert abs(excite(group, eta, theta)-(2*theta**2-q*theta**4)) <= 64*2/45*theta**6
gap = excite("unitary", 1, math.pi)-excite("rotation", 1, math.pi)
certified = gap-2*(2*math.pi**2*.001+.01)
assert certified > .255 and .9995**11000 < .005
for eta in np.linspace(.999, 1., 31):
    for theta in np.linspace(math.pi-.005, math.pi+.005, 31):
        assert excite("unitary", eta, theta)-excite("rotation", eta, theta) >= certified
print(f"ideal_gap={gap:.10f} guaranteed_gap={certified:.10f} bad_calibration={.9995**11000:.10f}")
