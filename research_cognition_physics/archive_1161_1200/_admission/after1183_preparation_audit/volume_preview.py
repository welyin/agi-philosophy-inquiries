"""Exact initial effective-volume preview, not an actual readout certificate."""
from fractions import Fraction as F
import json

epsilon, nu, g, omega, zt_squared = F(3,20), F(7,20), F(1,5), F(1,2), F(3,4)
denominator = 1-g*g*zt_squared
slope = 4*epsilon*omega**2 * (3*epsilon-nu*(1+2/denominator))
common = -8*omega**2*g*g*zt_squared/denominator
assert slope == -F(1809,19400)
assert common+slope/2 == -F(4209,38800)
assert common-slope/2 == -F(591,38800)
print(json.dumps(dict(kind='unnumbered_effective_dictionary_preview',
    theta_dot_difference=str(slope), theta_dot_plus=str(common+slope/2),
    theta_dot_minus=str(common-slope/2), initial_volume=str(1/denominator),
    actual_readout_certified=False, new_science_groups=0)))
