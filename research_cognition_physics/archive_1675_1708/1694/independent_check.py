"""One shared broken-generator current, the measured proton width, and a baryon-conserving control."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def phase_space(r):
    # Integrate over lepton pair mass q^2/M^2. y=ymax(1-u^2) regularizes the endpoint.
    nodes, weights = np.polynomial.legendre.leggauss(64)
    u = (nodes+1)/2
    ymax = (1-math.sqrt(r))**2
    y = ymax*(1-u**2)
    lam = ((1-math.sqrt(r))**2-y)*((1+math.sqrt(r))**2-y)
    shape = (1-r)**2+(1+r)*y-2*y**2
    integral = weights @ (2*np.sqrt(lam)*shape*ymax*u)
    exact = 1-8*r+8*r**3-r**4-(12*r*r*math.log(r) if r else 0)
    np.testing.assert_allclose(integral, exact, rtol=1e-10, atol=1e-13)
    return float(integral)


def compute():
    heavy = json.loads((HERE.parent/'1685/results.json').read_bytes())
    m = 1000*min(heavy['conditional_masses_TeV_descending'])
    # Super-K 2010.16098; two independent scales connected by ratio of Wilson coefficients.
    hbar = 6.582119569e-25
    light_speed = 299792458.
    year = 365.25*86400
    mp, mpi = .938272, .134977
    limit = 2.4e34*year
    floor = .069
    assert (.159-.015-.020-.025)/math.sqrt(2)>floor
    coefficient = math.sqrt(32*math.pi*hbar/(limit*mp))/(floor*(1-(mpi/mp)**2))
    rate_per_coeff = m**5/(1536*math.pi**3)
    time = hbar/(rate_per_coeff*coefficient**2)
    f = phase_space(.25)
    assert abs(phase_space(0)-1)<1e-12
    for r in (.01,.1,.4,.64,.9):
        assert 0<phase_space(r)<1
    scale = 1.e5
    alternate = hbar*light_speed/(rate_per_coeff*f/scale**4)
    # The operator J_5^dagger J_bar5 is invariant: use 6Y, 3B, both additive.
    for d, odd_doublet, dc, lepton in [(-2,3,2,-3),(1,0,-1,0)]:
        assert -d+odd_doublet-dc+lepton==0
    return dict(colored_mass_GeV=m, proton_matrix_element_floor_GeV2=floor,
                proton_single_source_coefficient_upper_GeV_minus2=float(f'{coefficient:.9e}'),
                massless_one_channel_width_per_C_squared_GeV5=float(f'{rate_per_coeff:.9e}'),
                unit_ratio_inclusive_lifetime_lower_seconds=float(f'{time:.9e}'),
                unit_ratio_proper_length_lower_m=float(f'{(time*light_speed):.9e}'),
                necessary_amplitude_norm_ratio_for_1mm=float(f'{math.sqrt(time*light_speed/.001):.9e}'),
                independent_mass_half_ratio_phase_space=round(f,10),
                baryon_preserving_contact_100TeV_ctau_m=float(f'{alternate:.9e}'))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check',action='store_true')
    output = compute()
    path = HERE/'independent_results.json'
    if p.parse_args().check:
        assert output == json.loads(path.read_bytes())
        print('PASS: inclusive phase space, proton coefficient, common-current lifetime and prompt control')
    else:
        with path.open('x',encoding='utf8') as fp:
            fp.write(json.dumps(output,indent=2)+'\n')
