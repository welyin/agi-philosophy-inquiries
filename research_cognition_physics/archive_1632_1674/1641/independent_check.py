"""Run 2 five-coefficient Gaussian cross-check of the Bs -> DsK phase."""
import argparse
import json
from math import cos, sin, radians, degrees
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def covariance(errors, upper):
    correlation = np.eye(5)
    correlation[np.triu_indices(5, 1)] = upper
    correlation += np.triu(correlation, 1).T
    return correlation * np.outer(errors, errors)


def evaluate():
    measured = np.array([.791, -.051, -.303, -.571, -.503])
    statistical = covariance([.061, .134, .125, .084, .084],
                             [.134, .130, .039, .022, .501, -.108, -.036, -.056, -.067, .006])
    systematic = covariance([.022, .058, .055, .023, .025],
                            [.008, .012, -.080, -.246, .878, .004, -.022, -.002, -.022, .085])
    precision = np.linalg.inv(statistical + systematic)

    def profile(angle):
        c, s = cos(radians(angle)), sin(radians(angle))
        design = np.array([[0, 0, 1], [-c, -s, 0], [-c, s, 0], [-s, c, 0], [-s, -c, 0]])
        gram = design.T @ precision @ design
        d, basis = np.linalg.eigh(gram)
        g = basis.T @ design.T @ precision @ measured
        low = np.nextafter(-d[0], np.inf)
        high = max(d[-1], np.linalg.norm(g))
        assert np.sum((g/(d+low))**2) > 1 and np.sum((g/(d+high))**2) < 1
        for _ in range(62):
            mid = (low+high)/2
            if np.sum((g/(d+mid))**2) > 1:
                low = mid
            else:
                high = mid
        u = basis @ (g/(d+high))
        residual = measured-design@u
        assert abs(u@u-1) < 1e-12 and np.linalg.norm(gram@u+high*u-design.T@precision@measured) < 1e-9
        return float(residual@precision@residual), u

    grid = np.array([profile(i)[0] for i in range(180)])
    i = int(grid.argmin())
    assert 15 < i < 165
    a, b = i-1., i+1.
    for _ in range(52):
        l, r = (2*a+b)/3, (a+2*b)/3
        if profile(l)[0] < profile(r)[0]:
            b = r
        else:
            a = l
    central = (a+b)/2
    minimum, u = profile(central)

    def interval(delta):
        inside = np.flatnonzero(grid < minimum+delta)
        assert np.all(np.diff(inside) == 1)
        endpoints = []
        for outer, inner in [(inside[0]-1, inside[0]), (inside[-1]+1, inside[-1])]:
            for _ in range(42):
                mid = (outer+inner)/2
                if profile(mid)[0] > minimum+delta:
                    outer = mid
                else:
                    inner = mid
            endpoints.append(round((outer+inner)/2, 5))
        return endpoints

    assert 65 < central < 83 and abs(profile(central+180)[0]-minimum)<1e-12
    fixed = json.loads((HERE/'../1635/results.json').read_bytes())
    angles = []
    for key in ['common_witness','high_mass_witness']:
        packed = np.asarray(fixed[key]['light_charged_K'])
        k = packed[:,:,0] + 1j*packed[:,:,1]
        z = k[:,1].conj()*k[:,2]
        angles.append(round(degrees(np.angle(z[1]/z[0])),6))
    assert max(abs(v-66.465) for v in angles)<.001
    lower = 68 - degrees(-.021)
    assert 2.73 < lower-66.4705 < 2.74
    return {'Run2_five_physical_coefficients': {'angle_deg': round(central,4), 'strong_phase_deg': round(degrees(np.arctan2(u[1],u[0])),3), 'r': round(float(np.sqrt((1-u[2])/(1+u[2]))),6), 'chi2_min': round(minimum,7), 'delta_1': interval(1), 'delta_4': interval(4)},
            'gamma_s_fixed_degrees': angles, 'combined_one_width_gap_deg_rounded': round(lower-66.4705,4)}


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--check', action='store_true')
    args = cli.parse_args()
    data = evaluate()
    path = HERE/'independent_results.json'
    if args.check:
        assert json.loads(path.read_bytes()) == data
        print('PASS: five-coefficient constrained time profile and fixed weak phases')
    else:
        path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
