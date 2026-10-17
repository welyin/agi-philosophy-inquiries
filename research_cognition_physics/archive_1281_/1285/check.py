"""1285: exact finite moments, source and FLRW consistency.
Standard library only. Default is read-only; --write creates results exclusively.
No finite-dimensional Hamiltonian compression and no propagation quadrature.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir() and (p / 'scripts').is_dir())
RESEARCH = ROOT / 'research_cognition_physics'


def add(*ps):
    out = {}
    for p in ps:
        for key, value in p.items():
            out[key] = out.get(key, F(0)) + value
    return {k: v for k, v in out.items() if v}


def scale(p, c):
    return {k: v*c for k, v in p.items() if v*c}


def mul(p, q):
    out = {}
    for k, v in p.items():
        for l, w in q.items():
            z = tuple(a+b for a, b in zip(k, l))
            out[z] = out.get(z, F(0)) + v*w
    return {k: v for k, v in out.items() if v}


def derivative(p, axis):
    out = {}
    for k, v in p.items():
        if k[axis]:
            z = list(k)
            z[axis] -= 1
            out[tuple(z)] = v*k[axis]
    return out


def radial_operator(p):
    one, x, y = {(0, 0): F(1)}, {(1, 0): F(1)}, {(0, 1): F(1)}
    dx, dy = derivative(p, 0), derivative(p, 1)
    v = add(scale(one, F(5, 4)), scale(add(x, y), F(-1, 8)),
            scale(mul(x, y), F(1, 2)), scale(mul(y, y), F(1, 16)))
    return add(scale(mul(x, derivative(dx, 0)), -2), mul(add(x, scale(one, -3)), dx),
               scale(mul(y, derivative(dy, 1)), -2), mul(add(y, scale(one, -2)), dy), mul(v, p))


def radial_moment(p):
    return sum((v*math.prod(range(3, 3+2*i, 2))*math.prod(range(2, 2+2*j, 2))
                for (i, j), v in p.items()), F(0))


def coordinate_moment(p):
    out = F(0)
    for k, value in p.items():
        if any(n % 2 for n in k):
            continue
        out += value*math.prod(math.prod(range(1, n, 2)) for n in k)
    return out


def coordinate_crosscheck():
    n = 5
    one = {(0,)*n: F(1)}
    coords = []
    for i in range(n):
        k = [0]*n
        k[i] = 1
        coords.append({tuple(k): F(1)})
    x = add(*(mul(z, z) for z in coords[:3]))
    y = add(*(mul(z, z) for z in coords[3:]))
    v = add(scale(mul(x, y), F(1, 2)), scale(mul(y, y), F(1, 16)))
    # Acting on P exp(-|q|^2/4), directly in the five original coordinates.
    def act(p):
        terms = [mul(v, p), scale(p, F(5, 4)), scale(mul(add(x, y), p), F(-1, 8))]
        for i, qi in enumerate(coords):
            terms += [scale(derivative(derivative(p, i), i), F(-1, 2)),
                      scale(mul(qi, derivative(p, i)), F(1, 2))]
        return add(*terms)
    hp = act(one)
    h2p = act(hp)
    for i, j in ((0, 1), (1, 2), (3, 4)):
        rotation = add(mul(coords[i], derivative(hp, j)),
                       scale(mul(coords[j], derivative(hp, i)), -1))
        assert not rotation
    return [coordinate_moment(hp), coordinate_moment(mul(hp, hp)),
            coordinate_moment(mul(hp, h2p)), coordinate_moment(mul(h2p, h2p))]


def fs(x):
    return str(x.numerator) + '/' + str(x.denominator)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def calculate():
    one = {(0, 0): F(1)}
    ps = [one]
    for _ in range(4):
        ps.append(radial_operator(ps[-1]))
    moments = [radial_moment(p) for p in ps]
    exact = [F(1), F(33, 8), F(2683, 64), F(495575, 512), F(161114457, 4096)]
    assert moments == exact
    assert coordinate_crosscheck() == moments[1:]
    assert radial_moment(mul(ps[1], ps[1])) == moments[2]
    assert radial_moment(mul(ps[1], ps[2])) == moments[3]
    assert radial_moment(mul(ps[2], ps[2])) == moments[4]
    energy = moments[1]
    variance = moments[2] - energy**2
    d4 = 2*moments[4] - 8*moments[3]*energy + 6*moments[2]**2
    duration = F(1, 100)
    lower = variance*duration**2 - d4*duration**4/24
    upper = variance*duration**2
    assert variance == F(797, 32) and d4 == F(14661753, 256)
    assert F(3, 1250) < lower <= upper < F(1, 400)
    assert lower == F(505192749, 204800000000) and upper == F(797, 320000)
    mp2, slope, volume = F(11, 2), F(1, 2), F(1)
    assert slope**2 == energy/(3*mp2*volume)
    rows = []
    for eta in (F(0), duration/2, duration):
        a = 1+slope*eta
        rho = energy/(volume*a**4)
        pressure = rho/3
        proper_adot = slope/a
        proper_addot = -slope**2/a**3
        friedmann = 3*mp2*(proper_adot/a)**2-rho
        acceleration = proper_addot/a+(rho+3*pressure)/(6*mp2)
        continuity_conformal = -4*energy*slope/(volume*a**5) + 3*slope*(rho+pressure)/a
        assert friedmann == acceleration == continuity_conformal == 0
        rows.append({'eta':fs(eta), 'a':fs(a), 'rho':fs(rho), 'pressure':fs(pressure),
                     'friedmann_residual':fs(friedmann), 'acceleration_residual':fs(acceleration),
                     'continuity_residual':fs(continuity_conformal)})
    source_ledger = {'gauge_kinetic':F(3, 8), 'scalar_kinetic':F(1, 4),
                     'gauge_scalar_interaction':F(3), 'scalar_quartic':F(1, 2)}
    assert sum(source_ledger.values()) == energy
    end_a = 1+slope*duration
    proper_time = duration+slope*duration**2/2
    assert end_a == F(201, 200) and proper_time == F(401, 40000)
    sources = [HERE/'proof.md', HERE/'check.py',
               RESEARCH/'archive_1009_1043/research_note_1031.md',
               RESEARCH/'archive_1219_1233/1221/proof.md',
               RESEARCH/'archive_585_628/research_note_625.md',
               RESEARCH/'archive_1159_1160/_shared/field_bridge_after_1160.md',
               RESEARCH/'10_认知操作公理与有效广义相对论_阶段论文.md']
    return {
        'round':1285, 'status':'PASS_WITH_SCOPE', 'author_groups':1, 'new_axioms':0,
        'method':'exact radial chi-square recurrence and independent five-coordinate Gaussian moment crosscheck',
        'source_bindings':{p.relative_to(ROOT).as_posix():sha(p) for p in sources},
        'parameters':{'g':'1/1', 'lambda':'1/4', 'coordinate_volume':fs(volume),
                      'gaussian_coordinate_variance':'1/1', 'M_P_squared':fs(mp2),
                      'initial_scale':'1/1', 'conformal_duration':fs(duration)},
        'energy_ledger':{k:fs(v) for k,v in source_ledger.items()},
        'spectral_moments':{str(i):fs(x) for i,x in enumerate(moments)},
        'radial_polynomial_term_counts':[len(p) for p in ps],
        'variance':fs(variance), 'independent_spectral_difference_fourth_moment':fs(d4),
        'finite_probability':{'nonreturn_lower':fs(lower), 'nonreturn_upper':fs(upper),
                              'lower_decimal':float(lower), 'upper_decimal':float(upper),
                              'strict_coarse_lower':'3/1250', 'strict_coarse_upper':'1/400'},
        'same_source_geometry':{'slope':fs(slope), 'end_scale':fs(end_a),
                                'proper_duration':fs(proper_time), 'exact_checks':rows},
        'mathematical_luders_poststate_energy_upper':fs(5*energy),
        'scope':{'declared_homogeneous_effective_quantum_model':True,
                 'normal_Gauss_rotation_invariant_initial_state':True,
                 'interaction_energy_included':True, 'full_homogeneous_mean_source_during_wait':True,
                 'classical_parent_reduction_with_boundary_dictionary':True,
                 'finite_hilbert_projection_used':False, 'quantum_full_4D_sector_invariance_proved':False,
                 'full_QFT_anomaly_or_omitted_modes_set_to_zero':False,
                 'actual_survival_instrument_or_source_proved':False,
                 'unknown_input_uniform_semiclassical_geometry':False,
                 'conditional_branch_geometry':False,
                 'thermal_reference_or_temperature_identity':False,
                 'Einstein_equation_newly_generated':False,
                 'stage_completion_claimed':False}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    out=calculate()
    path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(out,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('PASS_WITH_SCOPE: EXCLUSIVE_FIRST_WRITE')
    else:
        assert out == json.loads(path.read_text(encoding='utf-8')), 'Saved result mismatch'
        print('PASS_WITH_SCOPE: READ_ONLY_MATCH')
    print('nonreturn='+out['finite_probability']['nonreturn_lower']+' .. '+out['finite_probability']['nonreturn_upper'])
    print('scale='+out['same_source_geometry']['end_scale']+'; proper_duration='+out['same_source_geometry']['proper_duration'])
    print('results_sha256='+sha(path))


if __name__=='__main__':
    main()
