"""Finite witnesses for 1067; the dimension theorem is proved analytically.

Default mode recomputes and compares saved results. --write creates, never overwrites.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
I = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
AXES = np.eye(3)


def effect(v, center=0.5, kappa=0.4):
    return center * I + kappa * np.einsum('j,jab->ab', v, PAULI)


def state(v):
    return effect(v, kappa=0.5)


def probability(rho, e):
    p = np.trace(rho @ e)
    assert abs(p.imag) < 1e-12
    return float(p.real)


def rotation(axis, angle):
    a = np.array(axis, float)
    a /= np.linalg.norm(a)
    cross = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return math.cos(angle) * AXES + (1 - math.cos(angle)) * np.outer(a, a) + math.sin(angle) * cross


def unitary(axis, angle):
    a = np.array(axis, float)
    a /= np.linalg.norm(a)
    return math.cos(angle / 2) * I - 1j * math.sin(angle / 2) * np.einsum('j,jab->ab', a, PAULI)


def opnorm(m):
    return float(np.linalg.norm(m, 2))


def run():
    residuals = []
    a, b = rotation(AXES[2], math.pi / 2), rotation(AXES[0], math.pi / 2)
    ua, ub = unitary(AXES[2], math.pi / 2), unitary(AXES[0], math.pi / 2)
    x = AXES[2]
    ab, ba = a @ b @ x, b @ a @ x
    rho = state((AXES[0] + AXES[1]) / math.sqrt(2))
    p_ab, p_ba = probability(rho, effect(ab)), probability(rho, effect(ba))
    gap = abs(p_ab - p_ba)
    assert abs(gap - 0.8 / math.sqrt(2)) < 1e-12
    residuals += [float(np.linalg.norm(ab - AXES[0])), float(np.linalg.norm(ba + AXES[1]))]

    # Same family of transformations and effects, including intermediate settings.
    samples = [*AXES, *(-AXES), np.array([1, 2, 3]) / math.sqrt(14)]
    covariance_cases = 0
    for axis in AXES:
        for t in [0, 1/7, 1/3, 1/2, 1]:
            r, u = rotation(axis, t * math.pi/2), unitary(axis, t * math.pi/2)
            for v in samples:
                residuals.append(opnorm(effect(r @ v) - u @ effect(v) @ u.conj().T))
                covariance_cases += 1

    # Equation (4): two nonorthogonal vector images determine the rotation.
    v1, v2 = 0.8 * np.array([0.6, 0.8, 0]), 0.8 * np.array([0, 0.6, 0.8])
    frame = np.column_stack([v1, v2, np.cross(v1, v2)])
    rec = []
    for r in [a, b, a @ b, b @ a]:
        r1, r2 = r @ v1, r @ v2
        rebuilt = np.column_stack([r1, r2, np.cross(r1, r2)]) @ np.linalg.inv(frame)
        residuals.append(float(np.linalg.norm(rebuilt - r)))
        rec.append(rebuilt)
    residuals.append(float(np.linalg.norm(rec[0] @ rec[1] - rec[2])))

    orbit_points = np.concatenate([AXES, -AXES])
    contrast_rank = int(np.linalg.matrix_rank(orbit_points - orbit_points[0]))
    assert contrast_rank == 3

    # O(2) acts continuously and covariantly, but reflection is not in its identity component.
    flip = np.diag([1, -1, -1])
    uflip = PAULI[0]
    equator = [np.array([math.cos(t), math.sin(t), 0]) for t in np.linspace(0, 2*math.pi, 25)]
    for v in equator:
        residuals.append(opnorm(effect(flip @ v) - uflip @ effect(v) @ uflip.conj().T))
    eq_ab, eq_ba = a @ flip @ AXES[0], flip @ a @ AXES[0]
    eq_probs = [probability(state(AXES[1]), effect(v)) for v in [eq_ab, eq_ba]]
    assert abs(eq_probs[0] - eq_probs[1] - 0.8) < 1e-12

    # Möbius circle maps give identity paths and order sensitivity, without common unitary covariance.
    def mobius(z, r):
        return (z-r)/(1-r*z)
    def vec(z):
        return np.array([z.real, z.imag, 0])
    m_ab, m_ba = 1j * mobius(1+0j, .5), mobius(1j, .5)
    mob_probs = [probability(state(AXES[0]), effect(vec(z))) for z in [m_ab, m_ba]]
    assert abs(abs(mob_probs[0]-mob_probs[1]) - 0.32) < 1e-12
    gram_before = float(vec(1+0j) @ vec(1j))
    gram_after = float(vec(mobius(1+0j,.5)) @ vec(mobius(1j,.5)))
    assert abs(gram_after + .8) < 1e-12 and gram_before == 0
    for t in [0, .1, .3, .5]:
        for phi in np.linspace(0, 2*math.pi, 25):
            z = complex(math.cos(phi), math.sin(phi))
            residuals += [abs(abs(mobius(z,t))-1), abs(mobius(mobius(z,t),-t)-z)]

    # Nonminimal S^3 and the trace-warning S^3 share the same actual rotations.
    # The parameter t is an invariant fourth coordinate, not physical time.
    counter_points = []
    for t in np.linspace(-1, 1, 17):
        for n in samples:
            v = math.sqrt(max(0, 1-t*t)) * n
            counter_points.append((v, float(t)))
    gaps, tf_gaps, eigs, coords = [], [], [], []
    for v, t in counter_points:
        e = effect(v, center=.5+.1*t, kappa=.3)
        eo = effect(-v, center=.5-.1*t, kappa=.3)
        diff = e-eo
        g = opnorm(diff)
        expected = .2*abs(t)+.6*float(np.linalg.norm(v))
        residuals.append(abs(g-expected))
        gaps.append(g)
        tf_gaps.append(opnorm(diff-np.trace(diff)*I/2))
        eigs.extend(np.linalg.eigvalsh(e).tolist())
        coords.append([.1*t, *(.3*v)])
        residuals.append(opnorm(effect(a @ v, .5+.1*t, .3) - ua @ e @ ua.conj().T))
        residuals.append(opnorm(effect(b @ v, .5+.1*t, .3) - ub @ e @ ub.conj().T))
    coords = np.array(coords)
    full_rank = int(np.linalg.matrix_rank(coords - coords[0]))
    projected_rank = int(np.linalg.matrix_rank(coords[:,1:] - coords[0,1:]))
    assert min(eigs) >= 0 and max(eigs) <= 1
    assert abs(min(gaps)-.2) < 1e-12 and min(tf_gaps) == 0
    assert full_rank == 4 and projected_rank == 3
    fixed_poles = all(np.array_equal(r @ np.zeros(3), np.zeros(3)) for r in [a,b])
    assert fixed_poles
    max_res = float(max(residuals))
    assert max_res < 1e-12
    return {
        'round':1067,
        'scope':'Finite models verify formulas and distinguish hypotheses; no sampling proves topology or universal covariance.',
        's2_valid_contract':{'ab_endpoint':ab.tolist(),'ba_endpoint':ba.tolist(),
            'same_preparation_probabilities':[p_ab,p_ba],'order_gap':gap,
            'read_error_each':.01,'certified_gap':gap-.02,'contrast_rank':contrast_rank,
            'sampled_covariance_cases':covariance_cases},
        'without_identity_path':{'source':'S1','group':'O2','same_preparation_probabilities':eq_probs,
            'order_gap':float(eq_probs[0]-eq_probs[1]),'reflection_identity_reachable':False},
        'without_global_unitary_covariance':{'source':'S1','mobius_r':.5,
            'ab_endpoint':[m_ab.real,m_ab.imag],'ba_endpoint':[m_ba.real,m_ba.imag],
            'same_preparation_probabilities':mob_probs,'order_gap':abs(mob_probs[0]-mob_probs[1]),
            'gram_before':gram_before,'gram_after':gram_after},
        'nonminimal_contract':{'source':'S3','spatial_dimension':4,'fixed_poles':fixed_poles,
            'minimal_action':False,'order_gap_at_equator':.6/math.sqrt(2)},
        'trace_warning':{'source':'S3','spatial_dimension':4,'lambda':.1,'kappa':.3,
            'point_checks':len(counter_points),'minimum_full_antipodal_gap':min(gaps),
            'minimum_traceless_antipodal_gap':min(tf_gaps),
            'full_contrast_rank':full_rank,'traceless_contrast_rank':projected_rank,
            'sample_effect_eigenvalue_range':[float(min(eigs)),float(max(eigs))],
            'all_antipodal_full_effects_separated_analytic':True,
            'upper_bound_traceless_separation_satisfied':False},
        'max_formula_residual':max_res,
        'unconditional_spatial_dimension_generated':False,
    }


def compare(a,b,path='result'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,abs_tol=1e-11,rel_tol=1e-10),(path,a,b)
    else:
        assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run();p=HERE/'results.json'
    if args.write:
        with p.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(result,json.loads(p.read_text(encoding='utf8')))
    print(json.dumps({'round':1067,'passed':True,'write':args.write,
                      'max_formula_residual':result['max_formula_residual']},ensure_ascii=False))


if __name__=='__main__':main()
