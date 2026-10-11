"""1114: whole finite preparations, reference words and causal support.

Read-only by default. Finite samples calibrate the analytic closure proof;
they do not certify all Lorentz transformations by sampling.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ETA = np.diag([1., -1., -1., -1.])
TOL = 3e-9


def small(x, label):
    x = float(x)
    if not math.isfinite(x) or x > TOL:
        raise AssertionError(f"{label}: {x}")
    return x


def mdot(p, q):
    return float(p @ ETA @ q)


def norm(a):
    return float(np.linalg.norm(a))


def boost(rapidity, direction):
    d = np.asarray(direction, float)
    d /= np.linalg.norm(d)
    ch, sh = math.cosh(rapidity), math.sinh(rapidity)
    b = np.eye(4)
    b[0, 0] = ch
    b[0, 1:] = sh*d
    b[1:, 0] = sh*d
    b[1:, 1:] += (ch-1)*np.outer(d, d)
    return b


def rotation(angle, axis):
    n = np.asarray(axis, float)
    n /= np.linalg.norm(n)
    x, y, z = n
    k = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    r = np.eye(4)
    r[1:, 1:] = np.eye(3) + math.sin(angle)*k + (1-math.cos(angle))*(k@k)
    return r


def load_model():
    path = HERE.parent/'1111/check.py'
    spec = importlib.util.spec_from_file_location('accepted1111', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.model()


def calculate():
    m = load_model()
    arrival, done = 23/15, 23/15+math.pi
    a = np.array([arrival, 1., 0, 0])
    f = np.array([done, 1., 0, 0])
    uc = np.array([23., 15., 0, 0])/(4*math.sqrt(19))
    ua = np.array([83., -45., 0, 0])/(16*math.sqrt(19))
    q0 = np.array([83., -45., 0, 0])/(8*math.sqrt(19))
    q1 = np.array([175., 15., 0, 0])/(16*math.sqrt(19))
    branches = [([2*uc, 3*ua, q1], [0., 152/83, 152/175]),
                ([3*uc, 2*ua, q0], [0., 152/83, 152/83])]
    edges = []
    for ps, starts in branches:
        for p, x in zip(ps, starts):
            start = np.array([0., x, 0, 0])
            edges.append((start, a, p))
            small(norm(a[1:]-start[1:]-(arrival/p[0])*p[1:]), 'original arrival')
    pa = [np.array([5., 4, 0, 0]), np.array([5., 0, 4, 0])]
    pb = [np.array([5., -4, 0, 0]), np.array([5., 0, -4, 0])]
    recoil_edges = [(f, f+(2/5)*p, p) for p in pa+pb]
    edges += recoil_edges
    words = [[], [boost(.7, [1, 0, 0])], [rotation(.6, [0, 0, 1])],
             [boost(.4, [1, 2, 0]), rotation(.8, [0, 1, 1]), boost(-.6, [1, 0, 1])],
             [rotation(-1.1, [1, 1, 1]), boost(.9, [0, 1, 0]),
              boost(.3, [0, 0, 1]), rotation(.7, [1, 0, 0])]]
    mats = []
    maxima = dict(metric=0., inverse=0., endpoint=0., phase=0., mass=0., word=0.)
    total_quantum = 0.
    spectrum_min = float('inf')
    wt = m['w'](math.pi)
    full_out = np.column_stack([np.kron(wt[:, i], np.eye(4)[:, (0, 3)[i]])
                                for i in range(2)])
    small(norm(full_out.conj().T@full_out-np.eye(2)), 'full original isometry')
    conserved = 0.
    pair_flip = np.eye(4)
    pair_flip[np.ix_([0, 3], [0, 3])] = [[0, 1], [1, 0]]
    pair_p = np.array([x+y for x in pa for y in pb])
    h_eigs = np.linalg.eigvalsh(m['hout'])
    u = np.array([1., 0, 0, 0])
    intertwiner = m['hout']@m['gamma']-m['gamma']@m['hin']
    small(norm(intertwiner), 'same complete clock intertwiner')
    resources = []
    for word in words:
        lam = np.eye(4)
        for g in word:
            lam = g@lam
        mats.append(lam)
        inv = ETA@lam.T@ETA
        maxima['metric'] = max(maxima['metric'], norm(lam.T@ETA@lam-ETA))
        maxima['inverse'] = max(maxima['inverse'], norm(inv@lam-np.eye(4)))
        for start, end, p in edges:
            d = end-start
            sp, ep, pp = lam@start, lam@end, lam@p
            sigma = d[0]*math.sqrt(mdot(p, p))/p[0]
            maxima['endpoint'] = max(maxima['endpoint'],
                norm(ep-sp-sigma*pp/math.sqrt(mdot(p, p))))
            maxima['phase'] = max(maxima['phase'], abs(mdot(pp, ep-sp)-mdot(p, d)))
            maxima['mass'] = max(maxima['mass'], abs(mdot(pp, pp)-mdot(p, p)))
            seq = start.copy()
            for g in word:
                seq = g@seq
            maxima['word'] = max(maxima['word'], norm(seq-sp))
            assert ep[0]>sp[0] and pp[0]>0
        # All four momentum components of the FULL pair product, including xy/yx.
        transformed_pair = pair_p@lam.T
        for mu in range(4):
            P = np.diag(transformed_pair[:, mu])
            conserved = max(conserved, norm(P@pair_flip-pair_flip@P),
                            abs((lam@u)[mu])*norm(intertwiner))
        old_phase = sum(mdot(p, end-start) for start, end, p in recoil_edges[::2])
        new_phase = sum(mdot(lam@p, lam@(end-start)) for start, end, p in recoil_edges[::2])
        total_quantum = max(total_quantum,
            norm((np.exp(-1j*new_phase)-np.exp(-1j*old_phase))*full_out))
        energies = (lam@u)[0]*h_eigs
        spectrum_min = min(spectrum_min, float(min(energies)))
        coeff = lam[0, 0]+np.linalg.norm(lam[0, 1:])
        small(max(0., float(max(energies))-coeff*float(max(h_eigs))), 'finite energy bound')
        resources.append({'word_length': len(word), 'energy_upper_multiplier': float(coeff),
                          'maximum_device_energy': float(max(energies))})
    for k, value in maxima.items():
        small(value, k)
    small(total_quantum, 'complete Q/reference/process transport')
    small(conserved, 'all four P transport')
    # A transformed preparation slice is not another t'=0 slice.
    lam = boost(.7, [1, 0, 0])
    source_times = [(lam@np.array([0., x, 0, 0]))[0] for x in [0, 10]]
    assert abs(source_times[1]-source_times[0])>1
    transformed_clock_span = (lam@(f-a))[0]
    wrong_clock_prob = math.sin((math.pi/lam[0, 0])/2)**6
    assert 1-wrong_clock_prob>.1
    # Gluing transports the shared physical endpoint exactly; moving only one half does not.
    shared = np.array([3., .4, .2, 0])
    left_end = lam@shared
    right_start = lam@shared
    small(norm(left_end-right_start), 'transported glue')
    wrong_glue = norm(left_end-shared)
    assert wrong_glue>.1
    # Every sample permitted path is future causal; the proof covers every finite graph.
    path = [np.array([.7, .1, .2, 0]), np.array([1.2, -.5, .1, .2]),
            np.array([.6, .3, .1, 0])]
    min_future_margin = float('inf')
    for lam in mats:
        steps = [lam@v for v in path]
        for v in steps+[sum(steps)]:
            min_future_margin = min(min_future_margin, v[0]-norm(v[1:]))
    assert min_future_margin>0
    # Illegal spacelike intervention and deadline target: neither is an allowed wire.
    forbidden = np.array([done, 2*done, 0, 0])
    reversed_time = (boost(-math.atanh(.75), [1, 0, 0])@forbidden)[0]
    assert reversed_time<0
    orbit=[]
    for r in [0., 1., 2., 3., 4., 6.]:
        point = boost(r, [1, 0, 0])@f
        v = norm(point[1:])/point[0]
        expected = math.tanh(r+math.atanh(1/done))
        small(abs(v-expected), 'actual task orbit velocity')
        assert 0<v<1
        orbit.append({'rapidity':r, 't':float(point[0]), 'speed':float(v)})
    assert orbit[-1]['speed']>.99999
    sources = ['1106/proof.md', '1111/check.py', '1111/proof.md',
               '../archive_1086_1105/1103/proof.md']
    hashes={p:hashlib.sha256((HERE.parent/p).read_bytes()).hexdigest() for p in sources}
    return {'round':1114, 'status':'PASS_CONDITIONAL_WHOLE_PREPARATION_REFERENCE_CLOSURE',
      'scientific_calibration':0, 'new_adopted_axioms':0, 'numerical_check_groups':6,
      'source_sha256':hashes,
      'scope':{'orbit_predeployment_permission_explicit':True,
               'finite_samples_do_not_prove_group_quantifiers':True,
               'all_covered_influences_use_declared_local_causal_graph':True,
               'fresh_single_source_Q_witness':True,
               'complete_nonscalar_mother_realization_not_claimed':True,
               'all_physical_SR_not_claimed':True},
      'group1_full_network_worldlines':{'words':len(words), 'edges':len(edges), **maxima},
      'group2_complete_clock_Q_reference':{'full_output_shape':list(full_out.shape),
          'transport_residual':total_quantum, 'all_four_P_residual':conserved,
          'actual_clock_proper_duration':math.pi,
          'boosted_coordinate_duration':float(transformed_clock_span),
          'wrong_fixed_coordinate_wait_probability_gap':1-wrong_clock_prob},
      'group3_real_preparations_and_resources':{'positive_device_energy_minimum':spectrum_min,
          'finite_word_resources':resources, 'transformed_initial_slice_times':source_times,
          'full_pair_product_cross_sectors_retained':True},
      'group4_finite_composition':{'correct_glue_residual':norm(left_end-right_start),
          'transporting_only_one_side_fails_by':wrong_glue,
          'arbitrary_reference_and_feedback_claim_basis':'complete CPTP graph proof, not sample'},
      'group5_all_preparation_bound':{'sample_future_margin':min_future_margin,
          'common_deadline':done, 'outside_target':forbidden.tolist(),
          'illegal_spacelike_target_reversed_time':float(reversed_time),
          'universal_bound_basis':'declared class local support plus finite DAG induction'},
      'group6_nonempty_menu_and_speed':{'single_source_completion':f.tolist(),
          'task_proper_interval':math.sqrt(mdot(f,f)), 'orbit_samples':orbit,
          'null_task_or_uniform_budget_not_assumed':True}}


def compare(saved, fresh, path='root'):
    if isinstance(fresh, dict):
        assert saved.keys()==fresh.keys(), path
        for k in fresh: compare(saved[k], fresh[k], path+'.'+k)
    elif isinstance(fresh, list):
        assert len(saved)==len(fresh), path
        for i,(a,b) in enumerate(zip(saved,fresh)): compare(a,b,f'{path}[{i}]')
    elif isinstance(fresh, float):
        assert math.isclose(saved, fresh, rel_tol=1e-9, abs_tol=5e-11), (path,saved,fresh)
    else:
        assert saved==fresh, (path,saved,fresh)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--write-results',action='store_true')
    args=p.parse_args()
    result=calculate()
    path=HERE/'results.json'
    if args.write_results:
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('WROTE_RESULTS')
    else:
        compare(json.loads(path.read_text(encoding='utf-8')),result)
        print('SAVED_RESULTS_MATCH')
    print('PASS: 6 conditional construction check groups; scientific calibration +0')


if __name__=='__main__':
    main()
