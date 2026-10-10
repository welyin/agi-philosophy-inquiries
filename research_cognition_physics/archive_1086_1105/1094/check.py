"""Round 1094: finite pulse events and a conditional source-law classification.

No dynamical implementation or qualification of the full cognitive axioms is claimed.
Default is read-only recomputation; --write creates results once.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def F(alpha, c, q):
    return (c-q)*(c+(1-alpha)*q)/(c*(c-alpha*q))


def pulse(source, receiver, tick, alpha, c):
    xs, us, rs = source
    xr, ur, rr = receiver
    te = tick / rs
    n = 1 if xr + ur*te > xs + us*te else -1
    w = alpha*us + n*c
    tr = (xr-xs + (w-us)*te)/(w-ur)
    assert tr > te
    assert abs(xs+us*te+w*(tr-te)-(xr+ur*tr)) < 1e-11
    return dict(emission_tick=tick, emission_time=te, arrival_time=tr,
                arrival_tick=rr*tr, emission_x=xs+us*te, arrival_x=xr+ur*tr,
                velocity=w, direction=n)


def protocol(alpha, c, v, r):
    actors = {'A': (-10., 0., 1.), 'B': (0., v, r), 'C': (10., 0., 1.)}
    records, ratios, residual = [], {}, 0.
    for partner, n in [('A', 1), ('C', -1)]:
        intervals = {}
        for s, d in [(partner, 'B'), ('B', partner)]:
            events = [dict(source=s, receiver=d, mode_c=c,
                           **pulse(actors[s], actors[d], k, alpha, c))
                      for k in range(3)]
            records.extend(events)
            ds = np.diff([e['arrival_tick'] for e in events])
            assert np.max(np.abs(ds-ds[0])) < 1e-11
            intervals[s] = float(ds[0])
        q = n*v
        expected_s = (c+(1-alpha)*q)/(r*(c-alpha*q))
        expected_b = r*c/(c-q)
        residual = max(residual, abs(intervals['B']-expected_s),
                       abs(intervals[partner]-expected_b))
        ratios[partner] = intervals['B']/intervals[partner]
        residual = max(residual, abs(ratios[partner]-F(alpha,c,q)/r**2))
    end = max(e['arrival_time'] for e in records)
    assert v*end < 10  # B has not passed the right partner; one common window.
    return records, ratios, residual, end


def interval_certificate(logs, tolerances):
    lo = max(b-e for b,e in zip(logs,tolerances))
    hi = min(b+e for b,e in zip(logs,tolerances))
    return dict(lower=float(lo), upper=float(hi), feasible=bool(lo <= hi))


def calculate():
    groups = []
    def add(name, condition, **data):
        assert condition, name
        groups.append(dict(name=name, status='PASS', **data))

    max_error, total_events, largest_window = 0., 0, 0.
    for a in [0., .5, 1.]:
        for c in [1., 1.4, 2.]:
            for r in [.87, 1.]:
                events, _, err, end = protocol(a,c,.25,r)
                total_events += len(events)
                max_error = max(max_error,err)
                largest_window = max(largest_window,end)
    r_best = math.exp(.25*(math.log(F(.5,1.,.25))+math.log(F(.5,1.,-.25))))
    example, example_ratios, _, example_end = protocol(.5,1.,.25,r_best)
    add('actual_three_observer_pulse_intercepts', max_error < 1e-11,
        scenarios=18, actual_events=total_events, max_formula_residual=max_error,
        common_window_max=largest_window, right_partner_crossing_time=40.,
        example_clock_rate=r_best, example_end=example_end,
        example_ratios=example_ratios, example_event_records=example)

    diff_error = 0.
    for a in np.linspace(0.,1.,17):
        for c in [1.,1.4,2.]:
            for v in [.1,.25,.6]:
                theory = -2*a*(1-a)*v**3/(c*(c*c-a*a*v*v))
                diff_error = max(diff_error,abs(F(a,c,v)-F(a,c,-v)-theory))
    half, one, quarter = Fraction(1,2), Fraction(1), Fraction(1,4)
    fp, fm = F(half,one,quarter), F(half,one,-quarter)
    add('exact_factor_and_rational_witness', diff_error < 1e-12
        and fp==Fraction(27,28) and fm==Fraction(35,36)
        and fp/fm==Fraction(243,245), max_factor_residual=diff_error,
        F_plus=str(fp), F_minus=str(fm), clock_free_ratio=str(fp/fm))

    cases = [
        ('background_common', [(0.,1.),(0.,1.)], True),
        ('background_different', [(0.,1.),(0.,2.)], False),
        ('source_carried_different', [(1.,1.),(1.,2.)], True),
        ('mixed', [(0.,1.),(1.,2.)], False),
        ('intermediate', [(.5,1.)], False)]
    tests = []
    for name,modes,expected in cases:
        logs = [math.log(F(a,c,s*.25)) for a,c in modes for s in [1,-1]]
        cert = interval_certificate(logs,[1e-12]*len(logs))
        assert cert['feasible']==expected
        tests.append(dict(name=name,modes=modes,certificate=cert))
    add('multimode_shared_clock_branch_tests', True, tests=tests,
        finite_checks_are_not_exhaustive_proof=True)

    v = .25
    # Different clocks can conceal precisely the inconsistency being tested.
    separate_side = [math.sqrt(F(.5,1.,s*v)) for s in [1,-1]]
    side_res = max(abs(F(.5,1.,s*v)/r**2-1)
                   for s,r in zip([1,-1],separate_side))
    separate_mode = [math.sqrt(1-v*v/(c*c)) for c in [1.,2.]]
    mode_res = max(abs(F(0.,c,v)/r**2-1)
                   for c,r in zip([1.,2.],separate_mode))
    add('removing_shared_clock_hides_conflict', side_res<1e-12 and mode_res<1e-12
        and separate_side[0]!=separate_side[1] and separate_mode[0]!=separate_mode[1],
        side_dependent_clocks=separate_side, mode_dependent_clocks=separate_mode,
        largest_residual=max(side_res,mode_res))

    logs = [math.log(float(fp)),math.log(float(fm))]
    floor = .5*math.log(245/243)
    z = .5*sum(logs)
    residual = max(abs(b-z) for b in logs)
    fail = interval_certificate(logs,[.9*floor]*2)
    pass_ = interval_certificate(logs,[1.1*floor]*2)
    # Independently perturb the two observed log ratios within a declared bound.
    observed = np.array(logs)-z+np.array([1e-4,-2e-4])
    error = np.array([1e-4,2e-4])
    observed_contrast_lower = abs(observed[0]-observed[1])-sum(error)
    add('finite_error_interval_certificate', abs(floor-residual)<1e-14
        and not fail['feasible'] and pass_['feasible']
        and observed_contrast_lower > 2*.8*floor,
        min_uniform_log_tolerance=floor, infeasible=fail, feasible=pass_,
        synthetic_observed_log_ratios=observed.tolist(), declared_error=error.tolist(),
        contrast_lower=float(observed_contrast_lower),
        not_experimental_data=True, interval_theorem_reused_from_round963=True)

    # Compare the same actual ballistic emission/arrival events in three charts.
    events, ratios, _, _ = protocol(1.,1.4,.25,1.)
    frames = [(-10.,0.),(0.,.25),(10.,-.4)]
    matrices = []
    for a,u in frames:
        matrices.append(np.array([[1.,0.,0.],[-u,1.,-a],[0.,0.,1.]]))
    covariance_error = 0.
    source_velocities = {'A':0.,'B':.25,'C':0.}
    for e in events:
        for (a,u),m in zip(frames,matrices):
            emit = m@np.array([e['emission_time'],e['emission_x'],1.])
            recv = m@np.array([e['arrival_time'],e['arrival_x'],1.])
            dt = recv[0]-emit[0]
            assert dt>0
            measured = (recv[1]-emit[1])/dt
            expected = source_velocities[e['source']]-u+e['direction']*1.4
            covariance_error = max(covariance_error,abs(measured-expected))
    cocycle_error = 0.
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ij=matrices[i]@np.linalg.inv(matrices[j])
                jk=matrices[j]@np.linalg.inv(matrices[k])
                ik=matrices[i]@np.linalg.inv(matrices[k])
                cocycle_error=max(cocycle_error,float(np.max(np.abs(ij@jk-ik))))
    add('ballistic_same_events_three_chart_consistency', covariance_error<1e-11
        and cocycle_error<1e-11 and max(abs(q-1) for q in ratios.values())<1e-11,
        pulse_covariance_residual=covariance_error, chart_cocycle_residual=cocycle_error,
        one_signal_velocities_in_three_frames=[1.4-u for _,u in frames],
        full_cognitive_joint_model_certified=False)

    arrivals = np.array([e['arrival_tick'] for e in example
                         if e['source']=='B' and e['receiver']=='A'])
    constant = arrivals+.375
    varying = arrivals+np.array([0.,.001,.004])
    constant_error=float(np.max(np.abs(np.diff(constant)-np.diff(arrivals))))
    varying_error=float(np.max(np.abs(np.diff(varying)-np.diff(arrivals))))
    # Physical reflection changes both v and n, leaving q=n*v unchanged.
    add('calibration_and_reflection_controls', constant_error<1e-12
        and varying_error>1e-3 and (-1)*(-.25)==1*.25,
        constant_delay_interval_residual=constant_error,
        varying_delay_interval_change=varying_error,
        reflection_preserves_radial_q=True,
        actual_right_and_left_tasks_required=True)

    return dict(schema='round1094_same_clock_source_law_v1',round=1094,status='PASS',
        numpy_version=np.__version__,groups=groups,
        scope=dict(theorem_conditional_on_affine_emission_and_actual_reciprocity=True,
            all_six_protocols_implemented=False, whole_conjecture_decided=False,
            Lorentz_derived=False,new_adopted_axioms=0,scientific_count_increment=0,
            scientific_count_total=3860))


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a: compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert isinstance(b,(list,tuple))
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11),(a,b)
    else: assert a==b,(a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
    else: compare(result,json.loads(path.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1094,status='PASS',groups=len(result['groups']),
        saved_result_matches=True,result_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        whole_conjecture_decided=False),ensure_ascii=False))


if __name__=='__main__': main()
