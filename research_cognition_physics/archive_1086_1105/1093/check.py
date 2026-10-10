"""Finite diagnostics of the conditional bridge, not a physical simulation.

Default: recompute and compare the saved results without changing any file.
--write: save the initial result; refuses to replace existing results.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
SEED = 109320261010


def tau(r, c0=1.0, T0=1.0):
    # Rationalized inverse avoids cancellation for small positive distances.
    return 2.0*r/(c0*(1.0+np.sqrt(1.0+4.0*r/(c0*T0))))


def fan(r, L):
    k = math.ceil(L/r)
    v = np.zeros((k, 3))
    if k == 1:
        v[0, 0] = r
    else:
        v[:-2, 0] = r
        a = (L-(k-2)*r)/2.0
        b = math.sqrt(max(0.0, r*r-a*a))
        v[-2] = [a, b, 0.0]
        v[-1] = [a, -b, 0.0]
    return v


def pair_marginal(psi, dims, site, reference):
    axes = [site, reference]+[i for i in range(len(dims)) if i not in (site, reference)]
    a = psi.reshape(dims).transpose(axes).reshape(dims[site]*dims[reference], -1)
    return a@a.conj().T


def calculate():
    rng = np.random.default_rng(SEED)
    groups = []
    def add(name, passed, **data):
        assert passed, (name, data)
        groups.append(dict(name=name, passed=bool(passed), **data))

    # This checks the geometric construction; physical copy and relay
    # permissions remain hypotheses and are not inferred from these vectors.
    radii = np.unique(np.r_[np.geomspace(1e-4, 1.0, 65), [0.5, 0.75, 0.999, 1/3]])
    endpoint_error = length_error = factor_excess = 0.0
    for r in radii:
        v = fan(float(r), 1.0)
        endpoint_error = max(endpoint_error, float(np.linalg.norm(v.sum(axis=0)-[1,0,0])))
        length_error = max(length_error, float(np.max(np.abs(np.linalg.norm(v, axis=1)-r))))
        factor_excess = max(factor_excess, len(v)*float(r)-2.0)
    add('finite_baseline_rotated_copies', endpoint_error < 1e-10 and length_error < 1e-12
        and factor_excess == 0.0, cases=len(radii), max_copies=max(math.ceil(1/r) for r in radii),
        endpoint_residual=endpoint_error, segment_length_residual=length_error)

    # Finite certificates of the derived inequality, conditional on D4.
    L, tauL, A = 1.0, 0.7, 1.6
    local_bound = 2*L/tauL
    maximum_ratio = 0.0
    for _ in range(64):
        v = rng.normal(size=(12,3))
        v *= rng.uniform(0.01, L, size=(12,1))/np.linalg.norm(v, axis=1)[:,None]
        durations = np.linalg.norm(v, axis=1)/local_bound+rng.uniform(0,0.1,size=12)
        parent_time = float(durations.sum()/A)
        maximum_ratio = max(maximum_ratio, float(np.linalg.norm(v.sum(axis=0))/parent_time))
    add('conditional_refinement_inequality', maximum_ratio <= A*local_bound,
        cases=64, speed_envelope=A*local_bound, maximum_sampled_ratio=maximum_ratio,
        physical_refinement_permission_verified=False)

    c0, T0 = 2.0, 3.0
    def f(t):
        return c0*t+c0*t*t/T0
    max_violation = 0.0
    for _ in range(256):
        t = np.exp(rng.uniform(-8,5,size=5))
        vectors = rng.normal(size=(5,3))
        vectors *= (f(t)*rng.uniform(0,1,size=5)/np.linalg.norm(vectors,axis=1))[:,None]
        displacement = float(np.linalg.norm(vectors.sum(axis=0)))
        max_violation = max(max_violation, displacement-f(t.sum()),
                            float(tau(displacement,c0,T0)-t.sum()))
    times = np.geomspace(1e-8,1e5,100)
    inverse_residual = float(np.max(np.abs(tau(f(times),c0,T0)-times)/times))
    add('quadratic_permission_serial_closure', max_violation == 0.0 and inverse_residual<1e-12,
        cases=256, relative_inverse_residual=inverse_residual,
        autonomous_physical_realization=False)

    metric_violation = 0.0
    for _ in range(256):
        a,b = rng.normal(size=(2,3))*rng.uniform(0,100)
        metric_violation = max(metric_violation, float(tau(np.linalg.norm(a+b),c0,T0)
            -tau(np.linalg.norm(a),c0,T0)-tau(np.linalg.norm(b),c0,T0)))
    # Every subinterval of a path with |v|<=c0 is allowed.
    dt = np.geomspace(1e-8,1e3,101)
    add('positive_time_metric_and_local_signal_subfamily', metric_violation == 0.0
        and np.all(tau(dt,c0,T0)>0) and np.all(c0*dt<=f(dt)), triangle_samples=256,
        continuous_subfamily_max_speed=c0,
        long_edges_have_free_interception_permission=False)

    L = 4.0
    cL = L/tau(L,c0,T0)
    long_tasks = []
    for t in [1.0,10.0,100.0,1000.0,10000.0]:
        r = f(t)
        long_tasks.append(dict(time=t, length=r, average_speed=r/t,
            required_refinement_factor_lower_bound=r/(cL*t)))
    add('failure_of_uniform_refinement', all(b['average_speed']>a['average_speed']
        and b['required_refinement_factor_lower_bound']>a['required_refinement_factor_lower_bound']
        for a,b in zip(long_tasks,long_tasks[1:])), L=L, local_optimal_speed=cL,
        long_task_certificates=long_tasks,
        unboundedness_proved_analytically_not_by_samples=True)

    # Pure global random inputs are entangled across ALL sites plus R.
    # Adjacent swaps preserve all unknown registers: no blank or cloning.
    transfer_error = recovery_error = 0.0
    for sites in [2,3,5,7]:
        dims = [2]*sites+[3]
        psi = rng.normal(size=math.prod(dims))+1j*rng.normal(size=math.prod(dims))
        psi /= np.linalg.norm(psi)
        initial = psi.copy()
        rhoQR = pair_marginal(psi,dims,0,sites)
        for j in range(sites-1):
            psi = psi.reshape(dims).swapaxes(j,j+1).reshape(-1)
        transfer_error = max(transfer_error,float(np.linalg.norm(
            pair_marginal(psi,dims,sites-1,sites)-rhoQR)))
        for j in reversed(range(sites-1)):
            psi = psi.reshape(dims).swapaxes(j,j+1).reshape(-1)
        recovery_error = max(recovery_error,float(np.linalg.norm(psi-initial)))
    add('full_unknown_payload_reference_and_private_registers', transfer_error<1e-12
        and recovery_error<1e-12, site_counts=[2,3,5,7], reference_dimension=3,
        reduced_reference_residual=transfer_error, full_recovery_residual=recovery_error,
        physical_duration_or_reversal_permission_inferred=False)

    # Dephasing channels give an exactly computable composition error.
    # Half diamond distance to identity is p; Bell state attains it.
    eps, delta = 0.02, 0.01
    bell = np.array([1,0,0,1],complex)/math.sqrt(2)
    phase_bell = np.array([1,0,0,-1],complex)/math.sqrt(2)
    rho = np.outer(bell,bell.conj())
    flipped = np.outer(phase_bell,phase_bell.conj())
    budgets = []
    error_residual = 0.0
    for k in [2,8,64]:
        fixed_p = (1-(1-2*delta)**k)/2
        allocated_p = (1-(1-2*eps/k)**k)/2
        output = (1-allocated_p)*rho+allocated_p*flipped
        trace_distance = float(np.abs(np.linalg.eigvalsh(output-rho)).sum()/2)
        error_residual=max(error_residual,abs(trace_distance-allocated_p))
        budgets.append(dict(copies=k, fixed_per_hop_error=delta, fixed_total_error=fixed_p,
            allocated_per_hop_error=eps/k, allocated_total_error=allocated_p))
    add('composable_error_budget_and_fixed_error_failure', error_residual<1e-12
        and all(b['allocated_total_error']<=eps+1e-14 for b in budgets)
        and budgets[-1]['fixed_total_error']>eps, target_error=eps, certificates=budgets,
        trace_distance_residual=error_residual, uniform_precision_time_contract_proved=False)

    k, r, L, Delta, c = 7, 0.16, 1.0, 0.02, 2.0
    each_error, each_location_error = eps/k, Delta/k
    total_time_upper = k*(r/c+0.001)
    annular_time_lower = (L-Delta)/c  # Only for the separately stated linear example.
    add('finite_annular_time_certificate', k*each_error<=eps+1e-14
        and k*each_location_error<=Delta+1e-14 and annular_time_lower<=total_time_upper,
        target_error=eps, location_tolerance=Delta, time_lower=annular_time_lower,
        certified_time_upper=total_time_upper,
        source_of_lower_bound='explicit_linear_menu_example_not_derived_from_nine_axioms',
        total_time_upper_is_not_actual_speed_upper_bound=True)

    t, x, u = 1.5, 4.0, 3.0
    add('order_does_not_fix_signal_speed', t>0 and x/t != (x-u*t)/t,
        original_dt=t, transformed_dt=t, original_velocity=x/t,
        transformed_velocity=(x-u*t)/t, physical_boost_equivalence_proved=False)
    return dict(schema='round1093_nondegenerate_refinement_bridge_v1', round=1093,
        status='PASS', seed=SEED, numpy_version=np.__version__, groups=groups,
        scope=dict(new_adopted_axioms=0, scientific_count_increment=0,
            scientific_count_total=3860, theorem_conditional_on_D1_D4=True,
            joint_FUCP_CO_counterexample=False, entire_conjecture_decided=False,
            Lorentz_derived=False, universal_weak_influence_front_derived=False))


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for key in a: compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-11), (a,b)
    else:
        assert a==b, (a,b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args=parser.parse_args()
    output=calculate()
    result=HERE/'results.json'
    if args.write:
        with result.open('x',encoding='utf8') as stream:
            json.dump(output,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
    else:
        compare(output,json.loads(result.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1093,status='PASS',groups=len(output['groups']),
        result_sha256=hashlib.sha256(result.read_bytes()).hexdigest(),
        saved_result_matches=True, entire_conjecture_decided=False),ensure_ascii=False))


if __name__=='__main__':
    main()
