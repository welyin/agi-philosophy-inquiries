"""1091 finite calibration: autonomous relative-contact clock and payload model.

Only stipulated finite typed-family primitives are tested. These computations do
not prove all axioms, discover natural laws, or exclude Lorentz reinterpretations.
Default recomputes and compares without writes; --write creates results once.
"""
from __future__ import annotations
import argparse
import functools
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SEED = 109120261010
EPS = .001
OMEGA = .2
EXACT_TOL = 3e-10
NS = (256, 512, 1024)
I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.diag([1.,-1.]).astype(complex)
HC = OMEGA/2 * np.kron(Z, I2)
CNOTY = np.kron((I2+Y)/2, X) + np.kron((I2-Y)/2, I2)
KY = np.pi/2 * (np.eye(4)-CNOTY)
SWAP = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
KS = np.pi/2 * (np.eye(4)-SWAP)
PLUS = np.array([1.,1.], dtype=complex)/np.sqrt(2)
ZERO = np.array([1.,0.], dtype=complex)
CD0 = np.kron(PLUS, ZERO)
D1 = np.kron(I2, np.diag([0.,1.]))


def maxabs(a):
    return float(np.max(np.abs(a))) if np.size(a) else 0.


def opnorm(a):
    return float(np.linalg.norm(a, ord=2))


def exp_minus_i(h, dt=1.):
    vals, vecs = np.linalg.eigh((h+h.conj().T)/2)
    return (vecs*np.exp(-1j*dt*vals)) @ vecs.conj().T


def relative(z):
    xs,w,xr,v = z
    return xs-xr, w-v


def flow(z,t):
    xs,w,xr,v = z
    return (xs+w*t,w,xr+v*t,v)


def boost(z,u,t):
    xs,w,xr,v = z
    return (xs-u*t,w-u,xr-u*t,v-u)


def key(z):
    return tuple(round(float(v),12) for v in z)


def pulse_window(z):
    r,s = relative(z)
    if s == 0:
        return None
    tc = -r/s
    width = EPS/abs(s)
    return tc-width,tc+width


def coupling(z,t):
    r,s = relative(z)
    r += s*t
    if abs(r) >= EPS or s == 0:
        return 0.
    return abs(s)/EPS * .75*(1-(r/EPS)**2)


def area_cdf(v):
    v = min(1.,max(-1.,v))
    return .75*(v-v**3/3)+.5


def exact_area(z,duration):
    r,s = relative(z)
    if s == 0:
        return 0.
    tc = -r/s
    width = EPS/abs(s)
    return area_cdf((duration-tc)/width)-area_cdf(-tc/width)


def swap_gate(alpha):
    phase = np.pi*alpha/2
    return np.exp(-1j*phase)*(np.cos(phase)*np.eye(4)+1j*np.sin(phase)*SWAP)


@functools.lru_cache(maxsize=None)
def clock_gate(z,duration,n):
    window = pulse_window(z)
    if window is None:
        return exp_minus_i(HC,duration)
    lo,hi = max(0.,window[0]),min(duration,window[1])
    if hi <= lo:
        return exp_minus_i(HC,duration)
    u = exp_minus_i(HC,lo)
    dt = (hi-lo)/n
    for k in range(n):
        tmid = lo+(k+.5)*dt
        u = exp_minus_i(HC+coupling(z,tmid)*KY,dt) @ u
    return exp_minus_i(HC,duration-hi) @ u


def block_diag(blocks):
    dims = [b.shape[0] for b in blocks]
    out = np.zeros((sum(dims),sum(dims)), dtype=complex)
    at = 0
    for b in blocks:
        d = b.shape[0]
        out[at:at+d,at:at+d] = b
        at += d
    return out


def embedding(source,target,internal=1):
    targetkeys = [key(z) for z in target]
    j = np.zeros((len(target),len(source)), dtype=complex)
    for i,z in enumerate(source):
        j[targetkeys.index(key(z)),i] = 1
    return np.kron(j,np.eye(internal))


def typed_gate(configs,duration,n):
    images = [flow(z,duration) for z in configs]
    target = sorted(images,key=key)
    internal = [np.kron(clock_gate(tuple(z),duration,n),swap_gate(exact_area(z,duration))) for z in configs]
    return embedding(images,target,16) @ block_diag(internal),target


def typed_boost(configs,u,t):
    images = [boost(z,u,t) for z in configs]
    target = sorted(images,key=key)
    return embedding(images,target,16),target


def density(rng,d):
    a = rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    a = a @ a.conj().T
    return a/np.trace(a)


def partial_keep(rho,dims,keep):
    dims = list(dims)
    data = rho.reshape(dims+dims)
    for idx in reversed(range(len(dims))):
        if idx not in keep:
            data = np.trace(data,axis1=idx,axis2=idx+len(dims))
            dims.pop(idx)
    d = int(np.prod(dims))
    return data.reshape(d,d)


def p_record(gate):
    psi = gate @ CD0
    return float(np.vdot(psi,D1@psi).real)


def quadrature_bound(n):
    # ||Ky||=pi; sup |g'| <= 3 s^2/(2 EPS^2).
    # The integrated midpoint freezing error is <= 3/(2n).
    return 3*np.pi/(2*n)


def run():
    rng = np.random.default_rng(SEED)
    checks = []
    def checked(name,residual,bound=EXACT_TOL,detail=''):
        residual,bound = float(residual),float(bound)
        if not np.isfinite(residual) or residual > bound:
            raise AssertionError(f'{name}: {residual} > {bound}')
        checks.append({'name':name,'residual':residual,'allowed_bound':bound,'status':'PASS','detail':detail})

    checked('contact_generators_are_declared_gates',max(maxabs(exp_minus_i(KY)-CNOTY),maxabs(exp_minus_i(KS)-SWAP)))
    checked('pulse_payload_exponential_formula',max(maxabs(exp_minus_i(KS,a)-swap_gate(a)) for a in [0.,.125,.5,1.]))
    configs = [(0.,2.,1.,0.),(0.,2.,1.,1.)]
    duration = 1.25
    rows = []
    for z in configs:
        r,s = relative(z)
        tc = -r/s
        gates = [clock_gate(z,duration,n) for n in NS]
        pvals = [p_record(u) for u in gates]
        ideal = (1+np.sin(OMEGA*tc))/2
        protocol_bound = 2*OMEGA*EPS/abs(s)
        diff1,diff2 = opnorm(gates[1]-gates[0]),opnorm(gates[2]-gates[1])
        checked(f'clock_unitarity_s_{s:g}',max(maxabs(u.conj().T@u-np.eye(4)) for u in gates))
        checked(f'clock_midpoint_refinement_s_{s:g}',diff2,.35*diff1,'Convergence calibration only; not an error proof.')
        checked(f'clock_record_vs_ideal_with_independent_bounds_s_{s:g}',abs(pvals[-1]-ideal),protocol_bound+quadrature_bound(NS[-1]),'Finite contact-window bound plus rigorous Duhamel midpoint bound; no ideal instantaneous gate is substituted for the autonomous pulse.')
        checked(f'full_diameter_area_s_{s:g}',abs(exact_area(z,duration)-1))
        rows.append({'configuration':list(z),'relative_speed':s,'contact_center':tc,
                     'n_values':list(NS),'record_probabilities':pvals,'ideal_record_probability':float(ideal),
                     'operator_refinement_differences':[diff1,diff2],
                     'finite_contact_window_bound':protocol_bound,
                     'midpoint_Duhamel_bound_at_1024':quadrature_bound(NS[-1]),
                     'exact_payload_area':exact_area(z,duration)})

    # Independent direct 16-dimensional midpoint evolution checks tensor factorization.
    z = configs[0]
    lo,hi = pulse_window(z)
    n = NS[0]
    dt = (hi-lo)/n
    direct = exp_minus_i(np.kron(HC,np.eye(4)),lo)
    approximate_area = 0.
    for k in range(n):
        gm = coupling(z,lo+(k+.5)*dt)
        h = np.kron(HC,np.eye(4))+gm*(np.kron(KY,np.eye(4))+np.kron(np.eye(4),KS))
        direct = exp_minus_i(h,dt) @ direct
        approximate_area += gm*dt
    direct = exp_minus_i(np.kron(HC,np.eye(4)),duration-hi) @ direct
    factored = np.kron(clock_gate(z,duration,n),swap_gate(approximate_area))
    checked('direct_full_H_midpoint_factorization',maxabs(direct-factored),3e-9,
            'Direct numerical pulse area is retained here. Exact area one is used only for the analytically solved payload sector.')
    checked('midpoint_area_bias_is_not_hidden',abs(approximate_area-(1+1/(2*n*n))),2e-10)

    # Unknown coherent configuration and payload/reference input; C,D,B are prepared internally.
    initial_small = density(rng,8) # config, A, external reference
    insert = np.zeros((64,8),dtype=complex)
    for p in range(2):
        for a in range(2):
            for ref in range(2):
                col = (p*2+a)*2+ref
                for c in range(2):
                    row = (((((p*2+c)*2+0)*2+a)*2+0)*2+ref)
                    insert[row,col] = 1/np.sqrt(2)
    initial = insert @ initial_small @ insert.conj().T
    u,targets = typed_gate(configs,duration,NS[-1])
    ur = np.kron(u,I2)
    output = ur @ initial @ ur.conj().T
    expected_payload_ref = partial_keep(initial_small,[2,2,2],[1,2])
    actual_payload_ref = partial_keep(output,[2,2,2,2,2,2],[4,5])
    checked('unknown_payload_and_reference_reach_receiver',maxabs(actual_payload_ref-expected_payload_ref),3e-9,
            'The reduced coherent configuration may change when clock records are discarded; only the receiver payload/reference marginal is identified here.')
    bell = np.array([1.,0.,0.,1.],dtype=complex)/np.sqrt(2)
    abr = np.zeros(8,dtype=complex)
    abr[0]=abr[5]=1/np.sqrt(2)
    received = np.kron(SWAP,I2)@abr
    received_rho = np.outer(received,received.conj())
    checked('Bell_reference_payload_swap',maxabs(partial_keep(received_rho,[2,2,2],[1,2])-np.outer(bell,bell.conj())))

    # Full instrument/propagator naturality under a non-prefix finite directory embedding.
    bigconfigs = [configs[1],(0.,3.,1.,.5),configs[0]]
    ub,bigtargets = typed_gate(bigconfigs,duration,NS[-1])
    jin,jout = embedding(configs,bigconfigs,16),embedding(targets,bigtargets,16)
    checked('finite_directory_full_propagator_naturality',maxabs(ub@jin-jout@u))
    jinr,joutr = np.kron(jin,I2),np.kron(jout,I2)
    ubr = np.kron(ub,I2)
    larger_out = ubr @ (jinr@initial@jinr.conj().T) @ ubr.conj().T
    checked('naturality_with_coherent_configuration_and_external_reference',maxabs(larger_out-joutr@output@joutr.conj().T),3e-9)

    # Cocycle at cuts within each contact pulse, with actual changed output types.
    cuts = [.49975,.500125,.9995,1.00025]
    cocycles = []
    for cut in cuts:
        wt,midtypes = typed_gate(configs,cut,NS[-1])
        ws,finaltypes = typed_gate(midtypes,duration-cut,NS[-1])
        align = embedding(finaltypes,targets,16)
        residual = opnorm(align@ws@wt-u)
        allowance = 3*quadrature_bound(NS[-1])
        checked(f'autonomous_typed_cocycle_cut_{cut:g}',residual,allowance,
                'Comparison of time-ordered numerical propagators; three independent midpoint errors bound the discrepancy. A cut lies inside a nonzero pulse.')
        cocycles.append({'cut':cut,'operator_residual':residual,'conservative_bound':allowance})

    # Co-boost every velocity, and compare different source/output catalogue types.
    boosts = []
    for velocity in [.25,.5,-.75]:
        bin0,boostedtypes = typed_boost(configs,velocity,0.)
        wboost,boostedfinal = typed_gate(boostedtypes,duration,NS[-1])
        bout,boostedtargets = typed_boost(targets,velocity,duration)
        align = embedding(boostedfinal,boostedtargets,16)
        residual = maxabs(bout@u-align@wboost@bin0)
        relation = max(maxabs(np.array(boost(flow(z,duration),velocity,duration))-np.array(flow(boost(z,velocity,0.),duration))) for z in configs)
        checked(f'co_boost_full_typed_intertwining_{velocity:g}',max(residual,relation),3e-9,
                'This checks the declared Galilean time/synchronization contract, not a claim that every physical clock must follow it.')
        boosts.append({'velocity_shift':velocity,'residual':residual,'position_flow_residual':relation})

    # Finite moving clock arrays, separate from the earlier slow-signal example.
    array_velocity,array_T,array_L = .5,2.,1.
    array_configs = [(-1.,.5,0.,0.),(0.,.5,1.,0.)]
    read_T = 2.01
    arrays = [clock_gate(z,read_T,NS[-1]) for z in array_configs]
    array_ps = [p_record(a) for a in arrays]
    ideal_array_p = (1+np.sin(OMEGA*array_T))/2
    checked('matching_clock_arrays_identical_complete_contact_propagators',maxabs(arrays[0]-arrays[1]),EXACT_TOL,
            'Both have the identical relative history -1+0.5t. Each of the four clocks has the same H and same internal initial state; matching roles alone interact.')
    checked('matching_array_joint_record_distributions',maxabs(np.kron([1-array_ps[0],array_ps[0]],[1-array_ps[0],array_ps[0]])-np.kron([1-array_ps[1],array_ps[1]],[1-array_ps[1],array_ps[1]])))
    array_window_bound = 2*OMEGA*EPS/array_velocity
    prob_bound = array_window_bound+quadrature_bound(NS[-1])
    derivative_lower = OMEGA/2*np.cos(OMEGA*2.5) # monotone calibration interval [1.5,2.5]
    one_clock_time_error = prob_bound/derivative_lower
    if one_clock_time_error >= .5:
        raise AssertionError('calibration error leaves its monotonic interval')
    checked('array_record_matches_calibration_with_contact_and_numerical_margin',max(abs(p-ideal_array_p) for p in array_ps),prob_bound)
    candidates = []
    for c in [1.,2.,10.]:
        gamma = 1/np.sqrt(1-(array_velocity/c)**2)
        lorentz_dt = -gamma*array_velocity*array_L/c**2
        actual_dt = 0.
        margin = abs(lorentz_dt)-2*one_clock_time_error
        candidates.append({'c':c,'gamma':float(gamma),'standard_Lorentz_delta_t_prime':float(lorentz_dt),
                           'declared_array_delta_t_prime':actual_dt,'two_clock_conservative_time_error':float(2*one_clock_time_error),
                           'finite_numerical_margin':float(margin),'positive_margin':bool(margin>0)})
    if [c['positive_margin'] for c in candidates] != [True,True,False]:
        raise AssertionError('unexpected finite-margin classification')
    checked('finite_c_clock_array_comparison_scope',0.,detail='For c=1,2 this deterministic integration/contact bound leaves a margin; c=10 does not. This is not a uniform finite-error exclusion as c tends to infinity, and no experimental sampling confidence is claimed.')

    costs = []
    for w in [2.,3.,8.]:
        source_velocity = .25
        cost = math.ceil(1+abs(w-source_velocity)**2)
        shifted_cost = math.ceil(1+abs((w-.5)-(source_velocity-.5))**2)
        if cost != shifted_cost or not np.isfinite(cost):
            raise AssertionError('finite relative preparation cost failed')
        costs.append({'w':w,'emitter_velocity':source_velocity,'declared_preparation_cost':cost})
    checked('finite_internal_preparation_budget_samples',0.,detail='Budget is a stipulated finite resource count, not an energy law; no infinite velocity or infinite task is implemented.')

    return {'schema':'research_round1091_autonomous_contact_v1','round':1091,'status':'PASS','seed':SEED,
            'numpy_version':np.__version__,'parameters':{'epsilon':EPS,'omega':OMEGA,'n_values':list(NS)},
            'scope':{'finite_calibration_not_universal_axiom_proof':True,
                     'autonomous_relative_contact_rule_is_model_input':True,
                     'exact_payload_sector_solved_analytically':True,
                     'clock_sector_time_ordered_and_numerically_approximated':True,
                     'slow_signal_alone_is_not_a_relativistic_counterexample':True,
                     'clock_array_exclusion_uses_fixed_actual_synchronization_and_units':True,
                     'all_Lorentz_reinterpretations_excluded':False,
                     'uniform_finite_error_exclusion_for_c_to_infinity':False,
                     'microphysics_of_the_universe_claimed':False},
            'arrival_contact_calibrations':rows,'cocycle_checks':cocycles,'co_boost_checks':boosts,
            'clock_arrays':{'u':array_velocity,'T':array_T,'L':array_L,'event_delta_t':0.,'event_delta_x':array_L,
                            'matching_configurations':[list(z) for z in array_configs],'record_probabilities':array_ps,
                            'ideal_record_probability':float(ideal_array_p),'finite_window_probability_bound':array_window_bound,
                            'Duhamel_numerical_probability_bound':quadrature_bound(NS[-1]),
                            'calibration_interval':[1.5,2.5],'probability_derivative_lower_bound':float(derivative_lower),
                            'one_clock_time_error_bound':float(one_clock_time_error),'candidate_c_comparisons':candidates},
            'resource_samples':costs,'check_groups':len(checks),'checks':checks}


def compare(actual,saved,path='root'):
    if isinstance(actual,dict):
        if actual.keys()!=saved.keys():
            raise AssertionError(path+': keys differ')
        for key0 in actual:
            compare(actual[key0],saved[key0],path+'.'+key0)
    elif isinstance(actual,list):
        if len(actual)!=len(saved):
            raise AssertionError(path+': lengths differ')
        for i,(a,b) in enumerate(zip(actual,saved)):
            compare(a,b,f'{path}[{i}]')
    elif isinstance(actual,float):
        if not np.isfinite(saved) or abs(actual-saved)>2e-10:
            raise AssertionError(path+': numeric value differs')
    elif actual!=saved:
        raise AssertionError(path+': value differs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    target=ROOT/'results.json'
    if args.write and target.exists():
        raise SystemExit('Refusing to overwrite existing results.json')
    result=run()
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2,allow_nan=False)
            handle.write('\n')
        mode='created'
    else:
        if not target.exists():
            raise SystemExit('Missing results.json; use --write once for initial creation')
        compare(result,json.loads(target.read_text(encoding='utf-8')))
        mode='read_only_recompute_match'
    print(json.dumps({'round':1091,'status':'PASS','mode':mode,'check_groups':result['check_groups'],
                      'arrival_record_probabilities':[r['record_probabilities'][-1] for r in result['arrival_contact_calibrations']],
                      'array_record_probabilities':result['clock_arrays']['record_probabilities'],
                      'finite_c_positive_margins':[r['positive_margin'] for r in result['clock_arrays']['candidate_c_comparisons']],
                      'all_Lorentz_reinterpretations_excluded':False},ensure_ascii=False))


if __name__=='__main__':
    main()
