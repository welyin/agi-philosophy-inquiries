"""1251 independent exact finite-statistics and weak-residual certificate.

No author functions, old checker functions, ODE solver or wavefunction are
imported.  Frozen old JSON is used only as an inherited bound/certificate.
The source-dictionary example is a completion-symbol identity, not a new
physical trajectory. --write is exclusive; ordinary replay is read-only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULT = HERE / 'independent_results.json'
ETA = F(1, 10**8)
TAU = F(5, 10**9)
SYSTEM = F(4, 10**9)
BIN = F(1, 10**10)
DELTA = F(1, 2**40)
EPS_K = F(1, 10**4)


def pack(x):
    if isinstance(x, F):
        return {'exact': str(x), 'decimal': float(x)}
    if isinstance(x, dict):
        return {k: pack(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [pack(v) for v in x]
    return x


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def old_quantum():
    old = json.loads((ROOT/'1249/quantum_bounds_results.json').read_text(encoding='utf-8'))
    return old['group_constructive_finite_resource_bound']['ledger']


def finite_symbol_and_statistics():
    old = old_quantum()
    beta_norm = factorial(81)//factorial(40)**2
    B = 324*beta_norm*2**40
    Lcut = 56*B
    M = 2**98*(1+Lcut)**2
    Lsource = 4*(Lcut+14)
    assert B == int(old['B_decimal'])
    assert Lcut == int(old['Lcut_decimal'])
    # Product derivative rule for a_c=chi exp(Q), Q<=3 on support.
    # The exponential series gives e<3; hence exp(3)<27.
    # ||D^r(chi exp Q)|| <= 27 r! (Lcut+1)^r, including r=0.
    derivative_ratios = []
    for r in range(21):
        numerator = 27*factorial(r)*(Lcut+1)**r
        denominator = M*factorial(r)*Lsource**r
        assert numerator < denominator
        derivative_ratios.append(F(numerator, denominator))
    assert 48+1 <= 64

    # One finite coarse-grained positive coherent POVM returns a triple.
    # M and M/BIN are integers; clipping does not introduce out-of-range
    # rounded outcomes, and floor bins have deterministic error < BIN.
    imax = M*10**10
    def coarse(x):
        clipped = min(F(M), max(F(-M), x))
        return F((clipped/BIN).numerator//(clipped/BIN).denominator)*BIN
    samples = [-F(M), -F(M)+BIN/3, -BIN/3, F(0), BIN/3,
               F(M)-BIN/3, F(M)]
    for x in samples:
        y = coarse(x)
        assert -M <= y <= M and 0 <= x-y < BIN
    levels = 2*imax+1
    triple_outcomes = levels**3
    bits = (triple_outcomes-1).bit_length()
    assert 2**(bits-1) < triple_outcomes <= 2**bits

    batches, means = 2*129, 2*129*3
    repeats = 26*M*M*(2*10**8)**2
    # Hoeffding's bounded-variable exponent, independent repetitions only.
    # The three outputs within one shot may be correlated.
    exponent = F(repeats)*TAU**2/(2*M*M)
    assert exponent == 13
    # Independent rational majorant: lower Taylor polynomial for exp(13).
    exp13_lower = sum((F(13**k, factorial(k)) for k in range(21)), F(0))
    union_bound = F(2*means)/exp13_lower + 2*DELTA
    declared_bound = 2*means*F(3,8)**13 + 2*DELTA
    contract_union_bound = union_bound+2*EPS_K
    contract_declared_bound = declared_bound+2*EPS_K
    assert union_bound < declared_bound < F(1,100)
    assert contract_union_bound < contract_declared_bound < F(1,100)
    total_menu_error = TAU+SYSTEM+BIN+DELTA
    assert total_menu_error < ETA

    # The complete source dictionary is not replaced with H^2.  Here the
    # phase psi=96 lies halfway across the original 64..128 cutoff ramp.
    # Symmetry of its beta density makes chi=1/2 exactly; all other
    # coordinates are in the cutoff plateau.  This is not a trajectory.
    chi, a, pi_Q, pi_angle = F(1,2), F(2), F(-8), F(1,10)
    Frod = F(1)
    h_gravity = -pi_Q**2/(4*a**3)
    h_rod = F(3,20)*a*Frod
    h_material_without_J = h_rod+F(72,100)+pi_angle/10
    h = h_gravity+h_material_without_J
    Jc = -chi*h
    Mc = F(71,100)+Jc
    rho_full = F(3,20)*Frod/a**2 + (Mc+F(1,100)+pi_angle/10)/a**3
    H = -pi_Q/(2*a**3)
    difference = chi*rho_full-chi*H**2
    assert h == -F(97,100)
    assert difference == chi*(1-chi)*h/a**3 == -F(97,3200)
    # S(1/2)=1/2 follows by exact polynomial integration as an extra check.
    beta_half = beta_norm*sum((F((-1)**j*comb(40,j),40+j+1)
                              *F(1,2)**(40+j+1) for j in range(41)),F(0))
    assert beta_half == F(1,2)

    return {
        'status':'PASS',
        'symbol_bound': {'M':str(M),'Lsource':str(Lsource),
            'new_symbol':'a_c=chi*exp(Q)', 'orders':list(range(21)),
            'worst_derivative_ratio':max(derivative_ratios),
            'old_menu_count':48,'new_menu_count':49,'capacity':64},
        'bounded_finite_output': {'bin_width':BIN,'scalar_levels':str(levels),
            'triple_outcomes':str(triple_outcomes),'bits_per_triple':bits,
            'range_preserved':True,'floor_then_clip_cases':len(samples)},
        'statistics': {'batches':batches,'means':means,'repeats_per_batch':str(repeats),
            'hoeffding_exponent':exponent,'statistical_error':TAU,
            'unimplemented_system_mean_budget':SYSTEM,'old_quantum_bias':DELTA,
            'total_menu_error':total_menu_error,'nominal_menu_error':ETA,
            'independent_ideal_joint_failure_upper':union_bound,
            'published_simple_ideal_failure_upper':declared_bound,
            'unimplemented_K_probability_matching_budget_per_branch':EPS_K,
            'independent_contract_joint_failure_upper':contract_union_bound,
            'published_simple_contract_failure_upper':contract_declared_bound,
            'source_mean_bias_does_not_bound_K_classification':True,
            'separate_terminal_K_tasks':2,
            'requires_independent_repetitions_not_independent_triple_entries':True},
        'source_dictionary_deletion':{'chi':chi,'a':a,'psi':96,'h':h,
            'rho_full':rho_full,'H_squared':H**2,
            'chi_rho_minus_chi_H_squared':difference,
            'physical_trajectory_claimed':False},
        '_M':M,'_repeats':repeats,'_bits':bits,
    }


def integrate_coefficients(poly, lo, hi):
    return sum((c*(hi**(i+1)-lo**(i+1))/(i+1)
                for i,c in enumerate(poly)), F(0))


def linear_product(poly, intercept, slope):
    p = [F(0)]*(len(poly)+1)
    for i,c in enumerate(poly):
        p[i] += intercept*c
        p[i+1] += slope*c
    return p


def weak_and_resources(first):
    N, length = 128, F(1,4)
    h, ds = length/N, F(1,N)
    # Integrate hat functions against 30*s^2*(1-s)^2 and its second
    # derivative, rather than calling the inherited weight implementation.
    shape = [F(0),F(0),F(30),F(-60),F(30)]
    second = [F(60),F(-360),F(360)]
    wa, wb = [], []
    for j in range(N+1):
        mid = F(j,N)
        ca, cb = F(0), F(0)
        if j:
            lo = mid-ds
            ca += integrate_coefficients(linear_product(second,-lo/ds,1/ds),lo,mid)/length**2
            cb += integrate_coefficients(linear_product(shape,-lo/ds,1/ds),lo,mid)
        if j < N:
            hi = mid+ds
            ca += integrate_coefficients(linear_product(second,hi/ds,-1/ds),mid,hi)/length**2
            cb += integrate_coefficients(linear_product(shape,hi/ds,-1/ds),mid,hi)
        wa.append(ca); wb.append(cb)
    assert sum(wa)==0 and sum(wb)==1 and min(wb)>0
    assert sum(wa[j]*(2+F(j)*h) for j in range(N+1))==0
    l1=sum(map(abs,wa)); assert l1<384
    inherited=json.loads((ROOT/'1248/results.json').read_text(encoding='utf-8'))
    g=inherited['group_weights_and_error']
    assert wa==[F(x) for x in g['a_weights']]
    assert wb==[F(x) for x in g['b_weights']]
    eta_b=(3*4*ETA+F(3,10)*ETA+4*ETA**2)/2
    E=384*ETA+eta_b+h*h*(384+20)/8
    old_E=F(inherited['group_common_domain_and_separation']['errors']['total_full']['rational'])
    assert E==old_E
    # Baseline-only omission, no estimated K: r0(a)=1/(200 a^2).
    # a and interpolated noisy a >2 imply |r0'(a)|<1/800<1/500.
    assert F(1,800) < F(1,500)
    dropped_noise=ETA/500
    second_r0=F(1,100)*(F(11,100)**2+F(3,8)*F(3,25))
    assert second_r0<F(1,1000)
    Eminus=384*ETA+eta_b+dropped_noise+h*h*(384+21)/8
    gap=F(1,1800)-E-Eminus
    assert gap>F(16,100000)
    # The two records are obtained from the same a data; the exact decrement
    # is positive, and its total noise is bounded by dropped_noise.
    assert F(1,200)*F(1,3)**2==F(1,1800)
    assert F(1,1800)>E+Eminus

    repeats, bits = first['_repeats'], first['_bits']
    total_runs=258*repeats+2
    scalar_resource=258*repeats*bits+4 # two distinct 3-outcome K tasks
    oldN=old_quantum()['integer_n_power_exponent_decimal']
    assert int(oldN)>0 and len(oldN)==923
    return {
        'status':'PASS', 'nodes':N+1,'mesh':h,
        'weight_checks':{'a_sum':sum(wa),'b_sum':sum(wb),'derivative_hat_l1':l1,
            'inherited_all_weights_match':True},
        'budgets':{'eta_b':eta_b,'E_full':E,'baseline_omission_noise':dropped_noise,
            'baseline_second_derivative_bound':second_r0,'E_baseline_omitted':Eminus,
            'strict_omission_lower':F(1,1800),'strict_interval_gap':gap,
            'K_squared_retained_and_not_estimated':True},
        'finite_resources':{'all_source_runs':str(258*repeats),
            'separate_old_K_runs':2,'all_runs':str(total_runs),
            'output_storage_capacity_bits':str(scalar_resource),
            'n_formula':'ceil(max(1,hbar))*2^N in fixed original units',
            'N_decimal':oldN,
            'per_run_coordinate_volume':'1000*n^3',
            'sum_declared_coordinate_volumes':f'{total_runs} * 1000*n^3',
            'classical_initial_material_energy_sum':f'{total_runs} * 1110*n^3',
            'initial_quantum_energy':'only inherited matching budgets; not exactly the classical total',
            'physical_preparation_and_measurement_apparatus_source_realized':False,
            'simultaneous_noninteracting_copies_in_one_universe_claimed':False},
    }


def calculate():
    first=finite_symbol_and_statistics()
    second=weak_and_resources(first)
    for key in ('_M','_repeats','_bits'):
        del first[key]
    files={
        'independent_check.py':Path(__file__),
        '1248/proof.md':ROOT/'1248/proof.md',
        '1248/results.json':ROOT/'1248/results.json',
        '1249/canonical.md':ROOT/'1249/canonical.md',
        '1249/quantum_bridge.md':ROOT/'1249/quantum_bridge.md',
        '1249/quantum_bounds_results.json':ROOT/'1249/quantum_bounds_results.json',
    }
    return pack({'round':1251,'status':'PASS','scientific_groups':2,
        'group_1_symbol_statistics_and_complete_source':first,
        'group_2_weak_residual_and_finite_resources':second,
        'inherited_files_sha256':{name:sha(path) for name,path in files.items()},
        'new_actual_instrument_permissions_proved':False,'CT_adopted':False,
        'entire_joint_GR_delivery_completed':False,
        'method':'Independent exact algebra, beta/hat polynomial integration, bounded-statistics union certificate and inherited fixed resources; no author code imported.'})


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); data=calculate()
    if args.write:
        with RESULT.open('x',encoding='utf-8') as f:
            json.dump(data,f,ensure_ascii=False,indent=2); f.write('\n')
        print('INITIAL_SAVE PASS: 2 independent scientific groups')
    else:
        assert data==json.loads(RESULT.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH PASS: 2 independent scientific groups')
    print(json.dumps({'ideal_joint_failure_upper':data['group_1_symbol_statistics_and_complete_source']['statistics']['independent_ideal_joint_failure_upper'],
                     'contract_joint_failure_upper':data['group_1_symbol_statistics_and_complete_source']['statistics']['independent_contract_joint_failure_upper'],
                     'gap':data['group_2_weak_residual_and_finite_resources']['budgets']['strict_interval_gap']},ensure_ascii=False))
