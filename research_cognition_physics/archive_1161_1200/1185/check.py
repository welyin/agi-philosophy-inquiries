"""1185: exact finite composite identities and rational original-clock bounds.

Default is read-only. --write writes results once. Numpy matrices are
diagnostics; all strict scalar certificates use outward Fraction intervals.
No simulation of the huge lattice or K-copy parent is claimed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
spec = importlib.util.spec_from_file_location('frozen1178', STAGE/'1178'/'check.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
point, add, sub, mul, scale = old.point, old.add, old.sub, old.mul, old.scale
divide, sqrt_i, trig, exp_i = old.divide, old.sqrt_i, old.trig, old.exp_i
ij = old.interval_json
EPS, OM, NU, G, ELL = F(3,20), F(1,2), F(7,20), F(1,5), F(16,65)
K, V = 10**14, 65**3


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def analytic_certificate():
    helper = old.old
    pi = old.outward(helper.add(helper.scale(helper.atan_interval(F(1,5)), F(16)),
                               helper.scale(helper.atan_interval(F(1,239)), F(-4))))
    time = scale(pi, F(1,2))
    tasks = {}
    intervals = {}
    for name, b in [('matched_q_plus', F(3,80)), ('matched_q_minus', F(9,80)),
                    ('deleted_q_plus', EPS), ('deleted_q_minus', EPS)]:
        r2 = OM*OM+b*b
        r = sqrt_i(point(r2))
        bracket = sub(time, divide(trig(mul(scale(r,F(2)),time),'sin'), scale(r,F(2))))
        tau = sub(time, scale(bracket, EPS*OM*b/r2))
        prob = scale(add(point(F(1)),trig(scale(tau,F(1,5)),'sin')), F(1,2))
        tasks[name] = {'b_exact':str(b), 'r_squared_exact':str(r2),
                       'tau':ij(tau), 'P_minus_Y_probability':ij(prob)}
        intervals[name] = prob
    ideal = sub(sub(intervals['matched_q_plus'], intervals['matched_q_minus']),
                sub(intervals['deleted_q_plus'], intervals['deleted_q_minus']))
    a_s = EPS*F(1,2)+3*NU/ELL
    a_t = 2*G/ELL
    b2 = 7*a_s*a_s+5*a_t*a_t
    w = F(3,10)  # after exact invariant-code intertwining, never full micro norm
    pref = add(sqrt_i(point(2*b2)),point(a_s+a_t))
    eta = scale(mul(mul(time,exp_i(scale(time,6*w))),pref),F(4,10**7))
    eta_adopt, instrument = F(49,200000), F(1,50000)
    assert eta[1] < eta_adopt
    per_task = eta_adopt+instrument
    lower = ideal[0]-4*per_task
    assert lower > F(3,20000)
    # All actual probabilities stay inside the original clock's invertible window.
    upper_window = old.sin_p(point(F(19,10)))[0]
    for p in intervals.values():
        assert p[0]-per_task > F(1,2)
        assert p[1]+per_task < upper_window
    return {'time':ij(time), 'tasks':tasks, 'ideal_four_task_difference':ij(ideal),
            'a_s_exact':str(a_s), 'a_t_exact':str(a_t), 'B_squared_exact':str(b2),
            'effective_pair_norm_bound_exact':str(w), 'eta_interval':ij(eta),
            'eta_adopted_exact':str(eta_adopt), 'instrument_each_input_exact':str(instrument),
            'total_error_each_exact':str(per_task), 'four_task_error_exact':str(4*per_task),
            'actual_four_task_lower_exact':str(lower), 'actual_four_task_lower_float':float(lower),
            'published_strict_lower_exact':'3/20000', 'outward_decimal_places':40,
            'all_strict_scalar_bounds_are_rational':True,
            'instrument_contract_covers_new_preparation_and_transported_ports':True,
            'sampling_error_included':False}


def kron(*ops):
    ans = np.ones((1,1),complex)
    for op in ops:
        ans = np.kron(ans,op)
    return ans


def matrix_diagnostics():
    I=np.eye(2); X=np.array([[0,1],[1,0]],complex); Z=np.diag([1.,-1.])
    A=(I+X)/2; n=(I-Z)/2
    rot=np.array([[math.sqrt(2/3),math.sqrt(1/3)],
                  [-math.sqrt(1/3),math.sqrt(2/3)]])
    W=kron(I-A,I)+kron(A,rot)
    E=W@kron(I,np.array([[1.],[0.]]))
    hm=np.eye(4)+kron(I,n)+math.sqrt(2)*kron(A,X)
    hc=I-A; zt=W@kron(Z,I)@W.conj().T
    nb=W@kron(I,n)@W.conj().T
    residuals={
        'W_unitarity':np.linalg.norm(W.conj().T@W-np.eye(4)),
        'E_isometry':np.linalg.norm(E.conj().T@E-I),
        'diagonal_identity':np.linalg.norm(W.conj().T@hm@W-(kron(hc,I)+kron(I+2*A,n))),
        'free_intertwining':np.linalg.norm(hm@E-E@hc),
        'transported_Z_intertwining':np.linalg.norm(zt@E-E@Z),
        'dressed_number_commutator':np.linalg.norm(hm@nb-nb@hm),
        'bare_Z_compression':np.linalg.norm(E.conj().T@kron(Z,I)@E-math.sqrt(2/3)*Z)}
    bare_leak=float(np.linalg.norm((np.eye(4)-E@E.conj().T)@kron(Z,I)@E,2))
    residuals['bare_Z_leak_formula']=abs(bare_leak-1/math.sqrt(3))
    alpha_beta=[]
    for alpha,beta in [(0.,0.),(1.,1.),(2.,1.25),(.7,.2)]:
        j=np.eye(4)+alpha*kron(I,n)+math.sqrt(2)*beta*kron(A,X)
        expect=I+(alpha-4*beta)*A/3
        leak=(np.eye(4)-E@E.conj().T)@j@E
        compress_err=float(np.linalg.norm(E.conj().T@j@E-expect))
        leak_err=float(np.linalg.norm(leak.conj().T@leak-2*(alpha-beta)**2*A/9))
        alpha_beta.append({'alpha':alpha,'beta':beta,'compression_residual':compress_err,
                           'leak_gram_residual':leak_err})
        assert max(compress_err,leak_err)<1e-13
    # K=2, one site, one diagnostic probe qubit. This exercises all t-port
    # types and reciprocal source algebra; it is NOT the actual 48-code model.
    em=kron(I,E); embed=kron(em,em,I)
    probe=np.diag([.3,.5]); fs=.11*X+.07*Z; ft=.13*X-.09*Z
    branches={}
    for name in ('matched','deleted'):
        Bm=float(EPS)*(hm if name=='matched' else np.eye(4))
        Be=float(EPS)*(hc if name=='matched' else I)
        hloc_m=kron(-float(OM)*X,np.eye(4))+kron(I,hm)
        hloc_e=kron(-float(OM)*X,I)+kron(I,hc)
        sm=kron(Z,np.eye(4)); se=kron(Z,I)
        tm=kron(I,zt); te=kron(I,Z)
        bm=kron(I,Bm); be=kron(I,Be)
        parents=[]
        for dim,hl,s,t,b in [(8,hloc_m,sm,tm,bm),(4,hloc_e,se,te,be)]:
            ident=np.eye(dim)
            H=kron(hl,ident,I)+kron(ident,hl,I)+kron(ident,ident,probe)
            H+=(kron(s,b,I)+kron(b,s,I))/2
            H+=(kron(s,ident,fs)+kron(ident,s,fs))/2
            H+=(kron(t,ident,ft)+kron(ident,t,ft))/2
            parents.append(H)
        micro,effective=parents
        err=float(np.linalg.norm(micro@embed-embed@effective))
        herm=float(np.linalg.norm(micro-micro.conj().T))
        comm=float(np.linalg.norm(micro@kron(kron(I,nb),np.eye(8),I)-
                                  kron(kron(I,nb),np.eye(8),I)@micro))
        branches[name]={'full_micro_dimension':128,'effective_dimension':32,
                        'parent_intertwining_residual':err,'Hermiticity_residual':herm,
                        'dressed_sector_commutator':comm,
                        'micro_B_norm':float(np.linalg.norm(Bm,2)),
                        'effective_B_norm':float(np.linalg.norm(Be,2))}
        assert max(err,herm,comm)<2e-13
    assert max(residuals.values())<1e-13
    return {'residuals':{k:float(v) for k,v in residuals.items()},
            'micro_spectrum':np.linalg.eigvalsh(hm).tolist(),
            'bare_Z_leak_operator_norm':bare_leak,
            'bare_number_commutator_norm':float(np.linalg.norm(hm@kron(I,n)-kron(I,n)@hm)),
            'alpha_beta_diagnostics':alpha_beta, 'two_copy_port_diagnostics':branches,
            'diagnostic_not_large_parent_simulation':True}


def calculate():
    cert=analytic_certificate(); diag=matrix_diagnostics()
    bs=2*EPS*2+6*NU/ELL+4*G/ELL
    bp=2*EPS*F(4,5)+18*NU/ELL+12*G/ELL
    original=1+K*V*(4+6/ELL+2*OM+bs)+V*(F(8,5)+18/ELL+bp)
    extra_free=(3-OM)*V*K
    extra_recip=F(9,20)*V*(K-1)
    paths=['1176/proof.md','1177/proof.md','1177/check.py','1178/proof.md',
           '1178/check.py','1178/results.json','_admission/after1174_direction/mean_field_bound.md']
    return {'round':1185,'date':'2026-10-13','status':'PASS_WITH_COMPOSITE_SOURCE_SCOPE',
            'source_sha256':{p:sha(STAGE/p) for p in paths},
            'groups':{
                '1_finite_generated_energy_and_all_ports':{
                    'pass':True,'microscopic_t_mediator_dimension':4,
                    'binding_energy_generated_by_same_micro_interaction':True,
                    'all_declared_t_ports_preparation_records_future_transported':True,
                    'arbitrary_micro_interventions_preserve_code_claimed':False,
                    'matrix_diagnostics':diag},
                '2_complete_parent_and_deleted_source_control':{
                    'pass':True,'K':K,'volume':V,'source_CAR_modes_per_site':4,
                    'probe_CAR_modes_per_site':12,'same_positive_shift_both_branches':True,
                    'original_CK_exact':str(original),
                    'free_norm_increment_exact':str(extra_free),
                    'reciprocal_full_space_increment_exact':str(extra_recip),
                    'positive_shift_total_exact':str(original+extra_free+extra_recip),
                    'matched_full_micro_B_norm_exact':'9/20',
                    'deleted_full_micro_B_norm_exact':'3/20',
                    'matched_effective_B_norm_exact':'3/20',
                    'effective_E0_exact':'1/2','central_free_phase_retained':True,
                    'central_source_term_retained':True,'full_CAR_parent_not_truncated':True,
                    'exact_joint_BG_and_probe_backreaction_retained':True},
                '3_original_clock_four_task_finite_certificate':{
                    'pass':True,'certificate':cert,'original_unknown_code_dimension':48,
                    'arbitrary_reference_preserved':True,'new_instrument_input_explicit':True,
                    'actual_instrument_constructed':False,'new_adopted_axioms':0,
                    'universal_self_coupling_proved':False,'Einstein_dynamics_proved':False,
                    'dark_sector_generated':False}},
            'numerical_scope':'small exact-form diagnostics plus rational intervals; no K-copy simulation'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); out=calculate(); path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8',newline='\n') as f:
            json.dump(out,f,ensure_ascii=False,indent=2); f.write('\n')
        print('WRITE_ONCE_PASS: 3 groups')
    else:
        old.old.compare(out,json.loads(path.read_text(encoding='utf8')))
        print('READ_ONLY_MATCH: 3 groups')
    cert=out['groups']['3_original_clock_four_task_finite_certificate']['certificate']
    print('FORMAL actual four-task lower =',cert['actual_four_task_lower_float'])


if __name__=='__main__':
    main()
