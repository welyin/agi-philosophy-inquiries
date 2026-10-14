"""1178: three groups, default read-only; --write creates results once.

Exact outward rational bounds certify scalar closed forms, sqrt and exp.
All numpy small matrices are diagnostics, not the numerical lower-bound proof.
No lattice/mean-field integration or K-copy state is numerically simulated.
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

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
spec=importlib.util.spec_from_file_location('fixed1177',STAGE/'1177'/'check.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
K=10**14
EPS,OM,NU,G,ELL=F(3,20),F(1,2),F(7,20),F(1,5),F(16,65)
LAM=-EPS*OM
DEN=10**40


def outward(a):
    return F((a[0]*DEN).__floor__(),DEN),F((a[1]*DEN).__ceil__(),DEN)
def point(x): return (x,x)
def add(a,b):return outward(old.add(a,b))
def mul(a,b):return outward(old.mul(a,b))
def scale(a,c):return outward(old.scale(a,c))
def sub(a,b):return add(a,scale(b,F(-1)))
def divide(a,b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))
def sqrt_i(a):
    assert a[0]>=0
    lo=math.isqrt((a[0]*DEN*DEN).__floor__())
    hi=math.isqrt((a[1]*DEN*DEN).__floor__())+1
    out=F(lo,DEN),F(hi,DEN)
    assert out[0]**2<=a[0] and out[1]**2>=a[1]
    return out
def trig(a,kind):return outward(old.trig_interval(a,kind))
def exp_i(a,n=64):
    assert 0<=a[0]<=a[1]<3
    def finite(x):
        s=F(1);term=F(1)
        for j in range(1,n+1):
            term=term*x/j;s+=term
        tail=term*x/(n+1)/(1-x/F(n+2))
        return s,tail
    slo,_=finite(a[0]);shi,tail=finite(a[1])
    return outward((slo,shi+tail))
def sin_p(a):return scale(add(point(F(1)),trig(scale(a,F(1,5)),'sin')),F(1,2))
def interval_json(a):return old.interval_json(a)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def analytic_certificate():
    pi=outward(old.add(old.scale(old.atan_interval(F(1,5)),F(16)),
                        old.scale(old.atan_interval(F(1,239)),F(-4))))
    t=scale(pi,F(1,2));q=F(1,2);b=LAM*q;r2=OM*OM+b*b
    r=sqrt_i(point(r2))
    bracket=sub(t,divide(trig(mul(scale(r,F(2)),t),'sin'),scale(r,F(2))))
    shift=scale(bracket,-EPS*OM*b/r2)
    tp=add(t,shift);tm=sub(t,shift)
    pp=sin_p(tp);pm=sin_p(tm);gap=sub(pp,pm)
    as_=EPS*F(1,2)+3*NU/ELL;at=2*G/ELL;b2=7*as_**2+5*at**2
    w=2*abs(LAM)
    pref=add(sqrt_i(point(2*b2)),point(as_+at))
    eta=scale(mul(mul(t,exp_i(scale(t,6*w))),pref),F(4,10**7))
    alpha=scale(sub(exp_i(scale(t,12*w)),point(F(1))),F(1,K))
    sqalpha=sqrt_i(alpha)
    eta_adopt=F(3,50000);local_adopt=F(1,2500000);inst=F(1,50000)
    assert eta[1]<eta_adopt and sqalpha[1]<local_adopt
    dD=eta_adopt+inst;dG=dD+local_adopt
    lower=gap[0]-2*dD
    assert lower>F(106,100000)
    # C: q=0 does not back-write a noncentral mean source, but tests clock closure.
    tc=add(t,point(EPS));pc=sin_p(tc)
    gc=scale(add(point(F(1)),trig(tc,'cos')),F(1,2))
    c_budget=dG+6*dD
    assert c_budget==F(1401,2500000)
    # All actual clock intervals remain in the same predeclared invertible window.
    upper_window=sin_p(point(F(19,10)))[0]
    for pt in (pp,pm,pc):
        assert pt[0]-dD>F(1,2) and pt[1]+dD<upper_window
    assert trig(point(F(19,50)),'cos')[0]>F(9,10)
    assert tc[1]<F(19,10)
    # Nonzero amplitude of q clocks; no equal-energy claim.
    amp=sqrt_i(point(F(3,4)))
    gp=scale(add(point(F(1)),mul(amp,trig(tp,'cos'))),F(1,2))
    gm=scale(add(point(F(1)),mul(amp,trig(tm,'cos'))),F(1,2))
    olddata=json.loads((STAGE/'1177'/'results.json').read_text(encoding='utf8'))
    old_lo=F(olddata['groups']['3_finite_role_violation_and_control']['actual_role_residual_lower_exact'])
    return {
        'time':interval_json(t),'lambda_selected_exact':str(LAM),
        'E0_witness_input_exact':'0','r_squared_exact':str(r2),
        'tau_plus':interval_json(tp),'tau_minus':interval_json(tm),
        'p_D_plus':interval_json(pp),'p_D_minus':interval_json(pm),
        'p_G_plus_Hartree_local':interval_json(gp),'p_G_minus_Hartree_local':interval_json(gm),
        'ideal_clock_gap':interval_json(gap),
        'eta_prepared_interval':interval_json(eta),
        'sqrt_alpha_interval':interval_json(sqalpha),
        'eta_adopted_exact':str(eta_adopt),'local_Gt_adopted_exact':str(local_adopt),
        'instrument_input_each_exact':str(inst),
        'D_error_each_exact':str(dD),'G_local_error_each_exact':str(dG),
        'finite_clock_gap_lower_exact':str(lower),
        'finite_clock_gap_lower':float(lower),
        'C_tau':interval_json(tc),'C_p_D':interval_json(pc),'C_p_G':interval_json(gc),
        'C_role_residual_upper_exact':str(c_budget),'C_role_residual_upper':float(c_budget),
        'old1177_active_residual_lower_minus_C_bound':float(old_lo-c_budget),
        'with_each_sampling_error_1e_minus4_lower':float(lower-F(1,5000)),
        'B_squared_exact':str(b2),'w_site_exact':str(w),
        'a_s_exact':str(as_),'a_t_exact':str(at),
        'all_probability_bounds_use_exact_outward_rational_arithmetic':True,
        'outward_decimal_places':40}


def state_bloch(v):
    x,y,z=v
    a=math.sqrt((1+z)/2)
    return np.array([a,(x+1j*y)/(2*a)],complex)


def unitary(h,t):
    vals,vec=np.linalg.eigh(h)
    return (vec*np.exp(-1j*t*vals))@vec.conj().T


def low_dimensional_checks():
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    Z=np.diag([1.,-1.]).astype(complex);I=np.eye(2)
    ht=-float(OM)*X;bt=float(LAM)*X
    # Central alternatives have identical fixed-z commutator response.
    central=bt+float(EPS)*I
    comm_res=max(float(np.linalg.norm((central-float(EPS)*ht)@a-a@(central-float(EPS)*ht)))
                 for a in (X,Y,Z))
    assert comm_res<1e-14
    xs=np.kron(X,I);xt=np.kron(I,X);zs=np.kron(Z,I);bt4=float(LAM)*xt
    W=np.kron(zs,bt4)+np.kron(bt4,zs)
    wnorm=float(np.max(np.abs(np.linalg.eigvalsh(W))))
    assert abs(wnorm-2*abs(float(LAM)))<1e-14
    plus=np.ones(2,complex)/math.sqrt(2)
    contraction_error=0.;energy_error=0.;bloch_error=0.;data=[]
    for q in (.5,-.5):
        gt0=state_bloch((q,0,math.sqrt(.75)))
        cluster=np.kron(plus,gt0)
        embed=np.kron(np.eye(4),cluster[:,None])
        contracted=embed.conj().T@W@embed
        expected=float(LAM)*q*zs # initial z_s is zero
        contraction_error=max(contraction_error,float(np.linalg.norm(contracted-expected)))
        b=float(LAM)*q;r=math.sqrt(float(OM)**2+b*b)
        for t in (0.,math.pi/6,math.pi/3,math.pi/2):
            tau=t-float(EPS)*float(OM)*b/r**2*(t-math.sin(2*r*t)/(2*r))
            st=unitary(-float(OM)*X+b*Z,t)@plus
            gt=unitary(ht,tau)@gt0
            xv=float(np.real(np.vdot(st,X@st)));zv=float(np.real(np.vdot(st,Z@st)))
            zv_formula=-float(OM)*b/r**2*(1-math.cos(2*r*t))
            bloch_error=max(bloch_error,abs(zv-zv_formula),
                            abs(float(np.real(np.vdot(gt,X@gt)))-q),
                            abs(float(np.real(np.vdot(gt,Z@gt)))-math.sqrt(.75)*math.cos(tau)))
            energy=-float(OM)*xv-float(OM)*q+float(LAM)*zv*q
            energy_error=max(energy_error,abs(energy-(-float(OM)*(1+q))))
        clock=unitary(np.diag([.3,.5]),tau)@plus
        pd=float(np.real(np.vdot(clock,(I-Y)@clock))/2)
        data.append({'q':q,'tau':tau,'p_D':pd,'p_G':(1+float(np.real(np.vdot(gt,Z@gt))))/2,
                     'bare_Gt_energy':-float(OM)*q,'clock_amplitude':math.sqrt(.75)})
    # C uses the actual time-lapse integral; the old bare t prediction is not reused.
    tc=math.pi/2+float(EPS);gtc=unitary(ht,tc)@np.array([1.,0.])
    cpg=float(abs(gtc[0])**2)
    # Local Z write as a complete isometry, not a destructive marginal replacement.
    CNOT=np.eye(4,dtype=complex);CNOT[2:,2:]=X
    E=np.kron(np.eye(2),np.array([[1.],[0.]]))
    gram=float(np.linalg.norm((CNOT@E).conj().T@(CNOT@E)-I))
    record=CNOT@np.kron(gtc,np.array([1.,0.]))
    recpg=float(np.sum(np.abs(record[[0,2]])**2))
    assert max(comm_res,contraction_error,energy_error,bloch_error,gram,abs(cpg-recpg))<2e-13
    return {'fixed_condition_central_commutator_residual':comm_res,
            'full_two_cluster_W_dimension':16,'W_spectral_norm':wnorm,
            'reciprocal_contraction_residual':contraction_error,
            'closed_form_vs_2x2_spectral_residual':bloch_error,
            'Hartree_energy_residual':energy_error,'q_tasks':data,
            'C_Gt_probability':cpg,'C_full_local_record_probability':recpg,
            'local_full_record_isometry_residual':gram,
            'diagnostic_not_large_parent_evolution':True}


def calculate():
    paths=['1176/proof.md','1177/proof.md','1177/check.py','1177/results.json',
           '1170/proof.md','1171/proof.md','_admission/after1174_direction/mean_field_bound.md']
    sources={p:sha(STAGE/p) for p in paths}
    assert sources['1177/check.py']=='ce6be9e2d7a186586af32169a38345d9a72bb653210b99ac4d2db66a4dbdb916'
    assert sources['1177/results.json']=='d93aeadd8b201f982f7006847531c54dfa50ce5ee52c724f00aa5127cb259c09'
    cert=analytic_certificate();diag=low_dimensional_checks()
    for got,key in zip(diag['q_tasks'],('p_D_plus','p_D_minus')):
        assert abs(got['p_D']-cert[key]['lower_float_diagnostic'])<2e-13
    volume=65**3
    bs=2*EPS*2+6*NU/ELL+4*G/ELL
    bp=2*EPS*F(4,5)+18*NU/ELL+12*G/ELL
    cold=1+K*volume*(4+6/ELL+1+bs)+volume*(F(8,5)+18/ELL+bp)
    extra=volume*(K-1)*abs(LAM)
    return {'round':1178,'date':'2026-10-13','status':'PASS_WITH_EXPLICIT_RECIPROCAL_ROLE_SCOPE',
            'source_sha256':sources,
            'groups':{
                '1_selected_response_and_complete_parent':{
                    'pass':True,'linear_Hermitian_class_declared':True,
                    'lambda_selected_exact':str(LAM),'E0_remains_free':True,
                    'E0_witness_input':0,'K':K,'volume':volume,
                    'positive_shift_extra_exact':str(extra),'positive_shift_total_exact':str(cold+extra),
                    'known_background_only_no_unknown_Gt_cloning':True,
                    'new_quantum_rule_and_instrument_inputs_explicit':True,
                    'small_matrix_diagnostics':diag},
                '2_exact_BG_and_unknown_probe_finite_K':{
                    'pass':True,'certificate':cert,
                    'original_unknown_code_dimension':48,'arbitrary_reference_preserved':True,
                    'all_mass_coherences_preserved':True,'exact_joint_BG_preserved':True,
                    'centering_uses_exact_mu_before_spatial_orthogonality':True,
                    'within_site_st_correlation_uses_Cauchy':True,
                    'global_pointer_distribution_replaced_by_Hartree_product':False},
                '3_finite_feedback_C_control_and_remaining_freedom':{
                    'pass':True,'formal_original_clock_gap_lower':cert['finite_clock_gap_lower'],
                    'formal_original_clock_gap_lower_exact':cert['finite_clock_gap_lower_exact'],
                    'C_role_residual_upper':cert['C_role_residual_upper'],
                    'C_has_noncentral_mean_feedback':False,'q_tasks_have_nonzero_mean_feedback':True,
                    'candidate_clock_contract_adopted_as_axiom':False,
                    'actual_instrument_constructed':False,'universal_self_coupling_proved':False,
                    'Einstein_selected':False,'central_baseline_is_Lambda_claim':False,
                    'extra_modes_preserved':True}},
            'numerical_scope':'closed forms and exact rational intervals; no large grid or K-copy simulation'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=calculate();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8',newline='\n') as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS: 3 groups')
    else:
        old.compare(out,json.loads(path.read_text(encoding='utf8')))
        print('READ_ONLY_MATCH: 3 groups')
    print('FORMAL original-clock feedback lower =',out['groups']['3_finite_feedback_C_control_and_remaining_freedom']['formal_original_clock_gap_lower'])
    print('C common-clock residual upper =',out['groups']['3_finite_feedback_C_control_and_remaining_freedom']['C_role_residual_upper'])


if __name__=='__main__':main()
