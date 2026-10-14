"""1177: default read-only recomputation; --write creates results.json once.

Three groups only. Exact rational intervals certify the probability bound.
Small matrices diagnose source vacuum, local norms and complete records;
no 65^3 lattice evolution or K-copy parent is numerically simulated.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
K = 10**14
ELL, EPS, NU, G = F(16,65), F(3,20), F(7,20), F(1,5)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(a,b): return (a[0]+b[0], a[1]+b[1])
def mul(a,b):
    v = [x*y for x in a for y in b]
    return min(v),max(v)
def scale(a,c): return mul(a,(c,c))
def power(a,n):
    out = (F(1),F(1))
    for _ in range(n): out=mul(out,a)
    return out


def atan_interval(x,terms=40):
    # Alternating positive decreasing series, exact rational endpoints.
    s = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    rem = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(s,s+rem),max(s,s+rem)


def trig_interval(x,kind,terms=20):
    out = (F(0),F(0))
    offset = 1 if kind=='sin' else 0
    for j in range(terms):
        n=2*j+offset
        out=add(out,scale(power(x,n),F((-1)**j,math.factorial(n))))
    n=2*terms+offset
    # Taylor remainder: omitted even/odd derivative is bounded by 1.
    rem=max(abs(x[0]),abs(x[1]))**n/F(math.factorial(n))
    return out[0]-rem,out[1]+rem


def interval_json(x):
    # Store short exact outward decimal-rational bounds, not thousands of digits.
    den=10**40
    lo=F((x[0]*den).__floor__(),den)
    hi=F((x[1]*den).__ceil__(),den)
    assert lo<=x[0]<=x[1]<=hi
    return {'lower_exact':str(lo),'upper_exact':str(hi),
            'outward_decimal_places':40,
            'lower_float_diagnostic':float(x[0]),'upper_float_diagnostic':float(x[1])}


def fermion_annihilator(n,j):
    a=np.zeros((2**n,2**n),complex)
    for state in range(2**n):
        if (state>>j)&1:
            a[state^(1<<j),state]=(-1)**((state&((1<<j)-1)).bit_count())
    return a


def local_material_check():
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    beta=np.kron(sz,np.eye(2))
    alpha=[np.kron(sx,p) for p in (sx,sy,sz)]
    # Central site plus +/- x,y,z: all 12 internal modes are retained.
    bm=np.kron(beta,np.diag([0.,.3,.5]))
    aa=[np.kron(a,np.eye(3)) for a in alpha]
    fs=np.zeros((84,84),complex); ft=np.zeros_like(fs)
    fs[:12,:12]=float(EPS)*bm
    for axis in range(3):
        for sign in (1,-1):
            site=1+2*axis+(sign<0)
            block=-1j*sign*aa[axis]/(4*float(ELL))
            sl=slice(12*site,12*(site+1))
            fs[:12,sl]=float(NU)*block
            fs[sl,:12]=float(NU)*block.conj().T
            if axis<2:
                c=float(G)*(1 if axis==0 else -1)
                ft[:12,sl]=c*block; ft[sl,:12]=c*block.conj().T
    # A representative full normal-ordered CAR operator; vacuum argument is general.
    car=[fermion_annihilator(4,j) for j in range(4)]
    q=np.diag([2.,2.,-2.,-2.]).astype(complex)
    q[0,3]=.3j; q[3,0]=-.3j
    fq=sum((q[j,k]*car[j].conj().T@car[k]
            for j in range(4) for k in range(4)),np.zeros((16,16),complex))
    vac=np.eye(16)[:,0]
    vacuum_res=float(np.linalg.norm(fq@vac))
    s_norm=float(np.max(np.abs(np.linalg.eigvalsh(fs))))
    t_norm=float(np.max(np.abs(np.linalg.eigvalsh(ft))))
    assert vacuum_res==0
    assert s_norm<float(EPS*F(1,2)+3*NU/ELL)
    assert t_norm<float(2*G/ELL)
    assert np.linalg.norm(ft[60:])==0 and np.linalg.norm(ft[:,60:])==0
    return {'normal_ordered_four_mode_vacuum_residual':vacuum_res,
            'seven_site_full_internal_dimension':84,
            's_local_norm_diagnostic':s_norm,'t_local_norm_diagnostic':t_norm,
            's_support_multiplicity':7,'t_support_multiplicity':5}


def record_check(t,tau):
    # Position(2), original mass clock(2), Gt_0,Gt_1,pointer_0,pointer_1.
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    plus=np.ones(2,complex)/math.sqrt(2)
    clock=np.diag(np.exp(-1j*np.array([.3,.5])*tau))@plus
    z0=np.array([1.,0.],complex)
    gt=(math.cos(t/2)*np.eye(2)+1j*math.sin(t/2)*sx)@z0
    anc=np.kron(np.kron(gt,gt),np.kron(z0,z0))
    prep=np.kron(np.eye(4),anc[:,None])
    perm=np.zeros((64,64),complex)
    for j in range(64):
        bits=[(j>>(5-k))&1 for k in range(6)]
        bits[4]^=bits[2];bits[5]^=bits[3]
        target=sum(v<<(5-k) for k,v in enumerate(bits))
        perm[target,j]=1
    out=perm@prep
    gram=float(np.linalg.norm(out.conj().T@out-np.eye(4),ord=2))
    psi=out@np.kron(plus,clock)
    clock_eff=np.kron(np.kron(np.eye(2),(np.eye(2)-sy)/2),np.eye(16))
    pclock=float(np.real(np.vdot(psi,clock_eff@psi)))
    # Read the pointer at the observed position, retaining every other record.
    selected=np.zeros(64)
    for j in range(64):
        bits=[(j>>(5-k))&1 for k in range(6)]
        selected[j]=1-bits[4+bits[0]]
    pg=float(np.sum(np.abs(psi)**2*selected))
    assert gram<1e-13
    assert abs(pclock-(1+math.sin(.2*tau))/2)<1e-13
    assert abs(pg-(1+math.cos(t))/2)<1e-13
    return {'full_two_site_pointer_isometry_residual':gram,
            'p_original_minus_y':pclock,'p_colocated_gt_record':pg,
            'all_pointer_and_position_records_retained':True,
            'dimension':64,'diagnostic_only_not_large_parent_simulation':True}


def calculate():
    paths=['1176/proof.md','1176/results.json','1176/NEXT.md']
    sources={p:sha(STAGE/p) for p in paths}
    # Bind the accepted parent, never silently import a changed physics object.
    assert sources['1176/proof.md']=='966a751121980b976ad86f459d06088c72a56e3eca49109828293020b86320f9'
    assert sources['1176/results.json']=='d9b815dec37cd1597202b61a683d768e4d5eb692170992fcaedfb0452f877873'
    pi=add(scale(atan_interval(F(1,5)),F(16)),scale(atan_interval(F(1,239)),F(-4)))
    t=scale(pi,F(1,2));tau=add(t,(EPS,EPS))
    as_=EPS*F(1,2)+3*NU/ELL;at=2*G/ELL
    b2=7*as_**2+5*at**2
    eta_up=F(1,500000)
    assert t[1]**2*b2/F(K)<eta_up**2
    gap=scale(trig_interval((EPS,EPS),'sin'),F(1,2))
    inst=F(1,50000);lip=F(6);e=eta_up+inst
    residual_lower=gap[0]-(1+lip)*e
    calibrated_separation_lower=gap[0]-2*(1+lip)*e
    assert residual_lower>F(7456,100000)
    assert gap[0]>EPS/F(2)-EPS**3/F(12)
    assert calibrated_separation_lower>F(297643,4000000)  # 0.07441075
    # Verify calibration window and Lipschitz certificate without libm.
    assert tau[1]<F(19,10)
    assert trig_interval((F(19,50),F(19,50)),'cos')[0]>F(9,10)
    assert F(5)/(F(9,10))<lip
    pd=scale(add((F(1),F(1)),trig_interval(scale(tau,F(1,5)),'sin')),F(1,2))
    role_pg=scale(add((F(1),F(1)),scale(trig_interval((EPS,EPS),'sin'),F(-1))),F(1,2))
    window_lo=F(1,2)
    window_hi=scale(add((F(1),F(1)),trig_interval((F(19,50),F(19,50)),'sin')),F(1,2))[0]
    assert pd[0]-e>window_lo and pd[1]+e<window_hi
    baseline_pd=scale(add((F(1),F(1)),trig_interval(scale(t,F(1,5)),'sin')),F(1,2))
    assert baseline_pd[0]-e>window_lo and baseline_pd[1]+e<window_hi
    volume=65**3
    bs=2*EPS*2+6*NU/ELL+4*G/ELL
    bp=2*EPS*F(4,5)+18*NU/ELL+12*G/ELL
    ck=1+K*volume*(4+6/ELL+1+bs)+volume*(F(8,5)+18/ELL+bp)
    local=local_material_check()
    records=record_check(math.pi/2,math.pi/2+float(EPS))
    baseline_records=record_check(math.pi/2,math.pi/2)
    return {
        'round':1177,'date':'2026-10-13','status':'PASS_WITH_EXPLICIT_ROLE_CONTRACT',
        'source_sha256':sources,
        'groups':{
            '1_same_parent_exact_vacuum':{
                'pass':True,'K':K,'volume':volume,'source_modes_per_copy':4*volume,
                'relation_qubits':2*volume*K,'new_blank_pointer_qubits':volume,
                'parent_positive_shift_exact':str(ck),
                'source_vacuum_exactly_invariant':True,'W_on_vacuum_exactly_zero':True,
                'exact_background_product_not_interacting_parent_product':True,
                'baseline_is_unit_lapse_not_claimed_flat_geometry':True,
                'local_algebra_diagnostics':local},
            '2_full_input_prepared_isometry':{
                'pass':True,'one_particle_support_bound_all_sites_internal_reference':True,
                'original_code_dimension':48,'all_mass_coherences_preserved':True,
                'a_s_exact':str(as_),'a_t_exact':str(at),'B_squared_exact':str(b2),
                'eta_formula_float_diagnostic':math.pi/2*math.sqrt(float(b2))/math.sqrt(K),
                'eta_adopted_exact':str(eta_up),'eta_adopted':float(eta_up),
                'instrument_error_input_exact':str(inst),'marginal_error_adopted_exact':str(e),
                'joint_effect_full_record_diagnostic':records},
            '3_finite_role_violation_and_control':{
                'pass':True,'pi_interval':interval_json(pi),'tau_interval':interval_json(tau),
                'clock_p_interval':interval_json(pd),'relation_role_p_interval':interval_json(role_pg),
                'reference_relation_p_exact':'1/2','ideal_residual_interval':interval_json(gap),
                'role_curve_Lipschitz_adopted':6,
                'combined_residual_probability_budget_exact':str(7*e),
                'actual_role_residual_lower_exact':str(residual_lower),
                'actual_role_residual_lower':float(residual_lower),
                'separation_after_baseline_error_lower_exact':str(calibrated_separation_lower),
                'separation_after_baseline_error_lower':float(calibrated_separation_lower),
                'with_statistical_each_1e_minus4_lower':float(residual_lower-7*F(1,10000)),
                'baseline_clock_p_interval':interval_json(baseline_pd),
                'baseline_joint_record_diagnostic':baseline_records,
                'candidate_contract_adopted_as_axiom':False,
                'actual_instrument_constructed':False,'Einstein_selected':False,
                'all_extra_modes_retained':True}},
        'numerical_scope':{'certificate':'exact rational Machin and Taylor intervals',
                           'matrix_checks':'small diagnostics only',
                           'large_lattice_or_K_parent_simulated':False}}


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for j,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{j}]')
    elif isinstance(a,float):
        assert abs(a-b)<=2e-12*(1+abs(b)),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=calculate();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS: 3 groups')
    else:
        compare(out,json.loads(path.read_text(encoding='utf-8')))
        print('READ_ONLY_MATCH: 3 groups')
    print('FORMAL certified paired residual lower =',out['groups']['3_finite_role_violation_and_control']['separation_after_baseline_error_lower'])
    print('single-branch residual diagnostic =',out['groups']['3_finite_role_violation_and_control']['actual_role_residual_lower'])


if __name__=='__main__':main()
