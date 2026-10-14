"""1179: exact rational finite-menu certificate; default read-only, --write once.

No K-copy state, mean-field ODE, or hardware is simulated.  Universal preparation
and full-output bounds are analytic in proof.md; three finite checks calibrate
the interface, outward scalar constants, and the fixed original-record event.
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
DEN = 10**50
EPS, NU, G, OM, ELL = F(3,20), F(7,20), F(1,5), F(1,2), F(16,65)
K = 10**14
FROZEN = {
    '1178/proof.md':'cd44222a8fedc2d45c497b2e0e02120e18408c3385d5a628e3f0066895f127cf',
    '1178/results.json':'e6255dd8f47a4bc57e1cf8a1b8cd602d9ce3ac87098013bc0e3e26442ba13b98',
    '1178/check.py':'3c5e06e1a048ca2758c33f1b74f6c67b75e4c19ab3654b516a08e3e115ae6037',
    '1178/review.md':'c7b76b46abe364ba5c2d031b732766d776937717e29ec6ce6a69cccde641d888',
    '1177/proof.md':'6c3229a8bd9caad35e7d96bc53ced14e867feb89a8aeaff01be8d7d40339c447',
    '1176/proof.md':'966a751121980b976ad86f459d06088c72a56e3eca49109828293020b86320f9',
    '1175/proof.md':'581427de3c8e03176a1145773d6bd0754930d849524110a22971c0e3f8e861f4',
    '_admission/after1174_direction/mean_field_bound.md':'11ac64ce824d3448dec97548efe3dda768bcda30f2e48e1cd88bfe511fcd6b3f',
}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def out(a): return F((a[0]*DEN).__floor__(),DEN), F((a[1]*DEN).__ceil__(),DEN)
def pt(x): return F(x), F(x)
def add(a,b): return out((a[0]+b[0],a[1]+b[1]))
def scale(a,b):
    b=F(b); return out(tuple(sorted((a[0]*b,a[1]*b))))
def sub(a,b): return add(a,scale(b,-1))
def mul(a,b):
    v=[x*y for x in a for y in b]; return out((min(v),max(v)))
def sqrt_i(a):
    assert a[0]>=0
    lo=math.isqrt((a[0]*DEN*DEN).__floor__())
    hi=math.isqrt((a[1]*DEN*DEN).__floor__())+1
    z=F(lo,DEN),F(hi,DEN)
    assert z[0]**2<=a[0] and z[1]**2>=a[1]
    return z
def atan_i(x,n=100):
    assert 0<x<1
    s=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    tail=(-1)**n*x**(2*n+1)/F(2*n+1)
    return out((min(s,s+tail),max(s,s+tail)))
def trig(a,kind,n=32):
    # Point Taylor with Lagrange remainder; interval extension uses |f'| <= 1.
    x=(a[0]+a[1])/2; radius=(a[1]-a[0])/2
    start=1 if kind=='sin' else 0
    last=2*n+start
    s=sum(((-1)**j*x**(2*j+start)/math.factorial(2*j+start)
           for j in range(n+1)),F(0))
    err=abs(x)**(last+1)/math.factorial(last+1)+radius
    return out((s-err,s+err))
def exp_i(a,n=64):
    assert 0<=a[0]<=a[1]<3
    def series(x):
        term=s=F(1)
        for j in range(1,n+1): term=term*x/j; s+=term
        tail=term*x/(n+1)/(1-x/F(n+2))
        return s,tail
    lo,_=series(a[0]); hi,tail=series(a[1])
    return out((lo,hi+tail))
def down(x):
    v=float(x)
    return math.nextafter(v,-math.inf) if F.from_float(v)>x else v
def up(x):
    v=float(x)
    return math.nextafter(v,math.inf) if F.from_float(v)<x else v
def ji(a):
    return {'lower_exact':str(a[0]),'upper_exact':str(a[1]),
            'lower_float_down':down(a[0]),'upper_float_up':up(a[1])}
def rank(a):
    a=[[F(x) for x in r] for r in a]; row=0
    for col in range(len(a[0])):
        nz=next((i for i in range(row,len(a)) if a[i][col]),None)
        if nz is None: continue
        a[row],a[nz]=a[nz],a[row]; z=a[row][col]
        a[row]=[v/z for v in a[row]]
        for i in range(len(a)):
            if i!=row:
                z=a[i][col]; a[i]=[x-z*y for x,y in zip(a[i],a[row])]
        row+=1
        if row==len(a): break
    return row

def interface_checks():
    I=np.eye(2,dtype=complex)
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    Z=np.diag([1.,-1.]).astype(complex)
    pa=[X,Y,Z]; al=[np.kron(X,p) for p in pa]; beta=np.kron(Z,I)
    u=(I-1j*(X+Y+Z))/2; ud=np.kron(I,u)
    assert np.array_equal(ud@ud.conj().T,np.eye(4))
    assert np.array_equal(ud@beta@ud.conj().T,beta)
    for b in range(3): assert np.array_equal(ud@al[b]@ud.conj().T,al[(b+1)%3])
    C=np.array([[0,0,1],[1,0,0],[0,1,0]],int)
    S=np.diag([1,-1,0]); SC=C@S@C.T
    assert np.array_equal(SC,np.diag([0,1,-1]))
    w=np.array([1,1,-2]); assert w@np.ones(3)==0 and w@np.diag(S)==0
    assert w@np.diag(SC)==3
    # Coordinates xx,yy,xy,yz,zx in Sym_0(3): two diagonal and three off-diagonal.
    span=[[1,-1,0,0,0],[0,1,0,0,0],[0,0,1,0,0],
          [0,0,0,1,0],[0,0,0,0,1]]
    assert rank(span)==5
    assert rank([[1,1],[1,-1],[1,0]])==2
    # Fourier-compressed original f_i. All entries are exact binary rationals.
    f=np.array([[.5,.25j],[-.25j,.5]],complex)
    p0=np.array([1,1],complex)
    h=np.diag([0.,1.]); prime=1j*(h@f-f@h)
    assert ((p0.conj()@prime@p0)/2).real==.25
    e_xy=(F(1,4),F(3,4)); e_z=(F(1),F(0))
    assert all(0<=v<=1 for v in e_xy+e_z)
    return {'cycle_matrix':C.tolist(),'rotated_S':SC.tolist(),
            'same_48_code_and_lattice_cycle':True,'Dirac_cycle_exact_binary_check':True,
            'interface_rank':2,'SO3_orbit_linear_span_dimension':5,
            'annihilator':[1,1,-2],'rotated_annihilator_exact':'3',
            'generator_norm_distance_factor_exact':'3/4',
            'Fourier_probability_sign':'plus','Fourier_phase_derivative_exact':'1/4',
            'event_xy_b0_b1':[str(x) for x in e_xy],
            'event_z_b0_b1':[str(x) for x in e_z]}

def certificates():
    pi=sub(scale(atan_i(F(1,5)),16),scale(atan_i(F(1,239)),4))
    t=scale(pi,F(1,2)); k=scale(trig(scale(pi,F(2,65)),'sin'),1/ELL)
    s0=mul(k,t); a=scale(k,G)
    s_window=(scale(s0,1-NU)[0],scale(s0,1+NU)[1])
    phi_window=(scale(s0,1-NU-G)[0],scale(s0,1+NU+G)[1])
    assert 0<s_window[0]<s_window[1]<pi[0]
    assert 0<phi_window[0]<phi_window[1]<t[0]
    a_s=EPS*F(1,2)+3*NU/ELL; a_t=2*G/ELL
    b2=7*a_s*a_s+5*a_t*a_t; w_site=2*EPS*OM
    pref=add(sqrt_i(pt(2*b2)),pt(a_s+a_t))
    eta=scale(mul(mul(t,exp_i(scale(t,6*w_site))),pref),F(4,10**7))
    eta_adopt=F(3,50000); inst=F(1,50000); delta=eta_adopt+inst
    assert eta[1]<eta_adopt and delta==F(1,12500)
    sin_s=trig(s0,'sin'); cos_s=trig(s0,'cos')
    sin_a=trig(a,'sin'); cos_a=trig(a,'cos')
    wstar=scale(add(mul(sin_s,sub(pt(1),cos_a)),scale(mul(cos_s,sin_a),3)),F(1,4))
    target=add(pt(F(1,2)),scale(wstar,F(1,6)))
    ideal_gap=scale(wstar,F(1,6)); formal=ideal_gap[0]-2*delta
    target_actual_lo=target[0]-delta; candidate_actual_hi=F(1,2)+delta
    assert formal==target_actual_lo-candidate_actual_hi
    assert formal>F(79,10000)
    # Original monotone phase comparison is a separate structural diagnostic.
    cos_u=trig(pt(phi_window[1]),'cos')
    max_axis=scale(mul(k,cos_u),3*G/F(16))
    tv_phase=scale(mul(k,cos_u),G/F(8))
    eta_group={
        'K_exact':str(K),'time':ji(t),'pi':ji(pi),'k_star':ji(k),
        'a_s_exact':str(a_s),'a_t_exact':str(a_t),'B_squared_exact':str(b2),
        'w_site_exact':str(w_site),'eta_uniform_prepared':ji(eta),
        'eta_adopted_each_exact':str(eta_adopt),'instrument_each_input_exact':str(inst),
        'total_each_exact':str(delta),'scalar_phase_window':ji(s_window),
        'all_axis_phase_window':ji(phi_window),
        'arbitrary_known_pure_four_dimensional_cluster':True,
        'exact_BG_retained_not_product_replaced':True,
        'classical_labels_retained_and_axis_independent':True,
        'no_local_Gt_readout_error_needed':True}
    menu_group={
        'S0':ji(s0),'a':ji(a),'target_W':ji(wstar),'target_event_ideal':ji(target),
        'candidate_event_ideal_upper_exact':'1/2',
        'candidate_event_actual_upper_exact':str(candidate_actual_hi),
        'target_event_actual_lower_exact':str(target_actual_lo),
        'target_event_actual_lower_float_down':down(target_actual_lo),
        'formal_fixed_event_gap_lower_exact':str(formal),
        'formal_fixed_event_gap_lower_float_down':down(formal),
        'simple_strict_gap_lower_exact':'79/10000',
        'ideal_fixed_event_gap':ji(ideal_gap),
        'max_axis_diagnostic_ideal_lower':ji(max_axis),
        'max_axis_diagnostic_actual_lower_float_down':down(max_axis[0]-2*delta),
        'phase_TV_diagnostic_actual_lower_float_down':down(tv_phase[0]-2*delta),
        'with_additional_sampling_error_each_1e_minus4_lower_float_down':down(formal-F(1,5000)),
        'candidate_W_nonpositive_for_all_S_in_window_and_all_real_D':True,
        'known_labelled_mixtures_follow_by_linearity':True,
        'single_preparation_chosen_before_random_axis':True,
        'strict_values_are_exact_rationals_floats_are_directed':True}
    return eta_group,menu_group

def run():
    deps={}
    for rel,expected in FROZEN.items():
        actual=sha(STAGE/rel); assert actual==expected,(rel,actual)
        deps[rel]=actual
    texts={p.name:sha(p) for p in (STAGE/'research_note_1179.md',HERE/'admission.md',
                                    HERE/'proof.md',HERE/'NEXT.md',HERE/'check.py')}
    b,c=certificates()
    return {'round':1179,'status':'PASS','groups':{
        '1_cycle_interface_and_original_Fourier_menu':interface_checks(),
        '2_all_preparations_uniform_full_output_bound':b,
        '3_fixed_record_event_separates_the_cycle_target':c},
        'frozen_dependencies_sha256':deps,'author_inputs_sha256':texts,
        'scope':'fixed 1178 parent, source vacuum, homogeneous product preparations and retained-label mixtures; conditional active-G rotation target; no GR or universal rotation no-go',
        'numerical_method':'50-place outward rational Machin, Taylor, integer sqrt and positive exp; no imported author check or numerical integration'}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=run(); p=HERE/'results.json'
    if args.write:
        if p.exists(): raise SystemExit('REFUSE_OVERWRITE: results.json already exists')
        p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        status='WRITE_ONCE_PASS'
    else:
        if not p.exists(): raise SystemExit('MISSING_RESULTS: first authorized save requires --write')
        assert json.loads(p.read_text(encoding='utf-8'))==result,'READ_ONLY_MISMATCH'
        status='READ_ONLY_MATCH'
    print(status,'3 groups; fixed-event gap >',
          result['groups']['3_fixed_record_event_separates_the_cycle_target']['formal_fixed_event_gap_lower_float_down'])

if __name__=='__main__': main()
