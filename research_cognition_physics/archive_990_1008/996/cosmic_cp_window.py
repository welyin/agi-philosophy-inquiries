"""996: finite expansion response with a shared Weinberg coupling.
Thermal source normalization is left symbolic. This is not a baryon abundance.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'cosmic_cp_window_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def exact(x):return dict(exact=str(x),value=float(x))
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)

def run():
    old=read(STAGE/'995/internal_cp_bridge_results.json')
    state=old['parity_even_initial_state']
    m=F(state['omega']['exact']);vol=F(state['volume']['exact']);scale=F(state['scale']['exact'])
    q2=F(state['canonical_Q2']['exact']);p2=F(state['canonical_P2']['exact'])
    # Phase averaging is an additional preparation. Keep only the stated quadratic order.
    E2=(p2+m*m*q2)/2
    n=E2/m-F(1,2)
    fi=n/(m*vol*scale*scale)
    assert n==F(9,16) and fi==F(9,800000)
    # Move the initial mode-vacuum reference into the constant coefficient too.
    # This calibrates reference covariance, not full thermal/curvature matching.
    fv=1/(2*m*vol*scale*scale)
    c0_old=np.diag([.1,.2,.3]).astype(complex)
    c2=np.array([[.2j,.05,0],[.05,.1+.1j,.04j],[0,.04j,-.05+.15j]])
    c0=c0_old+float(fv)*c2
    I=F(old['invariant']['exact'])
    r2=F(233,2500)
    A=F(7,50)+fv/F(100)+fv*fv*r2
    B=fi/F(100)+2*fi*fv*r2;D=fi*fi*r2;gi=A+B+D
    assert abs(np.trace(c0.conj().T@c0).real-float(A))<1e-15
    assert abs(np.trace(c0.conj().T@c2).real-float(F(1,200)+fv*r2))<1e-15
    assert abs(np.trace(c2.conj().T@c2).real-float(r2))<1e-15
    assert np.max(abs(c0+float(fi)*c2-(c0_old+float(fi+fv)*c2)))<1e-15
    wi=F(1,10) # A declared dimensionless washout ratio, not an inferred physical rate.
    af=F(2)
    def norm(a):return float(A)+float(B)*a**-3+float(D)*a**-6
    def optical(a,b):
        return float(wi/gi)*(float(A)*(a**-1-b**-1)
            +float(B)/4*(a**-4-b**-4)+float(D)/7*(a**-7-b**-7))
    W=wi/gi*(A*(1-1/af)+B*(1-af**-4)/4+D*(1-af**-7)/7)
    unwashed=(1-af**-5)/5
    samples=[]
    for a in (1.,1.25,1.5,2.):
        f=float(fi)/a**3;c=c0+f*c2
        ng=float(np.trace(c.conj().T@c).real)
        assert abs(ng-norm(a))<1e-15
        samples.append(dict(a=a,f_excitation=f,coupling_norm_squared=ng,
            washout_over_H=float(wi)/a*ng/float(gi),
            optical_depth_to_endpoint=optical(a,float(af))))
    def response(order):
        x,w=np.polynomial.legendre.leggauss(order)
        a=1+(float(af)-1)*(x+1)/2
        return float((float(af)-1)/2*np.dot(w,a**-6*np.exp(-optical(a,float(af)))))
    y64=response(64);y128=response(128)
    lower=math.exp(-float(W))*float(unwashed)
    assert 0<lower<y128<float(unwashed)
    assert abs(y64-y128)<1e-13
    # Independent finite-time ODE for the SAME Green response, using RK4.
    def slope(a,y):return a**-6-float(wi)*a**-2*norm(a)/float(gi)*y
    def integrate(steps):
        h=(float(af)-1)/steps;y=0.
        for j in range(steps):
            a=1+j*h;k1=slope(a,y);k2=slope(a+h/2,y+h*k1/2)
            k3=slope(a+h/2,y+h*k2/2);k4=slope(a+h,y+h*k3)
            y+=h*(k1+2*k2+2*k3+k4)/6
        return y
    yode=integrate(1024)
    assert abs(yode-y128)<1e-12
    # |Y| / (|alpha/beta| * H_i/T_i), alpha is NOT evaluated.
    suppression=3*float(abs(I)*fi/gi*wi)*y128
    assert 0<suppression<5e-7
    # No external switch is required for a decreasing rate in this finite window.
    assert samples[-1]['washout_over_H']<samples[0]['washout_over_H']
    # A CP factor alone does not force a nonzero first thermal memory moment.
    # Integral tau*[exp(-tau)-4 exp(-2tau)] d tau = 0 exactly.
    counter_kernel_first_moment=F(1)-4*F(1,4)
    assert counter_kernel_first_moment==0
    inputs=[Path(__file__).resolve(),HERE/'drafts/STATUS.md',HERE/'drafts/selection.md',
        STAGE/'995/internal_cp_bridge.py',STAGE/'995/internal_cp_bridge_results.json',
        STAGE/'995/internal_weinberg_adoption.md',STAGE/'993/common_candidate_v1.md',
        STAGE/'research_note_994.md',STAGE/'research_note_964.md',STAGE/'research_note_992.md',
        STAGE.parent/'archive_301_341/research_note_332.md',
        STAGE.parent/'archive_301_341/research_note_333.md',
        STAGE.parent/'archive_301_341/332/occupied_scalar_gravity_audit_results.json']
    return dict(round=996,all_scientific_checks_passed=True,
        kind='conditional_slow_expansion_source_and_survival_screen',
        phase_averaged_quadratic_preparation=dict(occupation=exact(n),f_excitation_initial=exact(fi),
            initial_vacuum_reference_shift=exact(fv),constant_coefficient_shifted_together=True,
            excitation_energy=exact(m*n),quartic_and_thermal_corrections_certified=False),
        shared_coupling=dict(norm_constant=exact(A),norm_a_minus3=exact(B),norm_a_minus6=exact(D),
            norm_initial=exact(gi),CP_invariant=exact(I)),
        finite_window=dict(a_initial=1,a_final=int(af),washout_over_H_initial=exact(wi),
            optical_depth=exact(W),survival_lower=math.exp(-float(W)),
            unwashed_shape=exact(unwashed),response_lower_bound=lower,
            normalized_response=y128,quadrature_comparison=abs(y64-y128),
            independent_ODE_residual=abs(yode-y128),
            yield_per_abs_alpha_over_beta_times_H_over_T=suppression,samples=samples),
        counter_kernel_first_moment=exact(counter_kernel_first_moment),
        thermal_source_coefficient_evaluated=False,finite_physical_net_charge_certified=False,
        actual_thermal_hierarchy_verified=False,cosmological_abundance_fitted=False,
        full_goal_completed=False,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
