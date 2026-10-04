"""748 entry: original scalar-edge sixth coefficient, numerical not certified.

Duhamel reduction is stated in the companion working proof. No full graph
bosonic time simulation is performed. The edge is not weakened.
"""
import json,sys,math
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
TARGET=HERE/'scalar_edge_response_results.json'
G=109/625
def p4(a,b):
    return G*(-76/3-8*b-17*a/3+31*a*(a+b)/9+2*a*(a+b)**2/3)
def points(n):
    # Generalized Laguerre alpha=1, from its Jacobi matrix; no SciPy needed.
    j=np.arange(n,dtype=float);d=2*j+2;k=np.arange(1,n,dtype=float)
    J=np.diag(d)+np.diag(np.sqrt(k*(k+1)),1)+np.diag(np.sqrt(k*(k+1)),-1)
    a,U=np.linalg.eigh(J);wa=U[0,:]**2
    s,ws=np.polynomial.hermite.hermgauss(n)
    aa,ss=np.meshgrid(a,s,indexing='ij')
    u=1+(aa+ss*ss)/6
    wt=.5*wa[:,None]*ws[None,:]/np.sqrt(u)
    return aa.ravel(),ss.ravel(),u.ravel(),wt.ravel()
def stable_ratio(q):
    # acosh(q)/sqrt(q²-1), with the smooth value 1 at q=1.
    h=np.maximum(q-1,0)
    den=np.sqrt(h*(h+2))
    return np.divide(np.arccosh(1+h),den,out=np.ones_like(h),where=h>1e-13)
def edge_potential(x,y,k=1.):
    u=1+x@x/6;v=1+y@y/6;q=math.sqrt(u*v)-x@y/6
    return 3*k*math.acosh(max(q,1.))**2
def direct_force(x,y):
    u=1+x@x/6;v=1+y@y/6;q=math.sqrt(u*v)-x@y/6
    ah=y[:4]@y[:4];cros=x[:4]@y[:4]
    return -math.cos(ah/v)/(2*v)*float(stable_ratio(np.array(q)))*(math.sqrt(u/v)*ah-cros)
def force_check():
    rng=np.random.default_rng(748);errors=[]
    for _ in range(12):
        x=rng.normal(size=5)*.6;y=rng.normal(size=5)*.6
        step=2e-5;gradient=np.empty(5)
        for k in range(5):
            h=np.zeros(5);h[k]=step
            gradient[k]=(edge_potential(x,y+h)-edge_potential(x,y-h))/(2*step)
        v=1+y@y/6;ah=y[:4]@y[:4]
        native_gradient=np.r_[y[:4]*math.cos(ah/v)/(2*v),0.]
        numerical=-native_gradient@gradient
        err=abs(numerical-direct_force(x,y));assert err<1e-9;errors.append(err)
    # Duhamel beta-integral: integral_0^t s^4*(t-s) ds / 4! = t^6/6!.
    factor=(Q(1,5)-Q(1,6))/math.factorial(4)
    assert factor==Q(1,math.factorial(6))
    return dict(samples=12,max_force_identity_error=max(errors),exact_time_order_factor=str(factor))
def integrate(n,m):
    a,s,u,w=points(n);Z=float(sum(w));r=np.sqrt(a);z=np.sqrt(u)
    source=w*p4(a,s*s)
    # Sum over normalized relative S3 directions and independent radial ready states.
    angles=np.arange(1,m+1)*np.pi/(m+1)
    cs=np.cos(angles);wc=2/(m+1)*np.sin(angles)**2
    integral=0.
    for first in range(0,len(a),128):
        ab=a[first:first+128][None,:];sb=s[first:first+128][None,:]
        ub=u[first:first+128][None,:];zb=z[first:first+128][None,:]
        rb=r[first:first+128][None,:];wb=w[first:first+128][None,:]
        pa=source[:,None]*wb
        pref=-np.cos(ab/ub)/(2*ub)
        for c,ww in zip(cs,wc):
            cross=r[:,None]*rb*c
            q=z[:,None]*zb-(cross+s[:,None]*sb)/6
            force=pref*stable_ratio(q)*(z[:,None]*ab/zb-cross)
            integral+=ww*float(np.sum(pa*force))
    scalar=integral/(Z*Z)
    cert=json.loads((ROOT/'round747_drafts/remote_hop_sign_certificate_results.json').read_text('utf8'))
    lo,hi=map(float,cert['unit_hop_squared_sixth_unnormalized_interval'])
    hop=(177/1000)*(.5*(lo+hi))/Z
    return dict(radial_order=n,relative_angle_order=m,gaussian_Z=Z,
        scalar_sixth_for_k_wA_wB_1=scalar,
        original_tau_hop_sixth_for_hbar_wB_1=hop,
        total_sixth_for_common_epsilon_psi_hbar_1=scalar+hop,
        estimated_only_no_certified_error_bound=True)
def run():
    checks=force_check()
    rows=[integrate(n,m) for n,m in ((10,12),(16,18),(24,24),(32,32))]
    return dict(entry_round=748,latest_completed_round=747,formal_test_count_unchanged=3432,
        checks=checks,quadrature_rows=rows,
        assumptions=['Original two-node graph with epsilon=psi_A=psi_B=hbar=1, hence w_A=w_B=k_AB=1.',
        'Original tau and Yukawa coefficients from747; radial Gaussian comparison, later compact core required.',
        'All original potentials and fields retained; this computes the first scalar-edge coefficient, not a changed scalar-only dynamics.'],
        proposed_exact_formula='a6 = <P4(x_A) * [-grad_K f_B dot grad_B V_AB]> / (w_A² w_B), with the Gaussian pair measure; k included in V.',
        total_signal_not_yet_certified=True,full_field_time_evolution_not_simulated=True,
        next='Audit the Duhamel word exclusions and obtain a rigorous bound on the original total coefficient before claiming remote distinguishability.')
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
