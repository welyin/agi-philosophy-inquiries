"""808: cutoff-affinity and state-transport calibration in original 64 coefficients.

The constant auxiliary Fourier block is a coefficient diagnostic, not the original
coupled boson solution or future PDE. Analytic existence is proved in the note.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as old
TARGET=HERE/'cutoff_response_contract_results.json'

def step(x):
    x=np.asarray(x);out=np.zeros_like(x,dtype=float);out[x>=1]=1
    m=(x>0)&(x<1);a=np.exp(-1/x[m]);b=np.exp(-1/(1-x[m]));out[m]=a/(a+b)
    return out
def eta(t):return 1-step((np.abs(t)-.25)/.75)
def bump(t,center,width):
    x=(np.asarray(t)-center)/width;out=np.zeros_like(x,dtype=float);m=abs(x)<1
    out[m]=np.exp(-1/(1-x[m]**2));return out
def run():
    k=np.array([0.,1.,0.]);H=old.hamiltonian(k);P=old.I-old.occupied(k)
    e,v=np.linalg.eigh(H);K=.5*old.GAMMA[2]
    center=.55;width=.07;idx=[30,31]
    def U(t):return (v*np.exp(-1j*e*t))@v.conj().T
    a=P@U(center).conj().T@K@U(center)@(old.I-P)
    ell=(a-a.conj().T)/1j
    val,vec=np.linalg.eigh(ell[np.ix_(idx,idx)])
    f=np.zeros(64,complex);f[idx]=vec[:,-1];u=P@f;z=(old.I-P)@f
    def density(t):
        phases=np.exp(-1j*np.outer(t,e))
        ut=(phases*(v.conj().T@u))@v.T
        zt=(phases*(v.conj().T@z))@v.T
        return 2*np.imag(np.einsum('ti,ij,tj->t',ut.conj(),K,zt))
    def integrate(order):
        x,w=np.polynomial.legendre.leggauss(order);ans=np.zeros(3);slope=0.
        breaks=sorted([-1.,-.25,.25,center-width,center+width,1.])
        for lo,hi in zip(breaks[:-1],breaks[1:]):
            t=(lo+hi)/2+(hi-lo)*x/2;weights=w*(hi-lo)/2
            j=density(t);q=bump(t,center,width);base=eta(t)
            slope+=float(np.dot(weights,q*j))
            for ix,a in enumerate((-.1,0.,.1)):ans[ix]+=np.dot(weights,(base+a*q)*j)
        return ans,slope
    ans,slope=integrate(96);fine,slope_fine=integrate(192)
    residual=float(max(np.max(abs(ans-fine)),abs(slope-slope_fine)))
    affinity=float(max(abs(ans[0]-ans[1]+.1*slope),abs(ans[2]-ans[1]-.1*slope)))
    assert residual<2e-11 and affinity<1e-13 and slope>1e-4
    ts=np.linspace(-1,1,2001);base=eta(ts);q=bump(ts,center,width)
    assert np.all(base-.1*q>=-1e-14) and np.all(base+.1*q<=1+1e-14)
    inside=np.abs(ts)<.25
    assert np.max(abs(q[inside]))==0 and np.max(abs(base[inside]-1))==0
    # Same-process comparison requires opposite first-order transport of state.
    # This is a separate algebraic diagnostic; do not identify P with this density.
    rho=np.eye(64)/128+np.outer(f,f.conj())/2
    A=np.outer(u,u.conj());B=K
    dA=1j*(B@A-A@B);drho=1j*(B@rho-rho@B)
    observable_only=float(np.trace(rho@dA).real)
    state_only=float(np.trace(drho@A).real)
    # For alpha(A)=exp(-iaB) A exp(iaB), d(alpha A)=-dA,
    # while rho_a=exp(-iaB) rho exp(iaB) has derivative -drho.
    # Their joint pairing is invariant: trace(drho A)+trace(rho dA)=0.
    cancellation=abs(observable_only+state_only)
    assert cancellation<1e-13 and abs(observable_only)>1e-5
    return dict(round=808,all_checks_passed=True,original_Nambu_dimension=64,
        diagnostic_uses_auxiliary_constant_Fourier_block=True,
        original_boson_solution_or_future_PDE_computed=False,
        original_continuous_response_values_computed=False,
        diagnostic_cutoff_parameters=[-.1,0.,.1],diagnostic_responses=ans.tolist(),
        diagnostic_variation_slope=slope,quadrature_difference=residual,
        affine_response_residual=affinity,old_plateau_unchanged=True,
        observable_only_variation=observable_only,state_only_variation=state_only,
        common_transport_cancellation=cancellation,
        claimed_continuous_result='Existence in the permitted compact-cutoff process family, not nonzero response for the previously fixed process.',
        new_numbered_test_groups=1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
