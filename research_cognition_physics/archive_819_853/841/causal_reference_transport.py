"""841: exact covariance transport and a containing-region entropy bound.

The periodic two-band finite CAR diagnostic is NOT the original T^3 PDE.
It has a fixed pure reference and two finite-range unitary layers; neither
an instantaneous vacuum replacement nor a Gibbs reference is used.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'causal_reference_transport_results.json'
sys.path.insert(0,str(HERE.parent/'839'))
from gaussian_reference_entropy_split import kron_all,state,relative,entropy

def spectral(a,f):
    w,u=np.linalg.eigh(a);return (u*f(w))@u.conj().T
def logm(a):
    w=np.linalg.eigvalsh(a);assert w.min()>0
    return spectral(a,np.log)
def logs(a):
    ident=np.eye(len(a));la=logm(a);lc=logm(ident-a)
    return lc-la,-la-lc
def covariance_relative(g,a):
    ident=np.eye(len(a))
    return float(np.trace(g@(logm(g)-logm(a))+
                         (ident-g)@(logm(ident-g)-logm(ident-a))).real)
def gaussian_data(n):
    # n_ij = <c_j^* c_i>; the conventional one-body Gibbs matrix is
    # log((1-n)/n), without another transpose in c_i^* h_ij c_j.
    m=len(n);z=np.diag([1.,-1.]);lower=np.array([[0,1],[0,0]],complex)
    cs=[kron_all([z]*j+[lower]+[np.eye(2)]*(m-j-1)) for j in range(m)]
    h,_=logs(n);k=np.zeros((2**m,2**m),complex)
    for i in range(m):
        for j in range(m):k+=h[i,j]*cs[i].conj().T@cs[j]
    rho=state(k)
    measured=np.array([[np.trace(rho@cs[j].conj().T@cs[i]) for j in range(m)] for i in range(m)])
    assert np.linalg.norm(measured-n)<2e-11
    # Recovering log(rho) from very small Fock eigenvalues loses precision.
    # Use the independently constructed quadratic K and its partition function.
    ew=np.linalg.eigvalsh(k);logz=float(-ew.min()+np.log(np.exp(-(ew-ew.min())).sum()))
    return rho,-k-logz*np.eye(len(k))
def gaussian_density(n):return gaussian_data(n)[0]
def layer(dim,sites,angle,seed):
    rng=np.random.default_rng(seed);ids=[2*s+i for s in sites for i in range(2)]
    raw=rng.normal(size=(len(ids),len(ids)))+1j*rng.normal(size=(len(ids),len(ids)))
    h=(raw+raw.conj().T)/2
    u=np.eye(dim,dtype=complex);u[np.ix_(ids,ids)]=spectral(h,lambda x:np.exp(-1j*angle*x))
    return u

def run():
    sites=8;m=2*sites;sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.])
    momenta=2*np.pi*np.arange(sites)/sites
    fourier=np.exp(1j*np.outer(np.arange(sites),momenta))/np.sqrt(sites)
    band=np.zeros((m,m),complex)
    for j,k in enumerate(momenta):
        h=np.sin(k)*sx+(.31+1-np.cos(k))*sz
        band[2*j:2*j+2,2*j:2*j+2]=(np.eye(2)-h/np.sqrt(np.sin(k)**2+(.31+1-np.cos(k))**2))/2
    ft=np.kron(fourier,np.eye(2));past=ft@band@ft.conj().T
    assert np.linalg.norm(past@past-past)<1e-13
    # U=U_01 U_12. Backward transport of site 0 reaches exactly sites 0,1,2.
    u=layer(m,(0,1),.24,841)@layer(m,(1,2),.19,8411)
    frame=np.eye(m)[:,:2];back=u.conj().T@frame
    support_residual=float(np.linalg.norm(back[6:,:]));assert support_residual<1e-14
    w=back[:6,:];assert np.linalg.norm(w.conj().T@w-np.eye(2))<1e-13
    now=u@past@u.conj().T;small=now[:2,:2];large=past[:6,:6]
    compression_residual=float(np.linalg.norm(w.conj().T@large@w-small))
    assert compression_residual<1e-13
    q=back@back.conj().T;r=np.eye(m)-q
    reset_past=q/2+r@past@r
    reset_now=u@reset_past@u.conj().T
    qnow=frame@frame.conj().T;rnow=np.eye(m)-qnow
    direct=qnow/2+rnow@now@rnow
    reset_residual=float(np.linalg.norm(reset_now-direct));assert reset_residual<1e-13
    gsmall=reset_now[:2,:2];glarge=reset_past[:6,:6]
    bsmall=covariance_relative(gsmall,small);blarge=covariance_relative(glarge,large)
    assert 0<=bsmall<=blarge+1e-12
    hs,ls=logs(small);hd,ld=logs(large)
    # Self-dual doubling of a number-conserving diagnostic: q=4, d=4.
    zsmall=2*float(np.trace(spectral(hs,np.abs)).real)
    zlarge=2*float(np.trace(w.conj().T@spectral(hd,np.abs)@w).real)
    loglarge=2*float(np.trace(w.conj().T@ld@w).real)
    jensen_min=float(np.linalg.eigvalsh(w.conj().T@ld@w-ls).min())
    assert jensen_min>-1e-11 and zsmall<=loglarge+1e-11
    assert loglarge<=zlarge+8*np.log(2)+1e-11
    # Independent full Fock calculation in the backward code + remainder basis.
    _,_,vh=np.linalg.svd(w.conj().T,full_matrices=True)
    basis=np.column_stack((w,vh.conj().T[:,2:]))
    a=basis.conj().T@large@basis;nu=gaussian_density(a[2:,2:])
    reference,log_reference=gaussian_data(a);g=np.kron(np.eye(4)/4,nu)
    bias=.6;sigma=(np.eye(4)+bias*np.diag([1.,-1.,-1.,1.]))/4
    rho=np.kron(sigma,nu);c=float(np.log(4)-entropy(sigma))
    dsmall=relative(sigma,gaussian_density(small))
    dlarge=float(-entropy(rho)-np.trace(rho@log_reference).real)
    fock_residual=abs(-entropy(g)-np.trace(g@log_reference).real-blarge)
    split_residual=max(abs(dsmall-bsmall-c),abs(dlarge-blarge-c))
    assert fock_residual<2e-10 and split_residual<2e-10 and dsmall<=dlarge+1e-10
    # Keeping the old local covariance instead of transporting the same state
    # gives a detectably different cost. A coordinate identification is not transport.
    unchanged_cost=covariance_relative(gsmall,past[:2,:2])
    assert abs(unchanged_cost-bsmall)>1e-4
    # Full past region is the pure global reference; any nontrivial reset fails
    # support containment. The finite small-region cost above does NOT diverge.
    leakage=float(np.trace((np.eye(m)-past)@reset_past).real)
    assert leakage>0
    return dict(round=841,all_checks_passed=True,fresh_test_groups=1,
        diagnostic='Eight-site two-band periodic pure CAR reference and two compact unitary layers; not the original continuum background.',
        backward_support_outside_containing_region=support_residual,
        transported_compression_residual=compression_residual,reset_transport_residual=reset_residual,
        reference_compression_minimum_eigenvalue=float(np.linalg.eigvalsh(large).min()),
        final_region_gaussian_cost=bsmall,past_containing_region_gaussian_cost=blarge,
        nonGaussian_record_cost=c,final_region_organized_cost=dsmall,past_containing_region_organized_cost=dlarge,
        finite_Fock_covariance_cost_residual=fock_residual,entropy_split_residual=split_residual,
        final_self_dual_modular_moment=zsmall,past_self_dual_code_modular_moment=zlarge,
        past_log_sum_code_bound=loglarge,compression_Jensen_minimum=jensen_min,
        wrong_untransported_reference_cost=unchanged_cost,
        whole_past_reference_cost='infinite',pure_reference_forbidden_sector_occupation=leakage,
        whole_region_failure_implies_small_region_failure=False,
        original_reference_transfer_is_analytic=True,original_T3_PDE_solved=False,
        original_past_local_modular_domain_verified=False,physical_history_inferred_from_auxiliary_past=False,
        area_response_or_unification_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
