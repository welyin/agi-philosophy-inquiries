"""816: instantaneous principal identity and fixed-mode short-time bridge.
Only the auxiliary Fourier integral is evaluated numerically. The original
constrained-background lift is an analytic argument in the round report.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
import shear_mode_leakage as early
old=early.original
TARGET=HERE/'instantaneous_leakage_bridge_results.json'

def principal_identity():
    rng=np.random.default_rng(816);eta=np.diag([-1.,1.,1.,1.])
    C=[np.eye(64),*old.GAMMA];error=0.
    for i in range(10):
        k=rng.normal(size=(4,4));k=(k+k.T)/2
        n=rng.normal(size=3);n/=np.linalg.norm(n)
        dc=[-.5*sum((eta[mu,mu]*k[mu,nu]*C[nu] for nu in range(4)),np.zeros((64,64),complex)) for mu in range(4)]
        H=sum((n[i]*old.GAMMA[i] for i in range(3)),np.zeros((64,64),complex))
        actual=(dc[0]@H+H@dc[0])/2-sum((n[i]*dc[i+1] for i in range(3)),np.zeros((64,64),complex))
        expected=np.dot(k[0,1:],n)*np.eye(64)+sum((((k[1:,1:]+k[0,0]*np.eye(3))@n)[i]*old.GAMMA[i]/2 for i in range(3)),np.zeros((64,64),complex))
        error=max(error,float(np.max(abs(actual-expected))))
    assert error<1e-12
    # Nonconformal does not imply invertibility on both spin polarizations.
    spin=[a[np.ix_([30,31],[30,31])] for a in old.GAMMA]
    matrix=.5*(np.eye(2)+spin[0])
    assert np.linalg.matrix_rank(matrix,tol=1e-12)==1
    return dict(full_Nambu_principal_error=error,
        nonconformal_rank_one_counterexample=dict(k00=0.,k0=[.5,0.,0.],
        spatial_k=[[1,0,0],[0,0,0],[0,0,0]],active_polarization_rank=1))

def hat_normalized(omega,tau,order):
    x,w=np.polynomial.legendre.leggauss(order);rho=np.exp(-1/(1-x*x));weight=w*rho
    return np.einsum('t,tab->ab',weight,np.cos(tau*x[:,None,None]*omega[None,:,:]))/sum(weight)

def short_time(tau,order,kappa=8.):
    e0,v0=np.linalg.eigh(old.hamiltonian(np.array([0.,0.,kappa])))
    gram=np.zeros((2,2),complex);error=0.;idx=[30,31]
    for sign in (-1,1):
        e1,v1=np.linalg.eigh(old.hamiltonian(np.array([float(sign),0.,kappa])))
        vertex=v1.conj().T@old.GAMMA[2]@v0
        omega=e1[:,None]-e0[None,:]-sign
        K=-kappa/4*v1@(vertex*hat_normalized(omega,tau,order))@v0.conj().T[:,idx]
        target=-kappa/4*old.GAMMA[2][:,idx]
        gram+=K.conj().T@K
        error=max(error,float(np.linalg.norm(K-target,2)/kappa))
    return dict(width=tau,fixed_kappa=kappa,normalized_operator_error=error,
        Gram_eigenvalues_over_kappa_squared=(np.linalg.eigvalsh(gram)/(kappa*kappa)).tolist())

def proper_subwindow_variance(order):
    x,w=np.polynomial.legendre.leggauss(order);chi=np.exp(-1/(1-x*x))
    q=np.zeros_like(x);inside=abs(x)<.45;q[inside]=np.exp(-1/(1-(x[inside]/.45)**2))
    weights=w*chi*chi;weights/=sum(weights)
    # Polarization is an eigenvector of an actual original Gamma_z.
    mean=float(np.dot(weights,q/2));square=float(np.dot(weights,q*q/4))
    variance=square-mean**2
    assert variance>1e-4
    return dict(mean=mean,projected_external_norm_squared=variance,
        two_disjoint_packet_and_C_partner_Gram=[variance]*4)

def run():
    rows=[]
    for t in (.2,.1,.05,.025,.0125):
        a=short_time(t,96);b=short_time(t,192)
        assert abs(a['normalized_operator_error']-b['normalized_operator_error'])<1e-12
        rows.append(b)
    assert rows[-1]['normalized_operator_error']<rows[0]['normalized_operator_error']/100
    assert abs(rows[-1]['Gram_eigenvalues_over_kappa_squared'][0]-.125)<1e-4
    # det((T0+aI)^*(T0+aI)) has only isolated zeros; no all-parameter claim.
    determinant={}
    for a in (0.,1.,1.5,2.,3.,4.,5.):
        T=np.diag([a-i for i in (1,2,3,4)])
        det=float(np.linalg.det(T.T@T))
        direct=float(np.prod([(a-i)**2 for i in (1,2,3,4)]))
        assert abs(det-direct)<1e-8
        determinant[str(a)]=det
    return dict(round=816,all_checks_passed=True,principal=principal_identity(),
        fixed_mode_short_time=rows,instantaneous_Gram_limit=.125,
        two_packet_variance=proper_subwindow_variance(192),
        finite_rank_pencil_determinants=determinant,
        scope='Coefficient and continuity diagnostics; original background and state are not numerically replaced.',
        original_future_PDE_or_W_B_computed=False,
        previously_fixed_process_classified=False,
        all_possible_record_architectures_excluded=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

