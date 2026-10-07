"""888: bare versus explicitly transported quantum-coordinate energy.
Original continuum EM family, full matter trace, and exact electron Fock audit.
Transport is a declared change of kinetic law, not a free frame redefinition.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'887'))
import conditional_thermal_kinetic_certificate as prior
last=prior.working.last;old=prior.old
TARGET=HERE/'reference_transport_energy_results.json'

def sech(x):
    e=np.exp(-np.abs(x))
    return 2*e/(1+e*e)

def metrics(N,alpha,q,m2):
    k,r,c=last.weight.radial_counts(N);spec=0.;rot=0.
    for x in k:
        px=x+alpha*q
        e2=px[:,None]**2+m2[:,None]+r
        e=np.sqrt(e2);z=np.exp(-old.BETA*e);v=z/(1+z)**2
        spec+=old.BETA**2/8*float(np.sum(c*q[:,None]**2*px[:,None]**2/e2*v))
        rot+=float(np.sum(c*q[:,None]**2/(16*e2)*(1-px[:,None]**2/e2)*(1-sech(old.BETA*e))))
    return dict(N=N,alpha=alpha,spectral_quarter_Fisher=spec,
        eigenvector_quarter_Fisher=rot,bare_quarter_Fisher=spec+rot,
        eigenvector_term_over_N=rot/N)

def leading_coefficient(q,order=64):
    x,w=np.polynomial.legendre.leggauss(order);y=x/2
    area=float(np.sum((w[:,None]*w[None,:]/4)/(.25+y[:,None]**2+y[None,:]**2)))
    return float(np.sum(q*q)/8*area)

def fock_audit(m,q,k,alpha):
    ids=[26,27,28,29]
    h=old.H(m,q,17,k,alpha,kind='continuum')[np.ix_(ids,ids)]
    g=(m.GAMMA[0]@np.diag(q))[np.ix_(ids,ids)]
    e,u=np.linalg.eigh(h);ge=u.conj().T@g@u;occ=(e<0).astype(float)
    dP=np.zeros((4,4),complex)
    for i in range(4):
        for j in range(4):
            if occ[i]!=occ[j]:dP[i,j]=(occ[i]-occ[j])*ge[i,j]/(e[i]-e[j])
    P=np.diag(occ);A=dP@P-P@dP;A=u@A@u.conj().T
    cs=[prior.annihilator(j) for j in range(4)]
    def lift(x):return sum(x[i,j]*(cs[i].conj().T@cs[j]) for i in range(4) for j in range(4))
    B=lift(h);G=lift(g);AF=lift(A)
    ev,V=np.linalg.eigh(B);gp=V.conj().T@G@V;ap=V.conj().T@AF@V
    p=np.exp(-old.BETA*(ev-ev.min()));p/=p.sum()
    mean=float(np.dot(p,np.diag(gp).real))
    dg=np.zeros_like(gp)
    for i in range(16):
        for j in range(16):
            if abs(ev[i]-ev[j])<1e-9:dg[i,j]=-old.BETA*p[i]*gp[i,j]
            else:dg[i,j]=(p[i]-p[j])/(ev[i]-ev[j])*gp[i,j]
        dg[i,i]+=old.BETA*p[i]*mean
    quarter=.5*float(np.sum(abs(dg)**2/(p[:,None]+p[None,:])))
    spec=.25*float(np.sum(abs(np.diag(dg))**2/p))
    transported=dg-(ap@np.diag(p)-np.diag(p)@ap)
    off=float(np.max(abs(transported-np.diag(np.diag(transported)))))
    mass2=abs(np.exp(old.XI)*m.MASS[26,28])**2
    px=k[0]-6*alpha;e2=px*px+k[1]**2+k[2]**2+mass2;en=np.sqrt(e2)
    z=np.exp(-old.BETA*en);v=z/(1+z)**2
    es=old.BETA**2*36*px*px/e2*v
    er=36/(2*e2)*(1-px*px/e2)*(1-sech(old.BETA*en))
    assert abs(quarter-es-er)<2e-10 and abs(spec-es)<2e-10 and off<2e-12
    # Projector identity checked in original particle basis.
    pmat=u@P@u.conj().T;dp=u@dP@u.conj().T
    perr=float(np.max(abs(dp-(A@pmat-pmat@A))))
    assert perr<2e-13 and np.max(abs(A+A.conj().T))<1e-13
    return dict(momentum=k,alpha=alpha,Fock_quarter_Fisher=quarter,
        spectral_part=spec,eigenvector_part=quarter-spec,
        analytic_spectral_part=es,analytic_eigenvector_part=er,
        transported_density_offdiagonal_error=off,projector_transport_error=perr)

def run():
    m=old.load();_,matter=old.charges(m);q=last.em_charges(m,matter)
    qq,m2=last.joint_spectrum(m,q)
    checks=[fock_audit(m,q,k,.02) for k in ([0,0,0],[1,2,-1],[0,3,1])]
    rows=[metrics(N,.02,qq,m2) for N in (9,17,33,65,129)]
    coeff=leading_coefficient(q)
    assert abs(coeff-leading_coefficient(q,32))<1e-10
    assert abs(rows[-1]['spectral_quarter_Fisher']-rows[-2]['spectral_quarter_Fisher'])<1e-10
    assert rows[-1]['eigenvector_term_over_N']<coeff
    return dict(round=888,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3673,
        argument_scope='The original continuum EM thermal family already has a linearly divergent bare-coordinate Bures kinetic lower bound from changing vacuum eigenvectors. Its eigenvalue-only part is finite. Finite-cutoff projector transport can remove the eigenvector cost only by an explicitly changed covariant kinetic operator, not by free relabeling. It cannot remove887 spectral cut obstruction.',
        beta=old.BETA,full_original_components=64,full_charge_square_trace=float(np.sum(q*q)),
        original_Fock_checks=checks,continuum_cutoff_rows=rows,
        analytic_linear_coefficient=coeff,
        original_neutral_Majorana_block_kept=True,
        bare_uniform_budget_invalid_as_automatic_target_property=True,
        finite_cutoff_covariant_normal_state_construction=True,
        transport_requires_changed_kinetic_input=True,
        mere_basis_change_cancels_physical_energy=False,
        full_dynamic_Gauss_Hamiltonian_equivalence_proved=False,
        infinite_Fock_transport_unitary_claimed=False,
        full_interacting_Q_E_bridge_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
