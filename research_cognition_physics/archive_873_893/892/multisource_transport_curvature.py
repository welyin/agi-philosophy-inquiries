"""892: test the actual888 horizontal transport in several EM directions.
Uses original64 mass/charge matrices and an original neutral even electron
pair record. Nonzero canonical curvature does not exclude flat completions.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'891'))
import continuous_thermal_records as prior
old=prior.old
TARGET=HERE/'multisource_transport_curvature_results.json'
def comm(a,b):return a@b-b@a

def family(m,q,a,k):
    return sum(m.GAMMA[i]@np.diag(k[i]+q*a[i]) for i in range(3))+np.exp(old.XI)*m.MASS

def jets(m,q,a,k):
    H=family(m,q,a,k);e,U=np.linalg.eigh(H)
    inv=(U*(1/abs(e)))@U.conj().T
    P=(np.eye(64)-H@inv)/2
    dP=[];A=[]
    for i in range(3):
        G=m.GAMMA[i]@np.diag(q)
        Eder=np.diag(q*(k[i]+q*a[i]))@inv
        dp=-.5*(G@inv-H@inv@inv@Eder)
        dP.append(dp);A.append(comm(dp,P))
    return H,P,dP,A

def electron_H(m,q,a):
    ids=[26,27,28,29]
    return family(m,q,a,[0,0,0])[np.ix_(ids,ids)]

def segment(h1,h0):
    e1=np.sqrt(np.trace(h1@h1).real/4);e0=np.sqrt(np.trace(h0@h0).real/4)
    S1=h1/e1;S0=h0/e0;c=np.trace(S1@S0).real/4
    return (np.eye(4)+S1@S0)/np.sqrt(2*(1+c))

def wedge(v,w):
    out=np.zeros(16,complex)
    for i in range(4):
        for j in range(i+1,4):
            out[(1<<i)|(1<<j)]=v[i]*w[j]-v[j]*w[i]
    return out

def loop_audit(m,q,F,side,cs):
    origin=np.zeros(3)
    corners=[origin,np.array([side,0,0]),np.array([side,side,0]),np.array([0,side,0]),origin]
    H0=electron_H(m,q,origin)
    hol=np.eye(4,dtype=complex);flat=np.eye(4,dtype=complex)
    for a,b in zip(corners[:-1],corners[1:]):
        ha=electron_H(m,q,a);hb=electron_H(m,q,b)
        hol=segment(hb,ha)@hol
        # Explicit radial-frame flat completion; its transport telescopes.
        Ua=segment(ha,H0);Ub=segment(hb,H0)
        flat=(Ub@Ua.conj().T)@flat
    ev,V=np.linalg.eigh(H0);neg=V[:,:2];pos=V[:,2:]
    eig,C=np.linalg.eigh(1j*pos.conj().T@F@pos)
    basis=pos@C
    psi_particle=(basis[:,0]+basis[:,1])/np.sqrt(2)
    psiF=wedge(neg[:,0],psi_particle);psiF/=np.linalg.norm(psiF)
    readvec=(basis[:,0]+1j*basis[:,1])/np.sqrt(2)
    annih=sum(np.conj(readvec[j])*cs[j] for j in range(4))
    read=annih.conj().T@annih
    WF=prior.exterior(hol)
    before=float(np.vdot(psiF,read@psiF).real)
    after=float(np.vdot(WF@psiF,read@(WF@psiF)).real)
    B=sum(H0[i,j]*cs[i].conj().T@cs[j] for i in range(4) for j in range(4))
    e,Z=np.linalg.eigh(B);prob=np.exp(-old.BETA*(e-e.min()));prob/=prob.sum()
    gamma=(Z*prob)@Z.conj().T
    thermal=float(np.max(abs(WF@gamma@WF.conj().T-gamma)))
    vac=complex(np.linalg.det(neg.conj().T@hol@neg))
    error=float(np.linalg.norm((hol-np.eye(4))/side**2+F)/np.linalg.norm(F))
    assert abs(before-.5)<1e-12 and abs(after-before)>side**2*100
    assert thermal<2e-12 and abs(vac-1)<2e-12
    assert np.max(abs(flat-np.eye(4)))<2e-12
    return dict(side=side,horizontal_holonomy_norm=float(np.linalg.norm(hol-np.eye(4))),
        small_loop_relative_curvature_error=error,
        flat_completion_loop_error=float(np.max(abs(flat-np.eye(4)))),
        neutral_even_pair_read_before=before,neutral_even_pair_read_after=after,
        signed_record_change_over_area=(after-before)/side**2,
        original_Fock_Gibbs_holonomy_change=thermal,
        filled_negative_band_vacuum_phase_error=abs(vac-1))

def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.prior.last.em_charges(m,matter)
    a=np.zeros(3);k=np.zeros(3)
    H,P,dp,A=jets(m,q,a,k);F=comm(dp[0],dp[1])
    vals,V=np.linalg.eigh(H);invsq=(V*(1/vals**2))@V.conj().T
    expected=m.GAMMA[0]@m.GAMMA[1]@np.diag(q*q)@invsq/2
    analytic_error=float(np.max(abs(F-expected)))
    assert analytic_error<1e-10
    preserved=float(np.max(abs(comm(F,P))))
    assert preserved<1e-10 and np.max(abs(F+F.conj().T))<1e-10
    ids=[26,27,28,29];FE=F[np.ix_(ids,ids)]
    electron_mass2=float(abs(np.exp(old.XI)*m.MASS[26,28])**2)
    coeff=36/(2*electron_mass2)
    norm=float(np.linalg.norm(FE,2))
    assert abs(norm-coeff)<1e-9
    derivative=[]
    for pos,momentum in ((np.zeros(3),np.zeros(3)),(np.array([.02,-.01,.015]),np.array([1.,2.,-1.]))):
        H1,P1,p1,A1=jets(m,q,pos,momentum);f1=comm(p1[0],p1[1]);step=1e-6
        aa=pos.copy();aa[0]+=step;bb=pos.copy();bb[0]-=step
        d1A2=(jets(m,q,aa,momentum)[3][1]-jets(m,q,bb,momentum)[3][1])/(2*step)
        aa=pos.copy();aa[1]+=step;bb=pos.copy();bb[1]-=step
        d2A1=(jets(m,q,aa,momentum)[3][0]-jets(m,q,bb,momentum)[3][0])/(2*step)
        curvature=-d1A2+d2A1+comm(A1[0],A1[1])
        err=float(np.linalg.norm(curvature-f1)/max(1,np.linalg.norm(f1)))
        assert err<3e-7
        derivative.append(dict(alpha=pos.tolist(),momentum=momentum.tolist(),relative_covariant_curvature_identity_error=err))
    cs=[prior.prior.prior.prior.annihilator(i) for i in range(4)]
    loops=[loop_audit(m,q,FE,s,cs) for s in (.001,.0005,.00025,.000125)]
    assert loops[-1]['small_loop_relative_curvature_error']<loops[0]['small_loop_relative_curvature_error']/6
    return dict(round=892,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3677,
        argument_scope='The actual888 projector-horizontal connection has nonzero original Dirac curvature for independent EM holonomy directions, so it cannot be removed by one path-independent U solving all horizontal equations. Static thermal/vacuum checks can miss its neutral-even-pair record holonomy. A flat completion exists on this gapped patch but changes the band-internal kinetic connection; general transport is not ruled out.',
        original_full_components=64,original_electron_mass_squared=electron_mass2,
        full_original_curvature_formula_error=analytic_error,
        curvature_commutes_with_energy_projection_error=preserved,
        electron_curvature_operator_norm=norm,analytic_electron_curvature=coeff,
        full_original_covariant_derivative_checks=derivative,original_neutral_pair_loop_checks=loops,
        canonical_multidirection_connection_is_flat=False,
        single_path_independent_horizontal_frame_exists=False,
        gapped_local_flat_completion_exists=True,
        flat_completion_preserves_same_horizontal_kinetic=False,
        original_neutral_pair_net_EM_charge=0,original_neutral_pair_even_parity=True,
        filled_vacuum_and_thermal_invariance_implies_full_process_equivalence=False,
        local_Gauss_or_interacting_Q_E_bridge_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
