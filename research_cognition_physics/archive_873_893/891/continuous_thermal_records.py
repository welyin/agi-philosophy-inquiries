"""891: common thermal process and original finite physical-CAR record dictionary.
The full-species trace convergence is analytic; numerical joint propagation is
an original four-electron-mode calibration with a shared quantum coordinate.
The transported kinetic and pointwise vacuum subtraction remain explicit inputs.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'890'))
import shared_coordinate_heat_bound as prev
old=prev.old;prior=prev.prior
TARGET=HERE/'continuous_thermal_records_results.json'

def exterior(U):
    n=len(U);F=np.zeros((1<<n,1<<n),complex)
    sets=[[i for i in range(n) if (b>>i)&1] for b in range(1<<n)]
    for b,ib in enumerate(sets):
        for c,ic in enumerate(sets):
            if len(ib)==len(ic):
                F[b,c]=1 if not ib else np.linalg.det(U[np.ix_(ib,ic)])
    return F

def spectral_rotation(h,h0):
    # All four original electron modes share one mass and a Clifford plane.
    e=np.sqrt(np.trace(h@h).real/4);e0=np.sqrt(np.trace(h0@h0).real/4)
    S=h/e;S0=h0/e0;c=np.trace(S@S0).real/4
    U=(np.eye(4)+S@S0)/np.sqrt(2*(1+c))
    assert np.max(abs(U.conj().T@U-np.eye(4)))<1e-12
    assert np.max(abs(U@S0@U.conj().T-S))<1e-12
    return U

def functions(H,beta,time):
    e,V=np.linalg.eigh(H)
    p=np.exp(-beta*(e-e.min()));p/=p.sum()
    return (V*p)@V.conj().T,(V*np.exp(-1j*time*e))@V.conj().T

def records(rho,U,A,B,prepare):
    rho=prepare@rho@prepare.conj().T
    I=np.eye(len(rho));out=[]
    for P in (A,I-A):
        post=U@P@rho@P@U.conj().T
        for Q in (B,I-B):
            out.append(float(np.trace(Q@post@Q).real))
    assert min(out)>-1e-12 and abs(sum(out)-1)<1e-11
    return np.array(out)

def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.last.em_charges(m,matter)
    # Every original64 block, including neutral Majorana, matches on the
    # common interior dictionary; these are original physical matrices.
    matrix_checks=[]
    for N in (9,17,33):
        err=0.
        for a in np.linspace(prior.CENTER-prior.HALF_INTERVAL,prior.CENTER+prior.HALF_INTERVAL,5):
            for k in ([0,1,0],[1,0,-1],[-2,1,1]):
                err=max(err,float(np.max(abs(old.H(m,q,N,k,a)-old.H(m,q,N,k,a,kind='continuum')))))
        assert err<2e-12
        matrix_checks.append(dict(N=N,original64_interior_error=err))
    ids=[26,27,28,29];k=[0,1,0]
    cs=[prior.prior.prior.annihilator(j) for j in range(4)]
    def lift(h):return sum(h[i,j]*(cs[i].conj().T@cs[j]) for i in range(4) for j in range(4))
    n=11;dim=16*n;dx=2*prior.HALF_INTERVAL/(n+1)
    alphas=prior.CENTER+np.linspace(-prior.HALF_INTERVAL+dx,prior.HALF_INTERVAL-dx,n)
    h0=old.H(m,q,9,k,prior.CENTER,kind='continuum')[np.ix_(ids,ids)]
    e0,V0=np.linalg.eigh(h0)
    T=np.diag(np.full(n,2*prior.KAPPA/dx**2))+np.diag(np.full(n-1,-prior.KAPPA/dx**2),1)+np.diag(np.full(n-1,-prior.KAPPA/dx**2),-1)
    W=np.zeros((dim,dim),complex);Bphys=np.zeros_like(W);Bdiag=np.zeros_like(W)
    transport_residual=0.;anti_residual=0.;eig_error=0.
    for j,a in enumerate(alphas):
        h=old.H(m,q,9,k,a,kind='continuum')[np.ix_(ids,ids)]
        rotation=spectral_rotation(h,h0);V=rotation@V0;WF=exterior(V)
        hdiag=V.conj().T@h@V
        eig_error=max(eig_error,float(np.max(abs(hdiag-np.diag(np.diag(hdiag))))))
        ev=np.linalg.eigvalsh(h);vac=float(ev[ev<0].sum())
        sl=slice(16*j,16*(j+1));W[sl,sl]=WF
        Bphys[sl,sl]=lift(h)-vac*np.eye(16)
        Bdiag[sl,sl]=lift(np.diag(np.diag(hdiag)))-vac*np.eye(16)
        eps=1e-6
        def rotate_at(a2):
            hh=old.H(m,q,9,k,a2,kind='continuum')[np.ix_(ids,ids)]
            return spectral_rotation(hh,h0)
        dU=(rotate_at(a+eps)-rotate_at(a-eps))/(2*eps)
        E=np.sqrt(np.trace(h@h).real/4)
        G=(m.GAMMA[0]@np.diag(q))[np.ix_(ids,ids)]
        dE=np.trace(h@G).real/(4*E)
        S=h/E;dS=G/E-h*dE/E**2
        A=(dS@S-S@dS)/4
        transport_residual=max(transport_residual,float(np.max(abs(dU-A@rotation))))
        anti_residual=max(anti_residual,float(np.max(abs(A+A.conj().T))))
    assert transport_residual<1e-8 and eig_error<2e-12
    Hocc=np.kron(T,np.eye(16))+Bdiag
    # Construct the physical covariant hopping independently, with one T.
    Hphys=Bphys.copy()
    for i in range(n):
        for j in range(n):
            si=slice(16*i,16*(i+1));sj=slice(16*j,16*(j+1))
            Hphys[si,sj]+=T[i,j]*W[si,si]@W[sj,sj].conj().T
    joint_error=float(np.max(abs(Hphys-W@Hocc@W.conj().T)))
    assert joint_error<1e-10
    rho,U=functions(Hphys,old.BETA,.173)
    sigma,V=functions(Hocc,old.BETA,.173)
    thermal_error=float(np.linalg.norm(rho-W@sigma@W.conj().T,ord='fro'))
    unitary_error=float(np.max(abs(U-W@V@W.conj().T)))
    assert thermal_error<1e-11 and unitary_error<1e-10
    A0=cs[0].conj().T@cs[0]
    f=(cs[0]+1j*cs[2])/np.sqrt(2);B0=f.conj().T@f
    quartic=(cs[0].conj().T@cs[0])@(cs[1].conj().T@cs[1])
    P0=np.eye(16)+(np.exp(.71j)-1)*quartic
    A=np.kron(np.eye(n),A0);B=np.kron(np.eye(n),B0);prep=np.kron(np.eye(n),P0)
    p=records(rho,U,A,B,prep)
    transformed=records(sigma,V,W.conj().T@A@W,W.conj().T@B@W,W.conj().T@prep@W)
    wrong=records(sigma,V,A,B,prep)
    unprepared=records(rho,U,A,B,np.eye(dim))
    record_error=float(np.max(abs(p-transformed)));wrong_error=float(np.max(abs(p-wrong)))
    assert record_error<1e-11 and wrong_error>1e-3
    # Positivity/normalization and actual effect of original quartic preparation.
    prep_effect=float(np.max(abs(p-unprepared)))
    assert prep_effect>1e-6
    return dict(round=891,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3676,
        argument_scope='In the explicitly transported, pointwise-vacuum-subtracted free branch, full original-species joint Gibbs operators converge in trace norm on a common excitation representation. Compatible finite original physical CAR dictionaries preserve bounded preparations and finite time-ordered records. No infinite Fock implementation of the original reference transformation, original bare Gauss dynamics, vacuum sources, or interacting Q/E bridge is asserted.',
        original64_checks=matrix_checks,
        coordinate_grid_points=n,original_electron_Fock_dimension=16,total_joint_dimension=dim,
        spectral_rotation_derivative_residual=transport_residual,antihermitian_residual=anti_residual,
        original_block_diagonalization_error=eig_error,
        independently_built_covariant_Hamiltonian_error=joint_error,
        joint_thermal_frame_error_Frobenius=thermal_error,real_time_frame_error=unitary_error,
        original_physical_two_record_probabilities=p.tolist(),
        transported_physical_two_record_probabilities=transformed.tolist(),
        incorrectly_untransported_record_probabilities=wrong.tolist(),
        record_dictionary_error=record_error,missing_transport_probability_error=wrong_error,
        original_quartic_preparation_changes_records_by=prep_effect,
        full_original_species_continuous_joint_trace_norm_bridge=True,
        compatible_finite_original_CAR_record_dictionary=True,
        original_neutral_Majorana_preserved=True,
        bounded_even_preparation_need_not_be_Gaussian=True,
        fixed_original_infinite_Fock_transport_unitary_claimed=False,
        full_local_Gauss_dynamics_and_vacuum_sources_matched=False,
        unbounded_source_derivative_convergence_proved=False,
        original_interacting_Q_E_equivalence_proved=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
