"""756 candidate: same Hadamard two-point data, actual records, joint-source noise.

Finite-mode CAR checks and original background coefficients. Continuum positivity,
Hadamard regularity, Ward and linear-response statements are analytic, not inferred
from these finite matrices. No full quantum gravity or autonomous detector.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_reference_constraint_strata as background
import joint_matter_ground_source as matter_source

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_correlated_hadamard_noise_results.json'


def car(n):
    operators=[]
    for j in range(n):
        c=np.zeros((2**n,2**n),complex)
        for m in range(2**n):
            if (m>>j)&1:c[m^(1<<j),m]=(-1)**((m&((1<<j)-1)).bit_count())
        operators.append(c)
    return operators


def tau(p,q):
    assert max(0,2*p-1)-1e-14<=q<=p+1e-14
    return np.diag([1-2*p+q,p-q,p-q,q]).astype(complex)


def mean(rho,A):return np.trace(rho@A)
def covariance(rho,A,B):
    return (mean(rho,(A@B+B@A)/2)-mean(rho,A)*mean(rho,B)).real


def finite_mode_check():
    c=car(2);n=[a.conj().T@a for a in c]
    pair=c[0].conj().T@c[1].conj().T
    hop=c[0].conj().T@c[1]
    basis=[np.eye(4),*n,hop+hop.conj().T,1j*(hop-hop.conj().T),
           pair+pair.conj().T,1j*(pair-pair.conj().T)]
    parity=np.diag([1,-1,-1,1]);rows=[];err=comm=0.
    for p in (.2,.5,.8):
        ref=tau(p,p*p)
        for q in (max(0,2*p-1),p*p,p):
            r=tau(p,q);k=q-p*p
            assert np.linalg.eigvalsh(r).min()>-1e-14
            for A in basis:
                err=max(err,abs(mean(r,A)-mean(ref,A)))
                for B in basis:
                    comm=max(comm,abs(mean(r-ref,A@B-B@A)))
                    predicted=k*np.trace(parity@(A@B+B@A)/2).real
                    err=max(err,abs(covariance(r,A,B)-covariance(ref,A,B)-predicted))
            assert abs(mean(r,n[0]@n[1])-q)<1e-14
            pair_difference=covariance(r,basis[-2],basis[-2])-covariance(ref,basis[-2],basis[-2])
            assert abs(pair_difference-2*k)<1e-14
            rows.append(dict(p=p,q=q,kappa=k,record_probabilities=list(np.diag(r).real),
                             pair_source_variance_difference=float(pair_difference)))
    assert max(err,comm)<2e-14
    # Noise differences need not be positive: number-difference and pair-source.
    difference=np.diag([-2.,2.])*.25
    return dict(rows=rows,two_point_and_noise_identity_error=float(err),
                retarded_commutator_difference=float(comm),
                example_noise_difference_eigenvalues=list(np.linalg.eigvalsh(difference)),
                no_independent_positive_additive_noise_claim=True)


def fock_quadratic(h,d,c):
    out=np.zeros_like(c[0])
    for i in range(len(c)):
        for j in range(len(c)):
            out+=h[i,j]*(c[i].conj().T@c[j])
            out+=.5*(d[i,j]*(c[i].conj().T@c[j].conj().T)
                      +d[i,j].conjugate()*(c[j]@c[i]))
    return out


def complement_diagonal(probs):
    return np.array([np.prod([p if (m>>j)&1 else 1-p for j,p in enumerate(probs)])
                     for m in range(2**len(probs))])


def block_A(B):
    # Original lepton bits 6,7 = sterile modes 30,31. This partial trace
    # extracts their even quadratic coefficients, up to an irrelevant constant.
    return np.einsum('arbr->ab',B.reshape(4,64,4,64))/64


def original_lepton_transport_check():
    data=background.base(8)[0];phi=data['phi'][0,2,1]
    (h,d),(hs,ds)=matter_source.radial(float(np.linalg.norm(phi[:4])),float(phi[4]))
    ids=np.arange(24,32);ix=np.ix_(ids,ids)
    c=car(8);H=fock_quadratic(h[ix],d[ix],c);S=fock_quadratic(hs[ix],ds[ix],c)
    n6=c[6].conj().T@c[6];n7=c[7].conj().T@c[7]
    p=.5;q=.5;k=q-p*p
    rest=np.diag(complement_diagonal([.17,.29,.37,.43,.61,.73]))
    rho=np.kron(tau(p,q),rest);ref=np.kron(tau(p,p*p),rest)
    ev,U=np.linalg.eigh(H);parity=np.diag([1,-1,-1,1]);rows=[]
    mean_error=kernel_error=comm_error=0.
    evolved=[]
    for t in (0.,.3,.9):
        V=(U*np.exp(-1j*t*ev))@U.conj().T
        St=V.conj().T@S@V;Nt=V.conj().T@n6@V
        evolved.append((St,Nt))
        me=abs(mean(rho-ref,St));mean_error=max(mean_error,float(me))
        diff=covariance(rho,St,St)-covariance(ref,St,St)
        predicted=k*np.trace(parity@block_A(St)@block_A(St)).real
        kernel_error=max(kernel_error,abs(diff-predicted))
        rows.append(dict(time=t,same_scalar_force_mean_error=float(me),
                         scalar_force_noise_difference=float(diff),finite_mode_prediction=float(predicted)))
    for A,NA in evolved:
        for B,NB in evolved:
            comm_error=max(comm_error,float(abs(mean(rho-ref,A@B-B@A))))
            pred=k*np.trace(parity@(block_A(A)@block_A(B)+block_A(B)@block_A(A))/2).real
            kernel_error=max(kernel_error,abs(covariance(rho,A,B)-covariance(ref,A,B)-pred))
    assert max(mean_error,kernel_error,comm_error)<2e-12
    assert abs(mean(rho-ref,n6@n7)-k)<2e-14
    return dict(original_lepton_modes=8,Fock_dimension=256,
                original_Dirac_and_Majorana_retained=True,rows=rows,
                same_mean_error=mean_error,smooth_kernel_formula_error=float(kernel_error),
                retarded_difference=comm_error,
                calibration_only_not_continuum_PDE_or_full_quantum_scalar_H=True)


def original_curved_source_profile():
    q,psi,tensor,info,color=background.completed(12,1.)
    phi=q['phi'];s=phi[...,4]
    F= matter_source.matter.original.F(phi)
    Ys=matter_source.matter.Y['s'];volume=(2*np.pi)**3
    dmass=Ys*(F**-.5+s*s/(6*F**1.5))
    mass=Ys*s*F**-.5
    # chi = constant half-density, norm1 on the declared periodic spin bundle.
    eta_s=dmass/volume;eta_E=mass/volume
    integrated_s=eta_s.mean()*volume;integrated_E=eta_E.mean()*volume
    p=.5;rows=[]
    for joint in (0.,.25,.5):
        k=joint-p*p
        matrix=2*k*np.real(np.outer([integrated_E,integrated_s],np.conj([integrated_E,integrated_s])))
        rows.append(dict(p=p,q=joint,kappa=k,relative_integrated_mass_pair_energy_scalar_kernel=matrix.tolist()))
    assert rows[-1]['relative_integrated_mass_pair_energy_scalar_kernel'][1][1]>0
    # Independent finite differences of the original full32 mass pairing.
    point=phi[0,3,2].copy();ds=1e-6
    plus=point.copy();minus=point.copy();plus[4]+=ds;minus[4]-=ds
    numerical=(matter_source.matter.mass_matrices(plus)[1][30,31]
               -matter_source.matter.mass_matrices(minus)[1][30,31])/(2*ds)
    derivative_error=float(abs(numerical-dmass[0,3,2]))
    assert derivative_error<2e-10
    return dict(grid=12,original_color_amplitude=1.,F_min=float(F.min()),
                psi_min=float(psi.min()),psi_max=float(psi.max()),
                canonical_density_not_physical_volume_density=True,
                integrated_scalar_pair_coefficient=[float(integrated_s.real),float(integrated_s.imag)],
                integrated_mass_pair_energy_coefficient=[float(integrated_E.real),float(integrated_E.imag)],
                energy_label_is_mass_pair_component_not_complete_stress=True,
                derivative_error=derivative_error,rows=rows,
                smooth_noise_difference_only=True,
                full_noise_requires_spacetime_smearing=True,
                no_quantum_Einstein_covariance_numerically_solved=True)


def run():
    a=finite_mode_check();b=original_lepton_transport_check();c=original_curved_source_profile()
    deps=('research_note_620.md','research_note_634.md','research_note_730.md','research_note_732.md',
          'research_note_735.md','research_note_754.md','research_note_755.md',
          'joint_reference_constraint_strata.py','joint_matter_ground_source.py','joint_fermion_gauss_completion.py')
    return dict(round=756,tests_run=3,failures=0,errors=0,
                finite_positive_family=a,original_Majorana_transport=b,common_curved_background=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Positive finite neutral-mode correlations admit same-two-point Hadamard background extensions analytically. One joint-record parameter fixes the smooth difference of complete quadratic-source noise, while all one-point and quadratic commutator data agree. Finite matrices and753 background coefficients calibrate that connection, not full QFT or quantum Einstein dynamics. No global Gauss lift, autonomous preparation, finite-graph continuum identity, or independent gravity emergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
