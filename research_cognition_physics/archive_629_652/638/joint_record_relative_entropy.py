"""638: same original record, reference support and relative-entropy matching.

Full Gibbs statements are analytic, not a small-matrix replacement. Numeric
matrices below are Gram representations of ACTUAL original Gauss packets.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_record_band_control as records

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_record_relative_entropy_results.json'
prior=records.prior
instrument=records.entry.instrument


def entropy(matrix):
    eig=np.linalg.eigvalsh((matrix+matrix.conj().T)/2)
    assert eig.min()>-3e-13
    eig=eig[eig>1e-15]
    return float(-np.sum(eig*np.log(eig)))


def relative(rho,sigma):
    ev,U=np.linalg.eigh((sigma+sigma.conj().T)/2)
    assert ev.min()>1e-14,ev
    return float(-entropy(rho)-np.sum(np.diag(U.conj().T@rho@U).real*np.log(ev)))


def reset(matrix,n):
    out=matrix[:n,:n].copy()
    out[0,0]+=np.trace(matrix[n:,n:])
    return out


def record_energy_check():
    rng=np.random.default_rng(638)
    defects=[];kernel_error=0.
    for _ in range(40):
        phi=rng.normal(size=5)*.3
        K=prior.model.original.inverse(phi)
        L,dL,_=instrument(phi[-1])
        psi=rng.normal()+1j*rng.normal()
        grad=rng.normal(size=5)+1j*rng.normal(size=5)
        before=np.vdot(grad,K@grad).real
        after=0.
        for r in range(2):
            new=L[r]*grad
            new[-1]+=dL[r]*psi
            after+=np.vdot(new,K@new).real
        predicted=K[-1,-1]*np.sum(dL*dL)*abs(psi)**2
        defects.append(float(abs(after-before-predicted)))
        s1,s2=rng.uniform(-1,1,2)
        l1=instrument(s1)[0];l2=instrument(s2)[0]
        theta1=np.arccos(l1[0]);theta2=np.arccos(l2[0])
        kernel_error=max(kernel_error,float(abs(l1@l2-np.cos(theta1-theta2))))
    assert max(defects)<8e-14 and kernel_error<1e-14
    rows=[]
    for n in (96,160):
        _,s,p,Kss=prior.packet_data(n)
        L,dL,_=instrument(s)
        Q=prior.HBAR**2/(2*prior.W)*float(np.sum(p*Kss*np.sum(dL*dL,axis=-1)))
        rows.append(dict(quadrature_order=n,actual_Gauss_packet_injected_energy=Q,
                         probabilities=np.sum(p[...,None]*L*L,axis=(0,1)).tolist()))
    delta=prior.HBAR**2*prior.model.original.M/(32*prior.W)
    assert abs(rows[1]['actual_Gauss_packet_injected_energy']-
               rows[0]['actual_Gauss_packet_injected_energy'])<2e-13
    assert 0<rows[1]['actual_Gauss_packet_injected_energy']<delta
    return dict(original_full_inverse_metric_Dirichlet_error=max(defects),
        original_random_unitary_kernel_error=kernel_error,packet_rows=rows,
        full_model_universal_injection_bound=delta,
        chosen_beta=2.,full_Gibbs_unconditional_relative_entropy_upper=2*delta,
        full_Gibbs_record_retained_relative_entropy_upper=2*delta+float(np.log(2)),
        actual_full_Gibbs_energy_and_relative_entropy_not_numerically_evaluated=True)


def phase_support(n):
    _,s,p,_=prior.packet_data(n)
    functions=np.exp(1j*s.ravel()[:,None]*np.array([-24,0,24])[None,:])
    vectors=np.sqrt(p.ravel())[:,None]*functions
    G=vectors.conj().T@vectors
    ev,U=np.linalg.eigh(G)
    assert ev.min()>.05
    orth=vectors@(U*(ev**-.5))@U.conj().T
    L=instrument(s)[0].reshape(-1,2)
    leak=0.;trace=0.
    for r in range(2):
        image=L[:,r,None]*orth
        compressed=orth.conj().T@image
        trace+=float(np.linalg.norm(image,'fro')**2/3)
        # Projection complement evaluated as a vector residual, avoiding
        # subtraction of two nearly equal total traces.
        residual=image-orth@compressed
        leak+=float(np.linalg.norm(residual,'fro')**2/3)
    assert abs(trace-1)<2e-14 and leak>1e-7
    return dict(quadrature_order=n,Gram_eigenvalues=ev.tolist(),
        density_trace_after_record=trace,strictly_positive_support_leakage=leak)


def support_check():
    coarse=phase_support(96);fine=phase_support(160)
    error=abs(coarse['strictly_positive_support_leakage']-
              fine['strictly_positive_support_leakage'])
    assert error<2e-13
    return dict(rows=[coarse,fine],quadrature_change=error,
        reference_uniform_on_three_original_Gauss_phase_packets=True,
        relative_entropy_after_untruncated_record_is_infinite_by_support=True,
        phase_packet_span_is_not_claimed_to_be_an_H_spectral_band=True,
        all_nonzero_finite_H_band_obstruction_proved_analytically=True)


def original_pair(n):
    _,s,p,_=prior.packet_data(n)
    L=instrument(s)[0].reshape(-1,2)
    functions=np.column_stack((np.ones(s.size),L))
    vectors=np.sqrt(p.ravel())[:,None]*functions
    # Columns represent |psi>, L+|psi>, L-|psi> in the original Gauss space.
    orth,R=np.linalg.qr(vectors,mode='reduced')
    phase=np.sign(np.diag(R));phase[phase==0]=1
    orth=orth*phase;R=phase[:,None]*R
    assert np.linalg.norm(vectors-orth@R)<1e-13
    a=R[:,0];v=R[:,1:]
    pure=np.outer(a,a.conj())
    tau=v@v.conj().T
    # Explicit original physical reference, NOT a Gibbs state.
    sigma=.4*pure+.6*tau
    p_out=np.sum(abs(v)**2,axis=0)
    branches=[np.outer(v[:,r],v[:,r].conj())/p_out[r] for r in range(2)]
    D=relative(tau,sigma)
    average=float(sum(p_out[r]*relative(branches[r],sigma) for r in range(2)))
    holevo=entropy(tau)-sum(p_out[r]*entropy(branches[r]) for r in range(2))
    assert abs(average-D-holevo)<2e-11
    assert np.linalg.norm(tau@sigma-sigma@tau)>1e-7
    ds=[relative(reset(tau,k),reset(sigma,k)) for k in (1,2,3)]
    assert abs(ds[0])<2e-13 and abs(ds[-1]-D)<1e-13
    assert min(np.diff(ds))>-2e-11
    composition=float(np.max(abs(reset(reset(tau,2),1)-reset(tau,1))))
    assert composition<1e-14
    return dict(quadrature_order=n,reference_minimum_eigenvalue=float(np.linalg.eigvalsh(sigma).min()),
        noncommutativity_norm=float(np.linalg.norm(tau@sigma-sigma@tau)),
        full_relative_entropy=D,average_conditional_relative_entropy=average,
        Holevo_record_information=holevo,classical_record_entropy=float(-p_out@np.log(p_out)),
        common_reset_relative_entropies=ds,reset_composition_error=composition,
        actual_original_H_spectral_matrices_not_computed=True)


def common_reference_check():
    rows=[original_pair(96),original_pair(160)]
    error=max(abs(rows[0][key]-rows[1][key]) for key in
              ('full_relative_entropy','average_conditional_relative_entropy',
               'Holevo_record_information'))
    assert error<2e-10
    return dict(rows=rows,quadrature_change=error,
        all_fixture_vectors_are_original_compact_Gauss_states=True,
        non_Gibbs_reference_declared=True,
        original_full_Gibbs_convergence_is_analytic=True)


def run():
    deps=('research_note_305.md','research_note_317.md','research_note_578.md',
          'research_note_579.md','research_note_592.md','research_note_598.md',
          'research_note_603.md','research_note_617.md','research_note_637.md',
          'joint_record_band_control.py','joint_geometry_work_noise.py',
          'round592_drafts/record_band_entry.py','round638_drafts/relative_entropy_entry.md')
    return dict(round=638,tests_run=3,failures=0,errors=0,
        record_energy=record_energy_check(),reference_support=support_check(),
        common_reference=common_reference_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(full_original_Gauss_Gibbs_reference_and_record_used_analytically=True,
            finite_record_retained_relative_entropy_bound=True,
            naive_finite_reference_untruncated_record_rejected=True,
            common_Gauss_spectral_reset_transports_relative_entropy=True,
            fixed_region_relative_entropy_limit_retained=True,
            numerical_witnesses_are_original_Gauss_packets_not_arbitrary_matrices=True,
            no_spatial_continuum_or_area_law_or_GR_or_autonomous_instrument_claim=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False))
