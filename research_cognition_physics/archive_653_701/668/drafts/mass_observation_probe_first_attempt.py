"""Actual668 entry: original mass in the retained Weyl observation algebra.

Frozen phi, free finite Euclidean box, declared mass pairing convention.
This is not identification with the full original interacting Hamiltonian.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_physical_auxiliary_state as prior
import joint_spinor_subgroup_mass as dictionary
import joint_gauss_fermion_influence as car
old=prior.old


def mass_pairing(count):
    phi=car.PHI[0].copy()
    b=dictionary.ph_matrix()@dictionary.bdg(phi)@dictionary.ph_matrix().conj().T
    j=np.kron(dictionary.dictionary(),np.eye(2))
    delta=j@b[:32,32:]@j.T
    # Exterior internal,spin order to physical spin,internal order.
    delta=delta.reshape(16,2,16,2).transpose(1,0,3,2).reshape(32,32)
    p=np.kron(np.eye(count),delta.conj())
    # Declared Euclidean exp(-mass H) convention: ww=Delta*, barbar=-Delta.
    return old.diag(p,-p.conj())


def run():
    rows=[]
    for nx,seed in ((1,66801),(2,66802)):
        for sample in range(2):
            rng=np.random.default_rng(seed+sample)
            e=rng.normal(size=(2*nx,10))*.17;e[:,0]+=1
            e/=np.linalg.norm(e,axis=1)[:,None]
            q=old.data(nx=nx,nt=2,e=e);f=prior.weyl(q);r=q['r']
            pair=mass_pairing(len(q['sites']))
            assert old.err(pair+pair.T)<1e-13
            reflect_error=old.err(pair+f['Theta'].T@pair.conj()@f['Theta'])
            assert reflect_error<1e-13
            jm=np.kron(np.kron(np.eye(2*nx),old.old.old.VM),np.eye(16))
            jp=np.kron(np.kron(np.eye(2*nx),old.old.old.VP),np.eye(16))
            bare=old.diag(jm.conj().T,jp.T)
            pull=f['map'].T@pair@f['map']
            wrong=bare.T@pair@bare
            ti=np.linalg.inv(q['T'])
            select=f['map']@q['S']@ti
            expected=q['expected']+.37*select.T@pair@select
            actual=q['N']+.37*pull
            congruence=old.err(ti.T@q['S'].T@actual@q['S']@ti-expected)
            assert congruence<3e-12
            idx=np.r_[np.arange(r,2*r),np.arange(3*r,4*r)]
            # Work in original c,e' variables; w=J_minus^dagger v c need not
            # be invertible, so do not introduce an inverse physical covariance.
            phys0=q['expected'][np.ix_(idx,idx)]
            pmass=(select.T@pair@select)[np.ix_(idx,idx)]
            phys=phys0+.37*pmass
            pf=old.old.old.pfaffian
            ratio=pf(actual)/pf(q['N'])
            target=pf(phys)/pf(phys0)
            ratio_error=float(abs(ratio-target));assert ratio_error<2e-9
            src=.5*np.trace(np.linalg.solve(actual,pull))
            wanted=.5*np.trace(np.linalg.solve(phys,pmass))
            source_error=float(abs(src-wanted));assert source_error<2e-10
            step=2e-5
            fd=(pf(q['N']+(.37+step)*pull)-pf(q['N']+(.37-step)*pull))/(2*step*pf(actual))
            derivative_error=float(abs(fd-src));assert derivative_error<2e-7
            correct_map=select[:,idx]
            g=correct_map@(-np.linalg.inv(phys))@correct_map.T
            mapped=f['map']@(-np.linalg.inv(actual))@f['map'].T
            obs_error=old.err(g-mapped);assert obs_error<3e-11
            wrong_ratio=pf(q['N']+.37*wrong)/pf(q['N'])
            mismatch=old.norm(pull-wrong)
            assert mismatch>.1 and abs(wrong_ratio-ratio)>1e-5
            gram=(f['Theta']@g)[np.ix_(f['positive'],f['positive'])]
            hermiticity=old.err(gram-gram.conj().T)
            mineig=float(min(np.linalg.eigvalsh((gram+gram.conj().T)/2)))
            # Diagnostic only: covariance positivity is not an all-polynomial proof.
            assert hermiticity<3e-12 and mineig>-3e-12
            rows.append(dict(nx=nx,nt=2,sample=sample,Nambu_dimension=len(actual),
                mass_scale=.37,reflection_pair_error=reflect_error,
                transformed_full_mass_congruence_error=congruence,
                original_Weyl_mass_weight_ratio=old.cpair(target),
                pulled_back_local_candidate_ratio=old.cpair(ratio),
                bare_local_mass_ratio=old.cpair(wrong_ratio),
                relative_weight_identity_error=ratio_error,
                common_mass_source=old.cpair(src),source_identity_error=source_error,
                source_finite_difference_error=derivative_error,
                physical_covariance_identity_error=obs_error,
                omitted_pairing_matrix_norm=mismatch,
                massive_reflection_Gram_minimum=mineig,
                reflection_Gram_hermiticity_error=hermiticity))
    for a,b in ((rows[0],rows[1]),(rows[2],rows[3])):
        assert np.max(abs(np.array(a['original_Weyl_mass_weight_ratio'])-b['original_Weyl_mass_weight_ratio']))<1e-11
    deps=('joint_physical_auxiliary_state.py','joint_local_mirror_process.py',
          'joint_spinor_subgroup_mass.py','joint_gauss_fermion_influence.py',
          'research_note_604.md','research_note_611.md','research_note_614.md',
          'research_note_660.md','research_note_661.md','research_note_667.md')
    return dict(date='2026-10-02',entry_for_round=668,formal_round_complete=False,
        rows=rows,original_mass_and_exterior_dictionary_used=True,
        Wick_quadratic_extension_not_new_CAR_to_Euclidean_equivalence=True,
        covariance_RP_diagnostic_not_full_RP_proof=True,
        no_general_gauge_or_continuum_or_quantum_GR_claim=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in deps},
        all_checks_passed=True)


if __name__=='__main__':
    result=run();target=HERE/'mass_observation_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_checks_passed')}))
