"""Executed666 entry: transport original contact AND vacuum ordering.

Uses the actual614 Bogoliubov map and Spin(10) carrier. Eight-lepton Fock
checks are conditional; the full internal commutator uses all16 species.
No enlarged gauge symmetry or completed chiral measure is asserted.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_spinor_subgroup_mass as dictionary
import joint_connection_matter_matching as contact
import joint_gauss_fermion_influence as car
import joint_fermion_gauss_completion as matter


def norm(x):return float(np.linalg.norm(x,'fro'))
def bdg(h):return np.block([[h,np.zeros_like(h)],[np.zeros_like(h),-h.T]])


def run():
    w=dictionary.ph_matrix();currents=contact.current_matrices();chi=currents[0]
    lam=-chi
    errors=[];constants=[]
    targets=[lam,-np.eye(32),*list(currents[1:])]
    for c,target in zip([np.eye(32),*list(currents)],targets):
        new=w@bdg(c)@w.conj().T
        errors.extend([norm(new[:32,:32]-target),norm(new[:32,32:])])
        constants.append(float((np.trace(c)-np.trace(new[:32,:32])).real/2))
    assert max(errors)<1e-12 and constants==[16.,16.,0.,0.,0.]
    j=dictionary.dictionary();li=j@lam[::2,::2]@j.conj().T
    generators,_=dictionary.clifford()
    comms=[g@li-li@g for g in generators]
    realmap=np.column_stack([np.r_[c.real.ravel(),c.imag.ravel()] for c in comms])
    rank=int(np.linalg.matrix_rank(realmap,tol=1e-12))
    assert rank>0
    gauge=dictionary.exterior(dictionary.carrier(
        matter.gauge.group_exp(np.array([.12,-.14,.08,.04,.09,-.06,.05,.11]),3),
        matter.gauge.group_exp(np.array([.17,-.13,.09]),2),np.exp(.21j)))
    subgroup_error=norm(gauge@li-li@gauge)
    assert subgroup_error<1e-12
    # Actual original eight-lepton sector, four right-handed modes transformed.
    lm=currents[:,24:,24:];qold,_,js=contact.contact_operator(lm)
    size=256;identity=np.eye(size);zero=np.zeros((8,8),complex)
    number=car.fock(np.eye(8),zero);relative=car.fock(lam[24:,24:],zero)
    spin2=sum((x/2)@(x/2) for x in js[1:]);r=4
    transported=4*spin2-(r*identity-number)@(r*identity-number)-2*(r*identity+relative)
    reset=4*spin2-number@number-2*number
    mismatch=(2*r+2)*number-2*relative-(r*r+2*r)*identity
    assert norm(transported-reset-mismatch)<1e-12
    phi=car.PHI[0];h,d=matter.mass_matrices(phi)
    oldh=car.fock(h[24:,24:],d[24:,24:])
    mass_new=w@car.bdg(h,d)@w.conj().T
    newh=car.fock(mass_new[24:32,24:32],mass_new[24:32,56:64])
    beta=.7;k=3/(16*np.exp(.78))
    ev_old=np.linalg.eigvalsh(oldh+k*qold)
    ev_new=np.linalg.eigvalsh(newh+k*transported)
    spectrum_error=float(np.max(abs(ev_old-ev_new)))
    assert spectrum_error<1e-11
    logz=lambda ev:float(np.log(np.exp(-beta*ev).sum()))
    reset_z=logz(np.linalg.eigvalsh(newh+k*reset))
    correct_z=logz(ev_new)
    assert abs(reset_z-correct_z)>.01
    # The original vacuum is the right-filled b state, not the new empty state.
    occupied=sum(1<<i for i in range(4,8))
    assert abs(transported[occupied,occupied])<1e-12
    assert abs(transported[0,0]+24)<1e-12
    deps=('joint_spinor_subgroup_mass.py','joint_connection_matter_matching.py',
          'joint_fermion_gauss_completion.py','research_note_614.md','research_note_653.md','research_note_665.md')
    return dict(date='2026-10-02',entry_for_round=666,formal_round_complete=False,
                full32_bilinear_dictionary_error=max(errors),ordered_constants=constants,
                full_Spin10_generators=len(generators),anisotropy_commutator_rank=rank,
                preserving_Lie_algebra_dimension=45-rank,
                maximal_anisotropy_commutator_norm=max(norm(c) for c in comms),
                original_SM_subgroup_error=subgroup_error,
                conditional_lepton_modes=8,transported_contact_spectrum_error=spectrum_error,
                original_logZ=logz(ev_old),correctly_transported_logZ=correct_z,
                vacuum_reset_logZ=reset_z,operator_mismatch_norm=norm(mismatch),
                right_filled_original_vacuum_contact=float(transported[occupied,occupied].real),
                new_empty_vacuum_contact=float(transported[0,0].real),
                dependency_hashes={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in deps},
                scope='Actual614 particle-hole dictionary transports the finite contact with a reference-dependent bilinear term. Re-normal-ordering in a new vacuum changes the model. Original SM preserved, entire Spin(10) is not assumed; no continuum anomaly or chiral auxiliary measure equivalence.',
                all_checks_passed=True)


if __name__=='__main__':
    result=run();target=HERE/'particle_hole_contact_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_checks_passed')}))
