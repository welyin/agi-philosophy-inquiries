"""Executed665 entry: original mass, shared contact and positive time slices.

Uses the original eight-lepton conditional sector. The charge Gaussian and
spin-Casimir factorization are checked independently of exponentiating Q.
This is not yet the whole-graph auxiliary measure or a completed round665.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
sys.path.insert(0,str(BASE))
import joint_connection_matter_matching as previous
import joint_gauss_fermion_influence as car
import joint_fermion_gauss_completion as matter


def exp_sym(h,a):
    ev,vec=np.linalg.eigh(h)
    return (vec*np.exp(a*ev))@vec.conj().T


def norm(x):return float(np.linalg.norm(x,'fro'))


def run():
    mats=previous.current_matrices()[:,24:,24:]
    q,raw,js=previous.contact_operator(mats)
    spin=[j/2 for j in js[1:]]
    s2=sum(s@s for s in spin)
    number=car.fock(np.eye(8),np.zeros((8,8)))
    chi=js[0]
    casimir_error=norm(q-(4*s2-chi@chi-2*number))
    h,d=matter.mass_matrices(car.PHI[0]);h0=car.fock(h[24:,24:],d[24:,24:])
    spin_comm=max(norm(s@h0-h0@s) for s in spin)
    chi_comm=norm(chi@h0-h0@chi)
    number_comm=norm(number@h0-h0@number)
    beta=.7;theta=.13;k=3/(16*np.exp(6*theta));v=k*q
    contact_comm=norm(v@h0-h0@v)
    assert casimir_error<1e-12 and spin_comm<1e-12
    assert min(chi_comm,number_comm,contact_comm)>1e-3
    # Normalized real Gaussian charge factor. Spectral diagonal Q_chi gives a
    # separate quadrature check, rather than using the contact eigendecomposition.
    a=beta*k
    z,w=np.polynomial.hermite.hermgauss(24)
    qchi=np.diag(chi).real;n=np.diag(number).real
    charge=sum(weight*np.exp(2*np.sqrt(a)*point*qchi) for point,weight in zip(z,w))/np.sqrt(np.pi)
    factor=np.diag(np.exp(2*a*n)*charge)@exp_sym(s2,-4*a)
    direct=exp_sym(v,-beta)
    factor_error=norm(factor-direct)
    assert factor_error<3e-11
    exact=exp_sym(h0+v,-beta);z_exact=np.trace(exact).real
    source_exact=float(6*beta*np.trace(v@exact).real/z_exact)
    rows=[]
    for steps in (1,2,4,8,16):
        dt=beta/steps;half=exp_sym(h0,-dt/2);contact=exp_sym(v,-dt)
        transfer=half@contact@half
        dt_theta=half@(6*dt*v@contact)@half
        ev,vec=np.linalg.eigh(transfer)
        assert ev.min()>0
        powers=ev**steps;zn=float(powers.sum())
        source=float(steps*np.sum(ev**(steps-1)*np.diag(vec.conj().T@dt_theta@vec).real)/zn)
        total=(vec*powers)@vec.conj().T
        err=norm(total-exact)
        rows.append(dict(steps=steps,positive_step_minimum=float(ev.min()),
                         matrix_error=err,logZ=float(np.log(zn)),
                         source=source,source_error=abs(source-source_exact)))
    for arow,brow in zip(rows,rows[1:]):
        assert 3.5<arow['matrix_error']/brow['matrix_error']<4.5
        assert 3.5<arow['source_error']/brow['source_error']<4.5
    deps=('joint_connection_matter_matching.py','joint_fermion_gauss_completion.py',
          'joint_gauss_fermion_influence.py','research_note_603.md','research_note_643.md',
          'research_note_655.md','research_note_664.md')
    return dict(date='2026-10-02',entry_for_round=665,formal_round_complete=False,
                original_conditional_lepton_modes=8,casimir_identity_error=casimir_error,
                spin_mass_commutator=spin_comm,chirality_mass_commutator=chi_comm,
                number_mass_commutator=number_comm,contact_mass_commutator=contact_comm,
                charge_gaussian_points=24,independent_contact_factor_error=factor_error,
                exact_logZ=float(np.log(z_exact)),exact_source=source_exact,rows=rows,
                dependency_hashes={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in deps},
                scope='Actual original eight-lepton onsite sector; total charge Gaussian and spin Casimir recover the same contact factor, but original mass needs ordered positive time slices. No whole-graph Gauss integral, sign-free path weight, continuum torsion determinant or completed quantum-gravity matching.',
                all_checks_passed=True)


if __name__=='__main__':
    result=run();target=HERE/'spin_charge_auxiliary_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_checks_passed')}))
