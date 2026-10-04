"""723 entry: same four marginal laws need not give the same actual history.

This is an inherited finite radial diagnostic, not a full-model no-go proof.
The algebraic real/conjugate mechanism is elementary and not a new round.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_squared_reference_readout as old
TARGET=HERE/'joint_reference_entry_results.json'


def run():
    t=old.fixture(.017,'T');s=old.fixture(.017,'s')
    I=t['I'];n=len(I);T=t['f']
    e,u=np.linalg.eigh(T);F=(u*(.5+.25*np.sin(e)))@u.T
    rootF=(u*np.sqrt(.5+.25*np.sin(e)))@u.T
    K=t['Ks'][0]
    effect=K.conj().T@F@K
    B=(effect-effect.conj())/2;bnorm=float(np.linalg.norm(B,2))
    assert bnorm>1e-4
    rho_plus=(I+B/(2*bnorm))/n;rho_minus=(I-B/(2*bnorm))/n
    assert np.linalg.eigvalsh(rho_plus).min()>.49/n
    assert np.linalg.eigvalsh(rho_minus).min()>.49/n
    marginals=[]
    for label,A in (('T',T),('s',s['f']),('Q_T',t['b']*I-t['R']),('Q_s',s['b']*I-s['R'])):
        assert np.max(abs(A.imag))<1e-13
        _,basis=np.linalg.eigh(A.real)
        diff=np.diag(basis.T@(rho_plus-rho_minus)@basis)
        error=float(np.max(abs(diff)));assert error<1e-13
        marginals.append(dict(observable=label,all_rank1_spectral_probabilities_error=error))
    probabilities=[float(np.trace(rho@effect).real) for rho in (rho_plus,rho_minus)]
    predicted=float(np.trace(B@B).real/(n*bnorm))
    difference=probabilities[0]-probabilities[1]
    assert abs(difference-predicted)<1e-14 and difference>1e-4
    energy_diff=float(abs(np.trace((rho_plus-rho_minus)@t['H'])))
    assert energy_diff<1e-13
    # Same marginal effects with another Kraus phase are another instrument.
    E0=K.conj().T@K;e,u=np.linalg.eigh(E0.real)
    L=(u*np.sqrt(e))@u.T
    luders=L@F@L
    reversed_effect=rootF@E0@rootF
    erased=float(abs(np.trace((rho_plus-rho_minus)@luders)))
    reversed_diff=float(abs(np.trace((rho_plus-rho_minus)@reversed_effect)))
    assert erased<1e-13 and reversed_diff<1e-13
    names=('joint_squared_reference_readout.py','research_note_722.md','research_note_554.md',
           'research_note_652.md','research_note_704.md','research_note_707.md','research_note_708.md')
    return dict(entry_round=723,formal_round=False,dimension=n,marginals=marginals,
                same_original_diagnostic_mean_energy_error=energy_diff,
                actual_two_record_probabilities=probabilities,history_probability_difference=difference,
                analytic_matrix_difference=predicted,erasing_Kraus_phase_difference=erased,
                reversed_order_difference=reversed_diff,
                dependencies={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in names},
                scope='Exact real/conjugate matrix mechanism in inherited64D radial diagnostic with conditioned neighbour W. Equal full one-observable spectral laws of T,s,Q_T,Q_s and equal mean H do not fix the actual ordered record. Full original Gauss-state witness and common process transport not completed by this entry.')


if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(r,ensure_ascii=False,indent=2))
