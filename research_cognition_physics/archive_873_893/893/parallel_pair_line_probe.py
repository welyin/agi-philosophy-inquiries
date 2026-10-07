"""Supplement to893: original massive Dirac parallel neutral-pair line."""
from pathlib import Path
import json,sys
import numpy as np
sys.dont_write_bytecode=True
import curved_connection_heat_bridge as main
prior=main.prior;old=main.old
TARGET=Path(__file__).with_name('parallel_pair_line_results.json')
def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.prior.prior.last.em_charges(m,matter)
    ids=[26,27,28,29]
    mass=np.exp(old.XI)*m.MASS[np.ix_(ids,ids)]
    mm=np.sqrt(np.trace(mass@mass).real/4)
    gamma=[g[np.ix_(ids,ids)] for g in m.GAMMA]
    J=gamma[0]@gamma[1]@gamma[2]@mass/mm
    cs=[prior.prior.prior.prior.prior.annihilator(j) for j in range(4)]
    _,V0=np.linalg.eigh(mass)
    def lift(X):return sum(X[i,j]*cs[i].conj().T@cs[j] for i in range(4) for j in range(4))
    points=[np.array(a) for a in ([0.,0.,0.],[.012,-.009,.016],[-.083,.01,-.01])]
    states=[];matrices=[];checks=[]
    for a in points:
        H,P,dp,A=prior.jets(m,q,a,np.zeros(3))
        h=H[np.ix_(ids,ids)];p=P[np.ix_(ids,ids)];E=np.sqrt(np.trace(h@h).real/4)
        W=prior.prior.exterior(prior.segment(h,mass)@V0);vac=W[:,3]
        pair=lift((np.eye(4)-p)@J@p)@vac/np.sqrt(2)
        B=lift(h)+2*E*np.eye(16)
        checks.append(dict(alpha=a.tolist(),pair_norm=float(np.linalg.norm(pair)),
            excitation_identity_error=float(np.linalg.norm(B@pair-2*E*pair)),
            chiral_anticommutator_error=float(np.max(abs(J@h+h@J))),
            horizontal_commutator_error=float(max(np.max(abs(J@z[np.ix_(ids,ids)]-z[np.ix_(ids,ids)]@J)) for z in A))))
        states.append(pair);matrices.append(h)
    errors=[]
    # Radial straight segments through arbitrary original momenta likewise
    # lie in a fixed Clifford plane and preserve this parallel line.
    for i in (1,2):
        U=prior.prior.exterior(prior.segment(matrices[i],matrices[0]))
        errors.append(float(np.linalg.norm(U@states[0]-states[i])))
    assert max(errors)<1e-11 and max(x['horizontal_commutator_error'] for x in checks)<1e-10
    assert all(abs(x['pair_norm']-1)<1e-12 and x['excitation_identity_error']<1e-11 for x in checks)
    return dict(round=893,supplement_not_new_group=True,original_parallel_pair_line_checks=checks,
        parallel_pair_transport_errors=errors,
        scalar_ground_equality_has_structural_explanation=True,
        full_pair_connection_is_flat=False)
if __name__=='__main__':
    r=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
