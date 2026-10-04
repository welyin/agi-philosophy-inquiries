"""721 entry: reuse652 reference flows; test actual bounded velocity instruments."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_quantum_reference_forms as old
TARGET=HERE/'material_velocity_entry_results.json'


def run():
    d=old.diagnostic();H=np.diag(d['E']);rho=np.diag(d['prob'])
    rows=[];t=.13
    for label,V in zip(('T','s'),d['V']):
        e,u=np.linalg.eigh(V);half=(u*np.exp(-.5j*t*e/old.HBAR))@u.conj().T
        sin=(u*np.sin(t*e/old.HBAR))@u.conj().T
        K=[(half+1j*half.conj().T)/2,(half-1j*half.conj().T)/2]
        states=[k@rho@k.conj().T for k in K];post=sum(states)
        expected=(half@rho@half.conj().T+half.conj().T@rho@half)/2
        errs=[np.linalg.norm(sum(k.conj().T@k for k in K)-np.eye(len(H))),
              np.linalg.norm(K[0].conj().T@K[0]-(np.eye(len(H))-sin)/2),
              np.linalg.norm(post-expected)]
        delta=np.trace(H@(post-rho)).real
        assert max(errs)<2e-12 and delta>0
        rows.append(dict(reference=label,dimension=len(H),parameter=t,
                         max_identity_error=float(max(errs)),thermal_energy_change=float(delta),
                         probabilities=[float(np.trace(s).real) for s in states]))
    rng=np.random.default_rng(721);errors=[]
    for _ in range(30):
        p=rng.normal(size=5);p*=np.sqrt(6*old.M)*rng.uniform(.1,.98)/np.linalg.norm(p)
        F=old.original.F(p)
        for label in ('T','s'):
            X=old.smooth_data(p,label)[0]
            got=-p@X/(3*F)
            expected=-F*(np.sum(p[:4]**2) if label=='T' else p[4])/(3*old.M)
            bound=old.M/2 if label=='T' else 2*np.sqrt(2*old.M)/9
            errors.append(abs(got-expected))
            assert abs(got)<=bound+1e-14
    assert max(errors)<1e-14
    names=('research_note_652.md','joint_quantum_reference_forms.py',
           'joint_quantum_reference_forms_results.json','research_note_720.md')
    return dict(entry_round=721,formal_round=False,rows=rows,
        logF_directional_identity_error=max(errors),
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names},
        scope='652 reused: legal bounded spectral instruments and conditional64D diagnostic. Full original energy-form preservation and continuum localization NOT proved.')


if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:
        f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(r,ensure_ascii=False,indent=2))

