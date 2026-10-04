"""710 entry: inherited index kernel modulo the actual colour centre.
Direct joint scope consequence of629/709; not a new formal research round.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_charge_quantum_coarse as old
TARGET=HERE/'residual_phase_entry_results.json'


def run():
    omega=np.exp(2j*np.pi/3)
    centre=old.matter.representation(omega*np.eye(3),np.eye(2),1.+0j)
    phase=np.diag(np.exp(2j*np.pi*old.QUARK/3))
    centre_error=float(np.linalg.norm(centre-phase));assert centre_error<1e-13
    # Original even Gauss mode witness udd nu_R, not a continuum instanton vertex.
    base=old.normalized(old.create(old.baryon_create({0:1.},0),30))
    assert len(base)==3 and all(m.bit_count()==4 for m in base)
    rng=np.random.default_rng(710)
    C=old.matter.gauge.group_exp(rng.normal(size=8),3)
    W=old.matter.gauge.group_exp(rng.normal(size=3),2)
    R=old.matter.representation(C,W,np.exp(.173j))
    inputs=[([i for i in range(32) if (m>>i)&1],a) for m,a in base.items()]
    ge=0.
    for u,dd in itertools.product((12,14,16),itertools.combinations((18,20,22),2)):
        out=tuple(sorted((u,*dd,30)));mask=sum(1<<i for i in out)
        amplitude=sum(a*np.linalg.det(R[np.ix_(out,ins)]) for ins,a in inputs)
        ge=max(ge,float(abs(amplitude-base.get(mask,0))))
    assert ge<1e-13
    rows=[]
    for ng in (1,3):
        state={0:1.+0j}
        for generation in range(ng):
            state={m|(b<<(32*generation)):a*c for m,a in state.items() for b,c in base.items()}
        assert abs(old.inner(state,state)-1)<1e-13
        charge=3*ng
        assert all(m.bit_count()==4*ng for m in state)
        vals=np.array([0,3,6,9])
        phases=np.exp(2j*np.pi*vals/(3*ng))
        rho=np.ones((4,4),complex)/4
        residual=sum((phases**k)[:,None]*rho*(phases.conj()**k)[None,:] for k in range(ng))/ng
        expected=rho*((vals[:,None]-vals[None,:])%(3*ng)==0)
        assert np.linalg.norm(residual-expected)<1e-13
        if ng==1:assert np.linalg.norm(residual-rho)<1e-13
        # Coherence between vacuum and the explicit n_g-fold Gauss witness.
        kept=np.exp(2j*np.pi*charge/(3*ng));assert abs(kept-1)<1e-13
        rows.append(dict(generations=ng,phase_kernel_order=3*ng,physical_quotient_order=ng,
            even_witness_CAR_particles=4*ng,even_witness_quark_charge=charge,
            even_witness_sparse_components=len(state),phase_error=float(abs(kept-1)),
            pinching_error=float(np.linalg.norm(residual-expected)),
            vacuum_charge9_coherence_real=float(residual[0,3].real),
            vacuum_charge3_coherence_real=float(residual[0,1].real)))
    deps=('research_note_629.md','research_note_709.md','joint_charge_quantum_coarse.py',
          'joint_charge_quantum_coarse_results.json')
    return dict(date='2026-10-03',entry_round=710,latest_formal_round=709,new_formal_round=False,
        original_gauge_centre_identity_error=centre_error,even_mode_witness_gauge_error=ge,rows=rows,
        no_new_generation_selection=True,no_full_discrete_gauging_anomaly_claim=True,
        no_continuum_transition_or_autonomous_operation_constructed=True,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r))
