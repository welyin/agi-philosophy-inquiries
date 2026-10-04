"""Reproduce684 entry, verify explicit analytic witness, preserve all683 evidence."""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import verify_interaction_rounds as core
spec=importlib.util.spec_from_file_location('entry684',HERE/'compensation_involution_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
result=probe.run();assert result==core.read(probe.TARGET)
# Independently evaluate the predeclared analytic negative vector (no eigenvector fitting).
base=probe.base;body=probe.body
_,e,phis=body.prior.fixture();mat=base.fixed_matrices(e,phis)
links=np.tile(np.eye(16,dtype=complex),(2,4,1,1));_,_,_,h,_=base.kernel(links)
g5=mat[1]@mat[1].T-mat[0]@mat[0].T
ell,_,_=probe.reflection.reflection_matrices(mat,np.eye(4)[[2,3,0,1]])
R=np.kron(np.array([[0.,1.],[1.,0.]]),g5@ell)
a=.2;k,_,_=body.blocks(g5@h,g5,a,1);C=np.linalg.inv(k.conj().T@k)
pos=np.r_[np.arange(128,256),np.arange(128,256)+256]
B=(R@C)[np.ix_(pos,pos)]
gram=np.block([[np.zeros_like(B),B],[-B,np.zeros_like(B)]])
u=np.zeros(128);u[[0,64]]=1/np.sqrt(2)
v=np.zeros(512);v[:128]=u/np.sqrt(2);v[384:]=-u/np.sqrt(2)
value=np.vdot(v,gram@v);expected=-a/(4+a*a)**2
assert abs(value-expected)<1e-14 and expected<0
history=dict(core.read(ROOT/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
for n in range(584,684):
    receipt=core.read(ROOT/f'research_round_{n}_checks.json')
    for key in ('new_file_hashes','preserved_draft_hashes'):
        for name,digest in receipt[key].items():
            assert name not in history or history[name]==digest
            history[name]=digest
history.update(core.read(ROOT/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
assert len(history)==2555
for name,digest in history.items():assert core.digest(ROOT/name)==digest,name
names=('compensation_involution_probe.py','compensation_involution_probe_results.json',
       'compensation_involution_entry.md')
receipt=dict(entry_round=684,latest_formal_round=683,entry_reproduced=True,
    analytic_predeclared_witness=base.old.cpair(value),expected_value=expected,
    protected_evidence_files=len(history),numbered_tests_not_increased=True,
    evidence_hashes={p:core.digest(HERE/p) for p in names},all_checks_passed=True)
target=HERE/'entry_checks.json'
if target.exists():assert receipt==core.read(target)
else:
    with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt))
