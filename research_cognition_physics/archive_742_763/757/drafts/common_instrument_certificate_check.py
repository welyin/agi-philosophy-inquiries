"""Check reused instrument identities; no new numbered scientific round."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_correlated_hadamard_noise as prior
import verify_interaction_rounds as core
TARGET=HERE/'common_instrument_certificate_checks.json'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def verify():
    history=dict(core.read(ROOT/'round584_drafts/historical_evidence_manifest.json')['evidence_hashes'])
    for n in range(584,757):
        receipt=core.read(ROOT/f'research_round_{n}_checks.json')
        for key in ('new_file_hashes','preserved_draft_hashes'):
            for name,sha in receipt[key].items():
                assert name not in history or history[name]==sha
                history[name]=sha
    history.update(core.read(ROOT/'cognitive_foundation_bridge_605_navigation.json')['supplementary_artifact_hashes'])
    assert len(history)==3574
    for name,sha in history.items():assert digest(ROOT/name)==sha,name
    c=prior.car(2);n=[a.conj().T@a for a in c]
    pair=c[0].conj().T@c[1].conj().T;hop=c[0].conj().T@c[1]
    basis=[*n,pair+pair.conj().T,1j*(pair-pair.conj().T),hop+hop.conj().T]
    complementary=np.diag([.17,.83]);projectors=[np.diag(np.eye(4)[r]) for r in range(4)]
    err=0.;rows=[]
    saved=core.read(prior.TARGET)
    eta=complex(*saved['common_curved_background']['integrated_scalar_pair_coefficient'])
    S=eta*pair+eta.conjugate()*pair.conj().T
    branch_mean=[prior.mean(P,S).real for P in projectors]
    branch_var=[prior.covariance(P,S,S) for P in projectors]
    assert max(abs(x) for x in branch_mean)<1e-15
    assert np.max(np.abs(np.array(branch_var)-abs(eta)**2*np.array([1,0,0,1])))<1e-14
    for p,q in ((.2,0),(.2,.04),(.2,.2),(.5,0),(.5,.25),(.5,.5),(.8,.6),(.8,.64),(.8,.8)):
        rho=prior.tau(p,q);pi=np.diag(rho).real;whole=np.kron(rho,complementary)
        combined=np.zeros_like(whole)
        for r,P in enumerate(projectors):
            L=np.kron(P,np.eye(2));out=L@whole@L
            expected=np.kron(P@rho@P,complementary)
            err=max(err,float(np.linalg.norm(out-expected)),float(abs(np.trace(out)-pi[r])))
            combined+=out
        err=max(err,float(np.linalg.norm(combined-whole)))
        for A in basis:
            for B in basis:
                am=np.array([prior.mean(P,A) for P in projectors]);bm=np.array([prior.mean(P,B) for P in projectors])
                within=sum(pi[r]*prior.covariance(P,A,B) for r,P in enumerate(projectors))
                between=sum(pi[r]*(am[r]-pi@am)*(bm[r]-pi@bm) for r in range(4)).real
                err=max(err,abs(prior.covariance(rho,A,B)-within-between))
        v=prior.covariance(rho,S,S)
        err=max(err,abs(v-(1-2*p+2*q)*abs(eta)**2))
        rows.append(dict(p=p,q=q,finite_pair_source_variance=float(v),record_mean_variance=0.))
    assert err<2e-14
    names=('research_note_634.md','research_note_730.md','research_note_743.md','research_note_750.md','research_note_756.md',
           'joint_correlated_hadamard_noise.py','joint_correlated_hadamard_noise_results.json')
    links=0
    for link in core.link_parser()((HERE/'common_instrument_certificate.md').read_text('utf8')):
        target=(HERE/link).resolve()
        assert target.exists() or target==TARGET.resolve(),link
        links+=1
    return dict(date='2026-10-04',type='unnumbered_common_object_review',latest_scientific_round=756,
        cumulative_numbered_tests=3455,numbered_scientific_files=1580,
        protected_historical_files_checked=len(history),all_historical_hashes_unchanged=True,
        reused_CP_and_covariance_checks_passed=True,identity_error=float(err),finite_pair_rows=rows,
        dependency_hashes={name:digest(ROOT/name) for name in names},
        working_file_hashes={name:digest(HERE/name) for name in ('common_instrument_certificate.md',Path(__file__).name)},
        local_links_checked=links,broken_links=0,no_new_scientific_round=True,
        autonomous_instrument_or_full_graph_mapping_proven=False,active_goal_unchanged=True,all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-checks',action='store_true');a=p.parse_args();r=verify()
    if a.write_checks:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==core.read(TARGET)
    print(json.dumps({k:r[k] for k in ('latest_scientific_round','no_new_scientific_round','all_checks_passed')}))
