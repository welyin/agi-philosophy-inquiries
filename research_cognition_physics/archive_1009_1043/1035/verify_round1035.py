"""1035 delivery and independent-algorithm checks, not independent review."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
from urllib.parse import unquote
import numpy as np
import joint_selection_audit as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
NOTE=HERE.parent/'research_note_1035.md'
OUT=HERE/'research_round_1035_checks.json'
OWN=['joint_selection_audit.py','joint_selection_audit_results.json',
     'bridge_ledger_v0_2.json','drafts/joint_selection_derivation.md',
     'selection_audit.md','input_dependency_update_v0_24.md','review.md',
     'NEXT.md','verify_round1035.py']
LINK=re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)')


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def matrix_trace_checks():
    su3=[]
    for i,j in ((0,1),(0,2),(1,2)):
        a=np.zeros((3,3),complex);a[i,j]=a[j,i]=.5
        b=np.zeros((3,3),complex);b[i,j]=-.5j;b[j,i]=.5j
        su3.extend([a,b])
    su3.extend([np.diag([1.,-1.,0.])/2,np.diag([1.,1.,-2.])/(2*np.sqrt(3))])
    su2=[np.array([[0,1],[1,0]])/2,np.array([[0,-1j],[1j,0]])/2,np.diag([1.,-1.])/2]
    identity3=np.eye(3);identity2=np.eye(2)
    worst=0.;evaluations=0
    for k in (-2,0,1,3):
        # Separate explicit representation matrices; no science charge helper.
        hyper=[1/6+k,-2/3-k,1/3-k,-1/2-3*k,1+3*k,3*k]
        dims=[6,3,3,2,1,1]
        traces=[0j]*4
        for index,(y,d) in enumerate(zip(hyper,dims)):
            Y=y*np.eye(d)
            traces[0]+=np.trace(Y);traces[1]+=np.trace(Y@Y@Y)
        for a in su3:
            T=[np.kron(a,identity2),-a.conj(),-a.conj(),np.zeros((2,2)),np.zeros((1,1)),np.zeros((1,1))]
            value=sum(np.trace(y*(t@t)) for y,t in zip(hyper,T))
            worst=max(worst,float(abs(value)));evaluations+=1
        for a in su2:
            T=[np.kron(identity3,a),np.zeros((3,3)),np.zeros((3,3)),a,np.zeros((1,1)),np.zeros((1,1))]
            value=sum(np.trace(y*(t@t)) for y,t in zip(hyper,T))
            worst=max(worst,float(abs(value)));evaluations+=1
        worst=max(worst,float(abs(traces[0])),float(abs(traces[1])))
        evaluations+=2
    assert worst<1e-10,worst
    # Actual color singlet in u^c conjugate tensor d^c.
    singlet=np.eye(3).ravel()/np.sqrt(3)
    tensor_error=max(float(np.linalg.norm((np.kron(a,identity3)-np.kron(identity3,a.conj()))@singlet)) for a in su3)
    assert tensor_error<1e-13
    return evaluations,worst,tensor_error


def local_links(s):
    s=re.sub(r'\$\$.*?\$\$|^~~~.*?^~~~\s*$','',s,flags=re.S|re.M)
    for m in LINK.finditer(s):
        v=m[1].strip();v=v[1:v.index('>')] if v.startswith('<') else v
        path=unquote(v.split('#',1)[0])
        if path and not re.match(r'^[a-zA-Z]+:',path):yield v,path


def verify(prospective=False):
    fresh=science.run();current=science.ledger()
    science.compare(fresh,json.loads(science.OUT.read_text('utf8')))
    science.compare(current,json.loads(science.LEDGER.read_text('utf8')))
    assert len(fresh['family'])==13 and fresh['representation_words']['cases']==2772
    assert all(v==['0']*4 for v in fresh['anomaly_polynomial_coefficients'].values())
    for row in fresh['family']:
        k=row['k'];fields={r['field']:r for r in row['charges']}
        assert fields['N']['q_6Y']==18*k
        assert fields['dc']['q_6Y']-fields['uc']['q_6Y']==6
        assert row['real_S_Majorana_charge_6Y']==36*k
        assert len(row['central_kernel'])==6
        assert all(f['Z4']==1 for f in fields.values())
    assert len(fresh['mass_cases'])==6
    for row in fresh['mass_cases']:
        s=np.array(row['singular_values'])
        assert np.max(np.abs(s[::2]-s[1::2]))<1e-13
        assert len(set(np.round(s,12)))==6
        assert row['Majorana_Ward_defect']==102*abs(row['k'])
    assert len(current['rows'])==27 and len({r['id'] for r in current['rows']})==27
    before=json.loads((BASE/science.OLD).read_text('utf8'))
    old={r['id']:r for r in before['rows']}
    changed=[]
    for row in current['rows']:
        if row['id'] in old and row['status']!=old[row['id']]['status']:
            changed.append(row['id'])
            assert row['previous_status']=='missing_bridge' and row['status']=='mapped_interface'
        assert not row['cognitive_origin_closed'] and not row['full_parent_process_certified']
        assert all((BASE/name).is_file() for name in row['source_paths'])
    assert set(changed)=={'scalar_source_to_parent','gravity_ideal_to_parent'}
    checks,trace_error,tensor_error=matrix_trace_checks()
    assets=[NOTE]+[HERE/name for name in OWN];links=0
    for path in assets:
        assert path.is_file(),path
        if path.suffix=='.py':ast.parse(path.read_text('utf8'))
        if path.suffix=='.md':
            s=path.read_text('utf8');assert s.count('$$')%2==0,path
            for raw,local in local_links(s):
                target=(path.parent/local).resolve()
                assert target.exists() or (prospective and target==OUT),(path,raw)
                links+=1
    assert all(f'## {i}.' in NOTE.read_text('utf8') for i in range(1,11))
    for name,digest in fresh['historical_source_sha256'].items():assert sha(BASE/name)==digest
    assert '独立代理签审未完成' in (HERE/'review.md').read_text('utf8')
    assert not fresh['full_mixed_anomalies_certified'] and not fresh['goal_completed']
    return dict(round=1035,date='2026-10-08',local_delivery_checks_passed=True,
                scientific_result_reproduced=True,external_independent_agent_review_completed=False,
                new_calibration_groups=1,cumulative_research_groups=3812,new_adopted_cognitive_axioms=0,
                joint_family_cases=13,representation_word_cases=2772,mass_cases=6,
                independent_matrix_trace_checks=checks,matrix_trace_residual=trace_error,
                independent_tensor_singlet_generators=8,tensor_singlet_residual=tensor_error,
                bridge_rows=27,old_limited_gaps_now_mapped=changed,
                full_mixed_global_anomaly_classification_completed=False,
                physical_implementation_certified=False,all_observed_physics_matched=False,
                goal_completed=False,visual_checks_performed=False,
                frozen_current_files=len(assets),historical_input_files=len(science.HISTORY),
                local_links_checked=links,live_navigation_frozen=False,
                historical_source_sha256=fresh['historical_source_sha256'],
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in assets})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-receipt',action='store_true')
    args=parser.parse_args();result=verify(prospective=args.write_receipt)
    if args.write_receipt:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:science.compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},ensure_ascii=False))
