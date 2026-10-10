"""Read-only 1051 replay, independent author algorithm, links and asset hashes.

--write creates the author receipt once. This is not external peer review.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from urllib.parse import unquote
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1051_checks.json'
OWN=('../research_note_1051.md','proof.md','effective_lie_interface.py',
     'effective_lie_interface_results.json','input_dependency_update.md',
     'NEXT.md','review.md','verify_round1051.py',
     '../_admission/effective_lie/selection.md')


def dense_second_algorithm(result):
    # Different storage and twirl algorithm: one H x R x full-record matrix,
    # and the exact finite sign group, rather than the author's cq branch dict.
    d,r,m=6,2,2
    def rotation(j,theta):
        out=np.eye(d,dtype=complex)
        out[np.ix_([0,j],[0,j])]=[[math.cos(theta),-math.sin(theta)],
                                  [math.sin(theta),math.cos(theta)]]
        return out
    def average(rho,c):
        out=np.zeros_like(rho)
        for signs in product((-1,1),repeat=d-m-1):
            s=np.array([1]*(m+1)+list(signs))
            diag=np.repeat(s,r*c)
            out+=diag[:,None]*rho*diag[None,:]
        return out/(2**(d-m-1))
    def distance(a,b):
        return float(np.linalg.norm(a-b,ord='nuc')/2)
    def step(rho,stage,c):
        if stage in (0,2):
            op=np.zeros((d*r*c,d*r*c),complex)
            for old in range(c):
                angle=.002 if stage==0 else (-.0015 if old%2==0 else .0025)
                label=np.zeros((c,c));label[old,old]=1
                op+=np.kron(rotation(3,angle),np.kron(np.eye(r),label))
            return op@rho@op.conj().T,c
        probs=.35+.2*np.arange(d)/(d-1)
        answer=np.zeros((d*r*2*c,d*r*2*c),complex)
        for outcome in (0,1):
            label=np.zeros((2*c,c))
            for old in range(c):label[2*old+outcome,old]=1
            diagonal=np.diag(np.sqrt(probs if outcome==0 else 1-probs))
            operation=rotation(1,.001 if outcome==0 else -.001)@diagonal
            a=np.kron(operation,np.kron(np.eye(r),label))
            answer+=a@rho@a.conj().T
        return answer,2*c
    psi=np.zeros(d*r,complex);psi[0]=math.sqrt(.999);psi[3]=math.sqrt(.001)
    actual=np.outer(psi,psi.conj());effective=average(actual,1);c=1
    costs=[];losses=[]
    for stage in range(5):
        k=np.kron(np.diag([0]+[2**j for j in range(1,d)]),np.eye(r*c))
        costs.append(float(np.trace(k@actual).real))
        losses.append(distance(actual,average(actual,c)))
        if stage<4:
            next_actual,next_c=step(actual,stage,c)
            next_effective,check_c=step(effective,stage,c)
            assert next_c==check_c
            actual=next_actual;effective=average(next_effective,next_c);c=next_c
    history=result['finite_history']
    final=distance(actual,effective)
    error=max(abs(final-history['full_final_joint_distance']),
        max(abs(x-y) for x,y in zip(costs,history['true_prefix_budgets'])),
        max(abs(x-y) for x,y in zip(losses,history['local_joint_losses'])))
    assert error<1e-11 and c==4
    assert abs(np.trace(actual)-1)<1e-12
    assert np.linalg.eigvalsh(effective).min()>-1e-12

    exact_source=[]
    delta=F(1,100)
    for row in result['source_failure_truncations']:
        length=row['tail_levels']
        cost=2*delta**2*(1-delta**2)*F(length)/(1-F(1,2)**length)
        assert abs(float(cost)-row['output_budget'])<1e-12
        exact_source.append(str(cost))
    # For M=m+1, the same p in the stored reference witness is p_M.
    # This checks the coefficient; the proof, not this list, covers every quotient.
    dimension_certificates=[]
    for row in result['reference_witnesses']:
        M=row['m']+1
        p=F(1,2)**(M+1)
        squared=p*(1-p)
        assert abs(float(squared)-row['joint_distance']**2)<1e-12
        dimension_certificates.append(dict(head_dimension=M,budget='1/2',
            strict_error_threshold_squared=str(squared),
            conclusion='Every Lie quotient below this error has dimension at least M'))
    return dict(dense_record_dimension=4,finite_sign_group_order=8,
        dense_full_history_distance=final,maximum_difference_from_science=error,
        exact_source_truncation_costs=exact_source,
        any_quotient_dimension_bound_coefficients=dimension_certificates,
        additional_science_groups=0)


def local_links():
    count=0
    for name in OWN:
        path=(HERE/name).resolve()
        if path.suffix!='.md':continue
        for target in re.findall(r'(?<!_)\[[^\]]*\]\(([^)]+)\)',path.read_text(encoding='utf8')):
            target=target.strip().strip('<>')
            if re.match(r'[A-Za-z][A-Za-z0-9+.-]*://',target) or target.startswith('#'):continue
            full=(path.parent/unquote(target.split('#',1)[0])).resolve()
            if full!=RECEIPT:assert full.exists(),(name,target)
            count+=1
    return count


def build():
    science=HERE/'effective_lie_interface.py'
    replay=subprocess.run([sys.executable,'-B','-X','utf8',str(science)],
        check=True,capture_output=True,text=True,encoding='utf8')
    assert json.loads(replay.stdout)['all_scientific_checks_passed']
    result=json.loads((HERE/'effective_lie_interface_results.json').read_text(encoding='utf8'))
    for name,sha in result['historical_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    return dict(round=1051,date='2026-10-08',all_checks_passed=True,
        author_status='complete_pending_external_independent_final_review',
        external_independent_review_certified_here=False,
        new_science_groups=1,new_cognitive_axioms=0,roadmap_complete=False,
        science_read_only_compare_passed=True,
        author_second_algorithm=dense_second_algorithm(result),
        local_links_checked=local_links(),historical_sha256=result['historical_sha256'],
        owned_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in OWN},
        freeze_scope='Only nine named author assets; excludes pre_review, future independent review, navigation and mainline acceptance.')


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert np.isclose(a,b,rtol=2e-9,atol=2e-11),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=build()
    if args.write:
        with RECEIPT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(result,json.loads(RECEIPT.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1051,all_checks_passed=True,
        mode='exclusive_write' if args.write else 'read_only_compare',
        author_files=len(OWN),historical_files=len(result['historical_sha256']),
        links=result['local_links_checked'],
        dense_history_error=result['author_second_algorithm']['maximum_difference_from_science'],
        new_science_groups=1,new_cognitive_axioms=0,roadmap_complete=False),ensure_ascii=False))


if __name__=='__main__':main()
