"""1045 author second algorithm, historical integrity and explicit freeze scope.

Default is read-only. --write exclusively creates the first receipt.
This is not a claim that an independent person has reviewed the round.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import re
import subprocess
import sys
from urllib.parse import unquote
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1045_checks.json'
OWN=('../research_note_1045.md','proof.md','postwrite_source_bridge.py',
     'postwrite_source_bridge_results.json','selection_audit.md',
     'input_dependency_update.md','NEXT.md','review.md','verify_round1045.py')


def independent_checks(result):
    pars=result['parameters']
    p=np.array(pars['populations']); energies=np.array(pars['energies'])
    eta=np.array([1,1,-1,-1]); C=float(p@eta)
    # Direct closed expression on four original thermal sectors; no matrix
    # exponentiation, eigenvectors or call to the author's output function.
    u=np.array(pars['lapses'])*pars['feedback_wait']
    g=pars['g']; mu=20+2+.2+energies
    z=-np.cos(.2); ys=[-np.sin(.2),np.sin(.2)]
    joint=[]
    for y in ys:
        out=[]
        for outcome_sign in [-1,1]:
            matrix=np.zeros((2,2),complex)
            for l in range(2):
                for m in range(2):
                    difference=u[l]-u[m]; total=u[l]+u[m]
                    phase=np.exp(-1j*mu*difference)
                    matrix[l,m]=(np.sum(p*phase)*np.cos(g*difference)
                        +outcome_sign*np.sum(p*eta*phase)
                         *(z*np.cos(g*total)+y*np.sin(g*total)))/4
            assert np.linalg.eigvalsh(matrix)[0]>-1e-12
            out.append(matrix)
        joint.append(out)
    expected=np.array(result['feedback_population_path_distributions'])
    actual=np.array([[np.diag(b).real for b in q] for q in joint])
    probability_error=float(np.max(abs(actual-expected)))
    assert probability_error<2e-10
    offdiag_difference=max(abs(joint[0][b][0,1]-joint[1][b][0,1]) for b in range(2))
    assert offdiag_difference<1e-12
    D=sum(float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2) for a,b in zip(*joint))
    assert abs(D-C*np.sin(.2))<1e-11
    assert abs(D-result['joint_output_half_trace_gap'])<2e-10
    for out in joint:
        assert np.linalg.norm(sum(out)-np.eye(2)/2,2)<1e-10
        assert all(abs(np.trace(b)-.5)<1e-10 for b in out)
    # Exact matrix-unit support in the explicitly Hadamard diagonalized
    # source basis. No numerical SVD is used for this classification check.
    units={(j,j) for j in range(12)}
    effect_by_pair=[[F(1,2),F(sign,2),F(sign,2),F(1,2)] for sign in (1,1,-1,-1)]
    assert all(all(v for v in entries) for entries in effect_by_pair)
    for i in range(4):
        units|={(a,b) for a in (4+2*i,5+2*i) for b in (4+2*i,5+2*i)}
    assert len(units)==20
    assert sum(a!=b for a,b in units)==8
    # All 66 pairwise source-energy separations from the old rational root
    # intervals, not the scientific program's adjacent floating eigenvalues.
    certificate=json.loads((ROOT/'archive_956_989/980/finite_thermal_records_results.json').read_text(encoding='utf-8-sig'))
    certificate=certificate['eigensystem_certificate']['coarse_analytic_certificate']
    roots=[[F(v) for v in pair] for pair in certificate['root_brackets']]
    roots.insert(2,[F(1),F(1)])
    intervals=[]
    for lo,hi in roots:
        intervals.extend([(lo+2,hi+2),
                          (lo+F(219,100),hi+F(219,100)),
                          (lo+F(221,100),hi+F(221,100))])
    separations=[]
    for i in range(12):
        for j in range(i):
            a,b=intervals[i];c,d=intervals[j]
            gap=max(a-d,c-b)
            assert gap>0
            separations.append(gap)
    assert min(separations)==F(4633,250000)
    q_upper=sum(F(str(row[1])) for row in certificate['thermal_intervals'][2:])
    assert q_upper<F(188903,1000000)
    ideal=(1-2*F(188903,1000000))*(F(1,5)-F(1,5)**3/6)
    reset=F(13,14000)+F(1,10000)
    robust=ideal-2*reset
    assert ideal==F(result['certified_ideal_joint_gap_rational_lower'])
    assert robust==F(result['certified_thermal_and_reset_gap_rational_lower'])
    assert robust>F(12,100)
    assert result['round']==1045 and result['new_science_groups']==1
    assert result['new_cognitive_axioms']==0 and not result['stage_goal_completed_here']
    return {'closed_formula_probability_max_error':probability_error,
            'closed_formula_offdiag_difference':float(offdiag_difference),
            'closed_formula_joint_distance':D,
            'exact_visible_matrix_units':len(units),
            'off_diagonal_visible_matrix_units':8,
            'rational_pairwise_spectral_checks':len(separations),
            'rational_minimum_spectral_gap':str(min(separations)),
            'rational_ideal_gap_lower':str(ideal),
            'rational_reset_gap_lower':str(robust)}


def links(allow_missing_receipt=False):
    count=0
    for name in OWN:
        path=(HERE/name).resolve()
        if path.suffix!='.md': continue
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',path.read_text(encoding='utf8')):
            target=target.strip().strip('<>')
            if re.match(r'[a-zA-Z][a-zA-Z0-9+.-]*://',target) or target.startswith('#'): continue
            full=(path.parent/unquote(target.split('#',1)[0])).resolve()
            if not(allow_missing_receipt and full==RECEIPT): assert full.exists(),(name,target)
            count+=1
    return count


def build_receipt():
    result=json.loads((HERE/'postwrite_source_bridge_results.json').read_text(encoding='utf8'))
    process=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'postwrite_source_bridge.py')],
                           check=True,capture_output=True,text=True,encoding='utf8')
    replay=json.loads(process.stdout)
    assert replay['all_scientific_checks_passed'] and replay['mode']=='read_only_compare'
    history=result['historical_sha256']
    for path,sha in history.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    return {'round':1045,'date':'2026-10-08','all_checks_passed':True,
            'review_status':'author_complete_pending_independent_review',
            'external_independent_review_certified_here':False,
            'same_round_science_groups':1,'additional_science_groups':0,
            'science_read_only_compare_passed':True,
            'author_second_algorithm_checks':independent_checks(result),
            'local_links_checked':links(allow_missing_receipt=True),
            'historical_sha256':history,
            'owned_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in OWN},
            'freeze_scope':'Only nine named author assets; excludes navigation, admission scratch and future independent review files.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=build_receipt()
    if args.write:
        with RECEIPT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert result==json.loads(RECEIPT.read_text(encoding='utf8'))
    print(json.dumps({'round':1045,'all_checks_passed':True,
                      'mode':'exclusive_write' if args.write else 'read_only_compare',
                      'owned_files':len(OWN),'historical_files':len(result['historical_sha256']),
                      'local_links':links(),'checks':result['author_second_algorithm_checks']},ensure_ascii=False))


if __name__=='__main__': main()
