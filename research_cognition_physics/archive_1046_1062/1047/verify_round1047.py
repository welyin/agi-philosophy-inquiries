"""1047 author second algorithm, historical hashes, links and freeze scope.

The default invocation is read-only. --write creates the receipt once.
This does not certify review by an independent person.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, pi, sqrt
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
RECEIPT=HERE/'research_round_1047_checks.json'
OWN=('../research_note_1047.md','proof.md','newton_tt_matching.py',
     'newton_tt_matching_results.json','selection_audit.md',
     'input_dependency_update.md','NEXT.md','review.md','verify_round1047.py')


def second_algorithm(result):
    # Integration by parts: gradient of s wave times gradient of d wave.
    # No second derivatives or calls to the science program are used here.
    coefficients=[]
    for n in (1,2):
        c=F(1,n)+F(1,3)
        derivative=(-F(1),) if n==1 else (-F(2),F(1,2))
        def ri(power):
            return sum(v*F(factorial(power+j))/c**(power+j+1)
                       for j,v in enumerate(derivative))
        kinetic=F(4,3)*ri(3)-F(4,45)*ri(4)
        expected=F(result['exact_radial_factors'][str(n)]['kinetic'])
        assert kinetic==expected
        coefficients.append(str(kinetic))
    # Direct time-domain composite Simpson, versus the author's Gaussian rule.
    intervals=20000
    u=np.linspace(0,2*pi,intervals+1); du=2*pi/intervals
    weights=np.ones(intervals+1);weights[1:-1:2]=4;weights[2:-1:2]=2
    weights*=du/3
    pulse=(1-np.cos(u))*np.cos(u)/2
    time_errors=[];prob_errors=[]
    for q,tau,row in zip((1,F(5,32)),(-1/(8*sqrt(2)),F(1024,15625)),result['pulse']['rows']):
        J=np.sum(weights*pulse*np.exp(-1j*float(q)*u))
        stored=complex(row['fourier_real'],row['fourier_imag'])
        time_errors.append(float(abs(J-stored)))
        C=abs(float(tau)*J/(4/9))**2
        prob_errors.append(float(abs(C-row['coefficient_full'])))
    assert max(time_errors)<1e-11 and max(prob_errors)<1e-12
    a,b=F(16,25),F(1,64);x=F(3200,1049)
    assert a*x-1==1-b*x==F(999,1049)
    assert result['pulse']['sharp_relative_response_error']=='999/1049'
    assert result['pulse']['same_pulse_for_both_outputs']
    assert result['probabilities_are_leading_coefficients']
    assert not result['finite_nonzero_amplitude_remainder_certified']
    return dict(gradient_radial_kinetic=coefficients,
                simpson_fourier_max_error=max(time_errors),
                simpson_probability_coefficient_max_error=max(prob_errors),
                exact_two_output_minimax='999/1049',
                no_additional_science_groups=True)


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


def build():
    result=json.loads((HERE/'newton_tt_matching_results.json').read_text(encoding='utf8'))
    replay=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'newton_tt_matching.py')],
        check=True,capture_output=True,text=True,encoding='utf8')
    assert json.loads(replay.stdout)['all_scientific_checks_passed']
    history=result['historical_sha256']
    for path,sha in history.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    return dict(round=1047,date='2026-10-08',all_checks_passed=True,
        review_status='author_complete_pending_independent_review',
        external_independent_review_certified_here=False,
        same_round_science_groups=1,additional_science_groups=0,
        science_read_only_compare_passed=True,
        author_second_algorithm_checks=second_algorithm(result),
        local_links_checked=links(True),historical_sha256=history,
        owned_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in OWN},
        freeze_scope='Only nine named author assets; excludes navigation, admission and future independent review files.')


def compare(actual,saved,path='root'):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys(),path
        for k in actual: compare(actual[k],saved[k],path+'.'+k)
    elif isinstance(actual,list):
        assert len(actual)==len(saved),path
        for i,(a,b) in enumerate(zip(actual,saved)):compare(a,b,path+f'[{i}]')
    elif isinstance(actual,float):
        assert np.isclose(actual,saved,rtol=2e-9,atol=2e-11),(path,actual,saved)
    else: assert actual==saved,(path,actual,saved)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=build()
    if args.write:
        with RECEIPT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else: compare(result,json.loads(RECEIPT.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1047,all_checks_passed=True,
        mode='exclusive_write' if args.write else 'read_only_compare',
        author_files=len(OWN),historical_files=len(result['historical_sha256']),
        links=links(),checks=result['author_second_algorithm_checks']),ensure_ascii=False))


if __name__=='__main__': main()
