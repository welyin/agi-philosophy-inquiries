"""Read-only author replay and narrow artifact verification for round 1052."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from statistics import NormalDist
from urllib.parse import unquote

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1052_checks.json'
OWN=('../research_note_1052.md','analysis.md','weak_public_inputs.json',
     'weak_two_scale_check.py','weak_two_scale_results.json','dependency_update.md',
     'NEXT.md','review.md','verify_round1052.py')


def build():
    run=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'weak_two_scale_check.py')],
        check=True,capture_output=True,text=True,encoding='utf8')
    result=json.loads(run.stdout);assert result['all_checks_passed']
    inputs=json.loads((HERE/'weak_public_inputs.json').read_text(encoding='utf8'))
    x=inputs['belle_ii']['published_effective_amplitude_ratio'];sx=.0019
    y=inputs['atlas']['published_R_W'];sy=math.sqrt(.0022**2+.0036**2+.0014**2)
    # Direct derivative bisection; no polynomial root routine.
    def d(a):return (a-x)/sx**2+2*a*(a*a-y)/sy**2
    lo,hi=.9,1.1;assert d(lo)<0<d(hi)
    for _ in range(80):
        mid=(lo+hi)/2
        if d(mid)>0:hi=mid
        else:lo=mid
    best=(lo+hi)/2
    second=result['secondary_independent_gaussian_diagnostics']
    difference=abs(best-second['common_amplitude_best_fit'])
    assert difference<2e-13
    q=NormalDist().inv_cdf(.9875)
    corners=[a/math.sqrt(b)-1 for a in (x-q*sx,x+q*sx)
              for b in (y-q*sy,y+q*sy)]
    primary=result['primary_marginal_gaussian_bonferroni']
    assert max(abs(a-b) for a,b in zip([min(corners),max(corners)],
        primary['relative_transport_residual_interval']))<1e-14
    for name,sha in result['historical_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    count=0
    for name in OWN:
        p=HERE/name
        if p.suffix!='.md':continue
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text(encoding='utf8')):
            if re.match(r'[A-Za-z][A-Za-z0-9+.-]*://',target) or target.startswith('#'):continue
            dest=(p.parent/unquote(target.split('#',1)[0].strip('<>'))).resolve()
            assert dest==RECEIPT or dest.exists(),(name,target)
            count+=1
    return dict(round=1052,date='2026-10-08',all_checks_passed=True,
        author_status='complete_pending_independent_review',new_project_empirical_groups=1,
        new_theorems_claimed=0,new_cognitive_axioms=0,roadmap_complete=False,
        historical_sha256=result['historical_sha256'],local_links_checked=count,
        author_second_method=dict(derivative_bisection_best=best,root_algorithm_difference=difference,
                                  four_corner_interval=[min(corners),max(corners)]),
        owned_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in OWN},
        freeze_scope='Nine author assets only; excludes independent reviews and mutable navigation')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=build()
    if args.write:
        with RECEIPT.open('x',encoding='utf8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        old=json.loads(RECEIPT.read_text(encoding='utf8'))
        # All hashes and discrete assertions must agree; floats arise from stdlib only.
        assert out==old,'Receipt changed'
    print(json.dumps(dict(round=1052,all_checks_passed=True,author_assets=len(OWN),
        history=len(out['historical_sha256']),links=out['local_links_checked'],
        new_project_empirical_groups=1,new_theorems_claimed=0,roadmap_complete=False)))


if __name__=='__main__':main()
