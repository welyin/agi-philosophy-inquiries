"""Reuse696 endpoints to screen a proposed699 lower-bound shortcut.
No new sphere/holonomy experiment; exact arithmetic on inherited formulas.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
TARGET=HERE/'free_domination_entry_results.json'


def run():
    previous=json.loads((ARCHIVE/'joint_static_haar_majorant_results.json').read_text('utf8'))['exact_integral_upper']
    b00=F(previous['inherited_B00']);target=F(previous['unproved_offdiagonal_sufficient_lower'])
    required=target/b00
    #696: at identity holonomy,E=F a weak unit vector, physical ratio is
    #(1+1/sqrt2)^16 and auxiliary ratio is[(1+1/sqrt2)/2]^16.
    a,c=F(1),F(0)
    for _ in range(32):a,c=a+c,a/2+c
    a/=2**16;c/=2**16
    lo,hi=F(1414213,10**6),F(1414214,10**6)
    assert lo*lo<2<hi*hi and a>0 and c>0
    ratio_lo=a+c*lo;ratio_hi=a+c*hi
    assert ratio_hi<500<3000<required
    deps=('research_note_696.md','joint_dynamic_auxiliary_integral_results.json',
        'research_note_697.md','research_note_698.md','joint_static_haar_majorant_results.json')
    return dict(date='2026-10-02',entry_round=699,new_formal_round=False,
        inherited_identity_seam_weak_aligned_fields=True,
        endpoint_ratio_Qsqrt2=[str(a),str(c)],ratio_bounds=[str(ratio_lo),str(ratio_hi)],
        endpoint_ratio_upper_decimal=float(ratio_hi),required_uniform_multiplier=str(required),
        required_uniform_multiplier_decimal=float(required),
        constant_pointwise_free_domination_at_required_level_impossible=True,
        full_integral_lower_bound_not_decided=True,
        next_action='Use integrated tensor/character coefficients or a nonconstant lower minorant; do not repeat fixed-seam scans.',
        scope='A shortcut screened by696 and697 continuity; not an integral counterexample or new699 round.',
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps},all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(entry_round=699,endpoint_upper=r['endpoint_ratio_upper_decimal'],
        required=r['required_uniform_multiplier_decimal'],new_round=False,all_checks_passed=True)))
