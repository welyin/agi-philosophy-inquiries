"""699: exact original offdiagonal integral and full-measure negative Gram.

The Clifford annihilator, two complete S9 integrations and original Weyl
constant term use exact integer/Q(sqrt(2)) arithmetic. Matrix comparisons
are implementation checks, not the proof of an integral sign.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round699_drafts'))
import clifford_moment_certificate as alg
import offdiagonal_exact_sphere as sphere
import offdiagonal_exact_group as group
TARGET=HERE/'joint_offdiagonal_haar_certificate_results.json'


def native(value):return json.loads(json.dumps(value))


def calibrations(data,average):
    # Independent exactly soluble free limit of the same formal algebra.
    free=alg.run(free=True)
    x,xp,y,yp=alg.VAR[2:6]
    product=alg.pmul({0:(1,0),x:(1,0),y:(1,0)},
                     {0:(1,0),xp:(1,0),yp:(1,0)})
    square=alg.pmul(product,product)
    expect={m:alg.scale(c,F(free['e4_denominator'],16)) for m,c in square.items()}
    assert {m:tuple(c) for m,c in free['elementary'][4]}==expect
    assert all(c[1]==0 for c in expect.values())
    # Full two-sphere result at h=I agrees with independent 696 integration.
    pair=[0,0]
    for (b,c),v in average['auxiliary_character_coefficients']:
        multiplicity=math.comb(b+2,2)*(c+1)
        for j in range(2):pair[j]+=multiplicity*v[j]
    exact=[F(x,average['denominator']) for x in pair]
    reference=F(489381062671,199561360441344)
    assert exact==[reference,F(0)]
    # Verify actual original Dirac/weak-volume dictionary on all generators.
    ts=sphere.original.b.internal.T
    zero=np.zeros((16,16),complex);eye=np.eye(32)
    gammas=[np.block([[zero,t],[t.conj().T,zero]]) for t in ts]
    signs=np.array([-1 if odd else 1 for _,odd in group.previous.moment.weights()])
    rc=np.diag(signs);r=np.block([[rc,zero],[zero,rc]])
    max_cliff=max(float(np.linalg.norm(a@b+b@a-(2*eye if i==j else 0),2))
                  for i,a in enumerate(gammas) for j,b in enumerate(gammas))
    max_reverse=max(float(np.linalg.norm(r@g-((-1 if i>=6 else 1)*g@r),2))
                    for i,g in enumerate(gammas))
    assert max_cliff<1e-12 and max_reverse<1e-12
    # Independent full-free Haar calibration through the same sphere algorithm.
    free_average=sphere.integrate(free)
    coeff={tuple(key):alg.scale(tuple(v),2**16)
           for key,v in free_average['auxiliary_character_coefficients']}
    b00=F(group.previous.blocks.old.numerical.reference())
    primes=group.prime_stream();modular_checks=[]
    for n in (20,24):
        p=next(primes)
        value,root,count=group.modular(p,n,coeff,free_average['denominator'])
        assert value==[group.previous.rational_mod(b00,p),0]
        modular_checks.append(dict(n=n,p=p,root=root,points=count))
    return dict(free_quartic_polynomial_exact=True,
        same_seam_two_sphere_exact=[str(v) for v in exact],
        matches_independent_round696=True,max_clifford_error=max_cliff,
        max_weak_volume_error=max_reverse,free_B00_modular_calibrations=modular_checks)


def run():
    data=alg.run()
    saved=json.loads((HERE/'round699_drafts/clifford_moment_probe_results.json').read_text('utf8'))
    assert native(data)==saved
    average=sphere.run()
    assert average==json.loads(sphere.TARGET.read_text('utf8'))
    calibration=calibrations(data,average)
    integral=group.run()
    assert integral==json.loads(group.TARGET.read_text('utf8'))
    assert integral['target_met'] and integral['strict_full_group_and_sphere_B_negative']
    assert max(average['max_total_character_degree']+8+d for d in (3,3,3,2))<20
    deps=('research_note_382.md','research_note_383.md','research_note_384.md',
        'research_note_386.md','research_note_425.md','research_note_522.md','research_note_523.md',
        'research_note_591.md','research_note_673.md','research_note_675.md','research_note_695.md',
        'research_note_696.md','research_note_697.md','research_note_698.md',
        'joint_full_holonomy_reduction.py','joint_static_haar_majorant.py',
        'joint_static_haar_majorant_results.json','round699_drafts/clifford_moment_probe.py',
        'round699_drafts/clifford_moment_probe_results.json','round699_drafts/clifford_moment_certificate.py',
        'round699_drafts/offdiagonal_exact_sphere.py','round699_drafts/offdiagonal_exact_sphere_results.json',
        'round699_drafts/offdiagonal_exact_group.py','round699_drafts/offdiagonal_exact_group_results.json')
    return dict(date='2026-10-02',round=699,tests_run=2,failures=0,errors=0,
        exact_clifford=dict(identity_zero=data['characteristic_identity_zero'],
            power_counts=data['power_monomial_counts'],e4_terms=len(data['elementary'][4]),
            e4_denominator=data['e4_denominator']),
        exact_sphere=average,independent_calibrations=calibration,exact_integral=integral,
        scope=dict(whole_original_group_and_two_S9=True,original_B10_exact=True,
            full_B_strict_negative=True,fixed_candidate_large_tau_RP_failure_by_round675=True,
            numerical_finite_tau_threshold_NOT_computed=True,all_tau_failure_NOT_proved=True,
            original_HF_identity_NOT_established=True,positive_HF_branch_NOT_refuted=True,
            cognition_axioms_or_unification_NOT_refuted=True,continuum_limit_NOT_refuted=True,
            inherited_space_interfaces_unchanged=True),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=699,B10=result['exact_integral']['B10_decimal'],
        negative_upper=result['exact_integral']['quadratic_upper_decimal'],all_checks_passed=True)))
