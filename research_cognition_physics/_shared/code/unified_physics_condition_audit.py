"""Audit the joint-condition ledger; no new scientific round is counted.

Checks include exact familiar hypercharge identities under explicit inputs.
They illustrate joint constraints, not a derivation of the Standard Model.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
TARGET=HERE/'unified_physics_condition_audit_results.json'


def add(*polynomials):
    result={}
    for poly in polynomials:
        for power,value in poly.items():
            result[power]=result.get(power,F(0))+value
    return {key:value for key,value in result.items() if value}


def scale(poly,number):
    return {key:F(number)*value for key,value in poly.items() if number*value}


def multiply(first,second):
    result={}
    for (a,b),x in first.items():
        for (c,d),y in second.items():
            key=(a+c,b+d)
            result[key]=result.get(key,F(0))+x*y
    return {key:value for key,value in result.items() if value}


def cube(poly):
    return multiply(multiply(poly,poly),poly)


def serialize(poly):
    return {f'q^{a} h^{b}':str(value) for (a,b),value in sorted(poly.items())}


def evaluate(poly,q,h):
    return sum(value*q**a*h**b for (a,b),value in poly.items())


def anomaly_example():
    q,h={(1,0):F(1)},{(0,1):F(1)}
    u,d=add(scale(q,-1),scale(h,-1)),add(scale(q,-1),h)
    l=scale(q,-3)
    e=add(scale(q,3),h)
    n=add(scale(q,3),scale(h,-1))
    color=add(scale(q,2),u,d)
    weak=add(scale(q,3),l)
    grav=add(scale(q,6),scale(u,3),scale(d,3),scale(l,2),e)
    abelian=add(scale(cube(q),6),scale(cube(u),3),scale(cube(d),3),scale(cube(l),2),cube(e))
    expected=add(h,scale(q,-3))
    assert color==weak=={}
    assert grav==expected and abelian==cube(expected)
    assert add(grav,n)=={} and add(abelian,cube(n))=={}
    for term in (add(q,h,u),add(q,scale(h,-1),d),add(l,scale(h,-1),e),add(l,h,n)):
        assert term=={}
    fields=dict(Q=q,uc=u,dc=d,L=l,ec=e,nuc=n,H=h)
    hypercharge={name:str(evaluate(poly,F(1,6),F(1,2))) for name,poly in fields.items()}
    baryon_minus_lepton={name:str(evaluate(poly,F(1,3),F(0))) for name,poly in fields.items()}
    assert hypercharge==dict(Q='1/6',uc='-2/3',dc='1/3',L='-1/2',ec='1',nuc='0',H='1/2')
    assert baryon_minus_lepton==dict(Q='1/3',uc='-1/3',dc='-1/3',L='-1',ec='1',nuc='1',H='0')
    # Ordinary SU(2) global anomaly: 3 colored doublets + 1 lepton doublet.
    assert (3+1)%2==0
    return dict(inputs=['D=4','su(3)+su(2)+u(1)','one SM-like chiral generation',
                        'one Higgs doublet','specified Yukawa terms'],
        mixed_gravitational_polynomial=serialize(grav),cubic_polynomial=serialize(abelian),
        without_nuc_nontrivial_solution='h=3q; one overall normalization remains',
        with_dirac_nuc='two free parameters q,h; local anomalies vanish identically',
        bare_gauge_invariant_majorana_condition='2n=0 restores h=3q',
        usual_hypercharges=hypercharge,B_minus_L=baryon_minus_lepton,
        exact_polynomial_checks_passed=True,
        gauge_group_dimension_generations_or_parameters_derived=False)


def run():
    ledger=HERE/'unified_physics_condition_ledger.md'
    body=ledger.read_text('utf8')
    ids=re.findall(r'^\| (C\d\d) \|',body,re.MULTILINE)
    assert ids==[f'C{i:02}' for i in range(1,28)]
    required={'quantum':['C01','C02','C04'],
        'spacetime':['C03','C05','C06','C07','C08','C09'],
        'gravity':['C10','C11','C12','C13'],
        'matter':['C14','C15','C16','C17','C18'],
        'states_and_scales':['C19','C20'],
        'records_and_resources':['C21','C22'],
        'macro_cosmology_and_UV':['C23','C24','C25','C26'],
        'observations':['C27']}
    assert set(sum(required.values(),[]))==set(ids)
    links=[]
    # Markdown local targets in this report/status have no whitespace or titles.
    for path in (ledger,HERE/'round531_drafts/STATUS.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text('utf8')):
            if re.match(r'https?://',target):
                continue
            assert (path.parent/target).resolve().exists() or (path.parent/target).resolve()==TARGET.resolve(),target
            links.append((path.name,target))
    tags=[int(n) for n in re.findall(r'\\tag\{(\d+)\}',body)]
    assert tags==list(range(1,6))
    assert len(re.findall(r'^\$\$\s*$',body,re.MULTILINE))==10
    baseline=json.loads((HERE/'round530_navigation_checks.json').read_text('utf8'))
    assert (baseline['latest_round'],baseline['cumulative_tests'])==(530,2622)
    result=anomaly_example()
    preserved=['conditional_anchor_reference.py','conditional_anchor_reference_results.json',
        'round531_drafts/scope_review.txt','round531_drafts/construction_review.txt']
    return dict(date='2026-09-30',kind='direction_and_joint_constraint_inventory',
        latest_completed_scientific_round=530,cumulative_scientific_tests_unchanged=2622,
        conditions_count=len(ids),coverage=required,local_links_checked=len(links),
        broken_local_links=0,display_formulas=5,exact_cross_constraint_example=result,
        preserved_unfinished_candidate_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                                for name in preserved},
        anchor_candidate_not_completed_or_counted=True,
        old_theorems_reused_not_recounted=True,
        all_document_and_algebra_checks_passed=True,
        whole_physics_unified_or_selected=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert json.loads(TARGET.read_text('utf8'))==result
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(conditions=27,latest_round=530,checks_passed=True,
                         scientific_test_count_unchanged=2622),ensure_ascii=False))
