"""1027 limited delivery checks; no live navigation or neighboring assets."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import charged_source_selection as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
NOTE=HERE.parent/'research_note_1027.md'
OUT=HERE/'research_round_1027_checks.json'
OWN=['charged_source_selection.py','charged_source_selection_results.json',
     'drafts/source_selection_derivation.md','review.md','verify_round1027.py',
     'selection_audit.md','input_dependency_update_v0_16.md','NEXT.md']
HISTORICAL={
    'archive_1009_/research_note_1026.md','archive_1009_/1026/input_dependency_update_v0_15.md',
    'archive_1009_/1026/NEXT.md',*[f'archive_1009_/research_note_{n}.md' for n in (1017,1018,1024)],
    'archive_301_341/research_note_326.md',
    *[f'archive_702_741/research_note_{n}.md' for n in (732,735)],
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_990_1008/993/common_candidate_v1.md'}
LINK=re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$',re.M)
EXCLUDED=re.compile(r'\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$',re.S|re.M)
EXPECTED={
    'three_generation_charged_only':(58,57,1),
    'hypercharge_force_only_diagnostic':(63,59,4),
    'B993_zero_portal':(63,61,2),'B993_nonzero_portal':(63,62,1),
    'one_generation_nonzero_portal':(23,22,1),
    'three_free_neutral_fermions_extra':(72,62,10)}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def local_links(text):
    for match in LINK.finditer(EXCLUDED.sub(lambda match:' '*len(match[0]),text)):
        value=(match[1] if match[1] is not None else match[2]).strip()
        target=value[1:value.index('>')] if value.startswith('<') and '>' in value else re.split(r'\s+[\"\']',value,1)[0]
        local=unquote(target.split('#',1)[0])
        if local and not re.match(r'^[a-zA-Z]+:',local) and not local.startswith('//'):
            yield target,local


def exact_rank(rows):
    a=[[F(v) for v in row] for row in rows]
    rank=0
    for column in range(len(a[0])):
        selected=next((i for i in range(rank,len(a)) if a[i][column]),None)
        if selected is None:continue
        a[rank],a[selected]=a[selected],a[rank]
        pivot=a[rank][column]
        for i in range(rank+1,len(a)):
            if a[i][column]:
                factor=a[i][column]/pivot
                a[i]=[x-factor*y for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank


def rank_mod(rows,prime=103):
    a=[[int(F(v).numerator)*pow(int(F(v).denominator),-1,prime)%prime for v in row] for row in rows]
    rank=0
    for column in range(len(a[0])):
        selected=next((i for i in range(rank,len(a)) if a[i][column]),None)
        if selected is None:continue
        a[rank],a[selected]=a[selected],a[rank]
        inverse=pow(a[rank][column],-1,prime)
        for i in range(rank+1,len(a)):
            if a[i][column]:
                factor=a[i][column]*inverse%prime
                a[i]=[(x-factor*y)%prime for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank


def complex_array(matrix):
    return np.array([[complex(float(v.r),float(v.i)) for v in row] for row in matrix])


def su_basis(n):
    out=[]
    for i in range(n):
        for j in range(i+1,n):
            x=np.zeros((n,n),complex);y=x.copy()
            x[i,j]=x[j,i]=.5;y[i,j]=-.5j;y[j,i]=.5j
            out.extend([x,y])
    for k in range(1,n):
        t=np.zeros((n,n),complex)
        for i in range(k):t[i,i]=.5
        t[k,k]=-k/2;out.append(t)
    return out


def independent_representation_check():
    su2,su3=su_basis(2),su_basis(3)
    assert all(np.array_equal(a,complex_array(b)) for a,b in zip(su2,science.SU2))
    assert all(np.array_equal(a,complex_array(b)) for a,b in zip(su3,science.SU3))
    expected=[('q',3,2,F(1,6),False),('uc',3,1,F(-2,3),True),
              ('dc',3,1,F(1,3),True),('l',1,2,F(-1,2),False),('ec',1,1,F(1),False)]
    assert science.MODULES==expected
    count=0
    for item in expected:
        name,nc,nw,y,anti=item
        actual=science.module_generators(item)
        wanted=[('Y',float(y)*np.eye(nc*nw))]
        if nc==3:
            wanted += [('3',np.kron(-t.conjugate() if anti else t,np.eye(nw))) for t in su3]
        if nw==2:wanted += [('2',np.kron(np.eye(nc),t)) for t in su2]
        assert len(actual)==len(wanted)
        for (ag,am),(bg,bm) in zip(actual,wanted):
            assert ag==bg and np.array_equal(complex_array(am),bm)
            assert np.array_equal(bm,bm.conj().T)
            count+=1
    assert count==35
    ty=np.zeros((5,5))
    ty[0,2]=ty[1,3]=-.5;ty[2,0]=ty[3,1]=.5
    actual=complex_array(science.real_generator(science.scale(science.eye(2),F(1,2)),True))
    assert np.array_equal(actual,ty)
    assert np.array_equal(ty@ty,np.diag([-.25]*4+[0]))
    return count


def independent_jet_check(saved):
    eta=[F(-1),F(1),F(1),F(1)]
    field=[[F(0) for _ in range(4)] for _ in range(4)]
    field[1][0]=F(2);field[0][1]=F(-2)
    current=[F(1),F(0),F(0),F(1)]
    derivative=[[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    derivative[1][1][0]=F(1);derivative[1][0][1]=F(-1)
    derivative[0][0][3]=F(1);derivative[0][3][0]=F(-1)
    div=[sum(eta[mu]*eta[nu]*derivative[mu][mu][nu] for mu in range(4)) for nu in range(4)]
    assert div==[-v for v in current]
    assert all(derivative[a][b][c]+derivative[b][c][a]+derivative[c][a][b]==0
               for a in range(4) for b in range(4) for c in range(4))
    force=[sum(field[nu][mu]*current[mu] for mu in range(4))/2 for nu in range(4)]
    assert force==[0,1,0,0]
    assert saved['F_lower']==[[str(v) for v in row] for row in field]
    assert saved['current_upper']==list(map(str,current))
    assert saved['partial_mu_F_upper_mu_nu']==list(map(str,div))
    assert saved['bad_source_divergence_covector']==list(map(str,force))
    # Direct product variation: delta(A*j)=k*A'*j+r*A*j'.
    # Modulo the derivative r*(A*j)', the coefficient is (k-r)*A'*j.
    k,r,aprime,j=F(1),F(3,2),F(2),F(1)
    coefficient=(k-r)*aprime*j
    assert coefficient==F(saved['integrated_vertex_variation_coefficient'])==-1
    return list(map(str,force))


def verify(prospective=False):
    fresh=science.run()
    science.compare(fresh,json.loads(science.OUT.read_text('utf8')))
    assert fresh['round']==1027
    assert fresh['new_calibration_groups']==1 and fresh['cumulative_test_groups']==3804
    assert fresh['new_cognitive_axioms']==0
    assert fresh['code_sha256']==sha(HERE/'charged_source_selection.py')
    cases=fresh['exact_source_kernels'];assert cases.keys()==EXPECTED.keys()
    rank_certificate={}
    for name,expected in EXPECTED.items():
        case=cases[name];labels=case['labels']
        assert tuple(case[key] for key in ('unknown_count','independent_equation_rank','kernel_dimension'))==expected
        rows=[[F(row.get(label,'0')) for label in labels] for row in case['exact_constraint_rows']]
        assert len(rows)==case['deduplicated_equation_count']
        assert exact_rank(rows)==rank_mod(rows)==expected[1]
        kernel=[[F(row.get(label,'0')) for label in labels] for row in case['kernel_basis']]
        assert exact_rank(kernel)==len(kernel)==expected[2]
        assert all(sum(x*y for x,y in zip(row,vec))==0 for row in rows for vec in kernel)
        assert len(labels)-len(kernel)==expected[1]
        rank_certificate[name]={'unknowns':len(labels),'rank':expected[1],'kernel_dimension':len(kernel)}
    generator_count=independent_representation_check()
    rep=fresh['representation_checks']
    assert rep['exact_lie_closure_pairs']==73 and rep['three_generation_left_weyl_components']==45
    assert rep['all_matter_hypercharges_nonzero'] and rep['higgs_real_hypercharge_square']=='-I_4/4'
    bad=independent_jet_check(fresh['ward_and_gauss_jet'])
    portal=fresh['portal_certificate']
    assert F(portal['mixed_hessian'])==F(12,5)
    assert F(portal['bad_commutator_entry'])==(F(1)-F(7,4))*F(12,5)==F(-9,5)
    assert portal['vanishes_for_equal_weights_or_zero_portal']
    transport=fresh['basis_transport']
    assert transport['fermion_generator_cases']==35 and 0<=transport['max_transport_residual']<1e-12
    assert set(fresh['historical_source_sha256'])==HISTORICAL
    for name,digest in fresh['historical_source_sha256'].items():assert sha(BASE/name)==digest,name
    assets=[NOTE]+[HERE/name for name in OWN]
    assert len(assets)==9 and not any(p.name in {'README.md','research_direction.md','RESEARCH_STATE.md'} for p in assets)
    links=0
    for path in assets:
        assert path.is_file(),path
        if path.suffix=='.md':
            text=path.read_text('utf-8-sig');assert text.count('$$')%2==0,path
            for target,local in local_links(text):
                resolved=(path.parent/local.replace('\\','/')).resolve()
                assert resolved.exists() or (prospective and resolved==OUT),(path,target)
                links+=1
    assert all(f'## {n}.' in NOTE.read_text('utf8') for n in range(1,11))
    assert '独立科学签审通过' in (HERE/'review.md').read_text('utf8')
    return dict(round=1027,date='2026-10-08',all_delivery_checks_passed=True,scientific_result_reproduced=True,
                new_calibration_groups=1,cumulative_research_groups=3804,new_cognitive_axioms=0,
                calibration_kind='leading-source exact representation kernel and compatible point-jet audit',
                exact_source_kernel_certificates=rank_certificate,independent_rank_prime=103,
                exact_elimination_and_independent_kernel_upper_bounds=True,
                actual_representation_generators_checked=generator_count,Lie_closure_pairs=73,
                bad_weight_divergence_coefficient=bad,Bianchi_components_checked=64,
                portal_bad_weight_commutator='-9/5',maximum_basis_transport_residual=transport['max_transport_residual'],
                Grassmann_coefficient_not_a_quantum_state=True,full_QFT_or_gravity_generated=False,
                all_higher_derivative_source_classes_classified=False,
                live_navigation_frozen=False,neighboring_round_frozen=False,goal_completed=False,visual_checks_performed=False,
                frozen_current_files=9,historical_input_files=len(HISTORICAL),local_links_checked=links,
                historical_source_sha256=fresh['historical_source_sha256'],
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in assets})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write','--write-receipt',dest='write',action='store_true')
    args=parser.parse_args();result=verify(prospective=args.write)
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    else:assert result==json.loads(OUT.read_text('utf8'))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},ensure_ascii=False,indent=2))
