"""1028: independent retained-jet checks and limited asset verification."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import re
from urllib.parse import unquote

import full_parent_gravity_bridge as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
NOTE=HERE.parent/'research_note_1028.md'
OUT=HERE/'research_round_1028_checks.json'
OWN=['full_parent_gravity_bridge.py','full_parent_gravity_bridge_results.json',
     'drafts/full_parent_bridge_derivation.md','review.md','verify_round1028.py',
     'selection_audit.md','input_dependency_update_v0_17.md','NEXT.md']
HISTORICAL={
    'archive_1009_/research_note_1018.md','archive_342_369/research_note_358.md',
    'archive_1009_/research_note_1024.md','archive_1009_/research_note_1027.md',
    'archive_1009_/1027/charged_source_selection.py',
    'archive_1009_/1027/charged_source_selection_results.json',
    'archive_1009_/1027/input_dependency_update_v0_16.md','archive_1009_/1027/NEXT.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_990_1008/993/common_candidate_v1.md'}
LINK=re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$',re.M)
EXCLUDED=re.compile(r'\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$',re.S|re.M)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def local_links(text):
    for match in LINK.finditer(EXCLUDED.sub(lambda match:' '*len(match[0]),text)):
        value=(match[1] if match[1] is not None else match[2]).strip()
        target=value[1:value.index('>')] if value.startswith('<') and '>' in value else re.split(r'\s+[\"\']',value,1)[0]
        local=unquote(target.split('#',1)[0])
        if local and not re.match(r'^[a-zA-Z]+:',local) and not local.startswith('//'):
            yield target,local


def convolution(a,b,degree):
    result=[F(0)]*(degree+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=degree:result[i+j]+=x*y
    return result


def derivative(a):return [F(i)*a[i] for i in range(1,len(a))]


# Independent rational-complex representation, without importing 1027 arithmetic.
ZERO=(F(0),F(0))


def z(value):
    if isinstance(value,dict):return F(value['real']),F(value['imag'])
    if isinstance(value,tuple):return value
    return F(value),F(0)


def za(a,b):return a[0]+b[0],a[1]+b[1]


def zm(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]


def zs(a,c):return a[0]*c,a[1]*c


def conj(a):return a[0],-a[1]


def dot(a,b):
    result=ZERO
    for x,y in zip(a,b):result=za(result,zm(conj(x),y))
    return result


def mv(matrix,vector):
    answer=[]
    for row in matrix:
        value=ZERO
        for x,y in zip(row,vector):value=za(value,zm(x,y))
        answer.append(value)
    return answer


def generators():
    half=F(1,2)
    return {'Y':[[z(half),ZERO],[ZERO,z(half)]],
            'SU2_1':[[ZERO,z(half)],[z(half),ZERO]],
            'SU2_2':[[ZERO,(F(0),-half)],[(F(0),half),ZERO]],
            'SU2_3':[[z(half),ZERO],[ZERO,z(-half)]]}


def current(r,dr,t):
    difference=za(dot(r,mv(t,dr)),zs(dot(dr,mv(t,r)),F(-1)))
    return zm((F(0),F(1,2)),difference)


def verify(prospective=False):
    fresh=science.run()
    science.compare(fresh,json.loads(science.OUT.read_text('utf8')))
    assert fresh['round']==1028 and fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups']==1 and fresh['cumulative_test_groups']==3805
    assert fresh['new_cognitive_axioms']==0 and fresh['goal_complete'] is False
    assert fresh['code_sha256']==sha(HERE/'full_parent_gravity_bridge.py')
    assert 'complete first-order minimal transformations and gravity bracket' in fresh['retained_inputs']
    assert 'standard internal YM/local-Lorentz action fixes O at every order' in fresh['retained_inputs']
    branch=fresh['full_higgs_branch']
    assert branch['taylor_degree']==12 and branch['equation_orders']==list(range(11))
    f=list(map(F,branch['f_coefficients']));assert len(f)==13
    assert f[:2]==[F(5,4),F(2,7)]
    lam,v2=F(branch['lambda_h']),F(branch['v_squared'])
    assert lam==F(2,5) and v2==F(9,4)
    cubic=convolution(convolution(f,f,10),f,10)
    fpp=derivative(derivative(f))
    residual=[fpp[k]-lam*(cubic[k]-v2*f[k]) for k in range(11)]
    assert residual==[F(0)]*11 and branch['radial_equation_residual']==['0']*11
    # Four real Higgs components: h=(0,f,0,0), with a radial quartic potential.
    hh=[[F(0)]*13,list(f),[F(0)]*13,[F(0)]*13]
    norm=convolution(f,f,10);norm[0]-=v2
    full_higgs_residual=[]
    for component in hh:
        gradient=[lam*x for x in convolution(norm,component,10)]
        second=derivative(derivative(component))
        full_higgs_residual.extend(x-y for x,y in zip(second,gradient))
    assert full_higgs_residual==[F(0)]*44
    fp=derivative(f)
    for label,t in generators().items():
        saved=branch['gauge_current_series'][label]
        assert len(saved)==12
        for k in range(12):
            value=ZERO
            for j in range(k+1):
                value=za(value,current([ZERO,z(f[j])],[ZERO,z(fp[k-j])],t))
            assert value==z(saved[k])==ZERO
    assert branch['gauge_current_orders']==list(range(12))
    assert branch['color_current_identically_zero'] and branch['zero_fermion_branch_exact_by_homogeneity']
    assert branch['fermion_equation_minimum_fermion_degree']=={'kinetic':1,'Yukawa':1,'Weinberg':1}
    assert F(branch['sigma_equation_at_zero'])==0 and F(branch['nonzero_portal_calibration'])==F(2,9)
    shifted=convolution(f,f,11);shifted[0]-=v2
    potential=[lam*x/4 for x in convolution(shifted,shifted,11)]
    kinetic=[x/2 for x in convolution(fp,fp,11)]
    tyy=[x-y for x,y in zip(kinetic,potential)]
    assert list(map(F,branch['stress_T_y_y_coefficients']))==tyy
    assert derivative(tyy)==[F(0)]*11
    assert branch['full_stress_divergence_coefficients']==['0']*11
    assert F(branch['observable_at_origin'])==f[0]**2/2==F(25,32)
    assert F(branch['observable_y_derivative'])==f[0]*f[1]==F(5,14)

    internal=fresh['internal_and_phase_checks']
    rr=list(map(z,internal['arbitrary_complex_doublet']))
    generator=[[z(v) for v in row] for row in internal['arbitrary_hermitian_internal_generator']]
    assert all(generator[i][j]==conj(generator[j][i]) for i in range(2) for j in range(2))
    dr=[zm((F(0),F(1)),x) for x in mv(generator,rr)]
    variation=zs(za(dot(dr,rr),dot(rr,dr)),F(1,2))
    assert variation==z(internal['internal_O_variation'])==ZERO
    a=F(internal['phase_only_A_zero_theta_prime']);assert a==F(3,5)
    radial=[ZERO,z(f[0])];phase_derivative=[ZERO,(f[1],a*f[0])]
    expected_currents={'Y':F(-15,32),'SU2_1':F(0),'SU2_2':F(0),'SU2_3':F(15,32)}
    for label,t in generators().items():
        value=current(radial,phase_derivative,t)
        assert value==z(expected_currents[label])==z(internal['bad_phase_current_over_coupling'][label])
        assert current(radial,[ZERO,z(f[1])],t)==z(internal['restored_current_over_coupling'][label])==ZERO
    assert F(internal['restored_gprime_B_y'])==-2*a
    assert [z(x) for x in internal['covariant_derivative_after_restoring_connection']]==[ZERO,z(f[1])]
    assert internal['bad_phase_A_zero_violates_gauge_equation']
    assert internal['field_dependent_internal_parameters_also_annihilate_O']

    killing=fresh['invariant_killing_bridge']
    xi=list(map(F,killing['xi']));eta=list(map(F,killing['zeta_at_origin']))
    deta=[list(map(F,row)) for row in killing['derivative_zeta']]
    assert xi==[0,1,0,0] and eta==[0,0,0,0]
    assert deta==[[0,0,0,0],[0,0,1,0],[0,-1,0,0],[0,0,0,0]]
    bracket=[sum(xi[mu]*deta[mu][nu] for mu in range(4)) for nu in range(4)]
    assert bracket==list(map(F,killing['bracket']))==[0,0,1,0]
    signs=[-1,1,1,1]
    assert all(signs[nu]*deta[mu][nu]+signs[mu]*deta[nu][mu]==0 for mu in range(4) for nu in range(4))
    assert all(F(v)==0 for matrix in killing['killing_symmetric_derivatives'] for row in matrix for v in row)
    # Preserve derivatives of eta before restricting to r=r(y); eta(0)=0 is insufficient.
    nested_first=[ZERO,z(sum(xi[mu]*deta[mu][2]*f[1] for mu in range(4)))]
    nested_second=[ZERO,ZERO]
    assert [[z(v) for v in row] for row in killing['nested_field_variations']]==[nested_first,nested_second]
    first_terms=[dot(nested_first,radial),ZERO,ZERO,dot(radial,nested_first)]
    second_terms=[ZERO]*4
    assert [[z(v) for v in row] for row in killing['four_product_rule_terms']]==[first_terms,second_terms]
    response=zs(za(first_terms[0],first_terms[3]),F(1,2))
    assert response==z(F(5,14))==z(killing['composite_commutator'])==z(killing['bracket_action_on_O'])
    assert [z(v) for v in killing['ordered_O_variations']]==[response,ZERO]
    assert killing['all_higher_Killing_derivatives_zero_by_affinity']
    examples=killing['examples'];assert len(examples)==2
    for name,case in examples.items():
        q=list(map(F,case['q']));kappa=list(map(F,case['kappa']))
        wanted=[[z((q[i]*q[j]-(kappa[i]*q[i] if i==j else F(0)))*F(5,14)) for j in range(2)] for i in range(2)]
        assert wanted==[[z(v) for v in row] for row in case['direct_composite_closure_residual']]
        assert any(v!=ZERO for row in wanted for v in row)==(name=='two_active_incompatible')
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
    return dict(round=1028,date='2026-10-08',all_delivery_checks_passed=True,scientific_result_reproduced=True,
                new_calibration_groups=1,cumulative_research_groups=3805,new_cognitive_axioms=0,
                calibration_kind='same-parent radial Higgs and invariant Killing bridge',
                exact_Taylor_degree=12,radial_EOM_coefficients_checked=11,real_Higgs_EOM_coefficients_checked=44,
                electroweak_current_coefficients_checked=48,stress_divergence_coefficients_checked=11,
                composite_product_terms_per_ordering=4,actual_gravity_examples=2,nonzero_invariant_derivative='5/14',
                forbidden_phase_currents={'Y':'-15/32','SU2_3':'15/32'},restored_gprime_B_y='-6/5',
                finite_Taylor_not_a_global_solution=True,recoil_existence_is_reused_analytic_lemma=True,
                complete_first_order_minimal_representative_remains_input=True,all_first_order_classes_classified=False,
                quantum_parent_realization_certified=False,live_navigation_frozen=False,neighboring_round_frozen=False,
                goal_completed=False,visual_checks_performed=False,frozen_current_files=9,
                historical_input_files=len(HISTORICAL),local_links_checked=links,
                historical_source_sha256=fresh['historical_source_sha256'],
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in assets})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write','--write-receipt',dest='write',action='store_true')
    args=parser.parse_args();result=verify(prospective=args.write)
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    else:assert result==json.loads(OUT.read_text('utf8'))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},ensure_ascii=False,indent=2))
