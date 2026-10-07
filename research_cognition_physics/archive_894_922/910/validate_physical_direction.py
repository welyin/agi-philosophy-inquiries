"""910 verify the actual constrained direction and its inherited future-test mapping."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib
import numpy as np
import physical_direction_dual_bridge as bridge
HERE=Path(__file__).resolve().parent;TARGET=HERE/'physical_direction_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(reproduce=True):
    data=json.loads(bridge.TARGET.read_text('utf-8'))
    for p,s in data['source_hashes'].items():assert sha(HERE.parent/p)==s,p
    a,b,c=F(31,100),F(27,100),F(21,100)
    da=-(a*a+b*b)*c*c/(4*a*(b*b+c*c/4))
    assert 2*a*(b*b+c*c/4)*da+(a*a+b*b)*c*c/2==0
    for row in data['rows']:
        assert abs(row['initial_direction']['Ax_color_T1']-float(da))<1e-15
        assert row['initial_linear_constraints']['H']<1e-14
        assert max(row['initial_linear_constraints'][k] for k in ('M','Gauss','harmonic'))<1e-14
        assert row['full_source_pairing']['jet_pairing_identity_residual']<1e-17
        assert row['conditional_future_read_coefficient']==.75*row['full_source_pairing']['total'][2]
    if reproduce:
        bg=bridge.mb.Background(16);src=bridge.loop.LoopSource(bg,anchor_order=2,segments=8);src.kernels()
        u=bridge.PhysicalDirection(16);value=src.response(u.as_variation)
        assert np.allclose(value['total'],data['rows'][0]['full_source_pairing']['total'],rtol=1e-8,atol=1e-15)
    low,high=data['rows'];comp=low['previous_finite_family_comparison']
    ratio=comp[0]['finite_family_source_minus_exact_tangent']/comp[1]['finite_family_source_minus_exact_tangent']
    assert 3.9<ratio<4.1
    difference=abs(high['conditional_future_read_coefficient']-low['conditional_future_read_coefficient'])
    assert difference<2e-12
    return dict(round=910,date='2026-10-06',cumulative_numbered_groups=3695,
        implementation_checks_passed=True,one_actual_tangent_source_pairing_reproduced=reproduce,
        all_N16_N24_cases_previously_executed=True,exact_initial_Ax_derivative=str(da),
        initial_magnetic_constraint_derivative_exact_zero=True,
        earlier_finite_difference_errors_ratio=ratio,N16_N24_conditional_read_difference=difference,
        conditional_future_read_coefficients=[row['conditional_future_read_coefficient'] for row in data['rows']],
        original_background_source_and_parameters_unchanged=True,
        inherited_863_Green_identity_not_counted_as_new_theorem=True,
        argument_scope='Actual863 constrained-family tangent transported by904 on902 and paired with the complete908 original first-jet source. The863 conserved future test and870 pure receiver leg give a conditional finite mixed response coefficient; no completed physical instrument, continuous error bound, or872 total mean feedback.',
        future_test_kernel_array_computed=False,actual_continuous_source_support_certified=False,
        actual_continuous_response_error_certified=False,complete872_response=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(HERE)):sha(p) for p in (Path(__file__),HERE/'physical_direction_dual_bridge.py',bridge.TARGET)})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    v=run()
    if a.write:TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in v.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
