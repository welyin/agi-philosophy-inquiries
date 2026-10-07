"""908: original863 constraint-preserving color family through the actual source.
Fixed material curve/profile; whole902 background is evolved for each family.
"""
from pathlib import Path
import argparse,json
import numpy as np
import material_background as bg
import relational_loop_source as loops
HERE=Path(__file__).resolve().parent;TARGET=HERE/'physical_family_source_checks_results.json'
def build(epsilon,N=16):
    original=bg.ev.Model
    class FamilyModel(original):
        def initial(self):
            y,stats=super().initial();b=.27;c=.21*(1+epsilon);S=.008868915
            a=np.sqrt((S-b*b*c*c/4)/(b*b+c*c/4))
            y['A'][...,0,0]=a;y['A'][...,2,3]=c
            return y,stats
    bg.ev.Model=FamilyModel
    try:return bg.Background(N)
    finally:bg.ev.Model=original

def run():
    base=build(0.);q=loops.LoopSource(base,2,16);q.kernels();rows=[]
    for eps in (.02,.01):
        plus=build(eps);minus=build(-eps)
        def variation(points,primitives):
            a=plus.jets(points);b=minus.jets(points)
            return {k:(a[k]-b[k])/(2*eps) for k in ('g','phi','dphi','A','dA')}
        response=q.response(variation)
        wp=loops.LoopSource(plus,2,16).records;wm=loops.LoopSource(minus,2,16).records
        finite=(wp-wm)/(2*eps);pred=np.array(response['total'])
        constraints={sign:value.initial_constraints for sign,value in (('plus',plus),('minus',minus))}
        for sign in constraints:
            for key in ('H_max','M_max','Gauss_max','harmonic_max'):
                assert abs(constraints[sign][key]-base.initial_constraints[key])<1e-12
        assert abs(pred[2])>1e-8 and np.max(abs(finite-pred))<2e-9
        rows.append(dict(epsilon=eps,full_source=response,independently_rebuilt_record_derivative=finite.tolist(),difference=(finite-pred).tolist(),initial_constraints=constraints))
        print('actual constrained family',eps,'source',pred,'difference',finite-pred,flush=True)
    return dict(round=908,date='2026-10-06',all_checks_passed=True,
        exact_family='c=.21(1+epsilon); a^2=(S-.27^2*c^2/4)/(.27^2+c^2/4); S=.008868915; original E/g/phi/pi unchanged initially',
        baseline_initial_constraints=base.initial_constraints,rows=rows,
        full902_equations_evolved_for_each_member=True,fixed_original_material_curve_and_profile=True,
        nonzero_original_physical_direction_numerically_realized=True,
        analytic_nonzero_physical_source_argument_inherited_from_863=True,
        full_Hadamard_variance_evaluated=False,quantum_backreaction_computed=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
