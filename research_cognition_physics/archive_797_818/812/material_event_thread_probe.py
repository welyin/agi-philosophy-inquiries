"""812 working: causal event threading using the original 753 color background.

No new source, clock species, action or future spacetime evolution is supplied.
The exact inherited bound is separate from the original initial-grid diagnostic.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
with ResearchRuntime(Layout()).installed():
    import joint_reference_constraint_strata as old
TARGET=HERE/'material_event_thread_probe_results.json'

def run():
    lower=(Q('1.8')*Q('.97')*Q('.14')/Q('.033'))**2
    upper=Q(3,2)**8
    assert lower>upper
    rows=[]
    with ResearchRuntime(Layout()).installed():
        for N in (16,24):
            q,psi,tensor,info,color=old.completed(N,1.)
            f=old.coordinates.fields(N);idx=(0,N//4,N//8)
            ps=float(psi[idx]);phi=q['phi'][idx];momentum=q['p'][idx]
            vel=ps**-6*q['F'][idx]*(momentum-phi*np.dot(phi,momentum)/12)
            vh,vs=float(vel[1]),float(vel[4])
            dh=np.r_[vh,f['dh'][idx]];ds=np.r_[vs,f['ds'][idx]]
            metric=np.diag([-1.,ps**4,ps**4,ps**4]);inverse=np.linalg.inv(metric)
            A=float(dh@inverse@dh);B=float(ds@inverse@ds);C=float(dh@inverse@ds)
            w=inverse@dh/A
            assert A<0 and B<0 and vh>0
            assert abs(dh@w-1)<1e-13 and w@metric@w<0
            assert np.linalg.norm(w[1:])<1e-11
            b=-float(f['ds'][idx][1]);assert b>0
            # dh(V)=1 and ds(V)=0 force V0 and Vy at the original point.
            V=np.array([1/vh,0.,vs/(b*vh),0.])
            lower_norm=float(V@metric@V)
            algebraic_lower=(-B)*ps**4/(b*b*vh*vh)
            assert lower_norm>0 and abs(lower_norm/algebraic_lower-1)<1e-12
            assert abs(dh@V-1)<1e-13 and abs(ds@V)<1e-13
            exact_slope=-(12-phi[4]**2)/(phi[1]*phi[4])
            drift=float(ds@w)
            assert abs(drift-exact_slope)<2e-13
            rows.append(dict(N=N,color_amplitude=1.,psi_at_point=ps,
                h_clock_gradient_squared=A,s_clock_gradient_squared=B,
                h_time_thread_norm_squared=float(w@metric@w),
                all_fixed_s_h_unit_threads_norm_squared_lower_bound=lower_norm,
                required_s_label_drift_per_h=drift,
                source_formula_s_drift=float(exact_slope),
                clock_unit_residual=float(abs(dh@w-1)),
                original_solver_residual_info=info))
    return dict(working_round=812,all_checks_passed=True,
        exact_source_ratio_lower=str(lower),psi_power_eight_upper=str(upper),
        inherited_strict_s_timelike_bound=True,original_nonzero_color_rows=rows,
        fixed_s_h_coordinate_thread_cannot_be_causal_at_original_point=True,
        same_material_gradient_gives_local_timelike_thread=True,
        future_thread_numerically_integrated=False,
        derivatives_of_other_two_labels_computed=False,
        quantum_event_record_or_autonomous_apparatus_constructed=False,
        formal_round_completed=False,new_numbered_test_groups=0)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

