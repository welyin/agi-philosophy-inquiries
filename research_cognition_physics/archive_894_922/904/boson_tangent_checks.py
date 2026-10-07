"""904: independent family derivatives, whole-field JVP, and full coupled tangent flow."""
from pathlib import Path
import argparse,json,time,hashlib,sys
import numpy as np
import coupled_boson_tangent as t
HERE=Path(__file__).resolve().parent;TARGET=HERE/'boson_tangent_checks_results.json'

def checks():
    with t.base.ResearchRuntime(t.base.Layout()).installed():
        import joint_reference_constraint_strata as old
        f=t.Family(8,old);y,u,stats=f.data(with_tangent=True)
        assert max(float(np.max(abs(y[k]-f.base_y[k]))) for k in y)<1e-13
        real=f.model.rhs(y);z=f.analytic.rhs({k:v.astype(complex) for k,v in y.items()})
        extension=max(float(np.max(abs(real[k]-z[k]))) for k in real)
        rng=np.random.default_rng(904);x=f.q['grid'];v={}
        for key,value in y.items():
            profile=np.sin(x[...,0]+2*x[...,1])+.3*np.cos(x[...,2])
            coeff=rng.normal(size=value.shape[3:])*.03
            if key in ('g','gd'):coeff=(coeff+coeff.T)/2
            v[key]=profile.reshape(profile.shape+(1,)*coeff.ndim)*coeff
        j=f.analytic.jvp(y,v);j2=f.analytic.jvp(y,v,1e-18)
        step_error=max(float(np.max(abs(j[k]-j2[k]))) for k in j)
        both=f.analytic.jvp(y,{k:.7*u[k]-.3*v[k] for k in u});ju=f.analytic.jvp(y,u)
        linear=max(float(np.max(abs(both[k]-(.7*ju[k]-.3*j[k])))) for k in j)
        errors=[]
        for h in (2e-3,1e-3,5e-4):
            plus=f.model.rhs(t.base.add(y,v,h));minus=f.model.rhs(t.base.add(y,v,-h))
            errors.append(max(float(np.max(abs((plus[k]-minus[k])/(2*h)-j[k]))) for k in j))
        initial_errors=[]
        for h in (2e-3,1e-3,5e-4):
            plus=f.data(.02+h);minus=f.data(.02-h)
            initial_errors.append(max(float(np.max(abs((plus[k]-minus[k])/(2*h)-u[k]))) for k in u))
        assert extension<3e-13 and step_error<1e-12 and linear<1e-12
        assert errors[-1]<errors[0]/8 and initial_errors[-1]<initial_errors[0]/8
        original=f.model.constraints(y);fields=t.constraints(f.analytic,{k:v.astype(complex) for k,v in y.items()})
        field_error=max(abs(float(np.max(abs(fields[a])))-original[b]) for a,b in (('H','H_max'),('M','M_max'),('Gauss','Gauss_max'),('harmonic','harmonic_max')))
        assert field_error<3e-12
        # Finite difference of independently solved family trajectories, not RHS differences only.
        row,final,tangent=t.evolve(f,8,.01);trajectory=[]
        for h in (2e-3,1e-3,5e-4):
            plus=f.data(.02+h);minus=f.data(.02-h)
            for _ in range(8):plus=t.base.rk4(f.model,plus,.01/8);minus=t.base.rk4(f.model,minus,.01/8)
            err={k:float(np.max(abs((plus[k]-minus[k])/(2*h)-tangent[k]))) for k in tangent}
            trajectory.append(dict(parameter_step=h,errors=err,max_error=max(err.values())))
        assert trajectory[-1]['max_error']<trajectory[0]['max_error']/8
        # p alone with frozen geometry fails even the linearized initial Hamiltonian constraint.
        bad={k:np.zeros_like(a) for k,a in u.items()};bad['phi'][...,5]=u['phi'][...,5]
        invalid=t.tangent_constraints(f.analytic,y,bad)
        assert invalid['H']>.02
        return dict(N=8,analytic_extension_vs_original=extension,complex_step_1e24_vs_1e18=step_error,
            linearity_error=linear,whole_field_central_RHS_errors=errors,initial_family_errors=initial_errors,
            independent_constraint_field_reconstruction_error=field_error,
            initial_tangent_solver=stats,independent_family_trajectory_errors=trajectory,
            probe_only_initial_constraint_defects=invalid,
            original902_sha256=hashlib.sha256(t.SOURCE.read_bytes()).hexdigest())

def run():
    out=dict(round=904,date='2026-10-06',checks=checks());print('Independent full-field and family checks passed.',flush=True)
    rows=[]
    with t.base.ResearchRuntime(t.base.Layout()).installed():
        import joint_reference_constraint_strata as old
        for N in (12,16,24):
            start=time.time();family=t.Family(N,old);row,y,u=t.evolve(family)
            rows.append(row);print('tangent',N,'constraints',row['final_linear_constraints'],'seconds',round(time.time()-start,2),flush=True)
            assert row['induced_H5']>1e-9 and row['final_tangent_max']['A']>1e-8 and row['final_tangent_max']['g']>.001
    for key in ('H','M','harmonic'):
        assert all(b['final_linear_constraints'][key]<a['final_linear_constraints'][key] for a,b in zip(rows,rows[1:])),key
    assert rows[-1]['final_linear_constraints']['Gauss']<1e-10
    out.update(rows=rows,all_checks_passed=True,formal_reports=904,fresh_numbered_groups=1,cumulative_numbered_groups=3689,
        argument_scope='Full104-independent-Cauchy-variable classical tangent evolution of the same902 constrained859 family, verified against independently varied trajectories. Nonzero probe/geometry/matter mixing is retained. No sourced physical Green inverse, receiver-source calculation, or rigorous continuum/time error enclosure yet.',
        all_classical_background_components_linearized=True,original_family_constraints_differentiated=True,
        fixed_material_rule_f_epsilon_differentiated=False,
        physical_retarded_Green_operator_constructed=False,receiver_modes_evolved=False,
        quantum_source_or_backreaction_computed=False,rigorous_time_error_enclosure=False,full_goal_completed=False)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--checks-only',action='store_true');a=p.parse_args()
    r=checks() if a.checks_only else run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif not a.checks_only:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
