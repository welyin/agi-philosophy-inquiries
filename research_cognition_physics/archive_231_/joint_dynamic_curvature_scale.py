"""Round 550: on-shell curvature bounds for the shared bare spectral branch.

The conformal trace identity is inherited mathematics. New checks join it to
the shared spectral moments, gauge normalization and round-549 local solution.
No assertion of a complete heat-kernel error bound or full spectral solution.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_reference_gravity_constraints as previous

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_dynamic_curvature_scale_results.json'


def exact_trace_check(phi, derivatives):
    """Contract the complete nonminimal source, including Hessian F."""
    eta=[Q(-1),Q(1),Q(1),Q(1)]
    m02,c0,v0,lh,p,ls=Q(3,2),Q(1,4),Q(2),Q(1,3),Q(2,7),Q(4,5)
    h,s=phi
    v=v0-c0*(h*h+s*s)/2+(lh*h**4+2*p*h*h*s*s+ls*s**4)/4
    dv=[-c0*h+lh*h**3+p*h*s*s,-c0*s+ls*s**3+p*s*h*h]
    radius2=h*h+s*s; f=m02-radius2/6
    r=(4*v-sum(phi[i]*dv[i] for i in range(2)))/m02
    boxes=[dv[i]+phi[i]*r/6 for i in range(2)]
    hess=[[[Q(0) for _ in range(4)] for _ in range(4)] for _ in range(2)]
    for i in range(2):
        hess[i][1][1]=Q(i+1,7)
        hess[i][3][3]=Q(i-2,9)
        hess[i][0][0]=sum(hess[i][k][k] for k in (1,2,3))-boxes[i]
    kinetic=sum(eta[k]*derivatives[i][k]**2 for i in range(2) for k in range(4))
    hf=[[-sum(derivatives[i][a]*derivatives[i][b]+phi[i]*hess[i][a][b]
              for i in range(2))/3 for b in range(4)] for a in range(4)]
    boxf=sum(eta[k]*hf[k][k] for k in range(4))
    rhs=[]
    for a in range(4):
        rhs.append(sum(derivatives[i][a]**2 for i in range(2))
                   -eta[a]*(kinetic/2+v)+hf[a][a]-eta[a]*boxf)
    trace=sum(eta[k]*rhs[k] for k in range(4))
    assert trace==-f*r
    assert m02*r==4*v0-c0*radius2
    return dict(phi=[str(a) for a in phi],kinetic=str(kinetic),R=str(r),
                full_source_trace=str(trace),minus_F_R=str(-f*r))


def curvature_from_coupled_jets(phi,eps,pars):
    data,grad,hess,jj,_,_=previous.origin_jets(phi,eps,pars)
    f,df=data['F'],data['dF']
    contract=grad@previous.ETA@grad.T
    boxes=np.einsum('mn,imn->i',previous.ETA,hess)
    dlog=df/f
    ddlog=-np.eye(2)/(3*f)-np.outer(df,df)/f**2
    boxlog=np.einsum('ij,ij',ddlog,contract)+dlog@boxes
    dlogsq=dlog@contract@dlog
    re=np.einsum('ij,ij',data['G'],contract)+4*data['U']
    rj=f*(re+3*boxlog-1.5*dlogsq)
    boxes_j=f*boxes-contract@df
    scalar_residual=boxes_j-data['dV']-np.asarray(phi)*rj/6
    boxf_j=-(np.asarray(phi)@boxes_j+f*np.trace(contract))/3
    metric_trace_residual=f*rj-f*np.trace(contract)-4*data['V']-3*boxf_j
    trace=(4*data['V']-np.asarray(phi)@data['dV'])/float(pars['M02'])
    return dict(R_J=float(rj),trace_R=float(trace),F=f,R_J_over_F=float(rj/f),
                scalar_residual=float(np.max(abs(scalar_residual))),
                metric_trace_residual=float(abs(metric_trace_residual)),
                reference_determinant=float(np.linalg.det(jj)))


def run():
    checks=[]
    phi=[Q(1,2),Q(2,3)]
    exact=[]
    for grad in (([Q(1),Q(0),Q(1,5),Q(0)],[Q(2),Q(1),Q(0),Q(1,3)]),
                 ([Q(0),Q(2),Q(1),Q(-1)],[Q(1,7),Q(0),Q(3),Q(0)])):
        exact.append(exact_trace_check(phi,grad))
    assert exact[0]['kinetic']!=exact[1]['kinetic'] and exact[0]['R']==exact[1]['R']
    checks.append('exact_full_nonminimal_tensor_trace_cancels_different_gradient_contractions')

    pars=previous.bare_parameters(); rows=[]; error=0.
    for state in ([1.,np.sqrt(1/3)],[.7,.4],[1.2,-.2]):
        for eps in (.5,.125,.03125):
            row=curvature_from_coupled_jets(state,eps,pars)
            row.update(fields=state,epsilon=eps)
            error=max(error,abs(row['R_J']-row['trace_R']),row['scalar_residual'],row['metric_trace_residual'])
            assert row['reference_determinant']!=0
            rows.append(row)
    assert error<1e-12
    assert abs(rows[0]['R_J']-2.75)<1e-14 and abs(rows[0]['R_J_over_F']-2.475)<1e-14
    checks.append('independent_Einstein_to_Jordan_curvature_and_both_Jordan_field_equations')

    certificates=[]
    for chi,delta,gw2,mu2 in ((Q(2),Q(5,6),Q(3,16),Q(1,16)),
                             (Q(1),Q(1,50),Q(2,5),Q(3,8)),
                             (Q(7,6),Q(3,4),Q(1,3),Q(9,4))):
        m02=4*mu2/gw2; c0=4*mu2
        # n/(24T)=1/gw2 from the same finite weak representation.
        v0=6*chi*mu2*m02
        radius2=6*m02*(1-delta); f=m02-radius2/6
        r=(4*v0-c0*radius2)/m02
        assert r/mu2==24*(chi-1+delta)
        assert r/f==6*gw2*(1+(chi-1)/delta) and r/f>=6*gw2
        certificates.append(dict(chi=str(chi),delta=str(delta),gw_squared=str(gw2),
                                 R_over_mu_squared=str(r/mu2),R_over_F=str(r/f)))
    checks.append('exact_joint_moment_and_gauge_curvature_bounds_including_near_F_boundary')

    vacua=[]
    for ng,r,q,s in ((1,Q(1),Q(2,5),Q(1,10)),(3,Q(7,3),Q(1,4),Q(3,28)),
                     (5,Q(4),Q(1,3),Q(1,6))):
        t=3*q+r*s; c0=Q(1,4); n=8*ng*(3+r)
        m02=n*c0/(24*t); h2=c0/q; s2=c0*r/t*(1-s/q)
        delta=1-(h2+s2)/(6*m02)
        assert delta==1-Q(1,2*ng)
        # These are field values at a flat-potential minimum, NOT static GR vacua.
        lower=24*delta
        assert lower>=12
        vacua.append(dict(generations=ng,delta=str(delta),sharp_moment_ratio_infimum=str(lower)))
    checks.append('flat_potential_minimum_field_values_fix_delta_without_claiming_constant_GR_vacuum')

    rescalings=[]
    f0,f2,f4,lam2=Q(1),Q(1,2),Q(1,2),Q(1,8)
    for b in (Q(1,3),Q(2),Q(10)):
        new2,new4,newlam=f2/b**2,f4/b**4,lam2*b**2
        assert new2*newlam/f0==f2*lam2/f0
        assert f0*new4/new2**2==f0*f4/f2**2
        rescalings.append(dict(argument_scale=str(b),cutoff_squared=str(newlam),
            moment_scale_squared=str(new2*newlam/f0),chi=str(f0*new4/new2**2)))
    for k in (Q(1,10),Q(10),Q(100)):
        m02,c0,radius2,v0=Q(4,3)*k*k,Q(1,4)*k*k,Q(4,3)*k*k,k**4
        r=(4*v0-c0*radius2)/m02; f=m02-radius2/6
        assert r/f==Q(99,40) and r/(c0/4)==44
    checks.append('cutoff_argument_reparametrization_and_common_mass_rescaling_do_not_change_bounds')

    limits=[]
    for p in (2,8,32,128,512):
        f2p=math.gamma(1/p)/(2*p); f4p=math.gamma(2/p)/(2*p)
        chi=f4p/f2p**2; delta=1/p; gw2=.25
        ratio=24*(chi-1+delta); planck_ratio=6*gw2*(1+(chi-1)/delta)
        assert chi>=1 and planck_ratio>=6*gw2
        limits.append(dict(smooth_cutoff_power=p,chi=chi,relative_F=delta,
            R_over_moment_scale_squared=ratio,R_over_F=planck_ratio))
    assert limits[-1]['R_over_moment_scale_squared']<.05
    assert limits[-1]['R_over_F']>1.5
    checks.append('smooth_sharp_cutoff_and_relative_F_limit_keeps_nonzero_fixed_gauge_bound')

    # Inherited 533 smooth positive NON-monotone example f=(1+10u^2)e^-u^2.
    chi=Q(42,121); nonmonotone=dict(pars)
    nonmonotone['V0']=pars['n']*chi*pars['C0']**2/(16*pars['T'])
    nonmonotone['cutoff_squared']=pars['C0']/22  # f2/f0=11/2; keep the same C0
    radius2=6*float(pars['M02'])*float(chi)
    state=[math.sqrt(radius2*.6),math.sqrt(radius2*.4)]
    nonmono=curvature_from_coupled_jets(state,.125,nonmonotone)
    assert nonmono['F']>0 and abs(nonmono['R_J'])<1e-12
    assert nonmono['reference_determinant']!=0
    assert float(chi)<1
    nonmono.update(cutoff='(1+10u^2)exp(-u^2)',chi=str(chi),field_values=state)
    checks.append('dropping_monotone_cutoff_allows_zero_scalar_curvature_with_positive_F_and_reference_rank')

    shifted=dict(pars); shifted['V0']-=Q(11,12)
    counterterm=curvature_from_coupled_jets([1.,math.sqrt(1/3)],.125,shifted)
    assert counterterm['F']>0 and abs(counterterm['R_J'])<1e-12
    assert counterterm['reference_determinant']!=0
    counterterm.update(delta_V0='-11/12',new_V0=str(shifted['V0']))
    checks.append('independent_vacuum_matching_evades_bound_and_is_an_explicit_changed_input')

    deps=('research_note_533.md','research_note_543.md','research_note_549.md',
          'joint_reference_gravity_constraints.py','joint_reference_gravity_constraints_results.json',
          'research_round_549_checks.json')
    return dict(round=550,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_tensor_trace_examples=exact,coupled_examples=rows,maximum_frame_and_equation_residual=error,
        exact_joint_bound_certificates=certificates,bare_minimum_field_value_certificates=vacua,
        equivalent_cutoff_parameterizations=rescalings,smooth_boundary_sequence=limits,
        nonmonotone_scope_counterexample=nonmono,independent_vacuum_scope_counterexample=counterterm,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(fixed_four_dimensional_classical_nonminimal_branch=True,
            trace_identity_and_cutoff_moment_inequality_inherited=True,
            new_on_shell_pointwise_joint_bounds=True,same_scale_bare_gauge_normalization=True,
            fixed_nonzero_gauge_strength_needed=True,constant_vacuum_not_silently_subtracted=True,
            R_over_F_is_Jordan_coefficient_ratio_not_measured_Newton_constant=True,
            moment_scale_not_arbitrary_sharp_spectral_threshold=True,
            weak_curvature_only_a_necessary_diagnostic_not_full_error_control=True,
            extra_Weyl_squared_preserves_necessary_trace_not_549_full_solutions=True,
            generic_higher_order_and_quantum_corrections_not_excluded=True,
            zero_scalar_curvature_not_zero_Riemann_tensor=True,
            no_general_spectral_action_no_go=True,completed_unification=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else: assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','maximum_frame_and_equation_residual')}))
