"""856: covariant cell subtraction with the existing composite references.

The point jet is off shell for the full Einstein-matter action. The exact
checks concern an added local functional and its complete chain-rule sources.
No continuum renormalization or finite-coupling stability is inferred.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'relational_subtraction_source_bridge_results.json'

def exact_point(y=F(-2),z=F(-3),c=F(3),k=F(1)):
    d=c*c-y*z;assert y<0 and z<0 and d>0
    G=np.zeros((4,4),object);g=G.copy()
    G[0,0]=y;G[1,1]=z;G[0,1]=G[1,0]=c;G[2,2]=G[3,3]=F(1)
    g[0,0]=-z/d;g[1,1]=-y/d;g[0,1]=g[1,0]=c/d;g[2,2]=g[3,3]=F(1)
    assert np.array_equal(G@g,np.eye(4,dtype=int))
    dg=np.zeros((4,4,4),object)
    dg[2,0,0]=-z*z/d**2;dg[2,0,1]=dg[2,1,0]=c*z/d**2;dg[2,1,1]=-c*c/d**2
    dg[3,0,0]=-c*c/d**2;dg[3,0,1]=dg[3,1,0]=c*y/d**2;dg[3,1,1]=-y*y/d**2
    rho=-k*d/y;e2=k*(2*c*c/y**2-z/y);e3=k
    fixed=rho*g;fixed[0,0]-=2*rho/y
    full=fixed.copy();full[0,0]-=2*e2;full[1,1]-=2*e3
    gamma=np.zeros((4,4,4),object)
    for a in range(4):
        for b in range(4):
            for v in range(4):
                gamma[a,b,v]=sum(G[a,j]*(dg[b,j,v]+dg[v,j,b]-dg[j,b,v]) for j in range(4))/F(2)
    def div(T):
        M=G@T
        partial=np.array([F(0),F(0),k*c*c/y**2,k],object)
        return [partial[v]+sum(gamma[m,m,l]*M[l,v]-gamma[l,m,v]*M[m,l] for m in range(4) for l in range(4)) for v in range(4)]
    df=div(fixed);dt=div(full)
    assert df==[0,0,e2,e3] and dt==[0]*4
    # u^mu=-g^{mu0}/sqrt(-y); its squared contraction is rational.
    energy_fixed=sum(G[a,0]*fixed[a,b]*G[b,0] for a in range(4) for b in range(4))/(-y)
    energy_full=sum(G[a,0]*full[a,b]*G[b,0] for a in range(4) for b in range(4))/(-y)
    assert energy_fixed==rho
    return dict(y=str(y),z=str(z),c=str(c),rho_at_fixed_independent_references=str(rho),
        reference_euler_sources=["0","0",str(e2),str(e3)],
        fixed_reference_stress=[[str(v) for v in row] for row in fixed],
        composite_reference_stress=[[str(v) for v in row] for row in full],
        divergence_fixed_reference=[str(v) for v in df],divergence_complete=[str(v) for v in dt],
        clock_frame_contraction_complete=str(energy_full),
        composite_h_and_s_euler_sources_vanish_in_this_jet=True,
        actual_753_solution_or_full_action_energy_positivity_claimed=False)

def quadrature(n):
    v,w=np.polynomial.legendre.leggauss(n)
    xi,ze=np.meshgrid(v,v,indexing='ij');weights=np.outer(w,w)*.2*.2
    y=-2+.2*xi;z=-3+.2*ze;c=3.;d=c*c-y*z
    f=(1-xi**2)**2*(1-ze**2)**2
    fy=-4*xi*(1-xi**2)*(1-ze**2)**2/.2
    fz=-4*ze*(1-ze**2)*(1-xi**2)**2/.2
    rho=-d/y;s=1/np.sqrt(d);e2=2*c*c/y**2-z/y;e3=np.ones_like(y)
    fixed00=rho*(-z/d-2/y);fixed11=rho*(-y/d)
    full00=fixed00-2*e2;full11=fixed11-2*e3
    rows=[]
    for component,deriv,Tf,Tfull,E in [("00",fy,fixed00,full00,e2),("11",fz,fixed11,full11,e3)]:
        def action(e,composite):
            yy=y+e*f if component=="00" else y
            zz=z+e*f if component=="11" else z
            jac=1+e*deriv if composite else np.ones_like(y)
            return float(np.sum(weights*jac**2*np.sqrt(c*c-yy*zz)/yy))
        expected=float(np.sum(weights*(-.5*s*Tfull*f)))
        frozen=float(np.sum(weights*(-.5*s*Tf*f)))
        missing=float(np.sum(weights*s*E*f));assert abs(expected-frozen-missing)<2e-15
        checks=[]
        for h in (.001,.0005,.00025):
            actual=(action(h,True)-action(-h,True))/(2*h)
            wrong=(action(h,False)-action(-h,False))/(2*h)
            checks.append(dict(step=h,full_derivative=actual,frozen_derivative=wrong,
                full_error=abs(actual-expected),frozen_error=abs(wrong-frozen)))
        assert checks[-1]['full_error']<checks[0]['full_error']/10
        assert checks[-1]['full_error']<2e-7 and abs(missing)>.001
        rows.append(dict(metric_component=component,analytic_complete_derivative=expected,
            analytic_frozen_derivative=frozen,missing_composite_source=missing,finite_differences=checks))
    return rows

def run():
    points=[exact_point(),exact_point(F(-9,5),F(-14,5)),exact_point(F(-11,5),F(-16,5))]
    # The material-coordinate, mass-lumped counterpart exactly matches the
    # inherited cell energy. This is not a proof of a spatial continuum limit.
    delta=F(1,2);roots=[F(1),F(2),F(3)];lapses=[F(1),F(6,5),F(4,5)]
    old=sum(N/(delta**3*r) for N,r in zip(lapses,roots))/3
    new=sum(delta**3*N/(r*delta**6) for N,r in zip(lapses,roots))/3
    assert old==new
    q32,q64=quadrature(32),quadrature(64)
    err=max(abs(a[k]-b[k]) for a,b in zip(q32,q64) for k in ('analytic_complete_derivative','analytic_frozen_derivative','missing_composite_source'))
    assert err<3e-14
    return dict(round=856,all_checks_passed=True,fresh_test_groups=1,
        cell_energy_match_exact=str(old),material_cell_spacing=str(delta),
        exact_source_jets=points,quadrature_order=64,quadrature_crosscheck_error=err,
        variational_checks=q64,
        scope='Local covariant completion of the inherited inverse-cell-volume subtraction. Existing references are (h,s,(dh)^2,(ds)^2), and their full metric/scalar variations are retained. The explicit jet is a functional diagnostic, not the original on-shell background.',
        residual_579_divergence_cancelled=False,
        old_graph_to_continuum_proved=False,
        finite_coupling_stability_or_observational_fit_proved=False,
        full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
