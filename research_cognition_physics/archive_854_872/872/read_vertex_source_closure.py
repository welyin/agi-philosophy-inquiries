"""872: actual read-vertex density derivatives and four-leg source closure.
Finite nilpotent equations calibrate the complete recurrence, not the original
curved PDE. Its conditional existence/source theorem is analytic in the note.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,random,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'846'))
import constrained_six_leg_recursion as alg
TARGET=HERE/'read_vertex_source_closure_results.json'

def density_check():
    old=json.loads((HERE.parent/'859/probe_reference_cone_compatibility_results.json').read_text('utf-8'))
    bg=[row for row in old['background_rows'] if row['N']==24 and row['probe_amplitude']==.02][0]
    psi=bg['point_psi'];g=np.diag([-1.,psi**4,psi**4,psi**4]);gi=np.linalg.inv(g)
    V=np.array([[.07,.02,0,0],[.02,.11,.01,0],[0,.01,-.05,0],[0,0,0,.08]])
    Z=np.array([[-.09,0,.02,0],[0,.12,0,.02],[.02,0,.04,0],[0,.02,0,.06]])
    c=.75;o=1.2;bv=.7;bz=-.3
    nuv=np.trace(gi@V)/2;nuz=np.trace(gi@Z)/2
    wanted=-c*o*(bv*nuz+bz*nuv)
    def action(s,t):
        gg=g+s*V+t*Z
        assert np.linalg.det(gg)<0
        return -c*o*np.sqrt(np.linalg.det(gg)/np.linalg.det(g))*(s*bv+t*bz)
    rows=[]
    for e in (.004,.002,.001):
        got=(action(e,e)-action(e,-e)-action(-e,e)+action(-e,-e))/(4*e*e)
        rows.append(dict(step=e,mixed_vertex=float(got),error=float(abs(got-wanted))))
    assert rows[-1]['error']<2e-8 and rows[-1]['error']<rows[0]['error']/10
    assert abs(wanted)>.01
    return dict(original859_point_psi=psi,analytic_mixed_density_vertex=float(wanted),
        finite_difference_rows=rows,wrong_frozen_density_vertex=0.,
        omitted_density_defect=float(abs(wanted)),
        scope='exact local derivative evaluated at the original metric value; directions are calibration inputs')

def rational_source_check():
    rng=random.Random(872)
    def skew(den):
        a=[[F(0) for _ in range(6)] for _ in range(6)]
        for i in range(6):
            for j in range(i+1,6):
                a[i][j]=F(rng.randint(-4,4),den);a[j][i]=-a[i][j]
        return a
    a0,b1,a20=skew(3),skew(5),skew(7)
    kappa=F(2,5)
    b2=[[2*kappa*x for x in row] for row in b1]
    theta=[{1<<i:F(1)} if i<4 else {} for i in range(6)]
    di=[[F(0) for _ in range(6)] for _ in range(6)]
    for i in range(3):di[2*i][2*i+1]=F(-1,i+1);di[2*i+1][2*i]=F(1,i+1)
    a=F(2);g=F(3)
    f=alg.bilinear
    x20=alg.scale(f(a0,theta,theta),-1/a)
    x21=alg.scale(f(b1,theta,theta),-1/a)
    psi31=alg.vs(alg.mv(di,alg.va(alg.pv(x21,alg.mv(a0,theta)),
                               alg.pv(x20,alg.mv(b1,theta)))),-1)
    psi32=alg.vs(alg.mv(di,alg.pv(x21,alg.mv(b1,theta))),-1)
    terms=[
        alg.scale(alg.mul(x21,x21),g/2),
        alg.scale(f(a0,theta,psi32),2),
        alg.scale(f(b1,theta,psi31),2),
        alg.mul(x21,f(b2,theta,theta))]
    wanted=alg.add(*terms)
    assert wanted and terms[-1]
    def solve(lam):
        a1=alg.matadd(a0,[[lam*x for x in row] for row in b1])
        a2=alg.matadd(a20,[[lam*x for x in row] for row in b2])
        xx={};pp=theta
        for _ in range(9):
            newx=alg.scale(alg.add(alg.scale(alg.mul(xx,xx),g/2),f(a1,pp,pp),
                                  alg.mul(xx,f(a2,pp,pp))),-1/a)
            newp=alg.va(theta,alg.vs(alg.mv(di,alg.va(alg.pv(xx,alg.mv(a1,pp)),
                      alg.pv(alg.scale(alg.mul(xx,xx),F(1,2)),alg.mv(a2,pp)))),-1))
            if newx==xx and newp==pp:break
            xx,pp=newx,newp
        else:raise AssertionError('nilpotent fixed point failed')
        eb=alg.add(alg.scale(xx,a),alg.scale(alg.mul(xx,xx),g/2),f(a1,pp,pp),alg.mul(xx,f(a2,pp,pp)))
        assert not eb
        q4={mask:v for mask,v in xx.items() if mask.bit_count()==4}
        return alg.scale(q4,-a)
    minus,zero,plus=solve(F(-1)),solve(F(0)),solve(F(1))
    extracted=alg.scale(alg.add(plus,minus,alg.scale(zero,-2)),F(1,2))
    assert wanted==extracted
    def poly(p):return {str(k):str(v) for k,v in sorted(p.items())}
    deletions=[]
    for i in range(4):
        truncated=alg.add(*(t for j,t in enumerate(terms) if i!=j))
        residual=alg.add(extracted,alg.scale(truncated,-1))
        deletions.append(poly(residual))
    assert all(deletions)
    return dict(active_exterior_generators=4,fermion_components_retained=6,
        read_density_parameter=str(kappa),
        source_parts=[poly(t) for t in terms],complete_source=poly(wanted),
        independent_fixed_point_lambda_squared_source=poly(extracted),
        all_exact_residuals_zero=True,omitted_part_residuals=deletions,
        scope='exact finite nilpotent Euler calibration; no original continuum source value is inferred')

def run():
    return dict(round=872,date='2026-10-06',formal_reports=872,
        cumulative_numbered_groups=3657,fresh_numbered_groups=1,all_checks_passed=True,
        original_metric_read_density_check=density_check(),
        four_leg_source_calibration=rational_source_check(),
        original_source_construction='full inherited E/F/D derivatives plus explicit read density and Dirac mass derivatives',
        conditional_relative_mean_scope='lambda_alpha^2 hbar^2, original regular local patch and common causal/Cauchy contract',
        extra_independent_second_read_vertex=False,
        full_projected_source_and_relative_mean_defined=True,
        original_continuous_total_source_numerically_evaluated=False,
        original_total_mean_difference_proved_nonzero=False,
        finite_coupling_or_full_nonlinear_semiclassical_geometry=False,
        full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
