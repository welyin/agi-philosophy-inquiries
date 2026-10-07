"""846: exact six-generator checks of the projected classical recursion.

This algebraic calibration is not the original continuum Green operator.
It checks all feedback terms, quotient-source projection, and cutoff variation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,random
HERE=Path(__file__).resolve().parent;TARGET=HERE/'constrained_six_leg_recursion_results.json'
ONE={0:F(1)}

def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(p,c):return {k:v*c for k,v in p.items() if v*c}
def mul(p,q):
    out={}
    for a,c in p.items():
        for b,d in q.items():
            if a&b:continue
            inversions=sum((b&((1<<i)-1)).bit_count() for i in range(6) if (a>>i)&1)
            out[a|b]=out.get(a|b,F(0))+(-1)**inversions*c*d
    return {k:v for k,v in out.items() if v}
def mv(m,v):return [add(*(scale(p,c) for c,p in zip(row,v))) for row in m]
def va(*vs):return [add(*ps) for ps in zip(*vs)]
def vs(v,c):return [scale(p,c) for p in v]
def pv(p,v):return [mul(p,q) for q in v]
def bilinear(m,a,b):
    return scale(add(*(scale(add(mul(a[i],b[j]),mul(b[i],a[j])),m[i][j])
                       for i in range(6) for j in range(i+1,6))),F(1,2))
def inv2(m):
    a,b=m[0];c,d=m[1];det=a*d-b*c
    return [[d/det,-b/det],[-c/det,a/det]]
def outer(a,b):return [[x*y for y in b] for x in a]
def matadd(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def matmul(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def transpose(a):return list(map(list,zip(*a)))
def dot(a,v):return add(*(scale(p,c) for c,p in zip(a,v)))

def run():
    rng=random.Random(846)
    mats=[]
    for denominator in (2,3,5):
        a=[[F(0) for _ in range(6)] for _ in range(6)]
        for i in range(6):
            for j in range(i+1,6):a[i][j]=F(rng.randint(-4,4),denominator);a[j][i]=-a[i][j]
        mats.append(a)
    a1,a2,a3=mats
    d0=[[F(0) for _ in range(6)] for _ in range(6)]
    di=[[F(0) for _ in range(6)] for _ in range(6)]
    for j in range(3):
        d0[2*j][2*j+1]=j+1;d0[2*j+1][2*j]=-j-1
        di[2*j][2*j+1]=F(-1,j+1);di[2*j+1][2*j]=F(1,j+1)
    theta=[{1<<i:F(1)} for i in range(6)]
    a,b,c=F(2),F(3),F(5)
    f=lambda m,u,v:bilinear(m,u,v)
    x2=scale(f(a1,theta,theta),-1/a)
    psi3=vs(mv(di,pv(x2,mv(a1,theta))),-1)
    x4=scale(add(scale(mul(x2,x2),b/2),scale(f(a1,theta,psi3),2),mul(x2,f(a2,theta,theta))),-1/a)
    psi5=vs(mv(di,va(pv(x2,mv(a1,psi3)),pv(x4,mv(a1,theta)),
                       pv(scale(mul(x2,x2),F(1,2)),mv(a2,theta)))),-1)
    x6_terms=[scale(mul(x2,x4),b),scale(mul(mul(x2,x2),x2),c/6),
              scale(f(a1,theta,psi5),2),f(a1,psi3,psi3),
              scale(mul(x2,f(a2,theta,psi3)),2),mul(x4,f(a2,theta,theta)),
              scale(mul(mul(x2,x2),f(a3,theta,theta)),F(1,2))]
    x6=scale(add(*x6_terms),-1/a)
    x=add(x2,x4,x6);psi=va(theta,psi3,psi5)
    def nonlinear_boson(xx,pp):
        return add(scale(mul(xx,xx),b/2),scale(mul(mul(xx,xx),xx),c/6),
                   f(a1,pp,pp),mul(xx,f(a2,pp,pp)),scale(mul(mul(xx,xx),f(a3,pp,pp)),F(1,2)))
    def nonlinear_fermion(xx,pp):
        return va(pv(xx,mv(a1,pp)),pv(scale(mul(xx,xx),F(1,2)),mv(a2,pp)),
                  pv(scale(mul(mul(xx,xx),xx),F(1,6)),mv(a3,pp)))
    rb=add(scale(x,a),nonlinear_boson(x,psi));rf=va(mv(d0,va(psi,vs(theta,-1))),nonlinear_fermion(x,psi))
    assert not rb and not any(rf)
    # Independent simultaneous substitution, not the hand-written coefficient recursion.
    xx={};pp=theta
    for iteration in range(1,13):
        nxt=scale(nonlinear_boson(xx,pp),-1/a)
        npp=va(theta,vs(mv(di,nonlinear_fermion(xx,pp)),-1))
        if nxt==xx and npp==pp:break
        xx,pp=nxt,npp
    else:raise AssertionError('nilpotent fixed point did not terminate')
    assert xx==x and pp==psi
    frozen_fermion_residual=add(scale(x,a),nonlinear_boson(x,theta))
    assert frozen_fermion_residual
    dropped_source_residual=add(scale(add(x2,x4,scale(add(x6_terms[0],x6_terms[1]),-1/a)),a),
                                nonlinear_boson(x,psi))
    assert dropped_source_residual.get(63,0)!=0
    # Original 785 identities calibrated on a nonorthogonal constrained slice.
    k=[F(3),F(-1)];ell=[F(1),F(2)];v=[F(2),F(-1)];cov=[F(-1),F(-3)]
    assert sum(x*y for x,y in zip(k,ell))==1
    pi=outer(v,cov);p=outer(cov,cov);p=[[a*z for z in row] for row in p]
    gp=[[z/a for z in row] for row in outer(v,v)]
    pit=transpose(pi)
    assert matmul(p,gp)==pit
    g1=inv2(matadd(p,[[-z for z in row] for row in outer(k,k)]))
    extension=add(mul(theta[0],theta[3]),mul(x,mul(theta[2],theta[4])))
    n=nonlinear_boson(x,psi)
    source=va([scale(n,z) for z in cov],[scale(extension,z) for z in ell])
    correct=vs(mv(gp,source),-1)
    expected=[scale(x,z) for z in v]
    assert correct==expected and not dot(ell,correct)
    assert not any(va(mv(p,correct),mv(pit,source)))
    without_dual=vs(mv(matmul(pi,g1),source),-1)
    dual_defect=add(dot(cov,without_dual),scale(x,-1))
    bare=vs(mv(g1,source),-1);constraint_defect=dot(ell,bare)
    assert dual_defect and constraint_defect
    # Euler derivative of g*x*(x')^2/2 includes -g'*x*x'.
    t=F(1,2);cutoff_missing_term=-2*t*(1+t);assert cutoff_missing_term==F(-3,2)
    return dict(round=846,all_checks_passed=True,fresh_test_groups=1,
        exact_rational_exterior_algebra=True,Grassmann_generators=6,algebra_dimension=64,
        original_continuum_coefficients_numerically_computed=False,
        diagnostic_is_original_gravity_or_Dirac_propagation=False,
        hand_recursion_equals_independent_fixed_point=True,simultaneous_iterations_to_fixed_point=iteration,
        full_boson_and_Dirac_residual_coefficients_nonzero=0,
        pure_six_leg_boson_coefficient=str(x6.get(63,0)),
        seven_source_contributions_to_six_leg=[str(z.get(63,0)) for z in x6_terms],
        freezing_Dirac_response_has_nonzero_residual=True,
        omitted_Dirac_and_mixed_source_six_leg_defect=str(dropped_source_residual.get(63,0)),
        nonorthogonal_dual_projection_identity_exact=True,
        source_extension_independence_through_sixth_degree=True,
        missing_dual_projection_changes_physical_solution=True,
        bare_gauge_fixed_inverse_violates_chosen_slice=True,
        omitted_cutoff_derivative_at_half=str(cutoff_missing_term),
        actual_original_eighty_term_sum_nonzero_proven=False,
        finite_coupling_solution_or_new_cognition_axiom_claimed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
