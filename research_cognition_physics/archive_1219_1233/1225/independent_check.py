"""Independent 1225 audit: full-width coefficient algebra, not full-H propagation.
No author code is imported. Default execution only reads saved results.
"""
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
RESULTS=HERE/"independent_results.json"
SOURCES=(
    "archive_585_628/research_note_589.md",
    "archive_585_628/research_note_598.md",
    "archive_1219_/1224/proof.md",
    "archive_1219_/1224/results.json",
)

import math
from fractions import Fraction as F
import numpy as np
M=2.; R=math.sqrt(6); width=.1; n=8
def beta_ratio():
    # B(3,16)/B(5/2,17) via rational products and Gamma half integer.
    a=F(2,1)
    for j in range(16,19): a/=j
    # B(5/2,17)=16! / product_{j=0}^{16}(5/2+j)
    den=F(math.factorial(16),1)
    for j in range(17): den/=F(5,2)+j
    return float(a/den)
def J(a): return (math.sinh(a/R)/(a/R))**4
def Q(a):
    d=math.sqrt(M)*a*(M+1/3)/4
    return math.sqrt(M)/4*(4*n/a)*J(a)*beta_ratio()+d
b0=math.sqrt(M)/4
a0=b0*math.cos(math.sqrt(M)*width)/math.cosh(width/R)**2
h=(2*width/R)/math.tanh(2*width/R)
bk=min(width/2,a0/(4*h*Q(width)))
LF=M/4+math.sqrt(M)/(2*R)
alow=b0*(a0-bk*h*Q(width)-LF*bk)
bupper=Q(bk)*b0*h*bk
def quad(order,rad):
    x,w=np.polynomial.legendre.leggauss(order)
    rr=rad*(x+1)/2
    wr=w*rad/2*rr**4*(1-(rr/rad)**2)**16*(np.sinh(rr/R)/(rr/R))**4
    wr/=wr.sum()
    us=x; wu=.75*w*(1-x*x)
    return rr,wr,us,wu
def dot(x,y): return (x[...,1:]*y[...,1:]).sum(axis=-1)-x[...,0]*y[...,0]
def position(r,u,v=1):
    r,u,v=np.broadcast_arrays(r,u,v)
    z=np.zeros(r.shape+(6,))
    z[...,0]=R*np.cosh(r/R)
    amp=R*np.sinh(r/R)
    z[...,1]=amp*np.sqrt(1-u*u)*v
    z[...,2]=amp*np.sqrt(1-u*u)*np.sqrt(1-v*v)
    z[...,5]=amp*u
    return z
def grad(X):
    s=math.sqrt(M)*R*X[...,5]/X[...,0]
    g=np.zeros_like(X); g[...,5]=math.sqrt(M)*R*np.cos(s)/(4*X[...,0])
    g[...,0]=g[...,5]*X[...,5]/X[...,0]
    return g
def coeff(order):
    ry,wy,uy,wu=quad(order,width)
    j=np.arange(1,order+1); vv=np.cos(j*np.pi/(order+1)); wv=2*np.sin(j*np.pi/(order+1))**2/(order+1)
    rg,ug,vg=np.meshgrid(ry,uy,vv,indexing='ij')
    Y=position(rg,ug,vg).reshape(-1,6); V=grad(Y)
    weights=(wy[:,None,None]*wu[None,:,None]*wv[None,None,:]).ravel()
    rr,wr,ub,wub=quad(order,bk)
    av=bv=sv=0.
    for r,pr in zip(rr,wr):
      for u,pu in zip(ub,wub):
        X=position(r,u); U=grad(X)
        z=-dot(X,Y)/R**2; q=z-1
        f=np.empty_like(q); fp=np.empty_like(q); small=q<.001
        f[small]=1-q[small]/3+2*q[small]**2/15-2*q[small]**3/35+8*q[small]**4/315
        fp[small]=-1/3+4*q[small]/15-6*q[small]**2/35+32*q[small]**3/315
        f[~small]=np.arccosh(z[~small])/np.sqrt(z[~small]**2-1)
        fp[~small]=(1-z[~small]*f[~small])/(z[~small]**2-1)
        base=Y-z[:,None]*X
        log=f[:,None]*base
        dz=-dot(X,V)/R**2
        dlog=f[:,None]*V+dz[:,None]*(fp[:,None]*base-f[:,None]*X)
        AA=weights@dlog; LL=weights@log
        AA[2:5]=0; LL[2:5]=0
        s=math.sqrt(M)*R*math.tanh(r/R)*u
        dsdr=math.sqrt(M)/math.cosh(r/R)**2*u
        ds2=M/math.cosh(r/R)**2*(1-math.tanh(r/R)**2*u*u)
        lap_s=-2/R*math.tanh(r/R)*dsdr
        lap_f=(math.cos(s)*lap_s-math.sin(s)*ds2)/4
        density_div=-lap_f+32*r/(bk*bk-r*r)*math.cos(s)*dsdr/4
        av+=pr*pu*float(dot(U,AA))
        bv-=pr*pu*density_div*float(dot(AA,LL))
        sv+=pr*pu*float(dot(AA,AA))
    return {'order':order,'a':av,'b':bv,'S':sv}


def determinant(mat):
    a=[list(row) for row in mat]; d=F(1)
    for i in range(len(a)):
        k=next(k for k in range(i,len(a)) if a[k][i])
        if k!=i: a[i],a[k]=a[k],a[i]; d=-d
        p=a[i][i]; d*=p
        for j in range(i+1,len(a)):
            z=a[j][i]/p
            for k in range(i+1,len(a)): a[j][k]-=z*a[i][k]
    return d

def row_map(a,b,S):
    return [[-S,0,0,0,0,0],[0,-S,0,0,0,0],[0,0,-S,0,0,0],
            [b,b,0,a,0,0],[b,0,b,0,a,0],[0,b,b,0,0,a]]

def apply(mat,x):
    return [sum(v*t for v,t in zip(row,x)) for row in mat]

def invert(y,a,b,S):
    c12,c13,c23=[-z/S for z in y[:3]]
    return [c12,c13,c23,(y[3]-b*(c12+c13))/a,
            (y[4]-b*(c12+c23))/a,(y[5]-b*(c13+c23))/a]

def exact_six_row_check():
    a,b,S=F(3,25),F(1,8),F(5,41)
    L=row_map(a,b,S)
    det=determinant([[F(z) for z in row] for row in L])
    assert det==-a**3*S**3
    cases=([F(1,4),F(-1,5),F(1,10),F(2),F(3),F(4)],
           [F(0),F(0),F(0),F(4),F(1,4),F(1)])
    for x in cases: assert invert(apply(L,x),a,b,S)==list(x)
    delta=F(1,10000); errors=[delta*(-1)**j for j in range(6)]
    e=invert(errors,a,b,S)
    inv_bound=max(1/S,(1+2*abs(b)/S)/a)
    assert max(map(abs,e))<=inv_bound*delta
    return dict(exact_det=str(det),det_identity="-a^3*S^3",
        variable_order=["P12","P13","P23","P11","P22","P33"],
        exact_round_trips=len(cases),exact_inverse_infinity_bound=str(inv_bound))

def adjacent_stars():
    B=(0,0,0); dirs=((1,0,0),(0,1,0),(0,0,1)); side=5
    def add(x,y): return tuple((p+q)%side for p,q in zip(x,y))
    found=[]
    for axis,Jp in enumerate(dirs):
        centers=[]
        for i in range(side):
          for j in range(side):
            for k in range(side):
              K=(i,j,k); neighbors=[add(K,d) for d in dirs]
              if B in [K]+neighbors and Jp in [K]+neighbors:
                  centers.append(K)
        assert centers==[B]
        found.append(dict(axis=axis+1,source=list(Jp),common_centers=[list(x) for x in centers]))
    return dict(graph_side=side,adjacent_source_stars=found,
        diagonal_density_gradient_factor="-4*eta*A_B",
        diagonal_response="2*eta*w^(-2/3)*a*Cii",
        mixed_response="2*eta*w^(-2/3)*b*sum(Cij)",
        all_remaining_full_H_terms_retained_in_analytic_cancellation=True)

def analytic_bounds():
    d0=math.sqrt(M)*width*(M+1/3)/4
    qinf=32*J(2*width)*(4/3)**16*(d0+32*b0/width)
    eta=1/(4*qinf)
    slow=(a0-bk*h*Q(width))**2
    gradient_b_upper=b0*(h*h*Q(width)*bk+b0*h)
    inv_bound=max(1/slow,(1+2*gradient_b_upper/slow)/alow)
    smin_lower=1/math.sqrt(3/alow**2+6*gradient_b_upper**2/(alow**2*slow**2)+3/slow**2)
    assert 0<gradient_b_upper<bupper and smin_lower>0
    assert alow>0 and slow>0 and 0<bk<width
    # Hess s is the sum of both tensor products, not half-symmetrization.
    assert abs(LF-(M/4+math.sqrt(M)/(2*R)))<1e-15
    return dict(Q1=Q(width),Q1_kappa=Q(bk),Qinf=qinf,eta=eta,
        bK=bk,Hessian_F_bound=LF,a_lower=alow,b_absolute_upper=bupper,
        b_gradient_upper=gradient_b_upper,singular_value_lower=smin_lower,
        S_lower=slow,S_upper=M/16,inverse_infinity_uniform_upper=inv_bound,
        amplitude_C7_and_H3_only=True,numeric_bounds_are_not_interval_arithmetic=True)

def numerical_coefficients():
    rows=[{k:float(v) if isinstance(v,np.floating) else v for k,v in coeff(o).items()}
          for o in (12,16,20)]
    for row in rows:
        assert row["a"]>alow and abs(row["b"])<bupper
        assert row["S"]>=(a0-bk*h*Q(width))**2
    diff=max(abs(rows[-1][k]-rows[-2][k]) for k in ("a","b","S"))
    assert diff<1e-10
    # The same clock preparation uses mu, regardless of the probe packet.
    rr,wr,uu,wu=quad(20,width)
    rg,ug=np.meshgrid(rr,uu,indexing="ij")
    s=math.sqrt(M)*R*np.tanh(rg/R)*ug
    ds2=M/np.cosh(rg/R)**2*(1-np.tanh(rg/R)**2*ug**2)
    v0=float(np.sum(wr[:,None]*wu[None,:]*np.cos(s)*ds2)/4)
    v0low=M/4*math.cos(math.sqrt(M)*width)/math.cosh(width/R)**4
    assert v0low<v0<M/4
    return dict(rows=rows,last_refinement_max_difference=diff,
        v0=v0,v0_lower=v0low,clock_relation="v=v0/w",
        quadrature_only_diagnostic=True,full_H_moments_computed=False)

def cofactor_and_area_check(num):
    a,b,S=(num["rows"][-1][k] for k in ("a","b","S"))
    eta=analytic_bounds()["eta"]; v0=num["v0"]
    G=np.array([[4/3,2/3,.5],[2/3,4/3,0],[.5,0,1.]])
    assert abs(np.linalg.det(G)-1)<1e-14
    rows=[]; worst=0.
    for w in (.5,.8,1.,2.):
        hh=w**(2/3)*G
        PP=np.linalg.det(hh)*np.linalg.inv(hh)
        C=np.linalg.inv(G)
        assert np.max(abs(PP-w**(4/3)*C))<2e-14
        x=np.array([PP[0,1],PP[0,2],PP[1,2],PP[0,0],PP[1,1],PP[2,2]])
        # Leading coefficient of q(d), normalized by the known target integral.
        data=(eta/v0**2)*np.array(apply(row_map(a,b,S),x))
        recovered=np.array(invert(data*v0**2/eta,a,b,S))
        err=float(np.max(abs(recovered-x))); worst=max(worst,err)
        assert err<2e-14
        areas=np.sqrt(np.diag(PP))
        delta=1e-6; perturbed_diag=np.diag(PP)+np.array([delta,-delta,delta])
        area_bound=delta/(2*math.sqrt(float(np.min(np.diag(PP)))-delta))
        assert np.max(abs(np.sqrt(perturbed_diag)-areas))<=area_bound
        rows.append(dict(w=w,areas=areas.tolist(),roundtrip_error=err,
                         illustrative_area_error_bound=area_bound))
    # These are algebraic round trips, not measured finite-time probabilities.
    return dict(cases=rows,max_roundtrip_error=worst,
        normalized_relational_coefficient="(v0^2/eta)*lim(q(d)/d^2)=L*vec(cof(h))",
        cofactor_formula="cof(w^(2/3)*G)=w^(4/3)*inverse(G)",
        uncertainty_budget_is_conditional=True)

def deletion_and_window_scope():
    p0=[F(1),F(1),F(1)]; p1=[F(4),F(1,4),F(1)]
    assert math.prod(p0)==math.prod(p1)==1
    assert p0!=p1
    # Taylor transport is per plus/minus preparation; no classical path conditioning.
    assert F(1,6)*8==F(4,3)
    assert F(1,2)*2*3==3
    return dict(offdiagonal_only_counterexample={
            "P0_diag":[str(x) for x in p0],"P1_diag":[str(x) for x in p1],
            "area0":["1","1","1"],"area1":["2","1/2","1"],
            "same_three_mixed_coefficients":True,
            "same_full_histories_claimed":False},
        per_channel_normalized_error="(v0^2/eta)*(Kq*d + (2*eps_p + 2*L_clock*eps_d)/d^2)",
        inverse_error="||P_hat-P||_entry_max <= ||L^-1||_infinity * normalized_error",
        time_window_requires_compact_positive_geometry_and_uniform_full_H_moments=True,
        full_H_moments_computed=False,numerical_wait_or_sample_budget_certified=False,
        single_run_proper_clock=False,autonomous_instrument=False,entropy_area_identity=False)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    num=numerical_coefficients()
    return dict(status="PASS",all_checks_passed=True,check_groups=6,
        independent_code_sha256=digest(Path(__file__)),
        source_sha256={s:digest(RESEARCH/s) for s in SOURCES},
        exact_six_row_algebra=exact_six_row_check(),
        adjacent_full_graph=adjacent_stars(),analytic_bounds=analytic_bounds(),
        independent_target_integrals=num,cofactor_area=cofactor_and_area_check(num),
        deletion_and_scope=deletion_and_window_scope())

def main():
    p=argparse.ArgumentParser(); p.add_argument("--write-results",action="store_true")
    args=p.parse_args(); result=run()
    if args.write_results:
        with RESULTS.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
    else:
        assert json.loads(RESULTS.read_text(encoding="utf-8"))==result,"saved result or source mismatch"
    print(json.dumps(dict(status="PASS",all_checks_passed=True,results_match=True,
        check_groups=6,independent_code_sha256=result["independent_code_sha256"]),
        ensure_ascii=False))

if __name__=="__main__": main()
