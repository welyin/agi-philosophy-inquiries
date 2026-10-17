"""1282: inverse quadratic edge classification with finite normal Gauss tests.

The default run is read-only; --write creates results.json exclusively.
Universal claims are proved in proof.md. No Monte Carlo or singular states.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
I4 = np.eye(4)


def left(q):
    a, b, c, d = q
    return np.array([[a,-b,-c,-d],[b,a,-d,c],[c,d,a,-b],[d,-c,b,a]], dtype=float)


def right(q):
    a, b, c, d = q
    return np.array([[a,-b,-c,-d],[b,a,d,-c],[c,-d,a,b],[d,c,-b,a]], dtype=float)


J = [right(np.eye(4)[i]) for i in (1, 2, 3)]
L = [left(np.eye(4)[i]) for i in (1, 2, 3)]


def norm(a):
    return float(np.linalg.norm(a, 2))


def block(a, d, b0, bv, u):
    b = b0*I4 + sum((float(c)*j for c, j in zip(bv,J)), start=np.zeros((4,4)))
    b = b @ u
    return np.block([[a*I4,-b],[-b.T,d*I4]])


def yz_to_x(u):
    return np.block([[I4,I4],[u.T,-u.T]])/np.sqrt(2)


def covariance(u, s, t):
    w = yz_to_x(u)
    return w @ np.diag([s]*4+[t]*4) @ w.T


def normal_expectation(m,u,s,t):
    return float(np.trace(m@covariance(u,s,t))/2)


def fraction_string(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def calculate():
    residuals=[]
    def close(a,b,tol=5e-11):
        err=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
        residuals.append(err)
        assert err<tol, (err,tol)

    for j in J:
        close(j.T,-j)
        close(j@j,-I4)
        for ell in L:
            close(j@ell,ell@j)
    for i,j in enumerate(J):
        for k,jk in enumerate(J):
            close(j@jk+jk@j,-2*I4 if i==k else np.zeros((4,4)))

    # Independently solve the real commutant and its symmetric part.
    basis=[]
    for i in range(4):
        for j in range(4):
            e=np.zeros((4,4)); e[i,j]=1; basis.append(e)
    linear=np.array([np.concatenate([(l@e-e@l).ravel() for l in L]) for e in basis]).T
    nullity=16-int(np.linalg.matrix_rank(linear,tol=1e-10))
    symmetric=[]
    for i in range(4):
        for j in range(i,4):
            e=np.zeros((4,4)); e[i,j]=1; e[j,i]=1; symmetric.append(e)
    symlinear=np.array([np.concatenate([(l@e-e@l).ravel() for l in L]) for e in symmetric]).T
    symnullity=10-int(np.linalg.matrix_rank(symlinear,tol=1e-10))
    assert (nullity,symnullity)==(4,1)

    rng=np.random.default_rng(1282)
    max_bound_ratio=0.0
    min_eigen=1.0
    gaussian_cases=[]
    for index in range(48):
        quat=rng.normal(size=4); quat/=np.linalg.norm(quat)
        q1=rng.normal(size=4); q1/=np.linalg.norm(q1)
        q2=rng.normal(size=4); q2/=np.linalg.norm(q2)
        u,g1,g2=left(quat),left(q1),left(q2)
        a,d=rng.uniform(.2,3,size=2)
        bfour=rng.normal(size=4); bfour*=np.sqrt(a*d)*rng.uniform(0,.999)/np.linalg.norm(bfour)
        b0,bv=bfour[0],bfour[1:]
        m=block(a,d,b0,bv,u)
        ev=np.linalg.eigvalsh(m)
        expected=np.sort(np.repeat([(a+d-np.sqrt((a-d)**2+4*np.dot(bfour,bfour)))/2,
                                   (a+d+np.sqrt((a-d)**2+4*np.dot(bfour,bfour)))/2],4))
        close(ev,expected)
        assert ev.min()>-1e-12
        min_eigen=min(min_eigen,float(ev.min()))
        transform=np.block([[g1,np.zeros((4,4))],[np.zeros((4,4)),g2]])
        u2=g1@u@g2.T
        close(transform.T@block(a,d,b0,bv,u2)@transform,m)
        q,r=a+d-2*b0,a+d+2*b0
        v=(a-d)*I4/2+sum((c*j for c,j in zip(bv,J)),start=np.zeros((4,4)))
        close(yz_to_x(u).T@m@yz_to_x(u),np.block([[q*I4/2,v],[v.T,r*I4/2]]))
        s,t=rng.uniform(.1,3,size=2)
        close(normal_expectation(m,u,s,t),s*q+t*r)
        close(normal_expectation(m,u,s+.7,t)-normal_expectation(m,u,s,t),.7*q)
        discr=(a-d)**2+4*np.dot(bv,bv)
        assert discr<=q*r+1e-11
        target=block(r/4,r/4,r/4,[0,0,0],u)
        bound=(q+np.sqrt(q*r))/2
        error=norm(m-target)
        assert error<=bound+1e-11
        max_bound_ratio=max(max_bound_ratio,error/bound if bound else 0)
        # Arbitrary normal covariance: test the finite second-moment bound.
        x=rng.normal(size=(8,8)); sigma=x@x.T/8
        moment=float(np.trace(sigma))
        actual=abs(float(np.trace((m-target)@sigma)/2))
        assert actual<=moment*bound/2+1e-10
        # Exact Gaussian fourth moment for the full group kinetic term.
        zmat=np.diag([s]*4+[t]*4)
        cmat=np.diag([.5]*4+[-.5]*4)
        p1=np.diag([1]*4+[0]*4)
        p2=np.diag([0]*4+[1]*4)
        # The determinant identity is most directly audited at x coordinates.
        sx=covariance(u,s,t)
        px=np.diag([1]*4+[0]*4); py=np.diag([0]*4+[1]*4)
        cx=np.block([[np.zeros((4,4)),u/2],[u.T/2,np.zeros((4,4))]])
        def fourth(a0,b0m):
            return 2*np.trace(a0@sx@b0m@sx)+np.trace(a0@sx)*np.trace(b0m@sx)
        det_mean=float(fourth(px,py)-fourth(cx,cx))
        close(det_mean,12*s*t)
        grad_log_coefficient=(1/t-1/s)/4
        close(grad_log_coefficient**2*det_mean/4,3*(s-t)**2/(16*s*t))
        if index<3:
            gaussian_cases.append(dict(s=float(s),t=float(t),q=float(q),r=float(r),
                                       expected_F=float(s*q+t*r)))

    # The chosen positive family has equal nonzero normal-state cost.
    kappa=F(3,2); t=F(1,4); s1=F(1); s2=F(3)
    selected=block(float(kappa),float(kappa),float(kappa),[0,0,0],I4)
    selected_mean=4*kappa*t
    close(normal_expectation(selected,I4,float(s1),float(t)),float(selected_mean))
    close(normal_expectation(selected,I4,float(s2),float(t)),float(selected_mean))

    # Exact finite-error fixture, with strictly positive candidate edge form.
    a=d=F(1); b0=F(99,100); b1=F(1,10)
    q=a+d-2*b0; r=a+d+2*b0
    assert a*d>b0*b0+b1*b1
    mean1=s1*q+t*r; mean2=s2*q+t*r
    eta1=eta2=F(1,1000)
    epsilon=(abs(mean2-mean1)+eta1+eta2)/(s2-s1)
    rmax=F(4); moment=F(13)
    assert q<=epsilon and r<=rmax and 4*b1*b1<=q*r
    finite_m=block(float(a),float(d),float(b0),[float(b1),0,0],I4)
    finite_target=block(float(r/4),float(r/4),float(r/4),[0,0,0],I4)
    coefficient_bound=(float(epsilon)+np.sqrt(float(epsilon*rmax)))/2
    coefficient_error=norm(finite_m-finite_target)
    assert coefficient_error<=coefficient_bound

    # Deletion 1: signed cancellation at normal-state level.
    indefinite=block(1,-1,0,[0,0,0],I4)
    for s in (1,3): close(normal_expectation(indefinite,I4,s,.25),0)
    assert np.linalg.eigvalsh(indefinite).min()<0
    # Deletion 2: anisotropy passes the Gaussian calibration without local covariance.
    aniso=np.diag([1.,2.,3.,4.])
    ma=np.block([[aniso,-aniso],[-aniso,aniso]])
    for s in (1,3): close(normal_expectation(ma,I4,s,.25),2.5)
    assert norm(L[0].T@aniso@L[0]-aniso)>0
    # Deletion 3: a positive quartic invariant has the same amplitude calibration.
    quartic_mean=4*F(1)*t+96*F(1)*t*t
    assert quartic_mean==7
    # Deletion 4: the zero edge remains allowed.
    close(block(0,0,0,[0,0,0],I4),np.zeros((8,8)))
    # Deletion 5: changing differential width can cancel q in the mean difference.
    nocross=block(1,1,0,[0,0,0],I4)
    close(normal_expectation(nocross,I4,1,2),6)
    close(normal_expectation(nocross,I4,2,1),6)

    budgets=[]
    for s in (s1,s2):
        coord=1/s+1/t
        group=3*(s-t)**2/(16*s*t)
        h=s+t
        # Original W expression at a declared algebra check point v=L11=L22=uh=us=1,L12=0.
        # These numbers test the moment formula only, not an old numerical parameter replay.
        w=(6*h*h-4*h+1+2)/2
        full=coord+F(1,2)+group+w+selected_mean
        budgets.append({"sigma_squared":fraction_string(s),"delta_squared":fraction_string(t),
                        "coordinate_kinetic_per_a_kin":fraction_string(coord),
                        "group_kinetic_per_b_el":fraction_string(group),
                        "edge_mean":fraction_string(selected_mean),
                        "W_algebra_fixture":fraction_string(w),
                        "full_algebra_fixture":fraction_string(full)})
    assert budgets[0]["edge_mean"]==budgets[1]["edge_mean"]
    assert budgets[0]["full_algebra_fixture"]!=budgets[1]["full_algebra_fixture"]

    return {
        "round":1282,"status":"PASS","new_axioms":0,
        "commutant":{"real_dimension":nullity,"symmetric_dimension":symnullity},
        "random_positive_cases":48,"max_matrix_residual":max(residuals),
        "minimum_sample_eigenvalue":min_eigen,"largest_coefficient_bound_ratio":max_bound_ratio,
        "finite_normal_examples":gaussian_cases,
        "selected_normal_cost":fraction_string(selected_mean),
        "finite_error":{"q":fraction_string(q),"r":fraction_string(r),
                        "mean1":fraction_string(mean1),"mean2":fraction_string(mean2),
                        "eta_each":fraction_string(eta1),"epsilon":fraction_string(epsilon),
                        "r_max":fraction_string(rmax),"second_moment_cap":fraction_string(moment),
                        "target_kappa":fraction_string(r/4),"coefficient_error":coefficient_error,
                        "coefficient_bound":coefficient_bound,
                        "expectation_bound":float(moment)*coefficient_bound/2},
        "deletion_tests":{"positivity":True,"full_local_covariance":True,"quadratic_class":True,
                          "nonzero_interaction_not_selected":True,"fixed_difference_width":True},
        "same_source_budget":budgets,
        "scope":{"normal_Gauss_tests":True,"all_fixed_polynomial_moments_finite":"analytic proof",
                 "actual_autonomous_preparation_certified":False,
                 "arbitrary_unknown_cost_readout_granted":False,
                 "potential_operator_norm_bound_claimed":False,
                 "gauge_group_or_Yang_Mills_dynamics_generated":False,
                 "old_probe_results_recounted":False}}


def equivalent(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equivalent(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return isinstance(b,(int,float)) and abs(a-b)<=1e-11*max(1,abs(a),abs(b))
    return a==b


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=calculate()
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2); stream.write("\n")
        print("PASS: results.json created exclusively")
    else:
        saved=json.loads(path.read_text(encoding="utf-8"))
        assert equivalent(result,saved), "Saved result mismatch"
        print("PASS: read-only replay MATCH")
    print("result_sha256="+hashlib.sha256(path.read_bytes()).hexdigest())
    print("max_matrix_residual="+str(result["max_matrix_residual"]))


if __name__=="__main__":
    main()
