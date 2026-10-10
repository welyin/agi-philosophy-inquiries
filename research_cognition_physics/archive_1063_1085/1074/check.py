"""Finite checks of explicit models, not a proof of the topology theorem."""
import json
import numpy as np

rng=np.random.default_rng(1074)
groups=0
max_eq=0.0
min_slack=float("inf")

def eq(x,y,tol=5e-12):
    global groups,max_eq
    err=float(np.max(np.abs(np.asarray(x)-np.asarray(y))))
    assert err<tol,(err,x,y)
    max_eq=max(max_eq,err)
    groups+=1

def nonneg(value,tol=5e-12):
    global groups,min_slack
    assert value>=-tol,value
    min_slack=min(min_slack,float(value))
    groups+=1

def mean(x,y,a):
    e=np.zeros_like(x);e[-1]=1
    return (x+y)/2+a*np.linalg.norm(x-y)*e

def F(v,a):
    e=np.zeros_like(v);e[-1]=1
    return v-2*a*np.linalg.norm(v)*e

def inverse_iterate(x,a):
    # Independent of the algebraic inverse used in the separate review.
    e=np.zeros_like(x);e[-1]=1
    v=x.copy()
    for _ in range(700):
        v=x+2*a*np.linalg.norm(v)*e
    return v

def solve_other(x,target,a):
    e=np.zeros_like(x);e[-1]=1
    y=2*target-x
    for _ in range(700):
        y=2*target-x-2*a*np.linalg.norm(x-y)*e
    return y

def rotation(n):
    z=rng.normal(size=(n,n))
    q,_=np.linalg.qr(z)
    if np.linalg.det(q)<0:q[:,0]*=-1
    return q

def unit(v):
    return v/np.linalg.norm(v)

for n in (1,2,3,4,7):
    for a in (.1,.3,.45):
        e=np.eye(n)[-1]
        for _ in range(8):
            x,y,z,v,w=rng.normal(size=(5,n))*.03
            eq(mean(x,y,a),mean(y,x,a))
            eq(mean(x,x,a),x)
            eq(mean(F(v,a),F(-v,a),a),np.zeros(n))
            eq(inverse_iterate(F(v,a),a),v)
            eq(solve_other(x,mean(x,y,a),a),y)
            ix=F(-inverse_iterate(x,a),a)
            eq(F(-inverse_iterate(ix,a),a),x)
            lower=np.linalg.norm(mean(x,y,a)-mean(x,z,a))-(.5-a)*np.linalg.norm(y-z)
            nonneg(lower)
            nonneg(np.linalg.norm(F(v,a)-F(w,a))-(1-2*a)*np.linalg.norm(v-w))
            h=np.eye(n)
            if n>1:h[:-1,:-1]=rotation(n-1)
            t=rng.normal(size=n)*.01
            eq(mean(h@x+t,h@y+t,a),h@mean(x,y,a)+t)
            eq(F(h@v,a),h@F(v,a))
        # Non-C1 diagnostic has a closed-form nonzero limit in the proof.
        hval=1e-6
        u=np.eye(n)[0]
        eq((mean(hval*u,np.zeros(n),a)+mean(-hval*u,np.zeros(n),a))/hval,2*a*e)

# Sphere: nonnormal isotropy, distinct dimensions of orbit/group.
o=np.array([0.,0.,1.])
def exp_o(v):
    r=np.linalg.norm(v)
    return o if r==0 else np.r_[np.sin(r)*v/r,np.cos(r)]
def midpoint(x,y):
    return unit(x+y)
def reverse(x):
    return 2*np.dot(o,x)*o-x
def Rz(t):
    return np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])
for _ in range(30):
    v,w=rng.normal(size=(2,2))*.04
    x,y=exp_o(v),exp_o(w)
    eq(midpoint(x,reverse(x)),o)
    eq(reverse(reverse(x)),x)
    eq(exp_o(-v),reverse(x))
    g=rotation(3)
    eq(midpoint(g@x,g@y),g@midpoint(x,y))
    h=Rz(float(rng.uniform(-np.pi,np.pi)))
    eq(reverse(h@x),h@reverse(x))
g=np.array([[0.,0.,1.],[0,1,0],[-1,0,0]])
conj=g@Rz(.4)@g.T
nonnormal_gap=float(np.linalg.norm(conj@o-o))
assert nonnormal_gap>.3
groups+=1

# Failure without cancellation: max is an endpoint member a=1/2.
for t in (.01,.1,1.):
    eq(mean(np.array([0.]),np.array([-t]),.5),[0])
    eq(F(np.array([t]),.5),[0])

# Failure of shared shell if a larger action is silently imposed.
a=.25;e=np.array([0.,0.,1.]);u=np.array([1.,0.,0.])
h=np.diag([1.,-1.,-1.])
covariance_gap=float(np.linalg.norm(mean(h@u,-h@u,a)-h@mean(u,-u,a)))
eq(covariance_gap,4*a)
opposite=F(-inverse_iterate(e,a),a)
eq(opposite,-(1+2*a)/(1-2*a)*e)

print(json.dumps({
    "round":1074,"passed":True,"assertion_groups":groups,
    "max_equality_residual":max_eq,"minimum_inequality_slack":min_slack,
    "tested_dimensions":[1,2,3,4,7],
    "author_inverse_method":"700-step contraction iteration; proof supplies all-domain bound",
    "nonnormal_isotropy_gap":nonnormal_gap,
    "full_rotation_covariance_failure_gap":covariance_gap,
    "old_sphere_opposite_radius":float(np.linalg.norm(opposite)),
    "topological_claim_proved_by_sampling":False,
    "actual_execution_generated":False,"three_dimensions_derived":False,
    "empirical_data":False
},ensure_ascii=False,indent=2))
