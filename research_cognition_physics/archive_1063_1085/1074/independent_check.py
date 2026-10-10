"""1074 independent audit: continuous meetings and local reversed shells.

NumPy only. No author computations are imported. This script only prints JSON.
The proof is read as bytes solely to bind the audited version by SHA-256.
Finite samples check explicit witnesses; they do not prove the general theorem.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np

TOL = 3e-11
RNG = np.random.default_rng(1074013)
groups = []

def norm(x):
    return float(np.linalg.norm(x))

def err(x,y):
    return float(np.max(np.abs(np.asarray(x)-np.asarray(y))))

def add_group(name, equalities=(), margins=(), diagnostics=None):
    equality = max([0.0]+[float(x) for x in equalities])
    margin = min([0.0]+[float(x) for x in margins])
    finite = np.isfinite(equality) and np.isfinite(margin)
    ok = bool(finite and equality <= TOL and margin >= -TOL)
    # Exactly one assertion group, irrespective of the number of sample points.
    assert ok, (name,equality,margin)
    groups.append({"name":name,"passed":ok,
                   "max_absolute_equality_residual":equality,
                   "minimum_inequality_margin":margin,
                   "diagnostics":diagnostics or {}})

def axis(d):
    v=np.zeros(d); v[-1]=1.0
    return v

def meeting(x,y,a):
    return (x+y)/2 + a*norm(x-y)*axis(len(x))

def forward(v,a):
    return v-2*a*norm(v)*axis(len(v))

def inverse(x,a):
    aa=1-4*a*a
    radial=(2*a*x[-1]+math.sqrt(x[-1]**2+aa*float(x[:-1]@x[:-1])))/aa
    return x+2*a*radial*axis(len(x))

def opposite(x,a):
    return forward(-inverse(x,a),a)

def so(k):
    if k==0:
        return np.zeros((0,0))
    if k==1:
        return np.eye(1)
    q,_=np.linalg.qr(RNG.normal(size=(k,k)))
    q[:,-1]*=np.linalg.det(q)
    return q

def isotropy(d):
    r=np.eye(d)
    if d>1:
        r[:-1,:-1]=so(d-1)
    return r

# Every dimension uses the same formulas; the proof covers arbitrary d>=1.
cases=[]
for d in [1,2,3,5,8]:
    for a in [1/7,0.23,0.4]:
        pts=[np.zeros(d),0.11*axis(d),-0.08*axis(d)]
        for _ in range(18):
            v=RNG.normal(size=d)
            v*=RNG.uniform(0.01,0.24)/norm(v)
            pts.append(v)
        cases.append((d,a,pts))

eq1=[]; eq2=[]; eq3=[]; mg3=[]; eq4=[]; eq5=[]; eq6=[]
case_report=[]; max_iterations=0
for d,a,pts in cases:
    e=axis(d)
    mr=0.0; ir=0.0; sr=0.0
    for j,x in enumerate(pts):
        y=pts[(j+5)%len(pts)]; z=pts[(j+11)%len(pts)]
        eq1.extend([err(meeting(x,y,a),meeting(y,x,a)),err(meeting(x,x,a),x)])
        fx=forward(x,a)
        cur=err(inverse(fx,a),x)
        eq2.extend([cur,err(forward(inverse(x,a),a),x),
                    err(opposite(opposite(x,a),a),x),
                    err(meeting(fx,forward(-x,a),a),np.zeros(d))])
        ir=max(ir,cur)
        mr=max(mr,err(meeting(fx,forward(-x,a),a),np.zeros(d)))
        dist=norm(meeting(x,y,a)-meeting(x,z,a))
        yz=norm(y-z)
        mg3.extend([dist-(0.5-a)*yz,(0.5+a)*yz-dist])
        fd=norm(forward(y,a)-forward(z,a))
        mg3.extend([fd-(1-2*a)*yz,(1+2*a)*yz-fd])
        # Independent fixed-point solver for an arbitrary target b. This is
        # numerical arithmetic on declared coordinates, NOT an execution right.
        b=0.07*e+0.1*y
        closed=b+opposite(x-b,a)
        iterate=np.zeros(d)
        for iteration in range(600):
            nxt=2*b-x-2*a*norm(x-iterate)*e
            if norm(nxt-iterate)<1e-14:
                iterate=nxt
                break
            iterate=nxt
        max_iterations=max(max_iterations,iteration+1)
        eq3.extend([err(meeting(x,closed,a),b),err(iterate,closed)])
        # Allowed full group: translations and rotations that fix e.
        r=isotropy(d); shift=0.19*z
        eq4.extend([err(r@r.T,np.eye(d)),abs(float(np.linalg.det(r))-1),
                    err(r@e,e),
                    err(meeting(r@x+shift,r@y+shift,a),r@meeting(x,y,a)+shift),
                    err(forward(r@x,a),r@fx),
                    err(opposite(r@x,a),r@opposite(x,a))])
        if norm(x)>0:
            radius=0.13
            v=radius*x/norm(x)
            f=forward(v,a)
            center=-2*a*radius*e
            shell=abs(norm(f-center)-radius)
            eq5.extend([shell,
                        abs(norm(opposite(f,a)-center)-radius),
                        err(opposite(f,a),2*center-f),
                        abs(norm(r@f-center)-radius)])
            sr=max(sr,shell)
    # The exact |t| expression, rather than numerical differentiation, gives
    # unequal one-sided derivatives at the diagonal.
    u=np.ones(d)/math.sqrt(d)
    for k in [3,8,15,24]:
        t=2.0**(-k)
        positive=meeting(t*u,-t*u,a)/t
        negative=meeting(-t*u,t*u,a)/(-t)
        eq6.extend([err(positive,2*a*e),err(negative,-2*a*e),
                    abs(norm(positive-negative)-4*a)])
    case_report.append({"dimension":d,"a":a,"sample_points":len(pts),
                        "explicit_inverse_residual":ir,
                        "meeting_at_anchor_residual":mr,
                        "complete_shell_formula_residual":sr})

add_group("01_symmetry_and_idempotence",eq1,diagnostics={"parameter_cases":len(cases)})
add_group("02_explicit_inverse_meeting_and_involution",eq2,diagnostics={"dimensions":[1,2,3,5,8]})
add_group("03_right_uniqueness_bounds_and_independent_solver",eq3,mg3,
          {"maximum_fixed_point_iterations":max_iterations,
           "scope":"Sample checks of bounds proved by triangle inequalities and contraction."})
add_group("04_full_allowed_group_covariance",eq4,
          diagnostics={"group":"R^d semidirect SO(d-1); d=1 uses the trivial stabilizer."})
add_group("05_complete_offset_shell_and_isotropy",eq5,
          diagnostics={"formula":"F(S_r) is the sphere centered at -2*a*r*e with radius r."})
add_group("06_non_C1_one_sided_derivatives",eq6,
          diagnostics={"derivative_jump":"4*a, exactly for every positive t."})

# A concrete nonnormal stabilizer in the affine homogeneous model.
def skew(v):
    x,y,z=v
    return np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])

def rot(v,angle):
    v=np.asarray(v,dtype=float); v=v/norm(v)
    k=skew(v)
    return np.eye(3)+math.sin(angle)*k+(1-math.cos(angle))*(k@k)

def affine(r,t):
    g=np.eye(4);g[:3,:3]=r;g[:3,3]=t
    return g

rz=rot([0,0,1],0.61)
translation=np.array([0.31,0.12,-0.03])
h=affine(rz,np.zeros(3)); g=affine(np.eye(3),translation)
conjugate=g@h@np.linalg.inv(g)
origin4=np.array([0.,0.,0.,1.])
motion=(conjugate@origin4-origin4)[:3]
add_group("07_affine_stabilizer_is_not_normal",
          [err(h@origin4,origin4),err(motion,(np.eye(3)-rz)@translation)],
          [norm(motion)-0.1],
          {"conjugated_stabilizer_moves_origin_by":norm(motion),
           "position_dimension":3,"group_dimension":4,"stabilizer_dimension":1})

# S^2 = SO(3)/SO(2), with the northern hemisphere as the local domain.
north=np.array([0.,0.,1.])

def sphere_exp(v):
    v=np.asarray(v,dtype=float)
    r=norm(v)
    if r==0:
        return north.copy()
    return np.array([math.sin(r)*v[0]/r,math.sin(r)*v[1]/r,math.cos(r)])

def sphere_log(p):
    radial=norm(p[:2])
    if radial==0:
        return np.zeros(2)
    angle=math.atan2(radial,float(p[2]))
    return angle*p[:2]/radial

def sphere_meet(p,q):
    total=p+q
    return total/norm(total)

def sphere_opposite(p):
    return 2*float(north@p)*north-p

def section(p):
    # Minimal rotation north -> p. Its local exponent lies in the reductive
    # complement of rotations around north.
    k=skew(np.cross(north,p))
    return np.eye(3)+k+(k@k)/(1+float(north@p))

tangent=[np.zeros(2)]
for _ in range(32):
    v=RNG.normal(size=2);v*=RNG.uniform(0.02,0.34)/norm(v)
    tangent.append(v)

eq8=[]; eq9=[]; eq10=[]
for j,v in enumerate(tangent):
    p=sphere_exp(v); q=sphere_exp(tangent[(j+7)%len(tangent)])
    r=rot(RNG.normal(size=3),0.17)
    hz=rot([0,0,1],0.27+j*0.1)
    eq8.extend([abs(norm(p)-1),err(sphere_meet(p,q),sphere_meet(q,p)),
                err(sphere_meet(p,p),p),
                err(sphere_meet(r@p,r@q),r@sphere_meet(p,q))])
    plus=sphere_exp(v); minus=sphere_exp(-v)
    center=sphere_meet(plus,minus)
    f=section(center).T@plus
    inv=sphere_opposite(f)
    eq9.extend([err(center,north),err(f,plus),
                err(sphere_log(f),v),err(inv,minus),
                err(sphere_meet(f,inv),north),
                err(sphere_opposite(inv),f),
                err(sphere_opposite(hz@f),hz@inv)])
    sp=section(p)
    t=sp.T@north
    eq10.extend([err(sp.T@sp,np.eye(3)),abs(float(np.linalg.det(sp))-1),
                 err(sp@north,p),err(t,sphere_opposite(p)),
                 err(section(t).T@north,p),
                 err(section(hz@p),hz@sp@hz.T)])
    if norm(v)>0:
        # A complete fixed-radius circle is parametrized explicitly.
        vv=0.2*v/norm(v)
        ff=sphere_exp(vv)
        eq9.extend([abs(ff[2]-math.cos(0.2)),
                    abs(sphere_opposite(ff)[2]-math.cos(0.2)),
                    abs((hz@ff)[2]-math.cos(0.2))])

add_group("08_S2_meeting_and_full_rotation_covariance",eq8,
          diagnostics={"sample_tangent_points":len(tangent),"domain":"northern cap"})
add_group("09_S2_explicit_reverse_and_complete_circle",eq9,
          diagnostics={"local_unique_solution":"y = 2*(north dot x)*north - x"})
add_group("10_S2_reductive_section_and_center_recovery",eq10)

h=rot([0,0,1],0.63)
g=rot([0,1,0],0.71)
changed=g@h@g.T@north
generators=[skew(np.eye(3)[j]) for j in range(3)]
infinitesimal=np.column_stack([a@north for a in generators])
rank=int(np.linalg.matrix_rank(infinitesimal,tol=1e-12))
rankg=int(np.linalg.matrix_rank(np.column_stack([a.reshape(-1) for a in generators]),tol=1e-12))
add_group("11_S2_stabilizer_non_normal_and_position_dimension",
          [err(h@north,north),abs(rank-2),abs(rankg-3)],
          [norm(changed-north)-0.1],
          {"tangent_rank":rank,"Lie_algebra_rank":rankg,"stabilizer_kernel_dimension":rankg-rank,
           "conjugated_stabilizer_motion":norm(changed-north),
           "scope":"This two-dimensional example satisfies the bridge; the bridge alone does not select three."})

# Delete fixed-output uniqueness while retaining the other algebraic properties.
# In one dimension, a=1/2 gives max(x,y).
grid=np.array([-0.21,-0.1,-0.03,0.,0.04,0.12,0.2])
eq12=[]; mg12=[]
for x in grid:
    for y in grid:
        form=(x+y)/2+abs(x-y)/2
        eq12.extend([abs(form-max(x,y)),abs(form-max(y,x)),
                     abs(max(x+0.17,y+0.17)-(max(x,y)+0.17))])
    eq12.append(abs(max(x,x)-x))
for y in [-0.19,-0.07,0]:
    eq12.append(abs(max(0,y)))
for x in [0.03,0.11,0.2]:
    eq12.append(abs(x-abs(x)))
    # max(x,y)>=x for EVERY y; these margins are implementation checks only.
    mg12.extend([max(x,y)-x for y in grid])
add_group("12_missing_uniqueness_collapses_half_neighborhood",eq12,mg12,
          {"distinct_compensations_for_zero":[-0.19,-0.07,0],
           "no_solution_for_positive_x":"max(x,y)>=x>0 for every real y",
           "F":"v-|v|"})

# Delete full-group covariance: keep the a=0.23 meeting but demand all SO(3).
# The meeting remains covariant under the smaller SO(2) stabilizer only.
a=0.23; e=axis(3); rr=rot([1,0,0],math.pi)
x=np.array([0.12,0.06,0.15]); y=np.array([-0.05,0.08,0.1])
covariance_defect=meeting(rr@x,rr@y,a)-rr@meeting(x,y,a)
expected=a*norm(x-y)*(e-rr@e)
k=(1+2*a)/(1-2*a)
eq13=[err(covariance_defect,expected),err(rr@e,-e)]
mg13=[norm(covariance_defect)-0.01,k-1]
start=1e-5
point=start*e
orbit=[]
for j in range(7):
    eq13.append(err(point,(k**j)*start*e))
    orbit.append(norm(point))
    point=rr@opposite(point,a)
eq13.extend([err(opposite(start*e,a),-k*start*e),
             abs(norm(opposite(start*e,a))/start-k)])
add_group("13_missing_full_covariance_precludes_common_compact_shell",eq13,mg13,
          {"covariance_defect":norm(covariance_defect),
           "radius_multiplier":k,"finite_growth_witness":orbit,
           "analytic_reason":"A full SO(d)-invariant nonzero compact set contains r*e; (R*i)^N gives k^N*r*e, contradicting compactness.",
           "scope":"d>=2; this explicit candidate loses the common rotation/reversal shell when the required covariance is deleted."})

here=Path(__file__).resolve().parent
out={
    "round":1074,
    "independent":True,
    "assertion_groups":len(groups),
    "all_passed":all(g["passed"] for g in groups),
    "maximum_equality_residual":max(g["max_absolute_equality_residual"] for g in groups),
    "minimum_inequality_margin":min(g["minimum_inequality_margin"] for g in groups),
    "numerical_tolerance":TOL,
    "groups":groups,
    "affine_case_diagnostics":case_report,
    "proof_sha256":hashlib.sha256((here/"proof.md").read_bytes()).hexdigest(),
    "independent_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "numpy_version":np.__version__,
    "scope":[
        "The 13 groups are distinct mathematical checks; sample counts are not new assertion groups.",
        "Explicit models do not derive the prior actual contact space or physical controls.",
        "The general continuous theorem relies on the analytic proof and invariance of domain, not sampling.",
        "Bi-Lipschitz bounds occur in this witness only, not as restored assumptions of the general theorem.",
        "The stronger missing-covariance boundary is stated independently; author proof.md is not edited.",
        "No author computation imports, images, dependency installations, or default file writes."
    ]
}
print(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False))

