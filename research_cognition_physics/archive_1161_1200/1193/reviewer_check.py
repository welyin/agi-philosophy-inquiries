"""1193 internal, non-blind independent review.
No import/execution of author check.py. Rational explicit/implicit Euler enclosure,
scalar Taylor recurrence, Hamiltonian-derived moment identities, deletion checks.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
ZERO = F(0)

def enc(q):
    return {"exact":str(q), "float":float(q)}

def author_results():
    return json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))

def comparison_budget():
    B,L,v,r,rf,sigma,eta = F(1,100),F(1,50),F(1,2),F(1,10000),F(1,10**7),F(1,10000),F(1,100000)
    S = 1+3*v
    a,b = (1+v)*L*S*S, B*(S*S+6*(1+v)*S)
    D,W = B*r+rf, (1+v)*(B*sigma+L*r+eta)
    expected = {"A_x":F(3,16),"A_v":F(23,80),"D":F(11,10**7),"W":F(39,2*10**6)}
    assert dict(zip(expected,(a,b,D,W))) == expected
    data = author_results()["finite_budget"]
    for key,value in expected.items():
        assert str(value)==data[key]["exact"]

    # For u'=M u+b with M,b nonnegative, each exact component is increasing.
    # Forward Euler is a lower endpoint bound. Integration bounded by final
    # endpoint gives backward Euler upper bound, since (I-hM)^-1 is nonnegative.
    steps=128
    h=F(1,steps)
    det=1-h*b-h*h*a
    assert det>0 and 1-h*b>0
    lower=[ZERO,ZERO]
    upper=[ZERO,ZERO]
    for _ in range(steps):
        x,y=lower
        lower=[x+h*(y+D), y+h*(a*x+b*y+W)]
        z0,z1=upper[0]+h*D,upper[1]+h*W
        upper=[((1-h*b)*z0+h*z1)/det,(h*a*z0+z1)/det]
    assert upper[0]<F(2,100000)
    assert upper[1]<F(3,100000)
    for i,name in enumerate(("position","velocity")):
        saved_low=F(data[name+"_lower"]["exact"])
        saved_high=F(data[name+"_upper"]["exact"])
        assert lower[i] <= saved_low <= saved_high <= upper[i]

    # Independently check proof's coarser 12-term enclosure through a scalar
    # second-order recurrence, not through author matrix powers.
    derivatives_x=[ZERO,D,W]
    derivatives_y=[ZERO,W,a*D+b*W]
    for n in range(2,13):
        derivatives_x.append(a*derivatives_x[n-1]+b*derivatives_x[n])
        derivatives_y.append(a*derivatives_y[n-1]+b*derivatives_y[n])
    tail=F(3)*max(D,W)/math.factorial(14)
    proof_high=[
        sum((derivatives_x[n]/math.factorial(n) for n in range(1,14)),ZERO)+tail,
        sum((derivatives_y[n]/math.factorial(n) for n in range(1,14)),ZERO)+tail]
    exact_proof=[
        F(30117024078568791167039466180833,2496186064167763968000000000000000000),
        F(359750129290938154233532896012643,15361145010263162880000000000000000000)]
    assert proof_high==exact_proof
    assert all(lower[i]<=proof_high[i]<=upper[i] for i in (0,1))
    return {"method":"exact forward/backward Euler monotone enclosures",
            "steps":steps, "position_interval":[enc(lower[0]),enc(upper[0])],
            "velocity_interval":[enc(lower[1]),enc(upper[1])],
            "proof_scalar_recurrence_upper":list(map(enc,proof_high)),
            "proof_tail":enc(tail), "fixed_thresholds_pass":True,
            "constants":{key:enc(value) for key,value in expected.items()}}

def full_lipschitz_check():
    # Independently bound the change of the full quadratic contraction.
    # sum |z| <= 1+3v and sum |z-ztilde| <= 3*max_i |delta z_i|.
    v=F(1,2); S=1+3*v; B=F(1,100); L=F(1,50)
    position=(1+v)*L*S*S
    quadratic_velocity=2*S*3*B
    outer_velocity=B*S*S
    velocity=(1+v)*quadratic_velocity+outer_velocity
    assert position==F(3,16) and velocity==F(23,80)
    # An explicit dense tensor makes the position bound an equality:
    # Gamma^i=-B and Gamma^0=B, all z spatial components =+v.
    assert (1+v)*S*S == F(75,8)
    return {"A_x":"3/16","A_v":"23/80",
            "spatial_variation_coefficient":"75/8",
            "comparison_requires_both_trajectories_and_velocity_domains":True}

def curved_moments():
    # Independent prescribed-background, single-time diagnostic.
    # Start with H=N(x)*sqrt(p^2+m^2). Obtain xdot,pdot,E_coord_dot
    # by differentiating H and E_coord=sqrt(p^2+m^2)/N.
    eps=F(1,100)
    xs=[F(-1,4),F(1,8),F(1,3)]
    masses=[F(1),F(2),F(3)]
    ps=[F(5,12)*masses[0],-F(5,12)*masses[1],ZERO]
    rests=[F(13,12)*masses[0],F(13,12)*masses[1],masses[2]]
    rows=[]
    for x,m,p,e0 in zip(xs,masses,ps,rests):
        assert e0*e0-p*p==m*m
        N=1+eps*x*x
        Np=2*eps*x
        velocity=N*p/e0
        pdot=-Np*e0
        e0dot=(p/e0)*pdot
        Ndot=Np*velocity
        E=e0/N
        Edot=e0dot/N-e0*Ndot/(N*N)
        vdot=Ndot*p/e0+N*pdot/e0-N*p*e0dot/(e0*e0)
        # Nonzero curvature at every point: N''=2eps, R00=N*N''.
        assert N*2*eps>0
        rows.append((x,velocity,E,Edot,vdot))
    total=sum((q[2] for q in rows),ZERO)
    X=sum((x*E for x,v,E,Edot,vdot in rows),ZERO)/total
    w=sum((v*E for x,v,E,Edot,vdot in rows),ZERO)/total
    Eprime=sum((q[3] for q in rows),ZERO)
    Pprime=sum((Edot*v+E*vdot for x,v,E,Edot,vdot in rows),ZERO)
    Xprime=(sum((Edot*x+E*v for x,v,E,Edot,vdot in rows),ZERO)-X*Eprime)/total
    wprime=(Pprime-w*Eprime)/total
    d=Xprime-w
    sigma=sum((E*(v-w)**2 for x,v,E,Edot,vdot in rows),ZERO)/total
    r=sum((E*abs(x-X)*(1+abs(v))**2 for x,v,E,Edot,vdot in rows),ZERO)/total
    N=1+eps*X*X; Np=2*eps*X
    Fcenter=-N*Np+2*(Np/N)*w*w
    e=wprime-Fcenter
    rhs_Eprime=ZERO; rhs_Pprime=ZERO; rhs_d=ZERO
    for x,v,E,Edot,vdot in rows:
        ni=1+eps*x*x; nip=2*eps*x
        source0=-2*(nip/ni)*E*v
        source1=-ni*nip*E
        rhs_Eprime+=source0
        rhs_Pprime+=source1
        rhs_d+=(x-X)*source0/total
    assert Eprime==rhs_Eprime and Pprime==rhs_Pprime and d==rhs_d
    B,L,vmax=F(1,80),F(21,1000),F(2,5)
    assert all(abs(x)<F(3,5) and abs(v)<vmax for x,v,E,Edot,vdot in rows)
    assert abs(X)<F(3,5) and abs(w)<vmax
    Db=B*r; Wb=(1+vmax)*(B*sigma+L*r)
    assert abs(d)<=Db and abs(e)<=Wb
    # Bounds valid on |x|<=.6: Gamma1_00=N*N', Gamma0_01=N'/N.
    radius=F(3,5)
    assert 2*eps*radius*(1+eps*radius*radius)<=B
    assert 2*eps+6*eps*eps*radius*radius<=L
    return {"method":"Hamiltonian differentiation, then exact finite stress moments",
            "stream_count":len(rows), "moment_identities_exact":True,
            "r_tau":enc(r),"Sigma_l1":enc(sigma),
            "centroid_slip":enc(d),"slip_bound":enc(Db),
            "momentum_defect":enc(e),"defect_bound":enc(Wb),
            "scope":"prescribed genuine-curvature diagnostic at one time; no coupled source or quantum instrument"}

def deletion_checks():
    # Smooth positive profile supported in [-radius,radius], two equal local
    # energy streams of coordinate velocities +/-v. Total tau01 vanishes
    # initially; tau11=v^2*tau00. Each stream can satisfy the local dust
    # transport equations, so this is a regular initial-slice diagnostic.
    alpha,v,radius=F(1,5),F(1,2),F(1,1000)
    assert (1+alpha*radius)*v<1
    sigma=v*v
    lower_accel=alpha*v*v/(1+alpha*radius)
    L=alpha*alpha/(1-alpha*radius)**2
    r_bound=radius*(1+v*v)
    false_bound=(1+v)*L*r_bound
    assert lower_accel>false_bound
    # Exactly coincident two-stream algebra, used by author separately.
    assert -alpha*v*v==F(-1,20)
    data=author_results()["deletion"]
    assert data["cold_closure"]["Sigma_11"]==str(sigma)
    assert data["cold_closure"]["centroid_acceleration"]==str(-alpha*v*v)
    # Constant-energy translated profile rho(x-X), X'=v, v'=a in Minkowski:
    # partial_t tau00+partial_x tau10=0;
    # partial_t tau01+partial_x tau11=a*rho.
    acceleration=F(1,100)
    assert data["exchange"]["required_eta"]==str(acceleration)
    return {"smooth_two_stream_initial_slice":{
                "support_radius":str(radius),"Sigma_11":str(sigma),
                "acceleration_magnitude_lower":enc(lower_accel),
                "invalid_bound_if_Sigma_discarded":enc(false_bound),
                "flat_nonlinear_coordinates":True,
                "scope":"initial-slice closure deletion, not spacetime curvature or an equivalence violation"},
            "exchange":{
                "required_eta":"1/100","cold_closure_zero":True,
                "environment_opposite_exchange_required":True,
                "scope":"prescribed smooth effective stress, not an implemented full quantum process"}}

def bindings():
    names=("admission.md","proof.md","sources.md","check.py","results.json","baseline.json","NEXT.md")
    paths=[HERE/name for name in names]+[HERE.parent/"research_note_1193.md"]
    out={}
    for p in paths:
        assert p.is_file(),str(p)
        out[p.relative_to(HERE.parent).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    return out

def run():
    return {"status":"PASS_WITH_SCOPE","review_type":"internal non-blind independent mathematical review",
            "comparison_budget":comparison_budget(),"lipschitz":full_lipschitz_check(),
            "curved_moments":curved_moments(),"deletions":deletion_checks(),
            "author_files_sha256":bindings(),
            "scope":{"conditional_centroid_theorem":True,
                     "actual_CO_instrument_certified":False,
                     "proper_clock_identification_certified":False,
                     "common_quantum_parent_contracts_certified":False,
                     "all_inputs_and_reference_channel_bound":False,
                     "Einstein_dynamics_generated":False,"gr_goal_complete":False}}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    args=ap.parse_args()
    data=run(); target=HERE/"review_checks.json"
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(data,f,ensure_ascii=False,indent=2); f.write("\n")
    else:
        assert json.loads(target.read_text(encoding="utf-8-sig"))==data
    print(json.dumps({"status":data["status"],
                      "position_upper":data["comparison_budget"]["position_interval"][1]["float"],
                      "velocity_upper":data["comparison_budget"]["velocity_interval"][1]["float"],
                      "bound_author_files":len(data["author_files_sha256"])},ensure_ascii=False))
