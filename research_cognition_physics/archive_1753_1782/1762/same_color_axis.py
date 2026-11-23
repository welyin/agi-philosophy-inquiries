"""Use the actual 10 PeV tagged down right frame at the split weak matching surface."""
import argparse
import json
import math
from pathlib import Path
import runpy
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent


def positive(mat):
    x,v=np.linalg.eigh(mat)
    return (v*np.sqrt(x))@v.conj().T


def endpoint(ub,step):
    old=runpy.run_path(str(HERE/"../1758/check.py"))
    save={}
    def load(path):
        mod=runpy.run_path(str(path))
        if Path(path).parent.name=="1755":
            matter=mod["matter_filter"]
            def capture(bg):
                save["bg"]=bg
                return matter(bg)
            mod["point"].__globals__["matter_filter"]=capture
        return mod
    old["endpoint"].__globals__["runpy"]=SimpleNamespace(run_path=load)
    _,data=old["endpoint"](ub,step)
    bg=save["bg"]
    dr=bg["tagged_light_columns"].conj()
    D=bg["left_quark_to_GJ"].conj().T@bg["yW"][1,:,:3]@bg["right_quark_to_GJ"]
    mass=np.linalg.svd(D,compute_uv=False)[::-1]
    np.testing.assert_allclose(D.conj().T@D,(dr*(mass**2))@dr.conj().T,atol=1e-15)
    np.testing.assert_allclose(bg["F"].T,(dr*[bg["cW"],1.,bg["fW"]])@dr.conj().T,atol=1e-14)
    np.testing.assert_allclose(data["Y"],bg["N"].conj(),atol=1e-15)
    B,L=data["B"],data["L"]
    root=positive(B.conj().T@B)
    U=B@np.linalg.inv(root)
    W=L@np.linalg.inv(root)
    np.testing.assert_allclose(U.conj().T@U,np.eye(2),atol=2e-14)
    ratio=math.sqrt(1+.0025/.0000754)
    roots=old["diagonal_clock_reactor"](L,ratio,.0851)
    # The last archived branch has the largest positive imaginary part and
    # is the primary endpoint's successful full oscillation candidate.
    chosen=max(roots,key=lambda r:r["physical_mu22_over_mu11"]["imag"])
    t=chosen["physical_mu22_over_mu11"]
    target=L@np.diag([1.,complex(t["real"],t["imag"])])@L.T
    return data,bg,mass,dr,U,W,target,chosen


def fit_fixed_axis(W,T,target,poles):
    A=W@T
    n=np.linalg.norm(A,axis=0)
    c=abs(A[0,:])/n
    p0=np.linalg.norm(np.linalg.qr(W,mode="reduced")[0][0])**2
    pmin,pmax=(1-np.sqrt(1-np.array([.0803,.0899])))/2
    p=np.clip(c**2,pmin,pmax)
    cap=(c*np.sqrt(p)/p0+np.sqrt((1-c**2/p0)*(1-p/p0)))**2
    rmin=math.sqrt(1+.00219/.0000773)
    required=rmin/(1+rmin)
    eig,sv,_=np.linalg.svd(target)
    eig=eig[:,:2]
    qdiag=np.diag(eig.conj().T@target@eig.conjugate())
    eig=eig*np.sqrt(qdiag/abs(qdiag))
    vals=sv[:2]/sv[0]
    q=abs(eig[:,0].conj()@A/n)**2
    y=1/(1+vals[1])
    out=dict(reactor_only_column_caps=cap,required_using_low_frequency=required,
             complete_window_excluded=bool(max(cap)<required),fixed_target_overlaps=q,
             fixed_target_required=y,allows_exact_all_oscillations=bool(min(q)<y<max(q)))
    if not out["allows_exact_all_oscillations"]:
        return out
    h2=n[1]**2*(y-q[1])/(n[0]**2*(q[0]-y))
    h=math.sqrt(h2)
    X=A@np.diag([h,1.])
    H=eig.conj().T@X@X.conj().T@eig
    np.testing.assert_allclose(H[0,0]/H[1,1],1/vals[1],rtol=1e-13)
    b=H[0,1]
    assert abs(b)>1e-8
    z=-b.conjugate()/b
    new=eig@np.diag([1.,vals[1]*z])@eig.T
    dual=np.linalg.pinv(X)
    P=dual@new@dual.T
    common=np.linalg.norm(P,ord="fro")/math.sqrt(2)
    P=P/common
    err=np.max(abs(P@P.conj().T-np.eye(2)))
    np.testing.assert_allclose(X@P@X.T*common,new,rtol=1e-10,atol=5e-14)
    np.testing.assert_allclose(new@new.conj().T,(target@target.conj().T)/sv[0]**2,atol=1e-12)
    assert err<1e-11
    gain=max(h,1/h)
    pole_ratio=(poles[1]/poles[0])**2
    stages=math.ceil(math.log(gain)/math.log(pole_ratio))
    stage=gain**(1/stages)
    mass=np.sqrt((poles[1]**2-stage*poles[0]**2)/(stage-1))
    out.update(required_f1_over_f2=h,min_positive_stages_per_leg=stages,
               common_mediator_mass_GeV=mass,fast_first_color_eigenvector=bool(h>1),
               rotation_relative_original_Majorana_deg=np.angle(z,deg=True),
               mbb_meV=abs(new[0,0])*np.sqrt(.0000754)*1000/vals[1],unitary_residual=err)
    return out


def clock_reactor_solutions(W,T,h,obs):
    """At fixed relative gain, leave all heavy-direction azimuth solutions."""
    A=W@T
    X=A@np.diag([h,1.])
    H=X@X.conj().T
    Q=np.linalg.qr(W,mode="reduced")[0]
    e=np.array([1,0,0])
    along=Q@(Q.conj().T@e)
    p0=np.vdot(along,along).real
    along/=np.sqrt(p0)
    # Inside the same rank-two plane, perpendicular to its electron projection.
    kernel=np.linalg.svd((along.conj()@Q)[None,:])[2].conj().T[:,-1]
    across=Q@kernel
    p=(1-np.sqrt(1-.0851))/2
    s,t=np.sqrt(p/p0),np.sqrt(1-p/p0)
    ratio=np.sqrt(1+.0025/.0000754)
    base=(s*s*np.vdot(along,H@along)+t*t*np.vdot(across,H@across)).real
    cross=s*t*np.vdot(along,H@across)
    co=(np.trace(H).real*ratio/(ratio+1)-base)/(2*abs(cross))
    if abs(co)>1:
        return []
    rows=[]
    for sig in [1,-1]:
        phase=np.exp(1j*(sig*np.arccos(co)-np.angle(cross)))
        m3=s*along+t*phase*across
        m2=-t*along+s*phase*across
        b=np.vdot(m3,H@m2)
        assert abs(b)>1e-12
        z=-b.conjugate()/b
        m=np.outer(m3,m3)+z*np.outer(m2,m2)/ratio
        dual=np.linalg.pinv(X)
        test=dual@m@dual.T
        U=test/(np.linalg.norm(test,"fro")/np.sqrt(2))
        np.testing.assert_allclose(U@U.conj().T,np.eye(2),atol=1e-12)
        pred=obs(m,np.zeros(3),0.)
        np.testing.assert_allclose(pred["angle_readings"][0],.0851,atol=1e-12)
        rows.append(dict(angles=pred["angle_readings"],delta_deg=pred["delta_deg"],
                         mbb_meV=abs(m[0,0])*ratio*np.sqrt(.0000754)*1000,
                         inside_angles_and_FC90=all(pred["angles_in_box"]) and pred["SK_T2K_FC68_90"][-1],
                         SK_T2K_68_90=pred["SK_T2K_FC68_90"]))
    return rows


def study(ub,step):
    d,bg,m,dr,U,W,target,chosen=endpoint(ub,step)
    c,f=bg["cW"],bg["fW"]
    weights={"down_mass_squared":m**2,"inherited_positive_gate":[c,1.,f],
             "physical_loss_probability":[1-c*c,0.,1-f*f],
             "down_only":[1.,0.,0.],"strange_only":[0.,1.,0.],"bottom_only":[0.,0.,1.]}
    out={}
    obs=runpy.run_path(str(HERE/"../1731/check.py"))["observables"]
    for name,w in weights.items():
        form=U.conj().T@(dr*w)@dr.conj().T@U
        val,T=np.linalg.eigh(form)
        T=T[:,::-1]
        fixed=fit_fixed_axis(W,T,target,d["M"])
        out[name]=dict(source_form_eigenvalues=val[::-1],first_axis_first_pole_probability=abs(T[0,0])**2,
                       high_to_low=fixed)
        if name=="inherited_positive_gate":
            # No new pole mass. Enumerate both assignments of the nondegenerate
            # poles; for each mass include the first integer reaching the
            # independently reconstructed exact-target gain. Smaller integers
            # are examined exhaustively to identify the first viable depth.
            tests=[]
            for key,pole in [("lighter_neutral",d["M"][0]),("upper_neutral_color",d["M"][1]),
                             ("weak_vectorlike",bg["Mweak"])]:
                stage=(d["M"][1]**2+pole**2)/(d["M"][0]**2+pole**2)
                max_depth=math.ceil(abs(math.log(fixed["required_f1_over_f2"]))/math.log(stage))
                for n in range(1,max_depth+1):
                    for sign in [1,-1]:
                        h=stage**(sign*n)
                        tests.append(dict(mass_choice=key,stages=n,largest_gate_fast=(sign==1),h=h,
                                          solutions=clock_reactor_solutions(W,T,h,obs)))
            out[name]["reusing_existing_poles"]=tests
            out[name]["unamplified_gate_mass_ratio_upper"]=np.linalg.cond(W)**2*(val[1]/val[0])**2
            natural=val[1]/val[0]
            max_depth=math.ceil(math.log(fixed["required_f1_over_f2"])/math.log(natural))
            first=None
            for depth in range(1,max_depth+1):
                response=clock_reactor_solutions(W,T,natural**depth,obs)
                passed=[a for a in response if a["inside_angles_and_FC90"]]
                if passed:
                    first=dict(first_passing_depth=depth,gain=natural**depth,solutions=passed)
                    break
            out[name]["repeat_compressed_color_gate"]=dict(
                single_stage_ratio=natural,
                minimum_depth_from_clocks_only=math.ceil(math.log(np.sqrt(
                    np.sqrt(1+.00219/.0000773))/np.linalg.cond(W))/math.log(natural)),
                ceiling_depth_for_exact_target_gain=max_depth,first_joint_example=first)
    return dict(weak_GeV=bg["Mweak"],color_GeV=bg["Mcolor"],down_Yukawa_at_weak_ascending=m,
                c=c,f_half=f,heavy_poles_GeV=d["M"],source_ub=ub,
                target_t=chosen["physical_mu22_over_mu11"],target_angles=chosen["angles"],
                target_passed_original_all_angle_and_FC=all(chosen["angles_in_box"]) and chosen["SK_T2K_68_90"][-1],
                projected_color_forms=out)


def tidy(v):
    if isinstance(v,dict):
        return {k:tidy(x) for k,x in v.items()}
    if isinstance(v,np.ndarray):
        return tidy(v.tolist())
    if isinstance(v,(list,tuple)):
        return [tidy(x) for x in v]
    if isinstance(v,(float,np.floating)):
        return float(f"{v:.10g}")
    if isinstance(v,np.bool_):
        return bool(v)
    return v


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    old=runpy.run_path(str(HERE/"../1758/check.py"))
    out=tidy({"inclusive_primary":study(old["WEAK"]["inclusive_ub"][0],.25),
              "two_clock_profile":study(old["PROFILE"]["Vub"],.25)})
    dest=HERE/"same_color_results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==out
        print("1762 same-color check OK: tagged SU(5) right down frame and full light-mixing reconstruction")
    else:
        dest.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf-8")
