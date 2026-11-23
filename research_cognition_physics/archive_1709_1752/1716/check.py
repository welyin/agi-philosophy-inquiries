"""Projected inverse of the rank-two Majorana form on the physical 1709 null plane."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent

def read(f):
    return json.loads((HERE/f).read_bytes())

CFG=read("inputs.json")
FLOW=runpy.run_path(str(HERE/CFG["transport"]))

def at_cluster(flow,rep,step):
    y=np.zeros((4,3,4),complex)
    D=flow["low_ckm"](rep["gamma_MZ_deg"],flow["angles"]("direct"))
    y[0,:,:3]=np.diag(flow["central"]("up"))
    y[1,:,:3]=D@np.diag(flow["central"]("down"))
    y[2,:,:3]=FLOW["unpack"](rep["charged_left_MZ"])@np.diag(flow["central"]("charged_leptons"))
    poles=flow["fields"]["neutral_poles_GeV"]
    for start,susy,dt,gg in flow["segments"](step):
        assert not susy
        for j,pole in enumerate(poles):
            if start==pole:
                y[3,:,j]=flow["physics"]["Y0"][:,j]
        if start==poles[2]:
            return y[3][:,[0,2]]
        for k in range(0,len(gg)-1,2):
            beta=flow["beta"]
            k1=beta(y,gg[k],False)
            k2=beta(y+dt*k1/2,gg[k+1],False)
            k3=beta(y+dt*k2/2,gg[k+1],False)
            k4=beta(y+dt*k3,gg[k+2],False)
            y+=dt*(k1+2*k2+2*k3+k4)/6
    raise ValueError("missing upper pole")

def original_target(step):
    old=runpy.run_path(str(HERE/"../1712/check.py"))
    api=old["MODEL"]["sim"]
    original=api["lower_tensors"]
    cache={}
    def record(*args):
        cache["mass"]=original(*args)
        return cache["mass"]
    api["lower_tensors"]=record
    params=old["CFG"]["points"]["interior_all_one_width"]
    out=old["run_point"](**params,step=step)
    return cache["mass"],out["quark_neutrino"]

def background(step):
    flow,rep=FLOW["charged_background"](step)
    fh,_=FLOW["congruence"](flow,rep,step)
    D=FLOW["unpack"](rep["down_texture_10PeV"])
    A=np.array(rep["up_texture_diag_10PeV"])
    N=FLOW["unpack"](rep["neutral_columns_10PeV"])
    U,S,_=np.linalg.svd(N)
    v=U[:,-1]
    assert S[2]<1e-12*S[1] and np.linalg.norm(N.conj().T@v)<1e-13
    labels=read(CFG["family_labels"])
    a,b=np.array([labels["source_a"],labels["source_b"]])
    vacuum=read(CFG["flavor_vacuum"])
    angle=np.pi*vacuum["scalar_angle_over_pi"]
    phase=np.exp(1j*np.array([0,2*angle,3*angle+np.angle(3+2*np.exp(4j*angle))]))
    dc=at_cluster(flow,rep,step).conj()
    ab=phase[:,None]*np.column_stack((a,b))
    coeff=np.linalg.lstsq(dc,ab,rcond=None)[0]
    np.testing.assert_allclose(dc@coeff,ab,atol=1e-10)
    low_ab=fh.T@N[:,[0,2]].conj()@coeff
    v0=np.cross(a,b)/phase
    v0/=np.linalg.norm(v0)
    n=np.linalg.solve(fh,v)
    n/=np.linalg.norm(n)
    Q=np.linalg.svd(n[:,None],full_matrices=True)[0][:,1:]
    W=D@fh
    X=Q.conj().T@np.linalg.solve(W,np.diag(np.sqrt(A)))
    pieces=np.einsum('ik,jk->kij',X,X)
    test=np.exp(1j*np.array([0,.32,-.71]))
    C=W.T@np.diag(1/(test*A))@W
    dec=C-np.outer(C@n,C@n)/(n@C@n)
    reduced=Q.T@dec@Q
    inv=np.tensordot(test,pieces,axes=(0,0))
    np.testing.assert_allclose(reduced@inv,np.eye(2),atol=2e-10)
    F,s,_=np.linalg.svd(dec)
    Ui,si,_=np.linalg.svd(inv)
    np.testing.assert_allclose(abs(Q@Ui)**2,abs(F[:,1::-1])**2,atol=5e-10)
    np.testing.assert_allclose(s[0]/s[1],si[0]/si[1],rtol=5e-10)
    mass,reference=original_target(step)
    e0=FLOW["unpack"](rep["charged_left_MZ"])
    goal=e0.T@mass@e0
    goal/=np.linalg.norm(goal,2)
    assert np.linalg.norm(goal@n)<1e-9
    np.testing.assert_allclose(np.linalg.svd(goal,compute_uv=False)[1],1/reference["mass_ratio_after_transport"],rtol=2e-9)
    return n,Q,pieces,Q.T@low_ab,Q.T@goal@Q,(W,A,low_ab,goal),{"cluster_to_10PeV_projective_sine":np.sqrt(max(0.,1-abs(np.vdot(v,v0))**2)),
                       "Dirac_lift_relative_residual": np.linalg.norm(dc@coeff-ab)/np.linalg.norm(ab),
                       "original_two_source_goal_null_residual": np.linalg.norm(goal@n),
                       "reference_angles":reference["angles"],
                       "reference_ratio":reference["mass_ratio_after_transport"],
                       "null_charged_probabilities":abs(n)**2,
                       "projected_rank_one_norms":np.linalg.norm(X,axis=0)**2}

def scan(Q,pieces,axis=96):
    pieces=pieces[np.argsort(-np.linalg.norm(pieces,axis=(1,2)))]
    phases=np.exp(2j*np.pi*(np.arange(axis)+.5)/axis)
    ph=np.stack(np.meshgrid(phases,phases,indexing="ij"),axis=-1)
    mat=pieces[0]+ph[:,:,0,None,None]*pieces[1]+ph[:,:,1,None,None]*pieces[2]
    Ui,sv,_=np.linalg.svd(mat)
    mix=abs(Q@Ui)**2
    R=sv[:,:,0]/sv[:,:,1]
    z=mix[:,:,0,1]
    ang=np.stack([4*z*(1-z),mix[:,:,0,0]/(1-z),mix[:,:,1,1]/(1-z)],axis=-1)
    return sv,mix,ang,R

def adj(m):
    return np.array([[m[1,1],-m[0,1]],[-m[1,0],m[0,0]]])

def vec(m):
    return np.array([m[0,0],np.sqrt(2)*m[0,1],m[1,1]])

def one_dyad(pieces,h,goal):
    basis=np.column_stack([vec(adj(m)) for m in pieces])
    target=np.linalg.solve(basis,vec(goal))
    h=h/np.linalg.norm(h)
    dyad=np.outer(h,h)
    direction=np.linalg.solve(basis,vec(dyad))
    def diff(x):
        return x[1:]-x[0]
    w=diff(target.conj()*direction)
    line=np.column_stack((-2*w.real,2*w.imag))
    constant,radial=np.linalg.solve(line,np.column_stack((-diff(abs(target)**2),-diff(abs(direction)**2)))).T
    poly=[radial@radial,2*radial@constant-1,constant@constant]
    roots=sorted(np.roots(poly),key=lambda number:-number.real)
    witnesses=[]
    for t in roots:
        if abs(t.imag)>1e-10 or t.real<0:
            continue
        xy=constant+radial*t.real
        lam=complex(*xy)
        z=target-lam*direction
        np.testing.assert_allclose(abs(z/z[0]),np.ones(3),atol=1e-9)
        raw=sum(zz*adj(m) for zz,m in zip(z,pieces))
        np.testing.assert_allclose(raw+lam*dyad,goal,atol=1e-10)
        inv=sum(zz*m for zz,m in zip(z/z[0],pieces))
        scale=z[0]*np.linalg.det(inv)
        np.testing.assert_allclose(scale*np.linalg.inv(inv),raw,atol=1e-10)
        witnesses.append({"rank_one_weight_in_units_target_norm":abs(lam),
                          "dimensionless_mass_coeff_complex":[lam.real,lam.imag],
                          "relative_phase_to_full_Schur_rad":np.angle(lam/scale),
                          "named_ray_coefficient_over_Schur_scale":[(lam/scale).real,(lam/scale).imag],
                          "two_inverse_row_phases_rad":np.angle(z[1:]/z[0]),
                          "relative_coefficient_to_projected_form_norm":abs(lam)/np.linalg.norm(raw,2)})
    disc=(poly[1]**2-4*poly[0]*poly[2])/(poly[1]**2+4*poly[0]*poly[2])
    return {"signed_discriminant_over_sum":disc,"witnesses":witnesses}

def bound(n,Q,pieces,count):
    s,prob,angles,R=scan(Q,pieces,count)
    w=np.sort(np.linalg.norm(pieces,axis=(1,2)))
    delta=2*np.sin(np.pi/(2*count))*sum(w[:-1])+1e-7
    eta=(2*s[:,:,0]+delta)*delta
    gap=s[:,:,0]**2-s[:,:,1]**2
    sine=eta/(gap-eta)
    assert np.max(sine)<1
    def row_interval(idx):
        a=np.sqrt(prob[:,:,idx,1]); b=np.sqrt(prob[:,:,idx,0])
        lo=np.maximum(0,a*np.sqrt(1-sine*sine)-b*sine)**2
        hi=np.minimum(1,(a+b*sine)**2)
        return lo,hi
    zl,zh=row_interval(0); ml,mh=row_interval(1)
    assert zh.max()<.5
    rlo=np.min((s[:,:,0]-delta)/(s[:,:,1]+delta))-1e-7
    window=read(CFG["frequencies"])
    def ends(key):
        value,d,u=window[key]
        return value-d,value+u
    solar,d31,d32=[ends(k) for k in ("juno_NO_dm21","juno_NO_dm31","accelerator_NO_dm32")]
    rcap=min(math.sqrt(d31[1]/solar[0]),math.sqrt(1+d32[1]/solar[0]))
    enclosure=np.array([[np.min(4*zl*(1-zl)),np.max(4*zh*(1-zh))],
                        [np.min(1-abs(n[0])**2/(1-zh)),np.max(1-abs(n[0])**2/(1-zl))],
                        [np.min(ml/(1-zl)),np.max(mh/(1-zh))]])+np.array([-1e-7,1e-7])
    obs=np.array(read(CFG["angles"])["joint_windows"])
    assert enclosure[0,1]<obs[0,0] and enclosure[1,0]>obs[1,1] and enclosure[2,1]<obs[2,0] and rlo>rcap
    return {"cells":count**2,"ratio_lower":rlo,"ratio_frequency_upper":rcap,
            "joint_three_angle_enclosures":enclosure,"grid_center_ratio_min":np.min(R),
            "raw_center_angle_ranges":np.stack([angles.min(axis=(0,1)),angles.max(axis=(0,1))],axis=1),
            "necessary_inverse_mass_relative_norm":(1-rcap/rlo)/(1+rcap),
            "two_plane_congruence_necessary_condition_number":math.sqrt(rlo/rcap)}

def weak_observables(matrix):
    u,s,_=np.linalg.svd(matrix)
    m=u[:,::-1].conj()
    prob=abs(m)**2
    z=prob[0,2]
    scale=math.sqrt(CFG["solar_normalization_eV2"])/s[1]
    return {"angles":[4*z*(1-z),prob[0,1]/(1-z),prob[1,2]/(1-z)],
            "J":(m[0,0]*m[1,1]*np.conj(m[0,1]*m[1,0])).imag,
            "mass_ratio":s[0]/s[1],"dm31_eV2":(s[0]*scale)**2,
            "dm32_eV2":(s[0]*scale)**2-(s[1]*scale)**2,"mbb_eV":abs(matrix[0,0])*scale}

def calculate(step,count):
    n,Q,p,h,target,data,sc=background(step)
    W,A,low_ab,original=data
    limits=bound(n,Q,p,count)
    before=weak_observables(original)
    tests0=[one_dyad(p,h[:,i],target) for i in range(2)]
    ug,sg,vhg=np.linalg.svd(target)
    major=[s*np.outer(ug[:,j],vhg[j]) for j,s in enumerate(sg)]
    def solutions(deg):
        form=major[0]+np.exp(1j*np.radians(deg))*major[1]
        return one_dyad(p,h[:,1],form)["witnesses"]
    result=[]
    for branch in (0,1):
        lo,hi=CFG["real_coupling_root_bracket_deg"]
        def phase_part(deg):
            return solutions(deg)[branch]["named_ray_coefficient_over_Schur_scale"][1]
        start,end=phase_part(lo),phase_part(hi)
        assert start*end<0
        for _ in range(44):
            mid=(lo+hi)/2
            if phase_part(mid)*start>0:
                lo=mid
            else:
                hi=mid
        root=(lo+hi)/2
        witness=solutions(root)[branch]
        real,imag=witness["named_ray_coefficient_over_Schur_scale"]
        assert abs(imag)<1e-10*abs(real)
        z=np.exp(1j*np.r_[0.,witness["two_inverse_row_phases_rad"]])
        C=W.T@np.diag(1/(z*A))@W
        dec=C-np.outer(C@n,C@n)/(n@C@n)
        # Enforce a purely real high-scale coefficient in direct 3x3 forward propagation.
        rho=real/np.linalg.norm(h[:,1])**2
        form=dec+rho*np.outer(low_ab[:,1],low_ab[:,1])
        actual=weak_observables(form)
        np.testing.assert_allclose(actual["angles"],before["angles"],atol=1e-10)
        np.testing.assert_allclose([actual["J"],actual["mass_ratio"]],[before["J"],before["mass_ratio"]],atol=1e-10)
        target3=Q.conj()@(major[0]+np.exp(1j*np.radians(root))*major[1])@Q.conj().T
        multiplier=np.vdot(form,target3)/np.vdot(form,form)
        assert np.linalg.norm(multiplier*form-target3,2)<1e-10
        f=read(CFG["frequencies"])
        for val,name in [(CFG["solar_normalization_eV2"],"juno_NO_dm21"),(actual["dm31_eV2"],"juno_NO_dm31"),(actual["dm32_eV2"],"accelerator_NO_dm32")]:
            cen,d,u=f[name]
            assert cen-d<val<cen+u
        result.append({"relative_Majorana_shift_deg":root,"relative_real_coefficient_on_high_lifted_b":rho,
                       "full_forward":actual,"two_inverse_row_phases_rad":witness["two_inverse_row_phases_rad"],
                       "relative_rank_one_norm":witness["relative_coefficient_to_projected_form_norm"],
                       "imaginary_endpoints_of_root_bracket":[start,end]})
    return {"background":sc,"no_additional_term":limits,"original_goal":before,
            "fixed_original_Majorana_phase_a_b_discriminants":[row["signed_discriminant_over_sum"] for row in tests0],
            "real_b_port_solutions":result}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate(CFG["step"],CFG["torus_axis_cells"])
    coarse=calculate(CFG["crosscheck_step"],CFG["torus_axis_cells"])
    errors=[abs(result["no_additional_term"]["ratio_lower"]-coarse["no_additional_term"]["ratio_lower"])]
    for a,b in zip(result["real_b_port_solutions"],coarse["real_b_port_solutions"]):
        errors.append(abs(a["relative_Majorana_shift_deg"]-b["relative_Majorana_shift_deg"]))
    assert max(errors)<1e-6
    result["half_step_max_ratio_or_angle_degrees_change"]=max(errors)
    result=FLOW["tidy"](result)
    dest=HERE/"results.json"
    if args.check:
        assert result==json.loads(dest.read_bytes())
        print("1716 projected two-mode exclusion and positive real b-port reconstruction verified")
    else:
        dest.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
