"""Match a named SU(5) five-bar flavor contraction to charged and neutral leptons."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())


CFG = read("inputs.json")
TAGS = runpy.run_path(str(HERE/CFG["label_transport"]))
I = np.eye(3, dtype=complex)


def physical_inputs(whole, step):
    model = runpy.run_path(str(HERE/CFG["parent"]))
    if whole:
        model["CFG"]["new_Dirac_triplet_beta"] = CFG["new_Dirac_complete_five_beta"]
        model["CFG"]["new_two_scalar_triplet_beta"] = CFG["new_two_scalar_complete_five_beta"]
    data = model["CFG"]["nearby_measured_inputs_in_marginal_widths"][CFG["charged_neutral_point"]]
    summary, pt, nu = model["charged_endpoint"](2,step,*(data[k] for k in ("weak","up","g3","top","charm")))
    flow = pt["flow"]
    rep = dict(pt["rep"],down_texture_10PeV=flow["pack"](pt["D"]))
    D, u, _ = TAGS["tagged_matrices"](flow,rep,step)
    left, vals, right = TAGS["svd"](D)
    c = abs(np.vdot(u,left[:,0]))
    fb = math.sqrt((1+c)/2)
    Fd = (right*np.array([c,1.,fb]))@right.conj().T
    poles = flow["fields"]["neutral_poles_GeV"]
    segments = [a for a in flow["segments"](step) if a[0]<poles[2]]
    y = model["initial_state"](pt)
    S = I.copy()
    for start,susy,dt,g in segments:
        assert not susy
        for j,pole in enumerate(poles):
            if start==pole:
                y[3,:,j]=flow["physics"]["Y0"][:,j]
        y,S=model["advance"](y,S,dt,g,flow)
    _, ds, vre= TAGS["svd"](y[1,:,:3])
    np.testing.assert_allclose(ds,vals,rtol=3e-12,atol=1e-15)
    # D and E share the right (5bar) GJ index. (E Fd)^T = Fd^T E^T.
    top = flow["basis"](y[0])[0][:,0]
    chk = abs(np.vdot(top,flow["basis"](y[1])[0][:,0]))
    assert abs(chk-c)<1e-11
    # Check that the registered lepton simulation really uses E.T's left basis,
    # modulo the column phases of the charged mass eigenstates.
    pack = rep["charged_left_10PeV"]
    target = np.array(pack["real"])+1j*np.array(pack["imag"])
    np.testing.assert_allclose(abs(target.conj().T@flow["basis"](pt["E"].T)[0]),I,atol=2e-9)
    return model,summary,pt,nu,segments,y,S,Fd.T,c,fb,vre


def evaluate(background, common, override=None):
    model,summary,pt,nu,segments,y0,S,Fknown,c,fb,vref = background
    flow = pt["flow"]
    poles = flow["fields"]["neutral_poles_GeV"]
    F=override if override is not None else (Fknown if common else I)
    y=y0.copy()
    y[1,:,:3] = y[1,:,:3]@((vref*[c,1.,fb])@vref.conj().T)
    y[2] = F@y[2]
    y[3] = F@y[3]
    Tparent=S@pt["low_E"]
    inv=np.linalg.inv(Tparent)
    Cparent=inv.T@pt["C"]@inv
    bparent=inv.T@pt["b"]
    F_source=F.T@Cparent@F
    b_source=F.T@bparent
    rho=nu["FROZEN"][0]
    Mparent=Cparent+rho*np.outer(bparent,bparent)
    _,_,v0=np.linalg.svd(Mparent)
    null=v0.conj().T[:,-1]
    nullnew=np.linalg.solve(F,null)
    assert np.linalg.norm((F_source+rho*np.outer(b_source,b_source))@nullnew)<1e-12
    # Raw Dirac columns: columns 0 and 2 represent actual couplings at the
    # threshold (lower active, upper at its insertion); these are distinct
    # from the low completed effective Majorana terms C and b b.T.
    N=np.column_stack([y0[3,:,0],flow["physics"]["Y0"][:,2]])
    Q=np.linalg.qr(N)[0][:,:2]
    gram=np.linalg.eigvalsh(Q.conj().T@F.conj().T@F@Q)
    ye0=flow["basis"](y0[2])[1]
    ye1=flow["basis"](y[2])[1]
    assert abs(np.prod(ye1/ye0)-abs(np.linalg.det(F)))<1e-12
    Tdown=I.copy()
    res=[]
    for start,_,dt,g in reversed(segments):
        y,Tdown=model["advance"](y,Tdown,-dt,g[::-1],flow)
        for j,pole in enumerate(poles):
            if start==pole:
                old=F@flow["physics"]["Y0"][:,j]
                if np.linalg.norm(old):
                    res.append(np.linalg.norm(y[3,:,j]-old)/np.linalg.norm(old))
                y[3,:,j]=0.
    q,leptons,gamma,L,ups=flow["extraction"](y)
    phase=np.diag(L.conj().T@pt["low_E"])
    L=L@np.diag(phase/abs(phase))
    U=np.linalg.solve(Tdown,L)
    C=U.T@F_source@U
    b=U.T@b_source
    H=pt["low_E"].conj().T@np.linalg.solve(S,F@U)
    np.testing.assert_allclose(C,H.T@pt["C"]@H,atol=1e-13)
    np.testing.assert_allclose(b,H.T@pt["b"],atol=1e-12)
    obs=nu["observables"](C,b,rho)
    ref_obs=nu["observables"](pt["C"],pt["b"],rho)
    old=model["JOINT"]["LEPTONS"]
    normalized=leptons/old[:,0]
    tau_fit=1/normalized[-1]
    frac=2*old[:,1]/old[:,0]
    ratio_window=np.column_stack(((1-frac[:2])/(1+frac[-1])-1,(1+frac[:2])/(1-frac[-1])-1))
    mu_e=(normalized[1]/normalized[0]-1)
    obs["dm21_with_mother_solar_gain_eV2"]=nu["SOL"]*(ref_obs["scale_eV"]/obs["scale_eV"])**2
    obs["dm31_with_mother_solar_gain_eV2"]=obs["dm31_eV2"]*obs["dm21_with_mother_solar_gain_eV2"]/nu["SOL"]
    mass=C+rho*np.outer(b,b)
    _,_,v=np.linalg.svd(mass)
    nullweight=abs(v.conj().T[:,-1])**2
    return dict(c=c,fb=fb,parent_tan_beta=summary["positive_soft_tan_beta"],parent_Ye_at_M=ye0,
                squared_filter_elements_in_charged_eigenframe=abs(flow["basis"](y0[2])[0].conj().T@F@flow["basis"](y0[2])[0])**2,
                threshold_lepton_factors=ye1/ye0,threshold_volume_factor=np.prod(ye1/ye0),
                threshold_neutral_gram_singular_gains=np.sqrt(gram),
                threshold_null_line_angle_deg=np.degrees(np.arccos(np.clip(abs(np.vdot(null,nullnew))/np.linalg.norm(nullnew),0,1))),
                lower_active_column_residual_relative_to_fixed_proxy=res,
                MZ_actual_leptons=leptons,MZ_lepton_reference=old[:,0],
                MZ_absolute_marginal_shifts=(leptons-old[:,0])/old[:,1],
                MZ_charged_product_relative=np.prod(normalized),MZ_common_factor_keeping_tau=tau_fit,
                MZ_relative_e_mu_after_tau_calibration=normalized[:2]*tau_fit-1,
                two_marginal_width_ratio_windows_e_mu_over_tau=ratio_window,
                MZ_relative_mu_over_e=mu_e,MZ_zero_line_probability=nullweight,neutrinos=obs,
                MZ_down=q[:3],MZ_gamma=gamma,MZ_up=ups)


def untagged_completion(background):
    """Fit an arbitrary five-bar identification after allowing two separate scalar gains."""
    _,_,pt,_,_,y,_,_,c,fb,_=background
    N=np.column_stack([y[3,:,0],pt["flow"]["physics"]["Y0"][:,2]])
    W=np.linalg.qr(N,mode="complete")[0]
    l,s,_=TAGS["svd"](y[2,:,:3])
    m=fb**2
    ell=c**2
    z=(c*fb)**(2/3)
    k=1+ell-m
    length=math.sqrt((m-ell)*(1-m))
    A=W.conj().T@(l*(s*s/s.dot(s)))@l.conj().T@W
    inv_s=1/s
    B=W.conj().T@(l*(inv_s**2/inv_s.dot(inv_s)))@l.conj().T@W
    v=A[:2,2]
    axis=np.r_[v.real,v.imag]
    axis/=np.linalg.norm(axis)
    level=(z-m-(k-m)*A[2,2].real)/(2*length*np.linalg.norm(v))
    assert abs(level)<1
    C=np.linalg.svd(axis[None,:],full_matrices=True)[2][1:,:].T
    def scalar(direction):
        vector=axis*level+C@direction*math.sqrt(1-level**2)
        w=length*(vector[:2]+1j*vector[2:])
        G=np.diag([m,m,k]).astype(complex)
        G[:2,2],G[2,:2]=w,w.conj()
        assert abs(np.trace(A@G)-z)<2e-10
        return np.trace(B@np.linalg.inv(G)).real-1/z,G
    trials=np.vstack((np.eye(3),-np.eye(3)))
    values=[scalar(t)[0] for t in trials]
    u,vv=trials[np.argmin(values)],trials[np.argmax(values)]
    assert min(values)<0<max(values)
    # If antipodes are selected, a great circle through a distinct axis is valid.
    angle=math.acos(np.clip(u@vv,-1,1))
    tangent=vv-u*(u@vv)
    if np.linalg.norm(tangent)<1e-10:
        tangent=np.eye(3)[np.argmin(abs(u))]
    tangent/=np.linalg.norm(tangent)
    lo,hi=0.,angle
    for _ in range(55):
        theta=(lo+hi)/2
        residual,G=scalar(u*math.cos(theta)+tangent*math.sin(theta))
        if residual<0: lo=theta
        else: hi=theta
    w,V=np.linalg.eigh(W@G@W.conj().T)
    F=(V*np.sqrt(w))@V.conj().T
    np.testing.assert_allclose(w,[c*c,fb*fb,1],atol=5e-13)
    _,new,_=TAGS["svd"](F@y[2,:,:3])
    np.testing.assert_allclose(new/s,np.sqrt(z),rtol=2e-9)
    np.testing.assert_allclose(W[:,:2].conj().T@F@F@W[:,:2],m*np.eye(2),atol=1e-13)
    return F,dict(forward_sphere_height=level,inverse_trace_sign_bracket=[min(values),max(values)],
                  inverse_trace_at_root=residual,threshold_charged_mass_ratios=new/s,
                  fitted_F_real=F.real,fitted_F_imag=F.imag)


def frequency_bound(full, control, common):
    model,_,pt,nu,segments,y,S,_,c,_,_=full
    Q=np.linalg.qr(np.column_stack([y[3,:,0],pt["flow"]["physics"]["Y0"][:,2]]))[0]
    inv=np.linalg.inv(S@pt["low_E"])
    M=inv.T@(pt["C"]+nu["FROZEN"][0]*np.outer(pt["b"],pt["b"]))@inv
    proj=Q@Q.conj().T
    leak=np.linalg.norm(M-proj.conj()@M@proj,2)/np.linalg.norm(M,2)
    assert leak<2e-13
    flow=pt["flow"]
    poles=flow["fields"]["neutral_poles_GeV"]
    total=math.log(poles[2]/flow["MZ"])/(16*math.pi**2)
    active=math.log(poles[2]/poles[0])/(16*math.pi**2)
    max_g=max(max(np.max(seg[3][:,:2],axis=0)) for seg in segments)
    assert max_g<.6
    # Dropping positive universal traces in the descending lepton beta gives
    # simple uniform norm bounds, including up to 1/c separate calibrations.
    e0=np.linalg.norm(y[2],2)/c
    n0=np.linalg.norm(y[3],2)/c
    assert e0*math.exp((4.5*.6+1.5*.25**2)*total)<.02
    assert n0*math.exp((2.7*.6+1.5*.02**2)*active)<.25
    L=1.5*.02**2*total+.5*.25**2*active
    split=control["neutrinos"]["dm31_eV2"]
    solar,serr,_=nu["FREQ"]["juno_NO_dm21"]
    atm,alo,_=nu["FREQ"]["juno_NO_dm31"]
    gain=common["threshold_neutral_gram_singular_gains"]
    anisotropy=gain[1]/gain[0]
    intervals={}
    for n in (1,2):
        # The extremely small projected residue is kept as an explicit numerical
        # allowance on the registered mass pencil (its m3/m2 is below 7).
        control_upper=(solar+n*serr)*split/solar*math.exp(8*L)*(1+2e-10)
        required=((atm-n*alo)/control_upper)**.25
        intervals[str(n)]=dict(max_dm31_eV2_for_isometric_source_plane=control_upper,
                   required_ratio_of_two_source_singular_gains=required,
                   max_dm31_eV2_for_named_bridge=control_upper*anisotropy**4,
                   lower_dm31_window_eV2=atm-n*alo)
    return dict(mass_outside_source_plane_relative=leak,
                charged_and_neutral_norms_at_threshold_before_gain=[e0*c,n0*c],
                max_g1g2_squared=max_g, log_condition_number_upper_each_leg=L,
                actual_named_source_plane_anisotropy=anisotropy,intervals=intervals)


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [tidy(t) for t in x]
    if isinstance(x,(np.floating,float)): return float(f"{x:.11g}")
    if isinstance(x,np.bool_): return bool(x)
    return x


def solve(step):
    color=physical_inputs(False,step)
    full=physical_inputs(True,step)
    out={"archived_color_split":evaluate(color,False),"complete_multiplets_split_couplings":evaluate(full,False),
         "complete_multiplets_common_5bar":evaluate(full,True)}
    fit,certificate=untagged_completion(full)
    out["unlabelled_bridge_matched_two_independent_scales"]=evaluate(full,True,fit)
    expected=read("../1739/results.json")["rows"][CFG["charged_neutral_point"]]
    np.testing.assert_allclose(out["archived_color_split"]["MZ_down"],expected["MZ_down"],rtol=1e-8)
    assert abs(out["archived_color_split"]["neutrinos"]["dm31_eV2"]-expected["neutrinos"]["dm31_eV2"])<3e-9
    c,fb=full[-3:-1]
    ae=(c*fb)**(1/3)
    shape_gap=math.log(fb/ae)
    bound=frequency_bound(full,out["complete_multiplets_split_couplings"],out["complete_multiplets_common_5bar"])
    return {"cases":out,"constructive_independent_gain_competitor":certificate,"oscillation_necessary_bound":bound,
            "all_unitary_bridges":{"volume":c*fb,"charged_scalar_if_all_masses_scaled":ae,
            "neutral_common_scalar_if_two_column_gram_scaled":fb,"log_scalar_separation":shape_gap,
            "half_log_incompatible_equal_gains":shape_gap/2}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    coarse,fine=(solve(v) for v in CFG["steps"])
    a,b=(v["cases"]["complete_multiplets_common_5bar"] for v in (coarse,fine))
    steps=dict(absolute_lepton_product=abs(a["MZ_charged_product_relative"]-b["MZ_charged_product_relative"]),
               dm31=abs(a["neutrinos"]["dm31_eV2"]-b["neutrinos"]["dm31_eV2"]))
    assert steps["absolute_lepton_product"]<1e-7 and steps["dm31"]<1e-8
    result=tidy(dict(**fine,step_check=steps))
    if args.check:
        assert result==read("results.json")
        print("1743 common 5bar charge and neutral thresholds verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")
