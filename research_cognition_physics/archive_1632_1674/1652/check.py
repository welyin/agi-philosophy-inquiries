"""Absolute dimuon parity projections on the seven physical radial responses."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
EDGES=((0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2))
C=np.array([[0,0,0,0,0,-1,1],[0,0,-1,1,-1,0,1]],float)


def read(p):
    return json.loads((HERE/p).read_text(encoding="utf-8-sig"))


def box(a):
    return a[0]-2*sum(a[1:]),a[0]+2*sum(a[1:])


def texture(cfg):
    src=read(cfg["texture"])["reconstruction"]
    r,d=src["two_real_coefficients"]
    e=np.asarray([*src["normalized_real_edges_A_B_C_E_F"],r,d*cmath.exp(1j*math.radians(src["eta_example_deg"]))])
    Y=np.zeros((3,3),complex)
    P=np.zeros((7,6))
    for n,(pos,v) in enumerate(zip(EDGES,e)):
        Y[pos]+=v
        P[n,[pos[0],3+pos[1]]]=1
    ul,m,vh=np.linalg.svd(Y)
    ix=np.argsort(m); l,rr,m=ul[:,ix],vh.conj().T[:,ix],m[ix]
    err=max(np.linalg.norm(Y@rr-l*m),np.linalg.norm(Y.conj().T@l-rr*m),
            np.linalg.norm(l.conj().T@l-np.eye(3)),np.linalg.norm(rr.conj().T@rr-np.eye(3)))
    gap=min(np.diff(m).min(),2*m[0])
    assert err<1e-14 and gap>.00188
    return Y,l,rr,m,e,P,err,gap


def readers(l,r,e):
    cols=[]
    for (left,right),c in zip(EDGES,e):
        g=c*np.outer(l[left].conjugate(),r[right])
        minus=g-g.conj().T; plus=g+g.conj().T
        a=np.array([minus[0,2],minus[1,2],plus[1,2]])
        cols.append(np.column_stack([a.real,a.imag]).flatten())
    return np.array(cols).T


def form_floor(cf,k,mratio):
    # Including the squared Higgs propagator gives an integral in GeV^3.
    d=k["kinematics_GeV"]; fit=k["HPQCD_table_V"]
    b,ka,mu=d["Bplus"],d["Kplus"],d["mu"]
    tplus=(fit["B_average_GeV"]+fit["K_average_GeV"])**2
    def z(t):
        a=math.sqrt(tplus-t);c=math.sqrt(tplus)
        return (a-c)/(a+c)
    def lam(t):
        return ((b+ka)**2-t)*((b-ka)**2-t)
    lo,hi=k["q2_bin_GeV2"]; h=(hi-lo)/k["integration_cells"]
    common=(b*b-ka*ka)**2/(128*math.pi**3*b**3)
    total=0.
    for i in range(k["integration_cells"]):
        a=lo+i*h; c=a+h
        # z<0, and the lower endpoint of a2 is negative in this data set.
        poly=box(fit["a0"])[0]+box(fit["a1"])[1]*z(c)+box(fit["a2"])[0]*z(c)**2
        assert poly>0 and box(fit["a2"])[0]<0
        f=box(fit["L"])[0]*poly/(1-a/box(fit["pole_mass_GeV"])[1]**2)
        total+=h*math.sqrt(lam(c))*a*(1-4*mu**2/a)**1.5*f**2/(cf["mh_GeV"]**2-a)**2
    return common*total/(1-mratio)**2*(1-1e-11)


def bounds(cfg,m):
    q,rare,s=[read(cfg[x]) for x in ("semileptonic","dimuon","strange")]
    v,mh=cfg["vev_GeV"],cfg["mh_GeV"]
    mm=rare["const_GeV"]["m_mu"]
    common=(mm/(2*v*v))**2/q["hbar_GeV_ps"]
    transfer=cfg["scalar_amp_transfer_floor"]**2
    coeff=[]
    for i,(mes,f,tau,w) in enumerate((
        ("M_Bd",q["work_fBd_floor_GeV"],q["work_tau_Bd_floor_ps"],q["work_CMS_efficiency_floor"]),
        ("M_Bs",box(s["fBs_GeV"])[0],box(rare["tau_ps"]["s_heavy"])[0],s["run3_working_light_response"]))):
        M=rare["const_GeV"][mes]; beta=math.sqrt(1-4*mm**2/M**2)
        coeff.append(transfer*w*common*f*f*M**5*beta**3*tau/
                     (8*math.pi*(mh**2-M*M)**2*(1+m[i]/m[2])**2))
    val=form_floor(cfg,q,m[1]/m[2])
    coeff.append(transfer*common*cfg["K_bin_rate_efficiency_floor"]*cfg["charge_B_lifetime_floor_ps"]*val)
    b3=rare["LHCb_run3_2024_25"]["s_in_1e9"]
    bd=q["CMS_Bd_rate_upper_95"];bs=(b3[0]+2*(b3[2]+b3[4]))*1e-9
    blo,bhi=q["q2_bin_GeV2"]
    bk=box(q["LHCb_K_differential_per_GeV2"])[1]*(bhi-blo)
    coeff=np.array(coeff)
    cap=np.array([bd,bs,bk])
    return coeff,cap,np.sqrt(cap/coeff),val


def sensitivity(dual,radii,contrasts):
    """Radial distance in the product of three amplitude disks, after profiling ports."""
    ncell=65536
    step=math.pi/ncell
    phi=(np.arange(ncell)+.5)*step
    v=np.column_stack([np.cos(phi),np.sin(phi)])
    blocks=dual.reshape(2,3,2)
    support=np.zeros(ncell)
    slope=0.
    for j,rad in enumerate(radii):
        block=blocks[:,j,:]
        support+=rad*np.linalg.norm(v@block,axis=1)
        slope+=rad*np.linalg.norm(block,2)
    numerator=abs(v@contrasts)
    upper=np.max((numerator+np.linalg.norm(contrasts)*step/2)/(support-slope*step/2))
    low=np.max(numerator/support)
    assert min(support)>slope*step/2 and 0.<low<=upper<1
    return [float(low)*(1-1e-11),float(upper)*(1+1e-11)]


def calculate():
    cfg=read("inputs.json")
    Y,l,r,m,e,P,err,gap=texture(cfg)
    T=readers(l,r,e)
    singular=np.linalg.svd(T,compute_uv=False)
    neutral=T[:4]
    a=np.linalg.svd(neutral@P,compute_uv=False)
    assert np.linalg.norm(C@P)<1e-12 and np.linalg.norm(T@np.ones(7))<1e-13
    # A 1e-12 dilation residual and a spectral isolation of 0.00188
    # yield <2e-7 operator norm uncertainty for T, even rounding all 42 entries.
    ep=2e-7
    assert a[-1]>np.linalg.norm(P,2)*ep and singular[-1]>ep
    phases=np.array([.15,-.74,.33])
    T2=readers(l*np.exp(1j*phases),r*np.exp(1j*phases),e)
    np.testing.assert_allclose(T.T@T,T2.T@T2,atol=3e-15)
    dual=C@np.linalg.pinv(T)
    assert np.max(abs(dual@T-C))<1e-9
    coeff,cap,b,ff=bounds(cfg,m)
    bfh=b.copy(); bfh[-1]*=math.sqrt(cfg["optional_FH_upper"])
    work={}
    for label,lim in (("rate_only",b),("with_angular_template",bfh)):
        dual_norms=np.linalg.norm(dual.reshape(2,3,2),axis=2)
        terminal=dual_norms@lim
        support_pad=np.linalg.norm(dual,axis=1)*ep/(singular[-1]-ep)*np.linalg.norm(lim)
        assert np.all(terminal>0)
        test=[]
        for i,row in enumerate(dual):
            points=(row.reshape(3,2).T*(lim/dual_norms[i])).T.flatten()
            trial=np.linalg.pinv(T)@points
            test.append(float((C@trial)[i]))
        np.testing.assert_allclose(test,terminal,rtol=1e-10)
        g0,dg=cfg["atlas_signal_strength"]
        kmin=math.sqrt((g0-2*dg)/cfg["optional_production_over_width_upper"])
        work[label]={"abs_kmu_times_projection_upper":lim.tolist(),
                "support_abs_kmu_contrasts_nominal":terminal.tolist(),
                "upper_with_numeric_enclosure":(terminal+support_pad).tolist(),
                "if_higgs_production_width_template_abs_contrasts_upper":((terminal+support_pad)/kmin).tolist(),
                "nominal_minimum_common_radius_fraction_for_contrasts_at_kmu_one":{
                    name:sensitivity(dual,lim,np.array(values)) for name,values in
                    (("phase_path",[2.,2.]),("real_path",[-2.,0.]))}}
    radial=read(cfg["radial"])["example_slopes"]; points={}
    for name in ("phase_path_dim6","real_path_dim6","universal_linear"):
        n=np.array(radial[name])
        # Project onto the affine fibre C n = constant and impose both B amplitudes zero.
        x=np.linalg.lstsq(neutral@P,-neutral@n,rcond=1e-13)[0]
        nx=n+P@x
        assert np.linalg.norm(neutral@nx)<2e-14 and np.max(abs(C@nx-C@n))<1e-12
        val=np.linalg.norm((T@n).reshape(3,2),axis=1)
        changed=np.linalg.norm((T@nx).reshape(3,2),axis=1)
        points[name]={"unit_kmu_absolute_projections":val.tolist(),
                      "minimum_transferred_scalar_BF_unit_kmu":(coeff*val**2).tolist(),
                      "neutral_blind_slopes":nx.tolist(),"neutral_blind_projections":changed.tolist(),
                      "fixed_contrasts":(C@n).tolist()}
    # No raw SVD phase is needed for the physical bound: each pair rotates together.
    return {"central_normalized_masses":m.tolist(),"spectrum_gap":gap,"matrix_residual":err,
            "neutral_endpoint_singular_values":a.tolist(),"six_projection_singular_values":singular.tolist(),
            "dual_group_norms":np.linalg.norm(dual.reshape(2,3,2),axis=2).tolist(),
            "form_factor_integral_with_propagator_and_mass_ratio_GeV3":ff,
            "rate_floors_per_unit_kmu_times_projection_squared":coeff.tolist(),"experimental_rate_upper":cap.tolist(),
            "limits":work,"physical_examples":points}


if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--check",action="store_true");p.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    def tidy(x):
        if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
        if isinstance(x,list):return list(map(tidy,x))
        if isinstance(x,float):return 0.0 if abs(x)<1e-18 else float(format(x,".10g"))
        return x
    out=tidy(calculate())
    if args.check:
        assert out==read("results.json")
    elif args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2);f.write("\n")
    else:
        print(json.dumps(out,indent=2))
    print("1652 PASS: charged parity projection, neutral nulls and finite global support")
