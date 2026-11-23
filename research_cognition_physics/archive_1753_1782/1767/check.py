"""Check finite weak-current matching of the actual pair of colour singlet-weak down quarks."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


CFG = read("inputs.json")
ROW = read(CFG["charged_points"])["rows"][CFG["main_point"]]
MASS = read(CFG["color_spectrum"])["pole_mass_GeV"][0]
SHORT = read(CFG["short_distance"])
GF = SHORT["const_GeV"]["GF"]
VEV = (math.sqrt(2)*GF)**-.5


def threshold():
    """Read the actual tensors just before the colour gate in the frozen integrator."""
    joint = runpy.run_path(str(HERE / CFG["model"]))
    point = joint["point"]
    old_filter = joint["matter_filter"]
    captured = {}

    def keep_filter(bg):
        F,z = old_filter(bg)
        captured.update(bg=bg,F=F)
        return F,z

    point.__globals__["matter_filter"] = keep_filter
    target = np.array([ROW["moduli"][k] for k in ("Vus","Vub","Vcb")])
    trial = target.copy()
    for _ in range(3):
        row = point(HERE,*trial,CFG["integrator_step"])
        physical = np.array([row["low_moduli"][k] for k in ("Vus","Vub","Vcb")])
        if max(abs(physical-target)) < 2e-12:
            break
        trial += target-physical
    assert max(abs(physical-target)) < 2e-12 and abs(row["gamma_deg"]-ROW["gamma"]) < 5e-7
    bg = captured["bg"]
    F = captured["F"]
    y = bg["yW"].copy()
    y[2] = row["tau_calibration"]*(F.conj().T@y[2])
    y[3] = F.conj().T@y[3]
    running = np.eye(3,dtype=complex)
    descending=iter(reversed(bg["segments"]))
    for start,susy,dt,gauge in descending:
        assert not susy
        y,running = bg["model"]["advance"](y,running,-dt,gauge[::-1],bg["pt"]["flow"])
        if start == MASS:
            break
    else:
        raise ValueError("Colour threshold missing")
    Uu,yu,_ = np.linalg.svd(y[0,:,:3])
    Ud,yd,Vhd = np.linalg.svd(y[1,:,:3])
    V = Uu[:,::-1].conj().T @ Ud[:,::-1]
    c = abs(V[0,0])
    assert abs(c-row["port_overlap_at_color"]) < 2e-12
    attenuation=np.array([c,1.,math.sqrt((1+c)/2)])[::-1]
    y[1,:,:3]=y[1,:,:3]@(Vhd.conj().T@np.diag(attenuation)@Vhd)
    poles=bg["pt"]["flow"]["fields"]["neutral_poles_GeV"]
    for j,m in enumerate(poles):
        if m==MASS:
            y[3,:,j]=0.
    for start,susy,dt,gauge in descending:
        y,running=bg["model"]["advance"](y,running,-dt,gauge[::-1],bg["pt"]["flow"])
        for j,m in enumerate(poles):
            if m==start:
                y[3,:,j]=0.
    up_low=np.linalg.svd(y[0,:,:3])[0][:,::-1]
    down_low=np.linalg.svd(y[1,:,:3])[0][:,::-1]
    Vlow=up_low.conj().T@down_low
    assert abs(abs(Vlow[0,1])-target[0])<2e-12
    return yu[::-1],yd[::-1],V,Ud[:,::-1].conj().T@down_low,Vlow


def pole(a,f):
    """Compute with positive excess over the unbroken mass to keep all significant digits."""
    M = MASS
    gap = math.hypot(M*M-a*a,2*M*a*math.sqrt(1-f*f))
    excess = 2*(M*a)**2*(1-f*f)/(gap+M*M-a*a)
    heavy = math.sqrt(M*M+excess)
    loss = excess/gap
    light = a*f*M/heavy
    mat = np.array([[a*f,a*math.sqrt(1-f*f)],[0.,M]])
    left,sv,_ = np.linalg.svd(mat)
    np.testing.assert_allclose(sv,[heavy,light],rtol=2e-12)
    np.testing.assert_allclose(abs(left[1,-1])**2,loss,rtol=2e-5,atol=1e-25)
    q = (light/M)**2
    np.testing.assert_allclose(loss,q*(1-f*f)/(f*f*(1-2*q)+q*q),rtol=2e-12)
    return dict(unprojected_mass_GeV=a,overlap=f,light_GeV=light,
                heavy_squared_excess_GeV2=excess,heavy_mass_increment_GeV=excess/(heavy+M),
                lost_doublet_probability=loss)


def charged_holdout():
    """Delta Ms / the tau width cancels the actual bottom-leg probability."""
    tau=read(CFG["charged_tau"])
    legacy=read(CFG["belle_and_frozen_endpoints"])
    points=read(CFG["comparison_points"])["points"]
    build=runpy.run_path(str(HERE/"../1763/check.py"))["matrix"]
    clocks=read("../1736/independent_results.json")["data"]
    models={"1765_Bs_calibrated":(ROW["moduli"]["Vub"],ROW["tbs_squared"])}
    for key in ("QCDSF_B_central","QCDSF_reprofiled_two_clocks"):
        inv=points[key]["CKM"]
        us=math.sqrt(1-inv["ud"]**2-inv["ub"]**2)
        V=build(us,inv["ub"],inv["cb"],inv["gamma"])
        models[key]=(inv["ub"],abs(V[2,1]*V[2,2])**2)
    mt=tau["mu_tau_masses_GeV"][1]
    mb=tau["meson_masses_GeV"]["B"]
    mw=SHORT["const_GeV"]["MW"]
    pref=3*math.pi/4*(mt/mw)**2*(1-(mt/mb)**2)**2*mb
    def limits(x,n):
        delta=n*math.hypot(*x[1:])
        return np.array([x[0]-delta,x[0]+delta])
    # ms is registered in MeV, lifetime in ps and the oscillation clock in ps^-1.
    assert clocks["mass_s"][0]>5000
    def unit(n):
        mass=limits(clocks["mass_s"],n)/1000
        f=limits(CFG["decay_ratio_fBs_over_fBplus"],n)
        bag=limits(SHORT["bag"]["s"],n)
        eta=limits(SHORT["short_distance"]["etaB"],n)
        S=limits(SHORT["short_distance"]["S0_published"],n)
        return (pref*limits(tau["B_lifetime_ps"],n)*limits(clocks["delta_ms"],n)
                /(mass*f**2*bag*eta*S)[::-1])
    results={}
    for key,(ub,prod2) in models.items():
        ratio=ub**2/prod2
        central=unit(0)[0]*ratio
        nominal=unit(1)*ratio
        rows={}
        for label,exp in (("Belle",legacy["BR_Belle_combined"]),("Belle_II",tau["b_branching"])):
            # For an independent multiplicative dispersive amplitude h_extra, B=requested/h_extra.
            rows[label]=dict(abs_h_extra_center=central/exp[0],
                             abs_h_extra_one_width=nominal/limits(exp,1)[::-1])
        results[key]=dict(ub2_over_tsb2=ratio,BR_if_h_extra_one=central,
                          BR_one_width_hadronic_clock=nominal,h_extra_from_data=rows)
    return results


def analyze():
    yu,yd,V,rotate,Vlow = threshold()
    c = abs(V[0,0])
    factors = np.array([c,1.,math.sqrt((1+c)/2)])
    a = yd*VEV/math.sqrt(2)
    d = pole(a[0],factors[0])
    b = pole(a[2],factors[2])
    eps = b["lost_doublet_probability"]
    losses=np.array([d["lost_doublet_probability"],0.,eps])
    weak_Z=rotate.conj().T@np.diag(losses)@rotate
    # Transport the finite tree projector in a fixed weak basis. No separate anomalous
    # dimension or finite matching of higher-dimensional coefficients is assumed here.
    half_loss=losses/(1+np.sqrt(1-losses))
    delta_W=-Vlow@(rotate.conj().T@np.diag(half_loss)@rotate)
    relative_top = ((1+delta_W[2,1]/Vlow[2,1]).conjugate()
                    *(1+delta_W[2,2]/Vlow[2,2]))**2
    relative_cb=abs(1+delta_W[1,2]/Vlow[1,2])**2
    W = V @ np.diag(np.sqrt(1-np.array([d["lost_doublet_probability"],0.,eps])))
    sm = V[:,1].conj()*V[:,2]
    with_gate = W[:,1].conj()*W[:,2]
    np.testing.assert_allclose(with_gate,math.sqrt(1-eps)*sm,rtol=3e-15)
    assert abs(sum(with_gate)) < 1e-14
    assert np.linalg.norm(W.conj().T@W-np.diag(np.diag(W.conj().T@W))) < 2e-15
    # All up masses and any common real box-kernel yield the same relative amplitude.
    box = np.array([[.00021,.002,.026],[.002,.052,.114],[.026,.114,2.313]])
    np.testing.assert_allclose((with_gate@box@with_gate)/(sm@box@sm),1-eps,atol=0,rtol=1e-14)
    np.testing.assert_allclose(abs(W[1,2]/V[1,2])**2,1-eps,atol=0,rtol=5e-15)

    C0,sigmaC = read(CFG["clock"])["mixing_coefficient_ps_inv"]
    needed = sigmaC/C0
    # Bound valid for arbitrary left/right flavour orientation of identical complete tree spectra.
    delta = d["heavy_squared_excess_GeV2"]+b["heavy_squared_excess_GeV2"]
    lower_gap = MASS**2+min(d["heavy_squared_excess_GeV2"],b["heavy_squared_excess_GeV2"])-b["light_GeV"]**2
    leakage = delta/lower_gap
    U_bound = leakage/2
    qweak = ROW["tbs_squared"]
    sd = SHORT["short_distance"]
    G = GF**2*SHORT["const_GeV"]["MW"]**2*sd["S0_published"][0]*sd["etaB"][0]*qweak/(4*math.pi**2)
    etaZ = max(CFG["RGI_tree_Z_transport_range"])
    threshold_Usb = math.sqrt(needed*G/(GF*etaZ/math.sqrt(2)))
    min_pole_sq_excess = 2*threshold_Usb*(MASS**2-b["light_GeV"]**2)
    # sqrt(M^2+x)+sqrt(M^2+y)-2M >= sqrt(M^2+x+y)-M for x,y>=0.
    min_sum_mass_increase = min_pole_sq_excess/(math.sqrt(MASS**2+min_pole_sq_excess)+MASS)

    # Flavour-covariant right-handed spurion, evaluated below the mass threshold.
    daughter_y = yd*factors
    Hup = V.conj().T @ np.diag(yu**2) @ V
    A_sb = daughter_y[1]*daughter_y[2]*Hup[1,2]
    rgi = needed*G
    net = math.sqrt(1-eps)
    return dict(tau_vs_strange_B_clock=charged_holdout(),
        held_tree_coefficient_CKM_running=dict(abs_b_to_s_basis_overlap=abs(rotate[2,1]),
                      abs_Zsb=abs(weak_Z[1,2]),relative_top_box_minus_one_real=relative_top.real-1,
                      relative_top_box_phase_rad=np.angle(relative_top),
                      relative_cb_decay_rate_minus_one=relative_cb-1,
                      relative_top_box_to_cb_rate_minus_one=abs(relative_top)/relative_cb-1),
        matching=dict(mass_GeV=MASS,vev_GeV=VEV,high_yu=yu,daughter_yd=daughter_y,
                      c=c,threshold_CKM_sb_abs=abs(sm[-1]),gate_d=d,gate_b=b),
        aligned=dict(hs_box_real=1-eps,one_minus_hs=eps,phase_rad=0.,
                     fraction_lost_B_semileptonic_rate=eps,tree_effective_Vcb=ROW["moduli"]["Vcb"]*net,
                     effective_Vcb_shift=-ROW["moduli"]["Vcb"]*eps/(1+net),
                     physical_bs_frequency_shift_per_ps=-17.766*eps,
                     relative_mixing_to_b_semileptonic_rate=1.,
                     target_1765_fraction=needed,larger_by_factor=needed/eps),
        same_spectrum_any_flavour=dict(total_leakage_max=leakage,abs_Zsb_projector_max=U_bound,
                     Z_tree_abs_fraction_max=(GF/math.sqrt(2)*etaZ/G)*U_bound**2,
                     required_projector_for_target_if_Z_alone=threshold_Usb,
                     sums_of_heavy_mass_increases_GeV=sum(x["heavy_mass_increment_GeV"] for x in [b,d]),
                     sum_mass_increases_GeV_needed_for_target_if_Z_alone=min_sum_mass_increase,
                     box_leg_operator_norm_max=leakage/(1+math.sqrt(1-leakage))),
        independent_RR=dict(RGI_SM_VLL_coefficient_abs_GeV_minus2=G,
                     RGI_new_VRR_required_abs_GeV_minus2=rgi,
                     unit_RGI_unsuppressed_scale_TeV=rgi**-.5/1000,
                     absolute_projected_spurion=abs(A_sb),
                     normalized_hadronic_Wilson_at_color_mass=rgi*MASS**2,
                     coefficient_if_raw_Yukawa_spurion_squared_used=rgi*MASS**2/abs(A_sb)**2,
                     max_fraction_if_spurion_coefficient_one=(abs(A_sb)/MASS)**2/G,
                     A_sb_squared_relative_to_top_phase_rad=float(np.angle((A_sb/sm[-1])**2))))


def clean(value):
    if isinstance(value,dict):
        return {k:clean(v) for k,v in value.items()}
    if isinstance(value,np.ndarray):
        return [clean(v) for v in value.tolist()]
    return float(f"{value:.11g}")


if __name__ == "__main__":
    cmd = argparse.ArgumentParser()
    cmd.add_argument("--check",action="store_true")
    args = cmd.parse_args()
    actual = clean(analyze())
    if args.check:
        if actual != read("results.json"):
            raise SystemExit("1767 saved result differs")
        print("1767 OK: aligned finite-v box vs tree, same-spectrum neutral-current bound and independent RR spurion.")
    else:
        (HERE/"results.json").write_text(json.dumps(actual,indent=2)+"\n",encoding="utf-8")

