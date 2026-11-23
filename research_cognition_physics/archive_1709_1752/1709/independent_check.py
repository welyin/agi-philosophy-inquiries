"""Separate the two neutral source strengths with the two flavour-summed heavy residues."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(rel):
    return json.loads((HERE/rel).read_bytes())


def archived(rel):
    return runpy.run_path(str(HERE/rel))


def source_basis():
    """Return Dirac matrices for the two real monomials x and 1, before any oscillation fit."""
    heavy = read("../1689/inputs.json")
    fields = read("../1685/inputs.json")
    pattern = read("../1688/inputs.json")
    m, q = heavy["m_N0H0"], heavy["coupling_to_x"]
    phase = np.pi*fields["scalar_angle_over_pi"]
    rows = np.exp(1j*phase*np.array(pattern["su5_generation_charges_5bar"]))
    A = np.array(pattern["source_a"])
    B = np.array(pattern["source_b"])
    return np.array([rows[:,None]*np.column_stack((B/m,0*A)),
                     rows[:,None]*np.column_stack((0*B,np.exp(2j*phase)*A/q))])


def residue_coefficients(parts):
    linear, const = parts
    return np.stack([np.sum(abs(const)**2,axis=0),
                     2*np.sum(const.conj()*linear,axis=0).real,
                     np.sum(abs(linear)**2,axis=0)])


def residues(coeff, weight):
    return np.array([1,weight,weight**2]) @ coeff


def finite_loop(x):
    # G_gamma(0)=0 implements GIM subtraction. All heavy arguments in this test exceed 4000.
    return -x*(2*x*x+5*x-1)/(4*(1-x)**3) - 1.5*x**3*np.log(x)/(1-x)**4


def corrected_gamma_windows():
    table=read("../1707/independent_results.json")
    points=table["frozen_direct_tree_points"]
    roots=table["all_central_physical_triangle_inversions"]
    root=min(roots,key=lambda x:x["Vub"])["gamma_deg"]
    sets={"time_integrated":([60.1,64.8],[57.7,67.1]),
          "charged_B_only":([59.7,64.7],[57.2,67.2])}
    out={}
    for label,(narrow,wide) in sets.items():
        passing=[p["label"] for p in points if narrow[0]<=p["gamma_deg"]<=narrow[1]
                 and p["S_inside_one_width"] and p["frequency_lattice_inside_one_width"]]
        assert passing==table["one_width_joint_pass"]
        assert not narrow[0]<=root<=narrow[1] and wide[0]<=root<=wide[1]
        out[label]=dict(window68=narrow,window95=wide,number_checked=len(points),joint_pass=passing)
    return dict(inverted_gamma_deg=root,components=out)


def evaluate(model, coeff, weak, x):
    mu=model["m"]
    u=model["U"]
    vev=model["vev"]
    C=model["C0"]
    unit=(weak[0]*x+weak[1])@C.T
    target=model["src"]["neutral_only_Higgs_mass_squared_GeV2"]
    gain=vev*math.sqrt(8*np.pi**2*target/np.dot(mu**2,residues(coeff,x)))
    model["graph"].__globals__["D0"]=unit*gain
    dirac, rep, track, _, _, _, common=model["graph"](1,-1,"unit_third",.20)
    expected=residues(coeff,x)*(gain*common["dimensionless_common_D_gain"])**2
    np.testing.assert_allclose(np.sum(abs(dirac)**2,axis=0),expected,rtol=1e-12)
    fl=model["fl"]
    physical=np.linalg.norm(fl["physics"]["Y0"][:,[0,2]],axis=0)**2
    high_ratio=expected[0]/expected[1]
    read_ratio=physical[0]/physical[1]
    polynomial=coeff[:,0]-high_ratio*coeff[:,1]
    roots=np.polynomial.polynomial.polyroots(polynomial)
    assert np.isreal(roots).all() and np.isclose(roots,x).any()
    other=roots[np.argmax(abs(roots-x))]
    endpoints=sorted(coeff[[0,2],0]/coeff[[0,2],1])
    if endpoints[0]<high_ratio<endpoints[1]:
        assert other<0
    # Use a second, single-scale broken-phase comparison for the radiative decay.
    # Adjust the down-quark input of the charged-frame evolution to the same branch.
    fr=fl.copy()
    def shifted(key):
        vals=fl["central"](key)
        if key=="down":
            vals[0]-=fl["low"]["yukawa"]["down"][0][1]
        return vals
    fr["central"]=shifted
    frame=archived("../1706/independent_check.py")["charged_at_threshold"](
        fr,rep,model["src"]["neutral_poles_GeV"][2])
    lepton_dirac=frame.T@dirac@u.T
    pair=model["lower_tensors"](model["pair_sources"](dirac).sum(axis=1),track)
    light_sv=np.linalg.svd(pair,compute_uv=False)
    solar=read("../1668/inputs.json")["juno_NO_dm21"][0]
    common_H=math.sqrt(solar)/light_sv[1]
    H=model["H0"]*common_H
    full=np.block([[np.zeros((3,3)),lepton_dirac,np.zeros((3,2))],
                   [lepton_dirac.T,np.zeros((2,2)),C],
                   [np.zeros((2,3)),C.T,H]])
    _,mass,vh=np.linalg.svd(full)
    eig=vh.conj().T[:3,:4]
    mass=mass[:4]
    mw,alpha=80.4,1/137.036
    kernel=finite_loop((mass/mw)**2)
    amp=(eig[0]*eig[1].conj())@kernel
    loss=np.sum(abs(eig)**2,axis=1)
    br=3*alpha/(2*np.pi)*abs(amp)**2/((1-loss[0])*(1-loss[1]))
    approx=(frame.T@dirac).conj()/mu
    approx_amp=(approx[0]*approx[1].conj())@finite_loop((mu/mw)**2)
    assert abs(approx_amp/amp-1)<.001
    return dict(absolute_source=x,high_scale_reduced_width_ratio=high_ratio,
                actual_two_scales_reduced_width_ratio=read_ratio,
                partial_width_ratio_low_over_high_massless_W=read_ratio*mu[0]/mu[1],
                other_algebraic_root=float(other),common_dirac_multiplier=gain*common["dimensionless_common_D_gain"],
                common_terminal_multiplier=common_H, gamma_quark_deg=rep["gamma_MZ_deg"],
                atmospheric_dm31_eV2=solar*(light_sv[0]/light_sv[1])**2,
                mu_e_gamma_W_only=br,mu_e_gamma_ratio_to_published_MEGII_limit=br/1.5e-13,
                relative_approximation_error_in_W_amplitude=abs(approx_amp/amp-1))


def compute():
    model=archived("../1706/check.py")
    cfg=model["fl"]["cfg"]["direct_semileptonic_moduli"]
    cfg["Vub"][0]-=np.linalg.norm(cfg["Vub"][1:])
    weak=source_basis()
    coeff=residue_coefficients(weak@model["C0"].T@model["U"].conj())
    us=cfg["Vus"][0]
    table={name:evaluate(model,coeff,weak,x) for name,x in {
        "one_Cabibbo_projection":us,"two_terminal_normalisation":math.sqrt(2)*us,
        "both_vectors_normalised":math.sqrt(2/11),"unit_source":1.}.items()}
    primary=read("results.json")
    matched={}
    for tag,row in [("central_K",primary["same_1708_point"]["two_port_Cabibbo"]),
                    ("higher_K",primary["additional_Vus_marginal_endpoints"]["1"])]:
        cfg["Vus"][0]=row["Vus_MZ"]
        point=evaluate(model,coeff,weak,math.sqrt(2)*row["Vus_at_upper_neutral_cluster"])
        np.testing.assert_allclose(point["atmospheric_dm31_eV2"],row["delta_m31_squared_eV2"],rtol=2e-7)
        matched[tag]=point
    # Keep all physical weak-basis light matrix elements and the Higgs moment fixed at the common scale.
    Z=weak[0]*math.sqrt(2)*us+weak[1]
    h=model["H0"]
    m=Z@h@Z.T
    old=np.dot(model["m"]**2,residues(coeff,math.sqrt(2)*us))
    alts={}
    for factor in (.5,2.):
        x=math.sqrt(2)*us*factor
        a=math.sqrt(old/np.dot(model["m"]**2,residues(coeff,x)))
        comp=a*(weak[0]*x+weak[1])
        hp=h@np.diag([1/(factor*a)**2,1/a**2])
        np.testing.assert_allclose(comp@hp@comp.T,m,rtol=3e-14,atol=1e-21)
        nr=residues(coeff,x)
        alts[str(factor)]=dict(H66_over_H00=hp[1,1].real/hp[0,0].real,
                               g_over_original=a,common_scale_ratio=nr[0]/nr[1])
    uniqueness=sorted(coeff[[0,2],0]/coeff[[0,2],1])
    return dict(positive_unique_high_scale_ratio_interval=uniqueness,cases=table,
                same_physical_matching_scale=matched,common_scale_compensated_deformations=alts,
                archived_quark_window_correction=corrected_gamma_windows())


def stable(val):
    if isinstance(val,dict):
        return {k:stable(v) for k,v in val.items()}
    if isinstance(val,list):
        return [stable(v) for v in val]
    if isinstance(val,(float,np.floating)):
        return float(f"{val:.10g}")
    return val


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    data=stable(compute())
    if ap.parse_args().check:
        assert data==read("independent_results.json")
        print("Verified source inversion, paired residues, W-mediated decay and gamma holdouts.")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as out:
            out.write(json.dumps(data,indent=2)+"\n")
