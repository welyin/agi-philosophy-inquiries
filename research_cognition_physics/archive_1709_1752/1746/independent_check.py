"""Existence certificate: a single positive common filter tested on the split physical spectrum."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
OLD=runpy.run_path(str(HERE/"../1744/check.py"))
# Columns correspond to c, the half-angle gain, and one in the tagged 1744 left basis.
Q_REAL=[[-.23145959739247457,-.340041488785216,.15423287730845284],
        [-.036422195366467396,-.14665390941498183,.21304445099939143],
        [-.6807387720132495,-.42498917735360925,-.34827818768257107]]
Q_IMAG=[[.517915675744904,-.6462487967521159,.34804641321373003],
        [-.44033485465051686,.23338794974514315,.8266882103489953],
        [.13986840720870009,.45841157653416775,-.07054702374375134]]


def measure_at_threshold(Y,M,F):
    # Diagonalize mass squares, independently of the daughter SVD.
    yc=F@Y
    _,left=np.linalg.eigh(yc@yc.conj().T)
    n=left.T@F.T@M@F@left
    m2,U=np.linalg.eigh(n.conj().T@n)
    p=abs(U)**2
    z=p[0,2]
    return dict(m3_over_m2=math.sqrt(m2[2]/m2[1]),
                angle_readings=[4*z*(1-z),p[0,1]/(1-z),p[1,2]/(1-z)])


def width(x,table):
    central,lo,hi=table
    return (x-central)/(lo if x<central else hi)


def examine(step):
    b=OLD["mother"](step)
    nu=b["nu"]
    old_lep=b["model"]["JOINT"]["LEPTONS"]
    Q=np.array(Q_REAL)+1j*np.array(Q_IMAG)
    np.testing.assert_allclose(Q.conj().T@Q,np.eye(3),atol=2e-14)
    # QR removes the insignificant residual of decimal representation without changing the axes.
    Q=np.linalg.qr(Q)[0]
    gains=np.array([b["cW"],b["fW"],1.])
    F=(Q*gains)@Q.conj().T
    np.testing.assert_allclose(np.linalg.eigvalsh(F),gains,rtol=2e-15)
    rho=nu["FROZEN"][0]
    pt=b["pt"]
    T=np.linalg.inv(b["S"]@pt["low_E"])
    M=T.T@(pt["C"]+rho*np.outer(pt["b"],pt["b"]))@T
    Y=b["yW"][2,:,:3]
    s0=np.linalg.svd(Y,compute_uv=False)
    ss=np.linalg.svd(F@Y,compute_uv=False)
    volume=math.prod(gains)**(1/3)
    np.testing.assert_allclose(ss/s0,volume*np.ones(3),rtol=4e-8)
    shape=measure_at_threshold(Y,M,F)
    basis=np.linalg.svd(Y)[0][:,::-1]
    in_charged=basis.conj().T@F@basis
    orient=np.imag(in_charged[0,1]*in_charged[1,2]*in_charged[2,0]) / np.prod(
        [gains[1]-gains[0],gains[2]-gains[0],gains[2]-gains[1]])
    assert abs(shape["m3_over_m2"]-5.78)<4e-6 and abs(shape["angle_readings"][1]-.3155)<5e-7
    def run(a):
        # One common physical finite charged-sector multiplier applied at the weak pole,
        # followed by the complete daughter flow with that same multiplier.
        state=dict(b)
        state["yW"]=b["yW"].copy()
        state["yW"][2]*=a
        return OLD["daughter"](state,F)
    raw=run(1)
    a=1.
    row=raw
    for _ in range(3):
        a*=old_lep[2,0]/row["MZ_leptons"][2]
        row=run(a)
    assert abs(row["MZ_leptons"][2]/old_lep[2,0]-1)<2e-13
    o=row["neutrinos"]
    assert all(o["angles_in_box"]) and o["SK_T2K_FC68_90"][1]
    clocks=[width(o["dm31_eV2"],nu["FREQ"]["juno_NO_dm31"]),
            width(o["dm32_eV2"],nu["FREQ"]["accelerator_NO_dm32"])]
    assert max(abs(np.array(clocks)))<2
    charged=(np.array(row["MZ_leptons"])-old_lep[:,0])/old_lep[:,1]
    assert max(abs(charged))<2
    assert row["absolute_four_ratio_radius"]<2
    limits=nu["ANG"]
    margins=np.min(np.stack([o["angle_readings"]-limits[:,0],limits[:,1]-o["angle_readings"]]),axis=0)
    return dict(weak_scale_GeV=b["Mweak"],F_real=F.real.tolist(),F_imag=F.imag.tolist(),eig=gains.tolist(),
                relative_source_gain_rho=rho,weak_charged_singular_gains=ss/s0,
                rephasing_invariant_imaginary_cycle_over_eigengaps=orient,
                weak_shape=shape,charged_common_magnitude_at_weak=a,
                MZ_leptons_before_charged_gain=raw["MZ_leptons"],MZ_leptons=row["MZ_leptons"],
                MZ_lepton_marginal_offsets=charged,MZ_tau_calibrated_e_mu_ratios=row["MZ_relative_e_mu_after_tau_calibration"],
                neutrino_scale_eV_to_solar=o["scale_eV"],solar_center_eV2=nu["SOL"],dm31=o["dm31_eV2"],
                dm32=o["dm32_eV2"],clock_marginal_offsets=clocks,angles=o["angle_readings"],
                angle_box_margins=margins,delta_deg=o["delta_deg"],SK_T2K_68_90=o["SK_T2K_FC68_90"],
                mbb_meV=o["mbb_meV"],raw_dm21_with_unchanged_parent_scale_eV2=o["dm21_with_mother_solar_gain"],
                actual_mass_zero_column=row["MZ_zero_column_probabilities"],
                source_relative_residuals=row["matching_residual_neutral_columns"],
                quark_profile_radius=row["absolute_four_ratio_radius"],
                CKM_gamma_deg=row["MZ_gamma"],golden_S=row["S_Bd"],epsilonK_pull=row["epsilon_K_pull"])


def tidy(o):
    if isinstance(o,dict):return {k:tidy(v) for k,v in o.items()}
    if isinstance(o,np.ndarray):return tidy(o.tolist())
    if isinstance(o,list):return [tidy(v) for v in o]
    if isinstance(o,float):return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    fine,coarse=examine(.125),examine(.25)
    assert np.max(abs(np.array(fine["angles"])-coarse["angles"]))<2e-7
    assert abs(fine["dm31"]-coarse["dm31"])<2e-9
    output=tidy(dict(fine=fine,step_halving=dict(max_angle_difference=np.max(abs(np.array(fine["angles"])-coarse["angles"])),
                                              dm31_difference=abs(fine["dm31"]-coarse["dm31"]))))
    path=HERE/"independent_results.json"
    if args.check:
        assert output==json.loads(path.read_bytes())
        print("1746 independent OK: one positive filter, common scales, both clocks and three angles.")
    else:
        path.write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")

