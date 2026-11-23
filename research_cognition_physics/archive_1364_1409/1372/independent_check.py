"""Dimensionless GW spectra and direct derivatives of the one-loop determinant."""
import argparse
import json
from math import factorial, pi, sqrt
from pathlib import Path

import numpy as np


V, VH = 248.401, .032
HOBS, HD = 125.11, .11
GAUGE = np.array([.461228,.65096])
GD = np.array([.000026,.00004])
YQ = np.array([7.04e-6,.00356,.967, 1.54e-5,.000306,.01630])
YQD = np.array([.15e-6,.00006,.004, .02e-5,.000004,.00009])
YL = np.array([2.77713e-6,.000585042,.0099378])
YLD = np.array([.00036e-6,.000000075,.0000014])


def charged_spectrum(g, q, e):
    masses = [g[1]/2, np.linalg.norm(g*np.sqrt([.6,1]))/2]
    masses.extend(q/sqrt(2)); masses.extend(e/sqrt(2))
    dof = [6,3]+[-12]*6+[-4]*3
    finite_parts = [5/6]*2+[3/2]*9
    return np.array(masses),np.array(dof),np.array(finite_parts)


def determinant_derivatives(m_over_f, dof, finite_parts):
    # All masses scale as r/f. Solve the stationary subtraction scale in the
    # original sum, then differentiate the unchanged sum by a Cauchy contour.
    fourth = dof*m_over_f**4
    log_squared_scale = -.5 - (fourth @ (np.log(m_over_f**2)-finite_parts))/fourth.sum()
    directions = np.exp(2j*pi*np.arange(96)/96)
    h = .13
    z = 1 + h*directions[:,None]
    m_z_squared = (z*m_over_f)**2
    veff = ((m_z_squared**2*(np.log(m_z_squared)+log_squared_scale-finite_parts)) @ dof)/(64*pi**2)
    derivatives = [float((factorial(k)*np.mean(veff*directions**(-k))/h**k).real)
                   for k in range(1,5)]
    m2 = fourth.sum()/(8*pi**2)
    assert np.allclose(derivatives,[0,m2,5*m2,11*m2],rtol=3e-9,atol=1e-11)
    return {"f2_curvature_ratio":m2, "fV3_div_V2":derivatives[2]/derivatives[1],
            "f2V4_div_V2":derivatives[3]/derivatives[1]}


def stored_sign(g, q, e):
    masses, signs, _ = charged_spectrum(g,q,e)
    return float(signs @ masses**4)


def singlet_trial(lam):
    ms, weights, finite_parts = charged_spectrum(GAUGE,YQ,YL)
    f = sqrt(2)*V
    # The two real Takagi blocks yield a negative algebraic light eigenvalue,
    # whose absolute value is its physical mass; compare using stable formulas.
    neutrino_light_eV = [sqrt(.0000749),sqrt(.002513)]
    heavy = V/2
    ratios = []
    for m_ev in neutrino_light_eV:
        m = m_ev*1e-9
        md = sqrt(m*heavy)
        eig = np.linalg.eigvalsh([[0.,md],[md,heavy-m]])
        assert np.allclose(eig,[-m,heavy],rtol=5e-13,atol=5e-21)
        ratios.extend([m/V,heavy/V])
    full_mass = np.concatenate([ms,[2*sqrt(lam)],ratios])/sqrt(2)
    full_weights = np.concatenate([weights,[1],[-2]*4])
    constants = np.concatenate([finite_parts,[1.5]*5])
    derivatives = determinant_derivatives(full_mass,full_weights,constants)
    curvature = derivatives["f2_curvature_ratio"]*f*f
    return {"lambda":lam,"transverse_GeV":2*sqrt(lam)*V,"signed_curvature_GeV2":curvature,
            "positive_mass_GeV":sqrt(curvature) if curvature>0 else None,"derivatives":derivatives}


def compute():
    masses, weights, finite_parts = charged_spectrum(GAUGE,YQ,YL)
    spectrum = stored_sign(GAUGE,YQ,YL)
    lo = stored_sign(GAUGE-GD,YQ+YQD,YL+YLD)
    hi = stored_sign(GAUGE+GD,YQ-YQD,YL-YLD)
    assert lo<spectrum<hi<0
    target = 8*pi**2*(HOBS/V)**2
    needed = target-spectrum
    required_box = [8*pi**2*((HOBS-HD)/(V+VH))**2-hi,
                    8*pi**2*((HOBS+HD)/(V-VH))**2-lo]
    # Illustrative portal potential V_int=(lambda_i/2)|H|^2 S_i^2, doublet-only vev.
    # m_i^2=lambda_i v^2/2 and the equation fixes sum lambda_i^2.
    return {"source":"https://arxiv.org/html/2510.01312v2", "basis":"masses in units of v(MZ)",
            "SM_signed_spectral_fourth":spectrum, "full_box_signed_spectral_fourth":[lo,hi],
            "curvature_GeV2":spectrum*V**2/(8*pi**2),
            "critical_top_yukawa_if_other_fermions_ignored":
                ((6*masses[0]**4+3*masses[1]**4)/3)**.25,
            "if_f_equals_v_125_target_excess_over_SM_v4":needed,
            "same_box_including_observed_Higgs_reference":required_box,
            "equal_single_boson_GeV_if_f_equals_v":V*needed**.25,
            "sum_portal_squares_if_portal_only_and_f_equals_v":4*needed,
            "SM_stationary_max_detector":determinant_derivatives(masses,weights,finite_parts),
            "two_scalar_trials":[singlet_trial(.25),singlet_trial(.5)]}


def rounded(obj):
    if isinstance(obj,dict):
        return {k:rounded(v) for k,v in obj.items()}
    if isinstance(obj,list):
        return [rounded(v) for v in obj]
    return float(format(obj,".9g")) if isinstance(obj,float) else obj


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true",help="verify frozen results (default)")
    parser.add_argument("--write",action="store_true")
    op=parser.parse_args()
    file=Path(__file__).with_name("independent_results.json")
    want=rounded(compute())
    if op.write:
        with file.open("x",encoding="utf-8") as out:
            json.dump(want,out,indent=2,ensure_ascii=False)
            out.write("\n")
    else:
        assert want==json.loads(file.read_text(encoding="utf-8"))
    print("1372 independent spectrum, neutral sources and determinant derivatives verified.")
