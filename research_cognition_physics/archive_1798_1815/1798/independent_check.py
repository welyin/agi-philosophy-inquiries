"""Test the pure triplet neutrino sum rule on spectra, arbitrary embeddings, and one symmetric witness."""

import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SOLAR, ATM = 7.49e-5, 2.513e-3
TAU_CEILINGS = {"normal": .746, "inverted": .690}
FIT_ANGLES = [.308, .02215, .470, 212.]
# One solution obtained by damped least-squares using only d,e spectra and oscillations.
# The first two entries are log(neutral normalization), log(m1/sqrt(Delta31)).
POINT = [-2.485055539922243, -2.2994258677385915, -.17711536828030575,
         -3.2598496934582273, 2.0087929047818878, 2.8769246386582576, -.9848412926051147]
OTHER_POINT = [-2.310077582517645, -3.0059147630422105, .2097549831636534,
               -2.6845344184721767, 1.9070409497178737, -.9380379900559167, -1.124164916053417]


def values(a):
    return np.linalg.svd(a, compute_uv=False)[::-1]


def vacuum_rotation(i, j, sine_squared, phase=0):
    r = np.eye(3, dtype=complex)
    v = math.sqrt(sine_squared)*np.exp(-1j*phase)
    r[i, i] = r[j, j] = math.sqrt(1-sine_squared)
    r[i, j], r[j, i] = v, -v.conjugate()
    return r


def positive_gap_cosine(gap, perturbation):
    ratio = 2*perturbation/gap
    if ratio >= 1:
        return None
    return math.cos(.5*math.asin(ratio))


def required_perturbation(gap, pole_upper):
    return gap*math.sin(2*min(math.pi/4, math.acos(pole_upper)))/2


def solve_witness(e, d_spectrum, transport, p):
    # E is diagonal and normalized by its own largest singular value.
    u = (vacuum_rotation(1, 2, FIT_ANGLES[2])
         @ vacuum_rotation(0, 2, FIT_ANGLES[1], math.radians(FIT_ANGLES[3]))
         @ vacuum_rotation(0, 1, FIT_ANGLES[0]))
    modes = np.sqrt(math.exp(2*p[1])+np.array([0., SOLAR/ATM, 1.]))
    majorana = np.exp(1j*np.array(p[2:5]))
    row_phases = np.diag(np.exp(1j*np.array([p[5], p[6], 0.])))
    x_low = math.exp(p[0])*row_phases @ (u.conj() * (modes*majorana)) @ u.conj().T @ row_phases
    ti = np.linalg.inv(transport)
    d = e - ti.T @ x_low @ ti
    h, f = (3*d+e)/4, (d-e)/4
    for tensor in (d, e, h, f):
        np.testing.assert_allclose(tensor.T, tensor, atol=1e-14)
    np.testing.assert_allclose(h+f, d, atol=1e-14)
    np.testing.assert_allclose(h-3*f, e, atol=1e-14)
    np.testing.assert_allclose(values(d), d_spectrum, rtol=3e-10)

    ux, s, _ = np.linalg.svd(transport.T @ (-4*f) @ transport)
    # SVD gives left rather than conjugated Takagi columns for our Majorana convention.
    pmns = ux[:, ::-1].conj()
    np.testing.assert_allclose(abs(pmns), abs(u), atol=5e-14)
    j = (pmns[0, 0]*pmns[1, 1]*pmns[0, 1].conjugate()*pmns[1, 0].conjugate()).imag
    ref_j = (u[0, 0]*u[1, 1]*u[0, 1].conjugate()*u[1, 0].conjugate()).imag
    assert abs(j-ref_j)<2e-14
    eV_gain = math.sqrt(ATM/(s[0]**2-s[-1]**2))
    mass = s[::-1]*eV_gain
    np.testing.assert_allclose(mass[1]**2-mass[0]**2, SOLAR, rtol=6e-13)
    return d, {"parameters": p, "NuFIT_normal_fit_angles_and_delta_deg": FIT_ANGLES,
               "neutrino_masses_eV": mass.tolist(), "pmns_moduli": abs(pmns).tolist(), "J_lepton":j,
               "m_beta_eV":np.linalg.norm(pmns[0]*mass), "light_m_beta_beta_eV":abs(x_low[0,0])*eV_gain,
               "relative_down_mass_defect":float(max(abs(values(d)/d_spectrum-1))),
               "symmetric_down_in_tau_units_real": d.real.tolist(), "symmetric_down_in_tau_units_imag": d.imag.tolist()}


def quark_holdout(d, e, moments):
    # Use the main round's exact two-coefficient conditional interval on this different embedding.
    main = runpy.run_path(str(HERE/"check.py"))
    l, sd, r = main["frames"](d)
    down = np.diag(sd)
    other = l.conj().T@e@r
    forms = []
    for i in range(3):
        rows = np.column_stack([down[i], other[i]])
        forms.append(rows.conj().T@rows)
    res = main["conditional"](forms, 1, 0, 2, moments[1])
    return {"strange_over_bottom_input":moments[1], "down_over_bottom_target":moments[0],
            "down_over_bottom_allowed_at_strange":res["output_range"]}


def run(step=.1):
    old = runpy.run_path(str(HERE.parent/"1715/check.py"))
    flow, rep = old["charged_background"](step)
    decode = old["unpack"]
    d = decode(rep["down_texture_10PeV"])
    e = decode(rep["charged_texture_10PeV"]).T
    sd, se = values(d), values(e)
    gap = se[-1]-se[-2]
    w, _ = old["congruence"](flow, rep, step)
    transport = decode(rep["charged_left_10PeV"]).conj().T @ w
    transport = transport @ np.diag(np.diag(transport).conj()/abs(np.diag(transport)))
    distortion = np.linalg.norm(transport-np.eye(3), 2)
    norm = np.linalg.norm(transport, 2)
    cond = np.linalg.cond(transport)
    extra = se[-1]*distortion*(2+distortion)

    branches = {}
    for label, z in {"standard":1., "A_first":-e[0, 1]/(3*d[0, 1]),
                     "A_third":-e[2, 2]/(3*d[2, 2])}.items():
        strength = abs(z)*sd[-1]
        low = strength*norm**2+extra
        t = values(e-z*d)
        # Allow any right unitary in the matching of the two ends before the congruence.
        upper_ratio = (cond**2*t[1]/t[2])**2
        branches[label] = {"abs_z":abs(z),"source_phase_deg":np.degrees(np.angle(z)),
                           "high_tau_overlap_floor_any_symmetric_embedding":positive_gap_cosine(gap,strength),
                           "low_tau_overlap_floor_any_symmetric_embedding":positive_gap_cosine(gap,low),
                           "minimum_extra_low_matrix_norm_in_tau_units_to_reach_tau_data": {
                               k:max(0.,(required_perturbation(gap,c)-low)/se[-1]) for k,c in TAU_CEILINGS.items()},
                           "fixed_port_high_solar_lower_pair_over_full_gap":(t[1]**2-t[0]**2)/(t[2]**2-t[0]**2),
                           "fixed_port_low_lower_pair_over_full_gap_upper":upper_ratio,
                           "fixed_port_low_upper_pair_over_full_gap_lower":1-upper_ratio}
    standard_e = np.diag(se/se[-1])
    new_d, witness = solve_witness(standard_e, sd/se[-1], transport, POINT)
    new_d2, witness2 = solve_witness(standard_e, sd/se[-1], transport, OTHER_POINT)
    ckm = abs(np.linalg.svd(d)[0][:, ::-1])
    moments = np.array(rep["up_texture_diag_10PeV"])**2 @ ckm**2
    witness["up_quark_holdout"] = quark_holdout(new_d, standard_e, moments/moments[-1])
    witness2["up_quark_holdout"] = quark_holdout(new_d2, standard_e, moments/moments[-1])
    check = witness["up_quark_holdout"]
    assert check["down_over_bottom_allowed_at_strange"][0] > check["down_over_bottom_target"]
    return {"matching_scale_GeV":1e7,"down_and_lepton_Yukawas":[sd.tolist(),se.tolist()],
            "transport_norm":norm,"transport_distance_after_rephasing":distortion,
            "transport_condition_number":cond,"high_abs_z_threshold_for_isolated_top":gap/(2*sd[-1]),
            "NuFIT_6p0_two_fits_lepton_tau_modulus_caps":TAU_CEILINGS,
            "branches":branches,"standard_z_symmetric_oscillation_calibrated_witness":witness,
            "another_standard_z_symmetric_oscillation_calibrated_witness":witness2}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):
        return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(format(x,".10g"))
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    ans=clean(run())
    raw=json.dumps(ans,indent=2,allow_nan=False)+"\n"
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==ans
    else:
        path.write_text(raw,encoding="utf-8")
    print("1798 independent",hashlib.sha256(raw.encode()).hexdigest()[:12])
