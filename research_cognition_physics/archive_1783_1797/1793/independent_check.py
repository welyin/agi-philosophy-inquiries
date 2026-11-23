"""Joint LFV consequences of a single common non-aligned lepton insertion."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
MH = 125.
POLE = np.array([0.000510999, 0.1056584, 1.77686])
TAU_SM_BR = 0.0627
TAU_STRENGTH = (0.82, 0.11)
LFV_LIMIT = {"emu": 0.000044, "etau": 0.0022, "mutau": 0.0015}
MODES = {"emu": (0, 1), "etau": (0, 2), "mutau": (1, 2)}


def unpack(x):
    return np.array(x["real"]) + 1j * np.array(x["imag"])


def matrices(step):
    old = runpy.run_path(str(HERE.parent / "1715/check.py"))
    flow, rep = old["charged_background"](step)
    d = unpack(rep["down_texture_10PeV"])
    e = unpack(rep["charged_texture_10PeV"])
    u, s, vh = np.linalg.svd(e)
    u, s, v = u[:,::-1], s[::-1], vh.conj().T[:,::-1]
    family = u.conj().T @ d.T @ v
    pole_gain = np.sqrt(POLE/s)
    diag = pole_gain[:,None] * family * pole_gain[None,:]
    high = flow["shoot"](rep["gamma_MZ_deg"], flow["angles"]("direct"), unpack(rep["charged_left_MZ"]), step)[0]
    high_u, high_d, high_e = [field[:,:3] for field in high[:3]]
    reright = np.linalg.solve(e, high_e)
    assert np.linalg.norm(reright.conj().T @ reright - np.eye(3), 2) < 2e-6
    em = np.array(flow["central"]("charged_leptons"))
    def transport(m):
        physical = np.diag(em) @ np.linalg.solve(high_e, m @ reright)
        gain = np.sqrt(POLE / em)
        return gain[:,None] * physical * gain[None,:]
    np.testing.assert_allclose(transport(e), np.diag(POLE), atol=4e-13)
    dl, ds, dr = np.linalg.svd(d)
    btop = (dl.conj().T @ e.T @ dr.conj().T)[0,0]/ds[0]
    up = np.linalg.svd(high_u)[0][:,::-1]
    dw = np.linalg.svd(high_d)[0][:,::-1]
    align = dl[:,-1]/(up.conj().T @ dw)[:,0]
    align = (align / abs(align))[:,None] * up.conj().T
    right_d = np.linalg.solve(d,align@high_d)
    low_b = np.linalg.solve(high_d,align.conj().T @ e.T @ right_d)[2,2]
    return {"singular_frame_gains": (diag,btop), "fixed_background_left": (transport(d.T),low_b)}, {
        "first_edge": e.T[0,1]/d[0,1], "third_edge": e.T[2,2]/d[2,2]}


def width_coeff(t, i, j):
    a, b = (POLE[i]/MH)**2, (POLE[j]/MH)**2
    k = np.sqrt((1-a-b)**2 - 4*a*b)
    s = (abs(t[i,j])**2+abs(t[j,i])**2)*(1-a-b)
    h = 4*math.sqrt(a*b)*abs(t[i,j]*t[j,i])
    tau_beta = math.sqrt(1 - 4*(POLE[2]/MH)**2)
    # Divide by |k_mu-k_tau|^2/|k_tau|^2 and enclose all common CP phases.
    dt = t[1,1]/POLE[1] - t[2,2]/POLE[2]
    unit = POLE[2]**2*abs(dt)**2
    return [k*(s-h)/(unit*tau_beta), k*(s+h)/(unit*tau_beta**3)]


def entire_pure_B_box(tau,mu,b):
    # Cover the relaxed measured annulus by a finite union of closed disks.
    # Every complex entry may also vary in an independent radius-1e-6 disk.
    edge, steps, err = 4., 320, .000001
    grid = np.linspace(-edge,edge,steps+1)
    mid = (grid[:-1]+grid[1:])/2
    z = mid[:,None]+1j*mid[None,:]
    radius = math.sqrt(2)*edge/steps+1e-12
    sel = (abs(z)+radius >= 1.65) & (abs(z)-radius <= 3.78)
    z = z[sel]
    dif = tau-mu
    eta = 1-b*tau
    err_eta = err*(abs(b)+abs(tau))+err**2
    err_bdif = err*(abs(dif)+2*abs(b))+2*err**2
    wz = abs(z)+radius
    first = abs(eta*z-b*dif)+abs(eta)*radius+err_eta*wz+err_bdif
    second = abs(dif+tau*z)+abs(tau)*radius+2*err+err*wz
    denom = 3*(abs(dif)-2*err)*np.maximum(1.65,abs(z)-radius)
    amp = first*second/denom+1e-10
    al = np.maximum(1e-10, abs(1+z)-radius-1e-12)
    common_max = np.minimum(1.81,3.11/al**2)+1e-10
    b_hi = np.max(common_max*amp**2/(1-4*(5/MH)**2))+1e-10
    assert b_hi < .22
    return dict(cover_annulus=[1.65,3.78], covered_boxes=len(z), boxes_total=steps**2,
                independent_complex_input_radius=err, normalized_bottom_over_ZZ_uniform_upper=float(b_hi),
                CMS_normalized_bottom_over_ZZ_lower=.30)


def analyze(t, phases, bottom):
    ranges = {k:width_coeff(t,*ij) for k,ij in MODES.items()}
    strength_floor = TAU_STRENGTH[0] - 2*TAU_STRENGTH[1]
    caps = {k: math.sqrt(limit/(ranges[k][0]*TAU_SM_BR*strength_floor))
            for k,limit in LFV_LIMIT.items()}
    reference = ranges["mutau"]
    pairs = {}
    for k,(i,j) in MODES.items():
        b = t[i,j]
        a = t[j,i]
        xi, xj = (POLE[i]/MH)**2, (POLE[j]/MH)**2
        kin = np.sqrt((1-xi-xj)**2-4*xi*xj)
        common = (1-xi-xj)*(abs(a)**2+abs(b)**2)
        mass_term = 4*math.sqrt(xi*xj)*abs(a*b)
        numerator = kin*(abs(b)**2-abs(a)**2)
        pol = sorted(numerator/(common+sgn*mass_term) for sgn in [-1,1])
        pairs[k] = dict(ratio_to_mutau_outer=[ranges[k][0]/reference[1],
                                              ranges[k][1]/reference[0]] if k != "mutau" else [1.,1.],
                        RH_column_heavier_chiral_fraction=abs(b)**2/(abs(b)**2+abs(a)**2),
                        polarization_heavier_negative_helicity_R_minus_L=pol,
                        forward_fraction_unit_pion_analyzer=[0.5+p/4 for p in pol] if j == 2 else None,
                        opposite_sign_effective_branch_at_lower_tau=[4*TAU_SM_BR*strength_floor*x for x in ranges[k]])
    tau_hi = TAU_STRENGTH[0]+2*TAU_STRENGTH[1]
    cdiff = [ranges['mutau'][0]-ranges['etau'][1],ranges['mutau'][1]-ranges['etau'][0]]
    lower = math.sqrt((.0025-2*.0010)/(tau_hi*TAU_SM_BR*cdiff[1]))
    upper = min(caps.values())
    equal_L = 2*math.degrees(math.asin(lower/2))
    atlas = dict(deltaB_measured=[.0025,.001], first_order_difference_coeff_inner_to_outer=cdiff,
                 separate_2width_and_CMS_LFV_cap_abs_kmu_over_ktau_minus_one=[lower,upper],
                 min_phase_deg_if_equal_magnitudes=equal_L,
                 equal_magnitude_opposite_sign_deltaB_work_range=[4*strength_floor*TAU_SM_BR*cdiff[0],4*tau_hi*TAU_SM_BR*cdiff[1]],
                 pure_B_z_branches_that_fail_only_this_work_screen=[])
    exact = {}
    for label,z in {"canonical":1., **phases,
                    **{"A_vanishes_"+k:-z/3 for k,z in phases.items()}}.items():
        q = .75*(1-z*np.diag(t)/POLE)
        k_ratio = q[1]/q[2]
        predicted = {k:[x*abs(k_ratio-1)**2 for x in bounds] for k,bounds in ranges.items()}
        exact[label] = {"B_only_abs_kmu_over_ktau_minus_one":abs(k_ratio-1),
                        "B_only_width_ratios_to_tautau":predicted}
        if not lower <= abs(k_ratio-1) <= upper:
            atlas['pure_B_z_branches_that_fail_only_this_work_screen'].append(label)
    assert 1.65<lower<upper<3.78
    allz = entire_pure_B_box(t[2,2]/POLE[2],t[1,1]/POLE[1],bottom)
    return dict(LFV_C_squared_difference_enclosure=ranges, external_caps_abs_kmu_over_ktau_minus_one=caps,
                pure_B_all_complex_z_after_conditional_two_width_screen=allz,
                comparison_with_ATLAS_symmetry_method=atlas,
                same_magnitude_opposite_sign_can_pass_necessary_screen={k:2<=v for k,v in caps.items()},
                modes=pairs,pure_B_contrasts=exact)


def round_data(x):
    if isinstance(x, dict):
        return {a:round_data(b) for a,b in x.items()}
    if isinstance(x,list):
        return [round_data(y) for y in x]
    return float(f"{x:.8g}") if isinstance(x, (float,np.floating)) else x


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    baseline, phases = matrices(0.1)
    data = {k:analyze(t,phases,b) for k,(t,b) in baseline.items()}
    output = json.dumps(round_data(data),ensure_ascii=False,indent=2)+"\n"
    saved = HERE / "independent_results.json"
    if args.check:
        assert json.loads(saved.read_bytes()) == json.loads(output)
        print("1793 independent OK",hashlib.sha256(output.encode()).hexdigest()[:12])
    else:
        saved.write_text(output,encoding="utf-8")


if __name__ == "__main__":
    main()
