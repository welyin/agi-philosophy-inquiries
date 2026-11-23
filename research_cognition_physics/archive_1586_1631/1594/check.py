"""Rank-two common left image and independent charged-mass/dipole measurements."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def direction_bounds(ordering, s12, s13, s23):
    a, b = s12[1:]
    c, d = s13[1:]
    u, v = s23[1:]
    if ordering == "IO":
        return (math.sqrt((1-v)*(1-d)), math.sqrt((1-v)/v),
                math.sqrt((1-u)/u), math.sqrt(d/((1-d)*u)))
    tlo, thi = math.sqrt(u/(1-u)), math.sqrt(v/(1-v))
    rhi = math.sqrt((1-a)*d/a)
    assert tlo > rhi and thi*rhi < 1
    nt_min = math.sqrt(a*u) - math.sqrt((1-a)*(1-u)*d)
    ratio_lo = (tlo-rhi)/(1+tlo*rhi)
    ratio_hi = (thi+rhi)/(1-thi*rhi)
    e_over_mu_hi = math.sqrt((1-a)*(1-c)/a)*math.sqrt(1+thi**2)/(1-thi*rhi)
    return nt_min, ratio_lo, ratio_hi, e_over_mu_hi


def amplitude_cap(light, heavy, life, branching, dip, efficiency):
    phase = (1-(light/heavy)**2)**3
    return math.sqrt(4*dip["hbar_GeV_s"]*branching /
                     (life*dip["alpha"]*heavy**3*phase*efficiency))


def factorization_example():
    rng = np.random.default_rng(1594)
    raw = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    basis = np.linalg.qr(raw)[0]
    n, plane = basis[:, 0], basis[:, 1:]
    charged = np.diag([0.00051099895, 0.1056583755, 1.77686])
    neutral = plane @ np.diag([0.00865, 0.0501]) @ plane.T
    gains = np.array([1.e-7, 2.e-7])
    for t in [n, n + 0.4j*plane[:, 0], n - (2+0.3j)*plane[:, 1]]:
        third = np.outer(t, n.conj() @ charged)
        remainder = charged-third
        rows = plane.conj().T @ remainder
        dipole = plane @ np.diag(gains) @ rows
        assert np.linalg.norm(plane @ rows - remainder) < 5e-14
        assert np.linalg.matrix_rank(remainder) == 2
        assert np.linalg.matrix_rank(third) == 1
        assert np.linalg.norm(neutral.conj().T @ n) < 1e-14
        assert np.linalg.norm(n.conj() @ dipole) < 1e-20
    bypass = np.outer(n, [0, 1e-9, 0])
    assert abs((n.conj() @ bypass)[1]-1e-9) < 1e-20
    return True


def split_difference(d):
    # f(1+d)-f(1-d), f(x)=x log(x)/(x-1). This also resolves tiny neutral splits.
    return 2*sum(d**k/(k*(k+1)) for k in range(1, 18, 2))


def finite_matter_example(order, inputs, angles, atmospheric, dip):
    def trig(value):
        return math.sqrt(value), math.sqrt(1-value)
    s12, c12 = trig(angles[order]["s12"][0])
    s13, c13 = trig(angles[order]["s13"][0])
    s23, c23 = trig(atmospheric[order][0])
    phase = -1j  # exp(i delta) for a finite illustration
    u = np.array([[c12*c13, s12*c13, s13*phase.conjugate()],
                  [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
                  [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13]])
    solar = angles[order]["solar"][0]
    atm = angles[order]["atm_abs"][0]
    if order == "NO":
        n, cols, ev = u[:, 0], u[:, 1:], np.sqrt([solar, atm])
    else:
        n, cols, ev = u[:, 2], u[:, :2], np.sqrt([atm-solar, atm])
    point = inputs["finite_example"]
    mass, v = point["heavy_GeV"], point["vev_GeV"]
    a = point["cubic_over_heavy"]*mass
    d = a*v/(math.sqrt(2)*mass**2)
    g = mass*split_difference(d)/(32*math.pi**2)
    b = max(ev)*1e-9
    # Splittings are +/- lambda5*v^2/2; omitted terms are bounded by O(d_n^3).
    dn = b*32*math.pi**2/mass
    lambda5 = 2*mass**2*dn/v**2
    b_exact = mass*split_difference(dn)/(32*math.pi**2)
    left = cols @ np.diag(np.sqrt(ev*1e-9 / b_exact))
    charged = np.diag([dip["me_GeV"], dip["mmu_GeV"], point["tau_GeV"]])
    t = np.array([0., 0., 1/n[2].conjugate()])
    third = np.outer(t, n.conj() @ charged)
    rows = np.linalg.pinv(left) @ (charged-third)/g
    reconstructed = g*left @ rows + third
    assert np.linalg.norm(reconstructed-charged) < 1e-14
    assert np.linalg.norm(b_exact*left @ left.T - (cols*ev) @ cols.T*1e-9) < 1e-24
    # Direct positive Q=1 Feynman parameter integral, Gauss points include both masses.
    nodes, weights = np.polynomial.legendre.leggauss(32)
    coords, w = (nodes+1)/2, weights/2
    s, tfrac = coords[:, None], coords[None, :]
    delta = mass**2*(1+tfrac*(2*s-1)*d)
    weight = np.outer(w, w)
    gain = np.sum(weight*tfrac**2*(1-tfrac)/delta**2) / np.sum(weight*tfrac/delta)
    c = (reconstructed-third)*gain
    em = abs(c[0,1])**2+abs(c[1,0])**2
    tm = abs(c[1,2])**2+abs(c[2,1])**2
    te = abs(c[0,2])**2+abs(c[2,0])**2

    def br(light, heavy, lifetime, squared):
        return dip["alpha"]*heavy**3*(1-(light/heavy)**2)**3*lifetime*squared/(4*dip["hbar_GeV_s"])

    me, mmu, mtau = np.diag(charged)
    assert mass**2-v**2/2-a**2/4 > 0 and 0 < lambda5 < 1
    return {"common_neutral_mass_GeV":mass, "charged_loop_coefficient_GeV":g,
            "neutral_quartic_fitted":lambda5, "max_abs_right_yukawa":float(np.abs(rows).max()),
            "max_abs_left_yukawa":float(np.abs(left).max()),
            "tree_mass_source_op_norm_GeV":float(np.linalg.norm(third,2)),
            "internal_flip_gain_GeV_minus2":float(gain),
            "leading_amu":float(2*mmu*c[1,1].real), "leading_ae":float(2*me*c[0,0].real),
            "leading_Br_mu_e":float(br(me,mmu,dip["mu_lifetime_s"],em)),
            "leading_Br_tau_mu":float(br(mmu,mtau,inputs["tau_lifetime_lower_s"],tm)),
            "leading_Br_tau_e":float(br(me,mtau,inputs["tau_lifetime_lower_s"],te))}


def calculate():
    inp = read("inputs.json")
    dip = read(inp["dipole_inputs"])
    ref = read(inp["muon_reference"])
    me, mmu = dip["me_GeV"], dip["mmu_GeV"]
    mtau_lo, mtau_hi = inp["m_tau_GeV"]
    angles = read(inp["angles_12_13"])
    atmospheric = read(inp["angles_23"])["theta23_sin2"]
    historical = (ref["experiment_world"][0] - ref["WP20_theory"][0] -
                  2*math.hypot(ref["experiment_world"][1], ref["WP20_theory"][1]))*1e-11
    modern = [(ref["WP25_difference_as_published"][0]+s*2*ref["WP25_difference_as_published"][1])*1e-11
              for s in (-1, 1)]
    outcomes = {}
    for ordering in ("NO", "IO"):
        nt, lo, hi, emu = direction_bounds(ordering, angles[ordering]["s12"],
                                            angles[ordering]["s13"], atmospheric[ordering])
        cases = []
        for eff in inp["rates"]:
            em = amplitude_cap(me, mmu, dip["mu_lifetime_s"], dip["meg_upper"], dip, eff)
            tm = amplitude_cap(mmu, mtau_lo, inp["tau_lifetime_lower_s"],
                               inp["tau_mu_gamma_bound"], dip, eff)
            cap = 2*mmu*(emu*em+hi*tm)
            required_tm = max(0, historical/(2*mmu)-emu*em)/hi
            cases.append({"efficiency_floor":eff, "abs_additional_amu_cap":cap,
                          "WP20_min_Br_tau_mu_gamma_in_same_response":inp["tau_mu_gamma_bound"]*(required_tm/tm)**2,
                          "shortfall_in_amu_if_positive":max(0, historical-cap)})
        outcomes[ordering] = {"tau_component_min":nt, "tau_over_mu_range":[lo,hi],
                              "e_over_mu_cap":emu, "extra_mass_source_norm_min_GeV":mtau_lo*nt,
                              "extra_right_mu_over_tau_range":[mmu/(mtau_hi*hi),mmu/(mtau_lo*lo)],
                              "dipole_cases":cases}
    finite = {o: finite_matter_example(o, inp, angles, atmospheric, dip) for o in outcomes}
    return {"common_plane_examples":factorization_example(), "finite_matter_examples":finite, "WP20_lower":historical,
            "WP25_comparison":modern, "orders":outcomes}


def compare(saved, actual):
    if isinstance(actual, dict):
        assert saved.keys() == actual.keys()
        for key in actual:
            compare(saved[key], actual[key])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for left, right in zip(saved, actual):
            compare(left, right)
    else:
        assert math.isclose(saved, actual, rel_tol=4e-12, abs_tol=1e-27)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    out = calculate()
    if args.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(out, indent=2, allow_nan=False)+"\n")
    elif args.check:
        compare(read("results.json"), out)
        print("1594: common-plane mass source, oscillation bounds and dipole tests verified")
    else:
        print(json.dumps(out, indent=2))
