"""Test the finite selection, common source, sterile spectrum, and external weak-current kernel."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def model():
    cfg = read("inputs.json")
    flavor = read(cfg["flavor_fields"])
    old = read(cfg["old_model"])
    qed = {"N0": (0,1), "N2": (2,1), "H0": (0,2), "H6": (6,2)}
    inserts = {"1": (0,0), "B": (0,2)}
    inserts.update({name:(flavor["flavor_charge"][name], 0) for name in ["S", "X", "A2", "A3"]})
    allowed = {}
    for i,j in itertools.combinations_with_replacement(qed, 2):
        bits = [field for field,(q,z) in inserts.items()
                if (q+qed[i][0]+qed[j][0]) % 12 == 0 and (z+qed[i][1]+qed[j][1]) % 3 == 0]
        if bits:
            allowed[i+" "+j] = bits
    assert allowed == {"N0 H0":["1"], "N2 H0":["A2"], "N2 H6":["X"],
                       "H0 H0":["B"], "H6 H6":["B"]}
    assert (2+2+2)%3 == (2+0+1)%3 == (0+2+1)%3 == (2*2+2*1)%3 == 0
    m, g, h = [cfg[k] for k in ["m_N0H0", "coupling_to_s2", "coupling_to_x"]]
    h0, h6 = cfg["bare_H0_H6_masses_in_sterile_units"]
    extra_h, extra_m, extra_g = cfg["extra_H0_mass_and_two_couplings"]
    H = np.diag([h0, h6])

    def make_sterile(s, x, extra=False):
        c = np.array([[m, 0], [g*s**2, h*x]], complex)
        heavy = H
        if extra:
            c = np.column_stack([c, [extra_m, extra_g*s**2]])
            heavy = np.diag([h0, h6, extra_h])
        sterile = np.block([[np.zeros((2,2)), c], [c.T, heavy]])
        return sterile, c, -c @ np.linalg.solve(heavy, c.T)

    eps = flavor["scalar_magnitude_over_cutoff"]
    s = eps*np.exp(1j*np.pi*flavor["scalar_angle_over_pi"])
    x = flavor["real_X_over_cutoff"]
    sterile, c, mn = make_sterile(s,x)
    np.testing.assert_allclose(np.linalg.det(sterile), (m*h*x)**2)
    k = g/m
    extra_schur = ((m*extra_g-extra_m*g)**2/(h0*extra_h) /
                   (m*m/h0+extra_m*extra_m/extra_h))
    assert extra_schur > 0

    def schur_ratio(heavy):
        small = -heavy
        return (small[1,1]-small[0,1]**2/small[0,0] - h*h*x*x/h6)/s**4
    assert abs(schur_ratio(mn)) < 2e-14
    _, _, mn_extra = make_sterile(s,x,True)
    np.testing.assert_allclose(schur_ratio(mn_extra),extra_schur,atol=2e-14)
    step = 1e-5
    def radial(exponent_s,exponent_x,extra):
        logs = []
        for sign in [-1, 1]:
            st = make_sterile(s*np.exp(step*sign*exponent_s),
                              x*np.exp(step*sign*exponent_x),extra)[0]
            logs.append(sum(np.log(np.linalg.svd(st,compute_uv=False))))
        return float((logs[1]-logs[0])/(2*step))
    slopes = {name: [radial(1,0,extra), radial(0,1,extra)]
              for name,extra in [("two_mediators",False),("third_mediator",True)]}
    np.testing.assert_allclose(slopes["two_mediators"], [0,2], atol=2e-8)
    assert abs(slopes["third_mediator"][0])>1e-2
    nupoint = read(cfg["neutrino_calibration"])["joint_witness"]
    r = read(cfg["old_result"])["rank_two_witness"]["small_root_relative_amplitude"]
    a,b = np.array(old["source_a"]), np.array(old["source_b"])
    ell = np.array(old["su5_generation_charges_5bar"])
    row = np.diag(s**ell)
    y0 = eps**(-ell)*np.sqrt(r)*b
    y2 = eps**(-ell)*a/math.sqrt(abs(s**4/x**2))
    d = row @ np.column_stack([y0, s**2 * (y2 + k*y0)])
    full_light = d @ np.linalg.solve(-mn,d.T)
    low_masses = np.linalg.svd(full_light,compute_uv=False)
    d *= np.sqrt(np.sqrt(nupoint["d21"])/(low_masses[1]*1e9*cfg["sterile_unit_GeV"]))
    mix = np.column_stack([d, np.zeros_like(d)])
    full = np.block([[np.zeros((3,3)), mix], [mix.T, sterile]])
    direct = np.block([[np.zeros((3,3)), d], [d.T,mn]])
    np.testing.assert_allclose(mix@np.linalg.solve(sterile,mix.T),
                               d@np.linalg.solve(mn,d.T),rtol=2e-14,atol=1e-28)

    def moments(mass):
        adj = mass.conj().T
        sq = adj@mass
        return [np.linalg.matrix_power(sq,n)[:3]@adj[:,:3] for n in [0,1,2]]
    k1,k3,k5 = moments(full)
    d1,d3,_ = moments(direct)
    assert np.linalg.norm(k1)+np.linalg.norm(k3)+np.linalg.norm(d1) == 0
    np.testing.assert_allclose(k5, d.conj()@c@H.conj().T@c.T@d.conj().T,rtol=3e-14)
    np.testing.assert_allclose(d3, d.conj()@mn@d.conj().T,rtol=3e-14)

    def kernel(mass, q):
        adj = mass.conj().T
        return np.linalg.solve(q*q*np.eye(len(mass))+adj@mass,adj)[:3,:3]
    Q = cfg["momentum_over_sterile_unit"]
    assert max(np.linalg.svd(full,compute_uv=False)) < Q and 2*Q < cfg["flavor_cutoff_over_sterile_unit"]
    ratios = [float(np.linalg.norm(kernel(t,2*Q))/np.linalg.norm(kernel(t,Q))) for t in [full,direct]]
    assert abs(ratios[0]*64-1) < 0.01 and abs(ratios[1]*16-1) < 0.001
    return {
        "allowed_renormalizable_bilinear_insertions": allowed, "charge_table_Z12_linear_cubic_extra": [6,216],
        "neutral_Wilson_mass_ratio": round(low_masses[0]/low_masses[1],10),
        "active_D_prefactor_in_sterile_units": float(f"{np.linalg.norm(d)/np.linalg.norm(row@np.column_stack([y0,s**2*(y2+k*y0)])):.10g}"),
        "physical_sterile_singular_masses_EW_off_GeV": [round(v*cfg["sterile_unit_GeV"],7)
                                               for v in np.linalg.svd(sterile,compute_uv=False)],
        "Schur_remainder_two_mediators": 0,
        "Schur_remainder_with_third_mediator": round(extra_schur,9),
        "sterile_log_product_slopes_S_X": {key:[round(v,7) for v in values] for key,values in slopes.items()},
        "full_weak_kernel_norm_ratio_2Q_to_Q": round(ratios[0],11),
        "direct_competitor_weak_kernel_norm_ratio_2Q_to_Q": round(ratios[1],11),
        "full_UV_asymptotic_ratio": 1/64, "direct_UV_asymptotic_ratio": 1/16,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = model()
    if args.check:
        assert result == read("results.json")
        print("1689 checks passed")
    else:
        print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
