"""Two right tags with both output rays fixed. All Yukawas are dimensionless."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def read(i):
    return json.loads((HERE.parent / str(i) / "inputs.json").read_text(encoding="utf8"))


def ckm(x):
    u, v, w = np.sin(x[:3])
    a, b, c = np.cos(x[:3])
    p = np.exp(1j*x[3])
    return np.array([[a*b, u*b, v/p], [-u*c-a*w*v*p, a*c-u*w*v*p, w*b],
                     [u*w-a*c*v*p, -a*w-u*c*v*p, c*b]])


def ckm_bounds(x):
    (u, v, w), (U, V, W) = [np.sin(x[:3, 0]+sign*2*x[:3, 1]) for sign in [-1, 1]]
    a, b, c = np.cos(x[:3, 0]+2*x[:3, 1])
    return np.array([[[a*b, 1], [u*b, U], [v, V]],
                     [[u*c-W*V, U+W*V], [a*c-U*W*V, 1], [w*b, W]],
                     [[u*w-V, U*W+V], [a*w-U*V, W+U*V], [c*b, 1]]])


def row_test(lo, hi, vbox, row, angles, y, p):
    down, ell = y[:3], y[3:]
    downlo, el = lo[:3], lo[3:]
    downhi, eh = hi[:3], hi[3:]
    qmin, qmax = vbox.T**2
    g2 = np.dot(qmin, downlo**2)
    bshare = 1/(1+np.dot(qmin[:2], downlo[:2]**2)/(qmax[2]*downhi[2]**2))
    dv2 = downhi[1]**2+(downhi[2]**2-downhi[1]**2)*bshare
    a, b = eh[0]**2/g2, dv2/el[2]**2
    disc = math.sqrt(p*p-4*a*b)
    gram = abs(ckm(angles)[row])**2
    h = ell**2
    def weight(i):
        j, k = [l for l in range(3) if l != i]
        return np.dot(gram, (down**2-h[j])*(down**2-h[k]))/((h[i]-h[j])*(h[i]-h[k]))
    denomlo = (el[2]**2-eh[0]**2)*(el[2]**2-eh[1]**2)
    denomhi = (eh[2]**2-el[0]**2)*(eh[2]**2-el[1]**2)
    tau_lo = (qmin[2]*(downlo[2]**2-eh[0]**2)*(downlo[2]**2-eh[1]**2)
              - np.dot(qmax[:2], downhi[:2]**2)*eh[1]**2)/denomhi
    tau_hi = qmax[2]*(downhi[2]**2-el[0]**2)*(downhi[2]**2-el[1]**2)/denomlo
    return {"coefficient_A": a, "coefficient_B": b, "p_mu_upper_if_k1": a+b,
            "IO_open_forbidden_gain": [math.sqrt(2*a/(p+disc)), math.sqrt((p+disc)/(2*b))],
            "required_ray_rotation_deg_at_k1": math.degrees(math.asin(math.sqrt(p))-math.asin(math.sqrt(a+b))),
            "ten_projection_central_e_mu_tau": [weight(i) for i in range(3)],
            "tau_projection_interval_k1": [tau_lo,tau_hi]}


def calculate():
    table = np.array(read(1402)["scales"]["1e16"])
    angles = np.array(read(1397)["high_scale"]["angles_theta12_theta13_theta23_delta"])
    io = read(1370)["IO"]
    p = (1-io["s13"][2])*read(1580)["theta23_sin2"]["IO"][1]
    y, se = table.T
    lo, hi = y-2*se, y+2*se
    bounds = ckm_bounds(angles)
    v = ckm(angles[:, 0])
    rows = {k: row_test(lo, hi, z, i, angles[:, 0], y, p)
            for i, (k,z) in enumerate(zip(["u","c","t"], bounds))}
    zeta = math.sqrt(1-bounds[2,2,0]**2)
    gram_gap = min(1-(hi[4]/lo[2])**2, (lo[5]/hi[2])**2-1) - zeta
    allowed = np.array([[l/(hi[2]*math.sqrt(1+zeta)), h/(lo[2]*math.sqrt(1-zeta))]
                        for l,h in zip(lo[3:], hi[3:])])
    down = v*y[:3]
    tags = v/y[:3]
    tags /= np.linalg.norm(tags, axis=1, keepdims=True)
    amp = np.array([.07915261452457506-.03707002843828199j,
                    -.8440208595069326+.03593099897859599j])
    mass = down.copy()
    mass[[0,2]] += y[5]*amp[:,None]*tags[[0,2]]
    left, singular, vh = np.linalg.svd(mass)
    np.testing.assert_allclose(singular[::-1], y[3:], rtol=2e-10)
    np.testing.assert_allclose((mass@mass.conj().T)[:,1], (down@down.conj().T)[:,1], atol=1e-18)
    np.testing.assert_allclose(abs(left[1,::-1])**2, rows["c"]["ten_projection_central_e_mu_tau"], atol=1e-12)
    assert rows["u"]["tau_projection_interval_k1"][1] < rows["c"]["tau_projection_interval_k1"][0]
    assert all(row["p_mu_upper_if_k1"] < p for row in rows.values())
    assert all(not l <= 1 <= h for l,h in allowed)
    return {"inputs": ["1402","1397","1370","1580"], "scale_GeV": 1e16, "IO_pmu_lower": p,
            "k1_by_unmodified_output": rows,
            "unused_top_charged_test": {"zeta": zeta,"minimum_gram_change_in_yb2": gram_gap,
                "additive_full_matrix_norm_min": math.sqrt(hi[5]**2+lo[2]**2*gram_gap)-hi[5],
                "outer_union_of_allowed_k": allowed.tolist()},
            "charged_only_example_unused_c": {"k": 1,"coefficients_over_ytau_re_im": [[c.real,c.imag] for c in amp],
                "singular_values_ascending": singular[::-1].tolist(),
                "ten_projection": (abs(left[1,::-1])**2).tolist(),
                "right_tag_u_lepton_e_mu_tau": (abs(vh[::-1]@tags[0].conj())**2).tolist(),
                "right_tag_t_lepton_e_mu_tau": (abs(vh[::-1]@tags[2].conj())**2).tolist()}}


def compare(a,b):
    if isinstance(b, dict):
        assert a.keys() == b.keys()
        for k in b:
            compare(a[k],b[k])
    elif isinstance(b, str) or (isinstance(b, list) and b and isinstance(b[0], str)):
        assert a == b
    else:
        np.testing.assert_allclose(a,b,rtol=1e-9,atol=1e-14)


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    opt=ap.parse_args()
    result=calculate()
    if opt.write:
        with OUT.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(result,json.loads(OUT.read_text(encoding="utf8")))
        print("1605 independent: passed")
