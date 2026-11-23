"""Common left isometry, finite-rank channel defect, and running three-family data."""
import argparse
import json
from math import cos, sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(filename):
    return json.loads((HERE / filename).read_text(encoding="utf-8"))


def mix(s, phase):
    s12, s13, s23 = s
    c12, c13, c23 = [sqrt(1-x*x) for x in s]
    z = np.exp(1j*phase)
    return np.array([[c12*c13, s12*c13, s13/z],
                     [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
                     [s12*s23-c12*c23*s13*z, -c12*s23-s12*c23*s13*z, c23*c13]])


def data():
    cfg = read("inputs.json")
    low = read(cfg["low_scale_yukawa_and_ckm"])
    high = cfg["high_scale"]
    old = read(cfg["shared_old_light_quark_calibrations"][1])
    for group, pair in [("MZ",low["yukawa"]["up"]), ("GUT_SM",high["up"])]:
        assert old[group]["up"] == pair[0] and old[group]["charm"] == pair[1]
    return (("MZ", low["yukawa"]["up"], low["yukawa"]["down"], low["quark_sines_and_phase"],
             read(cfg["low_scale_g2"])["g2"], False),
            ("1e16_GeV", high["up"], high["down"], high["angles_theta12_theta13_theta23_delta"],
             high["g2"], True))


def shape_bounds(up, down, angles, g, radian, n):
    """Analytic outer intervals. All angles and the CP phase stay in (0, pi/2)."""
    def ends(rows):
        return np.array([v-n*e for v,e in rows]), np.array([v+n*e for v,e in rows])
    ul, uh = ends(up)
    dl, dh = ends(down)
    al, ah = ends(angles)
    gl, gh = g[0]-n*g[1], g[0]+n*g[1]
    if radian:
        al[:3], ah[:3] = np.sin(al[:3]), np.sin(ah[:3])
    assert np.all(ul > 0) and np.all(ul[1:] > uh[:-1])
    assert np.all(dl > 0) and np.all(dl[1:] > dh[:-1])
    assert 0 < al[3] < np.pi/2 and ah[3] < np.pi/2
    cl = np.sqrt(1-ah[:3]**2)
    j_floor = al[0]*al[1]*al[2]*cl[0]*cl[1]**2*cl[2]*sin(al[3])
    c_floor = cl[1]*cl[2]
    leak = 1-c_floor*c_floor
    r_c_sq = (uh[1]/ul[2])**2
    r_s_sq = (dh[1]/dl[2])**2
    lower_light = 1-r_c_sq-r_s_sq-(1-r_s_sq)*leak
    u_trace = 1+r_c_sq+(uh[0]/ul[2])**2
    d_trace = 1+r_s_sq+(dh[0]/dl[2])**2
    trace_floor = 3-2*max(u_trace,d_trace)/(1+c_floor)
    # Ordered opposite-slope affine spectra have exactly equal ratios.
    up_ratio_hi = (uh[1]**2-ul[0]**2)/(ul[2]**2-uh[0]**2)
    down_reversed_lo = 1-(dh[1]**2-dl[0]**2)/(dl[2]**2-dh[0]**2)
    fixed_min_eigen = [1-(uh[2]**2+dh[2]**2)/gl**2, 1-ul[2]**2/gh**2]
    return {"squared_spectrum_up_gap_ratio_upper": float(up_ratio_hi),
            "down_reversed_gap_ratio_lower": float(down_reversed_lo),
            "J_lower": float(j_floor), "fixed_g2_D_least_eigenvalue_bracket": fixed_min_eigen,
            "free_gains_any_light_plane_D_lower": float(lower_light),
            "free_gains_trace_D_lower": float(trace_floor)}


def center(up, down, angles, g, radian):
    yu, yd, a = [np.array(rows)[:, 0] for rows in (up, down, angles)]
    s = np.sin(a[:3]) if radian else a[:3]
    v = mix(s, a[3])
    np.testing.assert_allclose(v.conj().T @ v, np.eye(3), atol=2e-15)
    hu, hd = np.diag(yu**2), (v*yd**2) @ v.conj().T
    j = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    triplet = (hd[0,1]*hd[1,2]*hd[2,0]).imag
    assert np.isclose(triplet, j*np.prod([yd[1]**2-yd[0]**2, yd[2]**2-yd[1]**2,
                                          yd[2]**2-yd[0]**2]), rtol=1e-8, atol=0.)
    spectrum = np.linalg.eigvalsh(hu+hd)
    gcrit = sqrt(spectrum[-1])
    frame = np.vstack([np.diag(yu)/gcrit, np.diag(yd) @ v.conj().T/gcrit])
    val, vec = np.linalg.eigh(np.eye(3)-frame.conj().T @ frame)
    assert abs(val[0])<1e-12 and val[1]>.9
    rest = np.diag(np.sqrt(val[1:])) @ vec[:, 1:].conj().T
    extended = np.vstack([frame,rest])
    np.testing.assert_allclose(extended.conj().T @ extended, np.eye(3), atol=1e-14)
    return {"up_ordered_spectrum_gap_ratio": (yu[1]**2-yu[0]**2)/(yu[2]**2-yu[0]**2),
            "down_reverse_gap_ratio": (yd[2]**2-yd[1]**2)/(yd[2]**2-yd[0]**2),
            "CKM_quartet": float(j), "Hd_odd_product": float(triplet),
            "equal_g_critical": gcrit, "measured_g2": g[0],
            "fixed_g2_D_eigenvalues": np.linalg.eigvalsh(np.eye(3)-(hu+hd)/g[0]**2).tolist(),
            "fitted_common_g_D_eigenvalues": val.tolist(),
            "rank_two_extra_rows": rest.shape[0],
            "bounds_at_center": shape_bounds(up,down,angles,g,radian,0)}


def calculate():
    answer = {}
    for name,u,d,a,g,radians in data():
        bounds = {str(n): shape_bounds(u,d,a,g,radians,n) for n in read("inputs.json")["widths"]}
        for b in bounds.values():
            assert b["squared_spectrum_up_gap_ratio_upper"] < b["down_reversed_gap_ratio_lower"]
            assert b["J_lower"] > 0 and b["free_gains_any_light_plane_D_lower"] > .996
            assert b["free_gains_trace_D_lower"] > 1.998
        lower, upper = bounds["2"]["fixed_g2_D_least_eigenvalue_bracket"]
        assert (upper < 0 if name=="MZ" else lower > 0)
        answer[name] = {"center": center(u,d,a,g,radians), "outer_bounds": bounds}
    return answer


def verify(a, b):
    if isinstance(a, dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k], b[k])
    elif isinstance(a, list):
        assert len(a)==len(b)
        for i,j in zip(a,b):
            verify(i,j)
    else:
        np.testing.assert_allclose(a,b,rtol=5e-9,atol=1e-30 if 0<abs(a)<1e-20 else 2e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write:
        with (HERE/"results.json").open("w", encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    elif args.check:
        verify(result,read("results.json"))
        print("1397 three-family common doublet and finite resource bounds verified")
    else:
        print(json.dumps(result,indent=2))
