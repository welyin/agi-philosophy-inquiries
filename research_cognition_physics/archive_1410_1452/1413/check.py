"""Two calibrated Householder phases between cubic clocks in a shared left flavour space."""
import argparse
import cmath
import json
import math
from itertools import permutations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EYE = np.eye(3, dtype=complex)
ROOTS = np.exp(2j*np.pi*np.arange(3)/3)


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def ckm(sines, delta):
    x, y, z = sines  # 12, 13, 23
    a,b,c = np.sqrt(1-np.array(sines)**2)
    phase = cmath.exp(1j*delta)
    return np.array([[a*b, x*b, y*phase.conjugate()],
                     [-x*c-a*z*y*phase, a*c-x*z*y*phase, z*b],
                     [x*z-a*c*y*phase, -a*z-x*c*y*phase, c*b]])


def gate(v, eigenphase):
    return EYE + (eigenphase-1)*np.outer(v, v.conj())


def filter_gate(t, j, cycle, global_factor):
    x = math.sqrt((4-t*t)/3)
    return (x+t)/2*EYE - (x-t)/2 * np.conj(global_factor*ROOTS[j]) * cycle


def yield_masses(sigmas, cycle, global_factor):
    a,b,c = sigmas
    y = c * filter_gate(b/c, 1, cycle, global_factor) @ filter_gate(a/c, 0, cycle, global_factor)
    found, basis = np.linalg.eigh(y @ y.conj().T)
    assert np.allclose(np.sqrt(np.maximum(0, found)), sigmas, rtol=2e-6, atol=1e-11)
    return basis


def run_scale(mix, y_up, y_down):
    assert np.allclose(mix.conj().T @ mix, EYE)
    t = np.diag(ROOTS)
    prob = abs(mix)**2
    branches = []
    for perm in permutations(range(3)):
        # The j-th target column carries root perm[j].
        down = (mix * ROOTS[list(perm)]) @ mix.conj().T
        w = t.conj().T @ down
        singular = np.linalg.svd(EYE-w, compute_uv=False)
        ordered = prob[:, np.argsort(perm)]
        defect = ordered[0,1] - ordered[1,0]
        assert np.allclose(np.linalg.det(EYE-w), -3j*math.sqrt(3)*defect, atol=2e-14)
        branches.append({"down_root_order": list(perm), "abs_mod_square_asymmetry": abs(defect),
                         "exact_operator_distance": float(singular[-1])})
    down = (mix * ROOTS) @ mix.conj().T
    w = t.conj().T @ down
    vals, vec = np.linalg.eig(w)
    assert np.allclose(vec.conj().T @ vec, EYE)
    selected = int(np.argmin(np.abs(vals-1)))
    factor = np.conj(vals[selected])/abs(vals[selected])
    phases = factor*vals
    free = EYE.copy()
    for k in range(3):
        if k != selected:
            free = free @ gate(vec[:,k], phases[k])
    rotation = t @ free
    assert np.allclose(rotation, factor*down)
    assert np.allclose(np.linalg.matrix_power(rotation, 3), factor**3 * EYE)
    up_basis = yield_masses(y_up, t, 1.)
    down_basis = yield_masses(y_down, rotation, factor)
    charged = up_basis.conj().T @ down_basis
    j = lambda v: float(np.imag(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()))
    assert np.allclose(abs(charged)**2, prob, atol=2e-9) and np.isclose(j(charged), j(mix), atol=2e-9)
    # Actual exact operator using three raw phases, each once.
    exact = EYE.copy()
    for k in range(3):
        exact = exact @ gate(vec[:,k], vals[k])
    assert np.allclose(t @ exact, down)
    return {"root_assignments": branches, "CKM_abs_J": abs(j(mix)),
            "best_distance_all_root_assignments": min(r["exact_operator_distance"] for r in branches),
            "mass_order_global_phase_rad": float(np.angle(factor)),
            "controlled_three_pass_distance_to_I": abs(factor**3-1),
            "remaining_two_gate_phases_rad": [float(np.angle(phases[k])) for k in range(3) if k != selected],
            "effective_dressed_weak_moduli_squared": (abs(charged)**2).tolist()}


def phase_free_windows(limits):
    """Simple interval certificate for every CP phase; only elementary amplitude bounds."""
    (x0,x1),(y0,y1),(z0,z1) = limits
    c12 = math.sqrt(1-x1*x1)
    c13 = math.sqrt(1-y1*y1)
    c23 = math.sqrt(1-z1*z1)
    tdlo, tdhi = max(0., x0*z0-y1), x1*z1+y1
    udlo, cdlo = c12*c13, x0*c23-z1*y1
    cslo = c12*c23-x1*z1*y1
    uslo, cblo = x0*c13, z0*c13
    # The six differences B0i-B1j with i != j; two use column/row unitarity
    # to replace near-cancellations with separated small moduli.
    gaps = [tdlo**2-y1**2, cblo**2-tdhi**2, udlo**2-z1**2,
            cdlo**2-y1**2, uslo**2-z1**2, cslo**2-y1**2]
    assert min(gaps) > 0
    return {"six_squared_modulus_gaps_lower": gaps,
            "operator_norm_distance_lower_all_CP_all_labels": 3*math.sqrt(3)*min(gaps)/4}


def tree_test(data, widths):
    vprod, ep = data["K_product"]
    form, ef = data["K_form"]
    vmin, vmax = (vprod-widths*ep)/(form+widths*ef), (vprod+widths*ep)/(form-widths*ef)
    ratio, *errors = data["ub_over_cb"]
    rhi = ratio + widths*sum(errors)
    fcb, *ecb = data["F_times_cb"]
    cbmax = (fcb+widths*sum(ecb))/data["F_lower"]
    ubmax = rhi*cbmax
    c13 = math.sqrt(1-ubmax**2)
    # |Vtd| <= |Vcb|*(vmax/c13^2 + rhi), |Vtd| >= |Vcb|*(vmin-rhi).
    td_to_cb_upper = vmax/c13**2+rhi
    cd_lower = vmin*math.sqrt(1-(cbmax/c13)**2) - (cbmax/c13)*ubmax
    ud_lower = math.sqrt(1-(vmax/c13)**2)*c13
    cs_lower = math.sqrt(1-(vmax/c13)**2)*math.sqrt(1-(cbmax/c13)**2) - vmax*cbmax*ubmax/c13**2
    assert vmin > 2*rhi and td_to_cb_upper < 1
    assert min(cd_lower,cs_lower) > ubmax and min(vmin,ud_lower) > cbmax
    return {"Vus_range": [vmin,vmax], "Vub_over_Vcb_upper": rhi, "Vcb_upper": cbmax,
            "td_over_cb_minus_ub_over_cb_lower": vmin-2*rhi,
            "td_over_cb_upper": td_to_cb_upper, "cd_lower": cd_lower, "cs_lower": cs_lower,
            "ud_lower": ud_lower, "all_six_separated": True}


def run():
    config = read("inputs.json")
    low, high = read(config["weak_scale"]), read(config["high_scale"])["high_scale"]
    data = [("MZ", low["quark_sines_and_phase"], low["yukawa"]),
            ("1e16", high["angles_theta12_theta13_theta23_delta"], high)]
    out = {}
    for name, entries, mass in data:
        angles = np.array(entries[:3])
        span = config["working_widths"]*angles[:,1]
        lo, hi = angles[:,0]-span, angles[:,0]+span
        values = angles[:,0]
        if name == "1e16":
            lo,hi,values = np.sin(lo),np.sin(hi),np.sin(values)
        mix = ckm(values, entries[3][0])
        out[name] = run_scale(mix, [p[0] for p in mass["up"]], [p[0] for p in mass["down"]])
        out[name]["two_width_all_delta_certificate"] = phase_free_windows(zip(lo, hi))
    out["historical_tree_two_width"] = tree_test(read(config["tree_ports"]), config["working_widths"])
    return out


def agrees(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(agrees(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    if isinstance(a,bool):
        return a is b
    return math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-11)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    result=run()
    if parser.parse_args().write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as file:
            json.dump(result,file,indent=2,allow_nan=False)
            file.write("\n")
    else:
        assert agrees(result,read("results.json"))
        print("1413 cubic clocks, raw/projective gates, tree hierarchy and shared weak spectra verified")
