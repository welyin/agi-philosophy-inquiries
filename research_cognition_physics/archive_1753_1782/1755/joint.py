"""Shoot one new weak modulus through the archived split colour/electroweak model."""
import math
import runpy
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def matter_filter(bg):
    # Charge ratios and the two physical thresholds are the four constraints of 1748.
    L, mass, _ = np.linalg.svd(bg["yW"][2, :, :3])
    L, mass = L[:, ::-1], mass[::-1]
    N = L.conj().T @ bg["N"]
    N /= np.linalg.norm(N, axis=0)
    t = 1/np.array([bg["cW"], bg["fW"]])**2 - 1
    gain = (bg["cW"]*bg["fW"])**(1/3)
    h = mass**2

    def evaluate(x, matrix=False):
        z = 1/(1+math.exp(-x[3]))
        g = np.exp(x[:3])*[1, -1, 1]
        X = (np.sqrt([z, 1-z])[:, None]*N.conj().T)*g
        G = X.conj().T @ X
        if matrix:
            u, V = np.linalg.eigh(G)
            return L @ ((V/np.sqrt(1+u)) @ V.conj().T) @ L.conj().T, z
        diag = np.diag(G).real
        # The signed off-diagonal entries drop out of both charged spectra.
        direct = np.diag(np.linalg.solve(np.eye(3)+G, G)).real
        return np.array([sum(diag)/sum(t)-1, np.linalg.det(X@X.conj().T).real/np.prod(t)-1,
                         (diag @ (1/h))/sum(1/h)/(gain**-2-1)-1,
                         (direct@h)/sum(h)/(1-gain**2)-1])

    x = np.array([math.log(.66), math.log(.25), math.log(.195), -.33])
    for _ in range(15):
        v = evaluate(x)
        if max(abs(v)) < 5e-13:
            break
        directions = np.eye(4)*2e-5
        jac = np.column_stack([(evaluate(x+d)-evaluate(x-d))/4e-5 for d in directions])
        x -= np.linalg.solve(jac, v)
    assert max(abs(evaluate(x))) < 2e-11
    F, z = evaluate(x, True)
    assert 0 < z < .60
    np.testing.assert_allclose(np.linalg.eigvalsh(F), [bg["cW"], bg["fW"], 1], atol=2e-12)
    return F, z


def point(here, us, ub, cb, step, mother_shifts=None):
    """Old API takes offsets. The reference means/widths cancel at the call boundary."""
    port = runpy.run_path(str(here/"../1748/check.py"))
    old = port["OLD"]

    def load(path):
        model = runpy.run_path(str(path))
        cfg = model["CFG"]
        scales = model["DIRECT"]
        offsets = {k: (v-scales["V"+k][0])/math.hypot(*scales["V"+k][1:])
                   for k, v in (("us", us), ("ub", ub), ("cb", cb))}
        shifts = dict(up=0, g3=0, top=0, charm=0)
        if mother_shifts is not None:
            if set(mother_shifts) != set(shifts):
                raise ValueError("Expected up, g3, top and charm marginal shifts")
            shifts.update(mother_shifts)
        cfg["nearby_measured_inputs_in_marginal_widths"][old["CFG"]["parent_point"]] = dict(
            weak=offsets, **shifts)
        return model

    old["mother"].__globals__["runpy"] = SimpleNamespace(run_path=load)
    bg = old["mother"](step)
    F, z = matter_filter(bg)
    scale, row = port["daughter"](bg, F)
    vu = row["MZ_moduli"]
    # Finite weak and colour matching only perturb the originally prescribed bare CKM.
    assert abs(us-vu["Vus"]) < 2e-7
    factors, _ = bg["model"]["MATCH"]["thresholds"](
        dict(moduli_MZ=vu, weak_offsets={k[1:]:(v-bg["model"]["DIRECT"][k][0])/math.hypot(*bg["model"]["DIRECT"][k][1:])
                                        for k,v in vu.items()}, gamma=row["MZ_gamma"]),
        bg["model"]["MATCH"]["source_thresholds"]())
    d,s,b = row["MZ_down"]
    obs = row["neutrinos"]
    return dict(trial_us=us, low_moduli=vu, gamma_deg=row["MZ_gamma"], z=z,
                down_over_s=d/s/factors[0], bottom_over_s=b/s/factors[2], tau_calibration=scale,
                sun_calibrated_dm31_eV2=obs["dm31_eV2"],sun_calibrated_dm32_eV2=obs["dm32_eV2"],
                angle_readings=obs["angle_readings"], SK_T2K_FC68_90=obs["SK_T2K_FC68_90"],
                effective_double_beta_meV=obs["mbb_meV"],port_overlap_at_color=row["color_scale_c"],
                delta_lepton_deg=obs["delta_deg"], signed_return=row["signed_return_before_and_after_actual_color_gate"])


def shoot(here, target, ub, cb, step, mother_shifts=None):
    # Initial abscissae use only the chosen strong number and the independent b tail.
    guess = math.sqrt(target + 52*ub**2)/(1+target)**.25
    x, y = guess-.0005, guess+.0005
    p, q = point(here,x,ub,cb,step,mother_shifts), point(here,y,ub,cb,step,mother_shifts)
    for _ in range(6):
        if abs(q["down_over_s"]-target) < 2e-11:
            break
        new = y-(q["down_over_s"]-target)*(y-x)/(q["down_over_s"]-p["down_over_s"])
        x,y,p,q = y,new,q,point(here,new,ub,cb,step,mother_shifts)
    assert abs(q["down_over_s"]-target) < 2e-10
    return q
