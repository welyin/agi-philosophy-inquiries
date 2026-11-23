"""Test two shared electroweak flavour responses against a same-resonance three-channel fit."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())


def pair(z):
    return [float(z.real), float(z.imag)]


def frames(y):
    left, mass, right_h = np.linalg.svd(y)
    return left[:, ::-1], mass[::-1], right_h.conj().T[:, ::-1]


def canonical_left(d, running_dn, running_up):
    """Identify both quark flags, using the d column to fix harmless eigenvector phases."""
    ud = frames(d)[0]
    uq = frames(running_up)[0]
    actual = uq.conj().T @ frames(running_dn)[0]
    phases = ud[:, 0] / actual[:, 0]
    assert max(abs(abs(phases) - 1)) < 2e-8
    return np.diag(phases / abs(phases)) @ uq.conj().T


def transported(texture, insertion, running, left, low_masses):
    # texture @ right = left @ running. The low-energy right basis is the unit matrix.
    right = np.linalg.solve(texture, left @ running)
    defect = np.linalg.norm(right.conj().T @ right - np.eye(3), 2)
    assert defect < 1e-6, defect
    k = np.linalg.solve(running, left.conj().T @ insertion @ right)
    return np.diag(low_masses) @ k, defect


def rectangle():
    n = CFG["marginal_multiplier"]
    return {k: (max(0, x - n*m), x + n*p)
            for k, (x, m, p) in CFG["CMS_B_over_ZZ_normalized"].items()}


def limits(q, envelope):
    """All pairwise constraints after eliminating the one common positive scale."""
    cuts = []
    interval = rectangle()
    for i, (li, ui) in interval.items():
        for j, (lj, uj) in interval.items():
            if i == j or lj == 0:
                continue
            # Width correction w between 1 and envelope. Use a conservative relaxed box.
            k = ui*envelope[j]/lj
            a, c = 1-q[i], q[i]
            b, d = 1-q[j], q[j]
            aa = abs(a)**2 - k*abs(b)**2
            bb = a*c.conjugate() - k*b*d.conjugate()
            cc = abs(c)**2 - k*abs(d)**2
            assert abs(aa) > 1e-10
            center = -bb.conjugate()/aa
            radius = math.sqrt(max(0, abs(center)**2-cc/aa))
            cuts.append((aa, bb, cc, center, radius, i+"/"+j))
    return cuts


def closest(cuts):
    def valid(u):
        return all(a*abs(u)**2 + 2*(b*u).real + c < 4e-11 for a,b,c,*_ in cuts)
    candidates = [0j, 1+0j]
    for index, (_, _, _, x, r, _) in enumerate(cuts):
        unit = x/abs(x) if x else 1
        candidates.extend([x-r*unit, x+r*unit])
        for _, _, _, y, s, _ in cuts[:index]:
            dist = abs(y-x)
            if dist < 1e-13 or dist > r+s+1e-12 or dist < abs(r-s)-1e-12:
                continue
            u = (y-x)/dist
            proj = (dist**2+r**2-s**2)/(2*dist)
            perp = math.sqrt(max(0, r*r-proj*proj))
            candidates.extend([x+u*(proj+1j*perp), x+u*(proj-1j*perp)])
    best = min((z for z in candidates if valid(z)), key=abs)
    active = [name for a,b,c,_,_,name in cuts if abs(a*abs(best)**2+2*(b*best).real+c)<1e-8]
    return {"min_abs_rA_over_rB":abs(best), "witness_complex":pair(best), "active":active}


def test(z, eigen_e, eigen_d, m_e, m_d):
    qe = np.diag(eigen_e)/m_e
    q = {"mu": qe[1], "tau": qe[2], "b": eigen_d[2,2]/m_d[2]}
    mh = CFG["h_mass_GeV"]
    masses = CFG["lepton_masses_GeV"]
    fm = {"mu": masses[1], "tau": masses[2], "b":CFG["max_b_kinematic_mass_GeV"]}
    env = {k:1/(1-4*(v/mh)**2) for k,v in fm.items()}
    k = rectangle()["mu"][1]/rectangle()["tau"][0]
    ratio = m_e[1]/m_e[2]*math.sqrt(k*env["tau"])
    repair = max(0,abs(eigen_e[1,1])-ratio*abs(eigen_e[2,2]))/(1+ratio)
    mt = m_e[2]
    # Off-diagonal prediction uses common observed pole masses of the charged leptons,
    # matching each open chiral matrix element by the geometric mean of the gains.
    to_pole = np.sqrt(np.array(masses)/m_e)
    pole = to_pole[:,None]*eigen_e*to_pole[None,:]
    mtphys, muphys = masses[2], masses[1]
    x, y = (mtphys/mh)**2, (muphys/mh)**2
    lam = math.sqrt((1-x-y)**2-4*x*y)
    common = (abs(pole[1,2])**2+abs(pole[2,1])**2)*(1-x-y)
    spread = 4*math.sqrt(x*y)*abs(pole[1,2]*pole[2,1])
    bta = math.sqrt(1-4*x)
    lfv_lower = lam*(common-spread)/(bta*abs(pole[2,2])**2)
    lfv_upper = lam*(common+spread)/(bta**3*abs(pole[2,2])**2)
    vfrac=CFG["pure_B_two_doublet_positive_witness_vB_over_v"]
    actual_rates={k:[abs(v)**2/vfrac**4,
                     (v.real**2+env[k]*v.imag**2)/vfrac**4] for k,v in q.items()}
    in_box=all(rectangle()[k][0]<=lo<=hi*(1+1e-14) and hi<=rectangle()[k][1]
               for k,(lo,hi) in actual_rates.items())
    return {"z":pair(z), "B_response_over_mass":{k:pair(v) for k,v in q.items()},
            "B_only_real_radial_tau_CP_deg":(np.degrees(np.angle(q["tau"]))+90)%180-90,
            "two_doublet_witness_ZZ_normalized_rates_envelope":actual_rates,
            "two_doublet_witness_in_all_three_windows":in_box,
            "B_only_mu_tau_SM_ratio":abs(q["mu"]/q["tau"])**2,
            "B_only_b_tau_SM_ratio":abs(q["b"]/q["tau"])**2,
            "closest_common_three_modes":closest(limits(q,{k:1 for k in q})),
            "necessary_with_CP_phase_space":closest(limits(q,env)),
            "minimal_repair_for_mu_tau_in_units_B_lepton_norm":repair/np.linalg.norm(eigen_e,2),
            "minimal_repair_for_mu_tau_in_units_tau_mass":repair/mt,
            "LFV_B_only_Gamma_mutau_over_tautau_range":[lfv_lower,lfv_upper],
            "largest_Btautau_to_respect_CMS_LFV_limit_at_SM_production":CFG["CMS_LFV_limit"]/lfv_lower}


def calculate():
    old = runpy.run_path(str(HERE/CFG["charged_background"]))
    fl, rep = old["charged_background"](CFG["step"])
    unpack = old["unpack"]
    d = unpack(rep["down_texture_10PeV"])
    e = unpack(rep["charged_texture_10PeV"]).T
    le, me, re = frames(e.T)
    ld, md, rd = frames(d)
    gamma = rep["gamma_MZ_deg"]
    y = fl["shoot"](gamma,fl["angles"]("direct"),unpack(rep["charged_left_MZ"]),CFG["step"])[0]
    yd, ye, yu = (y[i,:,:3] for i in (1,2,0))
    left = canonical_left(d,yd,yu)
    em, dm = np.asarray(fl["central"]("charged_leptons")), np.asarray(fl["central"]("down"))
    cases = {}
    z0, zt = e[0,1]/d[0,1], e[2,2]/d[2,2]
    options = {"canonical":1., "first_edge":z0, "third_edge":zt,
               "A_vanishes_first_edge":-z0/3, "A_vanishes_third_edge":-zt/3}
    for label, z in options.items():
        b = (d-e/z)/4
        a = d-b
        np.testing.assert_allclose(z*(a-3*b),e,atol=1e-14)
        el = (-3*z*b).T
        mel, dl = transported(e.T,el,ye,np.eye(3),em)
        mdn, dd = transported(d,b,yd,left,dm)
        el_high, b_high = le.conj().T @ el @ re, ld.conj().T @ b @ rd
        cases[label] = {"ranks_A_B":[int(np.linalg.matrix_rank(a,tol=1e-12)),int(np.linalg.matrix_rank(b,tol=1e-12))],
                        "diagonal_gain_matching": test(z,el_high,b_high,me,md),
                        "frozen_background_common_left":test(z,mel,mdn,em,dm)}
    assert cases["third_edge"]["frozen_background_common_left"]["B_only_mu_tau_SM_ratio"]>1e4
    assert cases["A_vanishes_third_edge"]["frozen_background_common_left"]["closest_common_three_modes"]["min_abs_rA_over_rB"] == 0.
    assert all(cases[k]["frozen_background_common_left"]["two_doublet_witness_in_all_three_windows"]
               for k in ("A_vanishes_first_edge","A_vanishes_third_edge"))
    phases = [cases[n]["frozen_background_common_left"]["B_only_real_radial_tau_CP_deg"]
              for n in ("A_vanishes_first_edge","A_vanishes_third_edge")]
    return {"inputs_high_scale_GeV":1e7, "right_unitarity_defects_e_d":[dl,dd],
            "separation_of_accepted_real_radial_CP_angles_deg": abs(phases[0]-phases[1]),
            "CMS_relative_windows":rectangle(), "CMS_mu_tau_max":rectangle()["mu"][1]/rectangle()["tau"][0],
            "cases":cases}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [clean(v) for v in x]
    if isinstance(x,(np.floating,float)):
        return float(f"{x:.10g}")
    return x


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    data=json.dumps(clean(calculate()),ensure_ascii=False,indent=2)+"\n"
    path=HERE/"results.json"
    if args.check:
        assert json.loads(data)==json.loads(path.read_bytes())
    else:
        path.write_text(data,encoding="utf-8")
    print("1793",hashlib.sha256(data.encode()).hexdigest()[:12])
