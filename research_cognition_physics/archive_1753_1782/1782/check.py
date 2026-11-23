"""Test the full one-singlet tree mass matrix against reactor frequencies and weights."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))


def symmetric_box(pair, n):
    val, error = pair
    return [val - n*error, val + n*error]


def positive_atmospheric_box(triple, n):
    return [triple[0] - n*triple[1], triple[0] + n*triple[2]]


def complement(interval):
    return [1 - interval[1], 1 - interval[0]]


def small_root(amplitude):
    return amplitude / (2 * (1 + math.sqrt(1 - amplitude)))


def detector_bounds(dsolar, z, x):
    """Enclose coherent P over all unresolved fast phases and experimental true L/E band."""
    probe = CFG["spectral_probe"]
    half = probe["relative_L_over_E_half_width"]
    central_dm = CFG["juno"]["d21_eV2"][0]
    phases = [math.pi/2*dsolar[0]/central_dm*(1-half),
              math.pi/2*dsolar[1]/central_dm*(1+half)]
    assert 0 < phases[0] < math.pi/2 < phases[1] < math.pi
    minimum_sine_squared = min(math.sin(t)**2 for t in phases)
    close_pair_oscillation_floor = 4*min(t*(1-t) for t in z)
    pair_max = math.sqrt(1 - close_pair_oscillation_floor*minimum_sine_squared)
    pair_min = max(0., 1-2*z[1], 2*z[0]-1)
    pmin = max(0., (1-x[1])*pair_min-x[1], x[0]-(1-x[0])*pair_max)**2
    pmax = (x[1]+(1-x[1])*pair_max)**2
    delta = (probe["Vcc_eV_per_rhoYe"]*probe["max_rhoYe_g_per_cm3"]*
             probe["max_L_km"]*probe["inverse_eV_per_km"]/2)
    return {"vacuum":[pmin, pmax], "arbitrary_standard_matter_in_budget":
            [max(0., pmin-delta),min(1.,pmax+delta)], "uniform_P_matter_budget":delta}


def disjoint_budget(left, right):
    gap = max(left[0]-right[1], right[0]-left[1])
    assert gap > 0
    return {"min_probability_gap_before_detector_tails":gap,
            "sufficient_outside_true_L_over_E_fraction_less_than":gap/(1+gap)}


def exact_matrix_check():
    rng = np.random.default_rng(1782)
    worst = {"squared_mass_relative_span":0., "electron_weights":0., "survival":0.,
             "sterile_probability":0.}
    for _ in range(120):
        d = rng.normal(size=3) + 1j*rng.normal(size=3)
        d *= rng.uniform(.4,1.8)/np.linalg.norm(d)
        soft = rng.uniform(.05,3)*np.exp(1j*rng.uniform(-math.pi,math.pi))
        mass = np.zeros((4,4),complex)
        mass[:3,3] = d
        mass[3,:3] = d
        mass[3,3] = soft
        msq, u = np.linalg.eigh(mass.conj().T @ mass)
        large = (math.hypot(abs(soft),2*np.linalg.norm(d)) + abs(soft))/2
        middle = np.vdot(d,d).real/large
        r = middle/large
        t = abs(d[0])**2/np.vdot(d,d).real
        w = np.array([1-t,t/(1+r),t*r/(1+r)])
        measured = abs(u[0])**2
        worst["squared_mass_relative_span"] = max(worst["squared_mass_relative_span"],
            np.max(np.abs(msq-[0,0,middle**2,large**2]))/large**2)
        worst["electron_weights"] = max(worst["electron_weights"],
            np.max(np.abs(w-[sum(measured[:2]), *measured[2:]])))
        factor = rng.uniform(-12.,12.)
        prop = (u*np.exp(-1j*factor*msq)) @ u.conj().T
        pee = abs(np.dot(w,np.exp(-1j*factor*np.array([0,middle**2,large**2]))))**2
        pes = 4*t*r/(1+r)**2*math.sin((large**2-middle**2)*factor/2)**2
        worst["survival"] = max(worst["survival"],abs(abs(prop[0,0])**2-pee))
        worst["sterile_probability"] = max(worst["sterile_probability"],abs(abs(prop[3,0])**2-pes))
    assert max(worst.values()) < 1e-12
    return {"complex_mass_matrices":120,"all_errors_below":1e-12}


def calculate():
    j, b = CFG["juno"], CFG["daya_bay"]
    widths = j["marginal_widths"]
    solar_widths = widths + j["extra_solar_width_for_ordering_and_external_constraint"]
    d = symmetric_box(j["d21_eV2"],solar_widths)
    s = symmetric_box(j["sin_squared_theta12"],solar_widths)
    amp = symmetric_box(b["sin_squared_two_theta13"],b["widths"])
    xsmall = [small_root(v) for v in amp]
    xlarge = complement(xsmall)
    target = detector_bounds(d,s,xsmall)
    out = {"round":1782, "experimental_rectangles":{"solar_frequency_eV2":d,
           "solar_share":s, "isolated_fast_amplitude":amp, "isolated_weight_roots":[xsmall,xlarge]},
           "target_three_pole_solar_dip":target, "assignments":{},
           "massless_duplicate_and_complex_basis_checks":exact_matrix_check()}
    for ordering in ("NO","IO"):
        a = positive_atmospheric_box(j[ordering+("_positive_d31_eV2" if ordering=="NO" else "_absolute_d31_eV2")],widths)
        if ordering == "NO":
            r = [math.sqrt(d[0]/a[1]),math.sqrt(d[1]/a[0])]
        else:
            r = [math.sqrt(a[0]/(a[0]+d[1])),math.sqrt(a[1]/(a[1]+d[0]))]
        o = {"absolute_atmospheric_frequency_eV2":a, "middle_over_heaviest_mass_required":r,
             "charged_current_ratio_branches":{}, "frequency_and_fast_amplitude_only":{}}
        if ordering == "NO":
            for flip, zz in [("direct",s),("solar_swap",complement(s))]:
                u = [r[0]*zz[0], r[1]*zz[1]]
                vals = [4*p/(1+p)**2 for p in u]
                assert min(vals) > amp[1]
                o.setdefault("fast_amplitude_from_frequencies_and_close_pair_share",{})[flip]=vals
        else:
            zz = [n/(1+n) for n in r]
            o["solar_share_from_frequencies_only"]=[zz,complement(zz)]
            assert s[1]<zz[0]
        for tag, xx in [("small_isolated",xsmall),("large_isolated",xlarge)]:
            for label, zz in [("solar",s),("reflected_solar",complement(s))]:
                if ordering == "NO":
                    ratios = [xx[0]/((1-xx[0])*zz[1]),xx[1]/((1-xx[1])*zz[0])]
                else:
                    ratios = [zz[0]/(1-zz[0]),zz[1]/(1-zz[1])]
                gap = max(r[0]-ratios[1],ratios[0]-r[1])
                assert gap > 0
                row = {"heavy_over_middle_electron_weight":ratios,"min_disjoint_gap":gap}
                if ratios[0]>1:
                    row["impossible_for_any_two_distinct_positive_masses"] = True
                elif ordering=="NO":
                    row["absolute_d31_from_solar_frequency_and_electron_weights_eV2"] = [
                        d[0]/ratios[1]**2,d[1]/ratios[0]**2]
                else:
                    row["absolute_d31_from_solar_frequency_and_electron_weights_eV2"] = [
                        d[0]*ratios[0]**2/(1-ratios[0]**2),
                        d[1]*ratios[1]**2/(1-ratios[1]**2)]
                o["charged_current_ratio_branches"][tag+"_"+label] = row
            if ordering=="NO":
                zz = [xx[0]/((1-xx[0])*r[1]),xx[1]/((1-xx[1])*r[0])]
                if zz[0]>1:
                    o["frequency_and_fast_amplitude_only"][tag]={"excluded_by_close_pair_fraction_min":zz[0]}
                    continue
            else:
                zz = [n/(1+n) for n in r]
            probe = detector_bounds(d,zz,xx)
            probe.update(disjoint_budget(target["arbitrary_standard_matter_in_budget"],
                                         probe["arbitrary_standard_matter_in_budget"]))
            probe["predicted_close_pair_fraction_or_complement"]=zz
            o["frequency_and_fast_amplitude_only"][tag]=probe
        out["assignments"][ordering]=o
    probe = CFG["spectral_probe"]
    out["nominal_dip_energy_MeV"]=1000*2*probe["oscillation_phase_per_dm2_km_over_GeV"]*j["d21_eV2"][0]*probe["nominal_L_km"]/math.pi
    return out


def close(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)


if __name__ == "__main__":
    cmd = argparse.ArgumentParser()
    cmd.add_argument("--check",action="store_true")
    target = HERE/"results.json"
    numbers = calculate()
    if cmd.parse_args().check:
        assert close(numbers,json.loads(target.read_text(encoding="utf-8")))
        print("1782 checked: all eight weight branches, exact 4x4 matrices and matter-safe solar dip")
    else:
        target.write_text(json.dumps(numbers,indent=2)+"\n",encoding="utf-8")
