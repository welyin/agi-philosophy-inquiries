"""Joint finite gauge and type-II bounds for an adjoint-split 15 with possibly hidden partners."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOW = 1e7
B = np.array([[3.6,3.2,.2],[4.,0.,3.],[0.,5.,2.]])
B0 = np.array([8.6,3.,-1.])
W = np.array([5.,-12.,7.])
PI = math.pi


def read(path):
    return json.loads((HERE / path).read_bytes())


def positive_root(function, left, right):
    assert function(left) < 0 < function(right)
    for _ in range(90):
        mid = (left+right)/2
        if function(mid) > 0:
            right = mid
        else:
            left = mid
    return (left+right)/2


def calculate():
    data = read("../1369/inputs.json")
    lattice = read("../1769/inputs.json")
    old = read("../1698/results.json")["shared_witness"]
    fits = read("../1802/independent_results.json")
    strong = lattice["strong_lattice_5f"]["ALPHA_2026"]
    top = lattice["top_MS_matching_GeV"]
    shifts = []
    for i,j,k,l in itertools.product((-2,2),repeat=4):
        row = []
        for key, width in (("g1_GUT",i),("g2",j)):
            coupling, sigma = data[key]
            row.append(4*PI*((coupling+width*sigma)**-2-coupling**-2))
        row.append(1/(strong[0]+k*strong[1])+math.log((top[0]+l*top[1])/data["mu_ref_GeV"])/(3*PI)
                   -4*PI/data["g3"][0]**2)
        shifts.append(row)
    shifts = np.array(shifts)
    invariant = 2*PI*shifts@W
    diff12 = 2*PI*(shifts[:,0]-shifts[:,1])
    ilo, ihi = min(invariant), max(invariant)
    clo, chi = math.exp(ilo/21), math.exp(ihi/21)
    # Capped ratios s=min(S,U)/T and z=min(Z,U)/T; also covers S or Z above U.
    slo = positive_root(lambda s: 2*chi*s**(17/7)+s-1,0.,1.)-1e-12
    shi = max(1.,math.exp(-ilo/30))
    zlo, zhi = clo*slo**(17/7),max(1.,shi)
    u0, a0 = old["U_GeV"],old["inverse_alpha_U"]
    # h = d12/5.6 - I/42 - 9 log(s)/14; correlated affine terms use the same corners.
    base = diff12/5.6-invariant/42
    hlo = min(base)-(9/14)*math.log(shi)
    hhi = max(base)-(9/14)*math.log(slo)
    ulo,uhi = u0*math.exp(hlo),u0*math.exp(hhi)
    # If S or Z opens before T, bound its effect on the weak couplings over that interval.
    # Inverse couplings decrease monotonically for a=1,2. T does not run below itself.
    sublogs = np.array([0.,max(0.,-math.log(slo)),max(0.,-math.log(zlo))])
    amin_below_T = a0+shifts.min(axis=0)-(B0*hhi+B@sublogs)/(2*PI)
    gmax = 4*PI/amin_below_T[:2]
    assert min(amin_below_T[:2])>1
    kappa_factor = math.exp(-2*(.6*gmax[0]+3*gmax[1])*math.log(uhi/LOW)/(16*PI**2))
    t_caps = {}
    for case,d in fits["five_cases"].items():
        k0 = d["old_absolute_matching"]["kappa_at_10PeV_GeV_inverse"]
        t_caps[case] = 4*PI/(k0*kappa_factor)
    thi = max(t_caps.values())*(1+1e-9)
    assert thi*max(shi,zhi)<ulo
    # No clipping can survive the mass-source upper bound. Tighten with the uncapped triangle.
    su = positive_root(lambda s: 2*clo*s**(17/7)-s-1,.4,3.)+1e-12
    zhi = min(zhi,(1+su)/2)
    hlo = min(base)-(9/14)*math.log(su)
    ulo = u0*math.exp(hlo)
    # Conservative rectangle in log variables gives a valid upper inverse coupling.
    au_upper = a0+max(shifts[:,0])-(8.6*hlo+7*math.log(ulo/thi)-3.2*math.log(su)-.2*math.log(zhi))/(2*PI)
    # Unit flux/charge normalization; C_X = g_U^2/(2 M_X^2) only for the named generator.
    at_cross = 2*PI/(au_upper*uhi**2)
    target = fits["two_marginal_gauge_box"]["triplet_plus_sextet_only_witnesses"]["ALPHA_central"]
    z_adj_max = (target["triplet_GeV"]+target["sextet_GeV"])/2
    # Construct at least one actual point for all equations, using the original interior solution.
    witness = read("../1800/independent_results.json")["threshold_cases"]["1698_two_pairs"]["experimental_working_boxes"]["0"]
    masses = np.array(witness["lowest_scale_positive_mass_witness_GeV_T_S_Z"])
    uu, aa = witness["witness_U_GeV_and_inverse_alpha_U"]
    np.testing.assert_allclose(2*masses[2],masses[:2].sum(),rtol=5e-9)
    assert max(masses)<uu and aa>1
    assert slo<masses[1]/masses[0]<su and zlo<masses[2]/masses[0]<zhi
    assert ulo<uu<uhi and aa<au_upper
    # Verify the signs and normalization against the full three-coupling reconstruction.
    da = aa-a0+(B0*math.log(uu/u0)+B@np.log(uu/masses))/(2*PI)
    input3 = 1/strong[0]+math.log(top[0]/data["mu_ref_GeV"])/(3*PI)-4*PI/data["g3"][0]**2
    np.testing.assert_allclose(da,[0.,0.,input3],atol=2e-8)
    np.testing.assert_allclose(W@B,[-30,51,-21])
    assert abs(W@B0)<1e-12
    return {"two_width_gauge_invariant_range":[ilo,ihi],
            "capped_mass_ratios_S_T":[slo,shi],"capped_mass_ratios_Z_T":[zlo,max(1.,shi)],
            "uncapped_ratios_outer_S_T":[slo,su],"uncapped_ratios_outer_Z_T":[zlo,zhi],
            "all_crossing_scales_outer_GeV":[ulo,uhi],"g1_g2_squared_upper_before_T":gmax.tolist(),
            "type_II_norm_at_T_over_norm_at_10PeV_lower":kappa_factor,
            "T_upper_by_mass_source_GeV":t_caps,
            "max_all_five_S_and_Z_GeV":[su*thi,zhi*thi],
            "weakly_coupled_inverse_alpha_U_upper":au_upper,
            "gU_squared_over_2U_squared_lower_GeV_minus2":at_cross,
            "Z_in_two_fragment_central_example_over_adjoint_max":target["new_crossing_GeV"]/z_adj_max,
            "compatible_gauge_mass_example":{"T_S_Z_GeV":masses.tolist(),"U_GeV":uu,"inverse_alpha_U":aa}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    out = json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    target = HERE/"adjoint_threshold_results.json"
    if parser.parse_args().check:
        assert json.loads(out)==json.loads(target.read_bytes())
    else:
        target.write_text(out,encoding="utf-8")
    print("1804 adjoint 15",hashlib.sha256(out.encode()).hexdigest()[:12])
