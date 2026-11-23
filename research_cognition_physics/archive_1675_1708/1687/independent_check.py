"""R1687: a competing type-II doublet direction, retaining the actual radial Yukawas."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def matrix(edges):
    a, b, c, e, f, r, d = edges
    return np.array([[0, a, 0], [b, c, r + d], [0, e, f]])


def form_factors(z):
    if z < 0.002:
        # asin(sqrt(z))**2 = sum_{n>=1} a_n z**n.
        a = [0] + [4**n / (2*n*n*math.comb(2*n, n)) for n in range(1, 10)]
        return (2 * sum((a[n]-a[n+1])*z**(n-1) for n in range(1, 9)),
                2 * sum(a[n]*z**(n-1) for n in range(1, 9)))
    f = (math.asin(math.sqrt(z))**2 if z <= 1 else
         -((2*math.acosh(math.sqrt(z)) - 1j*math.pi)**2)/4)
    return 2*(z + (z-1)*f)/z**2, 2*f/z


def bisect(function, lo, hi):
    low_value = function(lo)
    assert low_value * function(hi) < 0
    for _ in range(75):
        mid = (lo+hi)/2
        if function(mid)*low_value > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def calculate():
    inputs = read("inputs.json")
    powers = np.array(read(inputs["down_powers"])["powers_A_B_C_E_F_R_D"])
    edge = (np.array(read(inputs["down_paths"])["quark_path_amplitudes"]) *
            np.exp(1j * math.pi * inputs["phase_rad_over_pi"] * powers))
    u, m, vh = np.linalg.svd(matrix(edge))
    w = (np.diag(u.conj().T @ matrix(powers*edge) @ vh.conj().T)/m)[::-1]
    assert abs(w.sum()-4) < 1.e-11
    archived = runpy.run_path(str(HERE.parent/"1373"/"check.py"))
    quarks = archived["common_inputs"](archived["read"]("inputs.json"))[2][:6]
    masses = np.array([q[0] for q in quarks])
    mh = inputs["higgs_GeV"]
    forms, odd = np.array([form_factors((mh/(2*mf))**2) for mf in masses]).T
    heavy = read(inputs["chain"])
    heavyA = sum(slope*form_factors((mh/(2*mass*1000))**2)[0] for mass, slope in
                 zip(heavy["conditional_masses_TeV_descending"],
                     heavy["individual_log_mass_radial_slopes"]))
    beta2 = 1 - (2*masses[-1]/mh)**2
    fbb = 0.01
    chain = read("../1685/inputs.json")
    eps = chain["scalar_magnitude_over_cutoff"]
    cutoff = chain["cutoff_TeV"]*1000
    f = math.sqrt(2)*eps*cutoff*math.sqrt(1+4*eps**2+9*eps**4)
    v = archived["read"]("inputs.json")["v_MZ_GeV"]
    ratio_fv = f/v

    def predict(x, T):
        # x = (v/f)*e/c, l=d/c=-2*x*tan(beta) sets k_t/c to +1.
        tanbeta = math.sqrt(T/2)
        ell = -2*x*tanbeta
        c = 1/math.sqrt(1 + ell**2 + (ratio_fv*x)**2)
        k = np.r_[[1-2*x, 1-2*x, 1], 1+x*(T+w.real)]
        p = np.r_[[0, 0, 0], x*w.imag]
        loop = (abs(k@forms + x*heavyA)**2 + abs(p@odd)**2)/abs(sum(forms))**2
        group = (1-fbb)*loop + fbb*(k[-1]**2 + p[-1]**2)
        return dict(x=x, two_tan_beta_squared=T, tan_beta=tanbeta,
                    h_along_vev_orthogonal_doublet_radial=[c, c*ell, c*x*ratio_fv],
                    top_over_V=1., bottom_scalar_over_V=k[-1], bottom_odd_over_V=p[-1],
                    gg_loop_rate_over_V_squared=loop, grouped_ggF_bbH_over_VBF=group,
                    bottom_over_ZZ_branching_ratio_to_SM=k[-1]**2+p[-1]**2/beta2)

    x0 = -0.05
    single_T = bisect(lambda t: predict(x0, t)["grouped_ggF_bbH_over_VBF"] - 1, 30, 300)
    wide = predict(x0, single_T)

    def bottom_reflected_T(x):
        return (-math.sqrt(1-(x*w[-1].imag)**2/beta2)-1)/x - w[-1].real

    def difference(x):
        return predict(x, bottom_reflected_T(x))["grouped_ggF_bbH_over_VBF"]-1

    bracket = [-0.04, -0.003]
    solution = bisect(difference, *bracket)
    both = predict(solution, bottom_reflected_T(solution))
    for sample in (wide, both):
        assert abs(float(sample["grouped_ggF_bbH_over_VBF"])-1) < 5.e-13
        assert abs(np.linalg.norm(sample["h_along_vev_orthogonal_doublet_radial"])-1) < 1.e-13
    assert wide["bottom_over_ZZ_branching_ratio_to_SM"] > 15
    assert abs(both["bottom_over_ZZ_branching_ratio_to_SM"]-1) < 5.e-13
    assert both["bottom_scalar_over_V"] < 0
    # BR(Z*)+BR(bb) <= 1. Transfer the ratios after normalizing to the SM.
    # Lower reference fractions .55 and .025 conservatively round CERN YR3's
    # .575 and .0266.  ATLAS 1909.02845 Fig.7: .93[-.27,+.38], .87[-.12,+.14].
    ref = dict(bb_floor=.55, zz_floor=.025, bb_lower_width=.27, bb_central=.93,
               diphoton_central=.87, diphoton_upper_width=.14)
    bound = {}
    vector, low, _, sm_vector = inputs["atlas"]["vector"]
    for n in (1, 2):
        vlow = (vector-n*low)/sm_vector
        bmin = ref["bb_central"]-n*ref["bb_lower_width"]
        gmax = ref["diphoton_central"]+n*ref["diphoton_upper_width"]
        c2min = vlow*(ref["zz_floor"]+bmin*ref["bb_floor"])/gmax
        bound[str(n)] = dict(c_squared_lower=c2min,
                             abs_x_upper=math.sqrt(1/c2min-1)/ratio_fv)
    assert abs(solution) > bound["2"]["abs_x_upper"]
    small_x = -0.0005
    small = predict(small_x, bottom_reflected_T(small_x))
    small_c2 = small["h_along_vev_orthogonal_doublet_radial"][0]**2
    normalization = .98
    by_mode = dict(vector=normalization, top=normalization,
                   gluon=normalization*small["grouped_ggF_bbH_over_VBF"])
    for key, predicted in by_mode.items():
        central, minus, plus, standard = inputs["atlas"][key]
        assert central-minus <= standard*predicted <= central+plus
    # Explicit positive branching completion for the three displayed final states.
    # D_b=D_gamma=1 at SM reference branching fractions. The rest is unresolved.
    scale = normalization/small_c2
    shown = [scale*b for b in (.575, .0266, .00227)]
    assert 0 < sum(shown) < 1
    small.update(VBF_diphoton_rate_normalized=normalization, three_rates=by_mode,
                 branching_bb_ZZ_diphoton_other=[*shown, 1-sum(shown)])
    return dict(down_radial_coupling_order_dsb_re_im=np.array([w.real,w.imag]).T.tolist(),
                heavy_Dirac_radial_loop=heavyA, fixed_parameters=dict(f_over_v=ratio_fv, f_bbH=fbb),
                absolute_VBF_and_two_decay_ratios_inputs=ref,
                absolute_VBF_and_two_decay_ratios_bounds=bound,
                top_and_grouped_rate_matched=wide, wrong_sign_root_bracket=bracket,
                wrong_sign_bracket_residuals=[difference(v) for v in bracket],
                top_grouped_and_bottom_width_matched=both,
                wrong_sign_small_x_box_and_budgets_witness=small)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    check = parser.parse_args().check
    def rounded(data):
        if isinstance(data, dict):
            return {k: rounded(v) for k,v in data.items()}
        if isinstance(data, list):
            return [rounded(v) for v in data]
        return float(f"{data:.11g}")
    result = rounded(calculate())
    target = HERE/"independent_results.json"
    if check:
        assert result == read(target)
        print("R1687 independent OK")
    else:
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
