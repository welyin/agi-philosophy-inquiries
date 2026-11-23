"""Round 1545: exact representation traces, finite-window threshold and Higgs rate tests."""
import argparse
import cmath
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def beta_vector(data, n):
    scalars = data["scalar_dynkin_weight"].values()
    weyl = data["weyl_dynkin_weight_per_family"].values()
    return [(11*Fraction(c,3) - 2*Fraction(n,3)*sum(rep[i] for rep in weyl)
             - Fraction(1,3)*sum(rep[i] for rep in scalars))
            for i,c in enumerate(data["adjoint_Casimir"])]


def finite_window(data):
    rows = {}
    for n in data["family_multiples_to_display"]:
        b4, bl, br = beta_vector(data, n)
        assert b4-bl == Fraction(13,3) and bl-br == Fraction(22,3)
        rows[str(n)] = [str(b4),str(bl),str(br)]
    assert beta_vector(data,3)[1] == 3 and beta_vector(data,6)[1] == -1
    log_range = math.log(data["log_window_scale_ratio"])
    return {"b4_bL_bR":rows, "pairwise_differences":["13/3","22/3"],
            "three_to_six_change_in_each_inverse_alpha_increment":-4*log_range/(2*math.pi),
            "critical_common_threshold_over_lower_scale":math.exp(log_range/4)}


def heavy_loop(mh,m):
    tau = (mh/(2*m))**2
    if not 0 < tau <= 1:
        raise ValueError("The displayed comparison uses above-threshold heavy quarks")
    return 2*(tau+(tau-1)*math.asin(math.sqrt(tau))**2)/tau**2


def collider(data):
    old = read(HERE/data["data1377"])
    saved = read(HERE/data["frozen1377"])
    amplitude = complex(*saved["SM_quark_loop_real_imag"])
    atlas = old["atlas_four_leptons"]
    norm = atlas["central_normalization"]
    eta = norm["inclusive_bbh_pb"]*norm["br_zz"]/norm["stage_gg_plus_bb_times_br_pb"]
    g, v = atlas["gg_plus_bb"], atlas["vbf"]
    inflate = lambda line,direction: math.hypot(*line["width_"+direction])
    output = {}
    for widths in (1,2):
        upper = (g["mu"]+widths*inflate(g,"plus"))/(v["mu"]-widths*inflate(v,"minus"))
        pure_upper = (upper-eta)/(1-eta)
        radius = abs(amplitude)*math.sqrt(pure_upper)
        kappa_sum = .75*(math.sqrt(radius**2-amplitude.imag**2)-amplitude.real)
        assert abs(kappa_sum-saved["marginal_rectangles"][str(widths)]
                   ["upper_sum_nonnegative_heavy_kappa"]) < 1e-9
        window = {}
        for n in (6,9):
            shift = 8*(n-3)/3
            total = amplitude+shift
            pure_min = abs(total/amplitude)**2
            # A disk of arbitrary complex extra amplitudes, or a positive rate transfer.
            displacement = -(abs(total)-radius)*cmath.exp(1j*cmath.phase(total))
            assert abs(abs(total+displacement)-radius)<1e-10
            window[str(n)] = {
                "minimal_mixed_rate":(1-eta)*pure_min+eta,
                "minimum_extra_amplitude":abs(displacement),
                "minimum_fraction_of_heavy_limit_shift":abs(displacement)/shift,
                "maximum_gg_only_transfer":pure_upper/pure_min,
                "maximum_mean_nonnegative_heavy_response":kappa_sum/(2*(n-3))}
        output[str(widths)] = {"observed_upper":upper,"allowed_pure_rate":pure_upper,
                               "family_branches":window}
    form = heavy_loop(old["on_shell_h_GeV"],data["extra_quark_finite_witness_GeV"])
    assert 4/3 < form < 2
    x = 6*form
    assert abs(abs((amplitude+x-x)/amplitude)-1)<1e-12
    return {"ATLAS_central_ratio":g["mu"]/v["mu"],"SM_amplitude":[amplitude.real,amplitude.imag],
            "bb_fraction_ceil":eta,"working_rectangles":output,
            "Ng6_finite_witness":{"quark_GeV":data["extra_quark_finite_witness_GeV"],
              "six_heavy_shift":x,"real_counterterm_for_SM_rate":-x,
              "mixed_ratio_without_counterterm":(1-eta)*abs((amplitude+x)/amplitude)**2+eta}}


def rounded(v):
    if isinstance(v,dict):
        return {k:rounded(x) for k,x in v.items()}
    if isinstance(v,list):
        return [rounded(x) for x in v]
    return round(v,10) if isinstance(v,float) else v


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    cfg=read(HERE/"inputs.json")
    out=rounded({"running":finite_window(cfg),"Higgs":collider(cfg)})
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    elif args.check:
        assert read(target)==out, "Numerical record differs"
        print("1545: PS Dynkin traces, threshold blindness and two ATLAS rectangles verified.")
    else:
        print(json.dumps(out,indent=2))
