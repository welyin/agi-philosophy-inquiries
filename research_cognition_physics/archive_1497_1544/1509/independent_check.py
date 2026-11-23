"""Reconstruct the independently tagged OPAL momentum window and its shape contrast."""
import argparse
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent


def interval(central, stat, syst):
    half = 2*(stat+syst)
    return [central-half,central+half]


def audit():
    # Table 5: bin width, gluon dN/dpT (40.1 GeV), uds (45.6 GeV),
    # and the published g/uds differential ratio at matched 40.1 GeV.
    soft_bins = [(0.2,4.51,2.051,2.35),(0.4,2.07,1.031,2.18),
                 (0.6,0.91,0.3876,2.62),(1.0,0.16,0.1074,1.71)]
    g_soft = sum(w*g for w,g,q45,ratio in soft_bins)
    q_soft = sum(w*g/ratio for w,g,q45,ratio in soft_bins)
    q_unmatched = sum(w*q45 for w,g,q45,ratio in soft_bins)
    ratio_check = g_soft/q_soft
    assert abs(ratio_check-2.29)<.006
    g_full, r_full, r_soft = 14.28,1.514,2.29
    q_full = g_full/r_full
    fractions = [g_soft/g_full,q_soft/q_full]
    assert abs((fractions[0]/fractions[1])/(ratio_check/r_full)-1)<1e-13
    soft, full = interval(r_soft,.09,.15),interval(r_full,.019,.034)
    soft_lo, soft_hi = soft
    full_lo, full_hi = full
    assert soft_lo > full_hi
    floor = (soft_lo-full_hi)/(soft_lo+full_hi)
    saturating_x = (soft_lo+full_hi)/2
    assert abs(saturating_x*(1+floor)-soft_lo)<1e-13
    assert abs(saturating_x*(1-floor)-full_hi)<1e-13
    transfers = {}
    for name,x in (("SU3",2.25),("equal_casimirs",1.0)):
        transfers[name]={"x":x, "soft_K":[soft_lo/x,soft_hi/x],
                         "full_K":[full_lo/x,full_hi/x],
                         "minimum_common_fractional_error":
                         max(0.,soft_lo/x-1,1-soft_hi/x,full_lo/x-1,1-full_hi/x)}
    return {"OPAL":"https://arxiv.org/abs/hep-ex/9903027",
            "soft_gluon_mean":g_soft,"soft_uds_matched_mean":q_soft,
            "soft_uds_unmatched_mean":q_unmatched,"raw_unmatched_ratio":g_soft/q_unmatched,
            "soft_matched_ratio_from_table5":ratio_check,
            "soft_fractions_gluon_uds":fractions,"charged_complement_ratio":(g_full-g_soft)/(q_full-q_soft),
            "r_full_two_widths":full,"r_soft_two_widths":soft,
            "Ksoft_over_Kfull":[soft_lo/full_hi,soft_hi/full_lo],
            "Ksoft_over_Kfull_central":r_soft/r_full,
            "minimax_common_error":floor, "saturating_x":saturating_x,
            "minimax_central_error":(r_soft-r_full)/(r_soft+r_full),
            "fixed_x_transfers":transfers,
            "extra_multiplicative_raise_for_modified_jetset":soft_lo/1.31-1}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
