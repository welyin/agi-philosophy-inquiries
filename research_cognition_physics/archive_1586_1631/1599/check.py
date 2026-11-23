"""Single right-handed flavor projector in three two-body missing-energy lines."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def box(values):
    center, *errors = values
    width = 2 * sum(errors)
    return [center - width, center + width]


def divide(a, b):
    assert min(b) > 0
    return [a[0] / b[1], a[1] / b[0]]


def width_factor(m, daughter, form, lifetime, hbar):
    return lifetime * (m * (1 - (daughter / m) ** 2)) ** 3 * form ** 2 / (64 * math.pi * hbar)


def factor_range(m, daughter, form, clock, hbar):
    return [width_factor(box(m)[0], box(daughter)[1], box(form)[0], box(clock)[0], hbar),
            width_factor(box(m)[1], box(daughter)[0], box(form)[1], box(clock)[1], hbar)]


def ckm_modulus_box(tree):
    ulo, uhi = tree["Vus"]
    clo, chi = tree["Vcb"]
    rlo, rhi = tree["ratio_ub_cb"]
    zlo, zhi = rlo * clo, rhi * chi
    c13 = math.sqrt(1 - zhi ** 2)
    s12, s23 = uhi / c13, chi / c13
    c12, c23 = math.sqrt(1 - s12 ** 2), math.sqrt(1 - s23 ** 2)
    return {
        "u": [[c12*c13, 1.], [ulo, uhi], [zlo, zhi]],
        "c": [[ulo*c23-s23*zhi, s12+s23*zhi], [c12*c23-s12*s23*zhi, 1.], [clo, chi]],
        "t": [[ulo*clo-zhi, s12*s23+zhi], [c12*clo-s12*zhi, s23+s12*zhi], [c23*c13, 1.]],
    }


def calculate():
    inputs = read("inputs.json")
    tags = read(inputs["previous_tags"])
    ys = [box(v) for v in read(inputs["down_masses"])["yukawa"]["down"]]
    fK = read(inputs["kaon_trees"])["K_form"]
    clock = read(inputs["charged_B_clock"])["charged_tau_control"]
    mB = clock["mass_Bplus_GeV"]
    mK, mpi = inputs["charged_K_mass_GeV"], inputs["charged_pi_mass_GeV"]
    tauB = [x * 1e-12 for x in clock["tau_Bplus_ps"]]
    tauK = [x * 1e-9 for x in inputs["charged_K_lifetime_ns"]]
    hbar = inputs["hbar_GeV_s"]
    channels = {"Kpi": (mK, mpi, fK, tauK), "BK": (mB, mK, inputs["form_BK"], tauB),
                "Bpi": (mB, mpi, inputs["form_Bpi"], tauB)}
    nominal = {ch: width_factor(*(v[0] for v in data), hbar) for ch, data in channels.items()}
    boxes = {ch: factor_range(*data, hbar) for ch, data in channels.items()}
    kaon = inputs["kaon_90_inputs"]["NA62_public_sample_2016_24_dedicated_reanalysis"]
    rho = inputs["expanded_new_amplitude_fraction"]
    q = inputs["expanded_kaon_relative_response_min"]
    amp = (1 + rho) / (1 - rho)
    rows = {}
    bounds = ckm_modulus_box(tags["two_width_box"])
    bk_bands = []
    bpi_bands = []
    for j, (name, row) in enumerate(bounds.items()):
        yd, ys_quark, yb = ys
        rel = {}
        for ch, denominator in [("BK", 0), ("Bpi", 1)]:
            v = divide(row[2], row[denominator])
            y = divide(ys[denominator], yb)
            kin = divide(boxes[ch], boxes["Kpi"])
            rel[ch] = [k * (a * b)**2 for a, b, k in zip(v, y, kin)]
        bk_bands.append(rel["BK"])
        bpi_bands.append(rel["Bpi"])
        central = {}
        for branch, record in tags["central_gamma_offsets_degrees"].items():
            pd, ps, pb = record["tags"][name]["down_flavor_fractions"]
            Rk = nominal["BK"] / nominal["Kpi"] * pb / pd
            Rpi = nominal["Bpi"] / nominal["Kpi"] * pb / ps
            assert rel["BK"][0] <= Rk <= rel["BK"][1]
            assert rel["Bpi"][0] <= Rpi <= rel["Bpi"][1]
            central[branch] = {"BK_over_Kpi": Rk, "Bpi_over_Kpi": Rpi,
                "f_lower_GeV": math.sqrt(nominal["Kpi"]*pd*ps/kaon)}
            # Flavor rephasings let us display the two positive currents with real components.
            d, s, b = map(math.sqrt, (pd, ps, pb))
            orthogonal = [s/math.hypot(d,s), -d/math.hypot(d,s), 0.]
            assert abs(d*orthogonal[0] + s*orthogonal[1]) < 2e-15
            assert abs(d*s + (pd+ps)*orthogonal[0]*orthogonal[1]) < 2e-15
            assert abs(s*b) > 0 and abs(d*b) > 0

        p_d_min = tags["two_width_box"]["down_tag_fraction_lower_u_c_t"][j]
        ws_over_wd = divide(row[1], row[0])[0] * divide(yd, ys_quark)[0]
        f_min = p_d_min * ws_over_wd * math.sqrt(boxes["Kpi"][0] / kaon)
        caps = {key: {ch: limit * v[1] * amp**2 / q for ch, v in rel.items()}
                for key, limit in inputs["kaon_90_inputs"].items()}
        target = inputs["B_factory_benchmark_rate"]
        cancel_if_B = {ch: 1 - (1+rho)*math.sqrt(kaon*v[1]/(q*target))
                       for ch, v in rel.items()}
        rows[name] = {"central": central, "two_width_rate_ratio": rel,
                      "f_lower_GeV_in_full_box": f_min,
                      "f_lower_GeV_in_transfer_stress": f_min*(1-rho)*math.sqrt(q),
                      "B_rate_caps_with_transfer": caps,
                      "K_amplitude_cancel_fraction_needed_for_one_part_per_million": cancel_if_B}

    critical = []
    for one, two in zip(bk_bands, bk_bands[1:]):
        ratio = (two[0] / one[1]) ** .25
        critical.append((ratio-1)/(ratio+1))
    assert min(critical) > rho
    all_upper = {ch: max(rec["B_rate_caps_with_transfer"][
        "NA62_public_sample_2016_24_dedicated_reanalysis"][ch] for rec in rows.values())
        for ch in ["BK", "Bpi"]}
    cap_sum = sum(all_upper.values())
    # Bound applies even if both mesons in every Belle event are charged.
    nB = inputs["Belle_deliberate_overestimate_number_charged_mesons"]
    m_test = inputs["nearby_B_scan_mass_GeV"]
    b = mB[0]
    dp = {}
    for name, daughter in [("K",mK[0]),("pi",mpi[0])]:
        p0 = (b*b-daughter**2)/(2*b)
        p_test = math.sqrt((b*b-(daughter+m_test)**2)*(b*b-(daughter-m_test)**2))/(2*b)
        dp[name] = (p0-p_test)*1e6
    return {"branching_normalization_GeV_squared": nominal,"two_width_normalization": boxes,
        "rows": rows, "gap_closure_amplitude_fractions_for_adjacent_BK_bands": critical,
        "max_B_rate_over_all_tags_with_transfer": all_upper, "all_B_modes_mean_events_Belle_cap": nB*cap_sum,
        "charged_B_needed_for_one_expected_event_even_at_unit_efficiency": 1/cap_sum,
        "CLEO_limit_divided_by_upper_sum": inputs["CLEO_combined_charged_90"]/cap_sum,
        "B_rest_momentum_shift_0_to_1MeV_in_keV": dp}


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            compare(a[k], b[k])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    else:
        assert math.isclose(a, b, rel_tol=3e-10, abs_tol=1e-35), (a, b)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    actions = cli.add_mutually_exclusive_group(required=True)
    actions.add_argument("--check", action="store_true")
    actions.add_argument("--write-new", action="store_true")
    args = cli.parse_args()
    out = calculate()
    if args.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        compare(out, read("results.json"))
        print("1599: all three right-handed tags, two-body kaon and B windows, transfer and escape verified")
