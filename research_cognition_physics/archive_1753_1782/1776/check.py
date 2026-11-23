"""Gap-free charged-current spectral bounds for a net chiral family index."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE / path).read_bytes())


def window_overlap(delta_u, delta_d, cut_u, cut_d):
    su, sd = min(1., delta_u / cut_u), min(1., delta_d / cut_d)
    return max(0., math.sqrt((1-su*su)*(1-sd*sd))-su*sd)


def pair_rate(mass_u, mass_d, energy, left, right):
    """Contribution per colour, normalized to a massless left lepton current."""
    if mass_u + mass_d >= energy:
        return 0.
    x, y = (mass_u / energy)**2, (mass_d / energy)**2
    root = math.sqrt((1-((mass_u+mass_d)/energy)**2) *
                     (1-((mass_u-mass_d)/energy)**2))
    a = 1-(x+y)/2-(x-y)**2/2
    cross = 6*math.sqrt(x*y)*(left*right.conjugate()).real
    return root * (a*(abs(left)**2+abs(right)**2) + cross)


def optimized_floor(energy, delta, extras, colours):
    c = 8*(delta/energy)**2
    if c >= 1:
        return dict(nonzero_bound=False, relative_rate=0.)
    z = 8*c/(math.sqrt(c*c+48*c)+c)
    cut = energy*math.sqrt(z)/2
    f = window_overlap(delta, delta, cut, cut)
    kinematic_factor = (1-z)**1.5
    return dict(nonzero_bound=True, cut_GeV=cut, weak_strength=extras*f*f,
                relative_rate=colours*extras*f*f*kinematic_factor)


def coherent_texture(n, entry, bare, v):
    """One extra index; all entries of the direct Higgs block equal entry."""
    one = np.ones(n)/math.sqrt(n)
    # The remaining columns of Q are an orthonormal complement of one.
    q, _ = np.linalg.qr(np.column_stack([one, np.eye(n)[:, :-1]]))
    hard = bare*q[:, 1:]
    direct = np.full((n, n), v*entry/math.sqrt(2))
    m = np.block([[hard, direct], [np.zeros((n-1, n-1)), hard.T]])
    masses = np.sort(np.linalg.svd(m, compute_uv=False))
    expected = np.sort([n*v*entry/math.sqrt(2)]+[bare]*(2*n-2))
    np.testing.assert_allclose(masses, expected, rtol=1e-11)
    return dict(number_new_mass_lines=len(masses), lightest_GeV=masses[0],
                largest_Yukawa_entry=entry,
                operator_Yukawa=np.linalg.norm(direct, 2)*math.sqrt(2)/v,
                vector_partner_mass_GeV=bare)


def qr(rng, n):
    z = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    return np.linalg.qr(z)[0]


def rectangular(rng, m, n, singular):
    return (qr(rng, m)[:, :len(singular)]*singular) @ qr(rng, n)[:, :len(singular)].conj().T


def modal_sample(rng, n, m, maximum):
    k = min(n, m)
    vals = rng.uniform(.015, maximum, size=k)
    return rectangular(rng, n, m, vals)


def exercise(rng, i):
    index = 4 + i%2
    doublet_right = 1+i%3
    doublet_left = doublet_right + index
    hard_values = rng.choice([0., .03, .85, 1.8, 6., 32.], size=doublet_right)
    rotation = qr(rng, doublet_left)
    hard = (rotation[:, :doublet_right]*hard_values) @ qr(rng, doublet_right).conj().T
    # Choose any protected index-dimensional common left-doublet plane.
    k = rotation[:, doublet_right:]
    data = []
    for single_left in (1+(i//2)%3, 1+(i//5)%2):
        single_right = index+single_left
        hard_s = modal_sample(rng, single_left, single_right, rng.choice([.23, 3., 50.]))
        forward = modal_sample(rng, doublet_left, single_right, 2.5)
        reverse = modal_sample(rng, single_left, doublet_right, 5.)
        mat = np.block([[hard, forward], [reverse, hard_s]])
        u, descending, vh = np.linalg.svd(mat)
        masses = descending[::-1]
        u, v = u[:, ::-1], vh.conj().T[:, ::-1]
        restricted = k.conj().T @ forward
        delta = np.linalg.norm(restricted, 2)
        assert masses[index-1] <= delta + 2e-11
        assert np.sum(masses[:index]**2) <= np.linalg.norm(restricted, "fro")**2 + 1e-10
        cut = rng.uniform(1.30, 3.8)*delta
        mask = masses <= cut
        projected = u[:doublet_left, ~mask].conj().T @ k
        leakage = np.linalg.norm(projected, 2) if projected.size else 0.
        assert leakage <= delta/cut + 2e-11
        data.append((masses[mask], u[:doublet_left, mask], v[:doublet_right, mask], delta, cut))
    mu, ul, ur, du, eu = data[0]
    md, dl, dr, dd, ed = data[1]
    left, right = ul.conj().T @ dl, ur.conj().T @ dr
    f = window_overlap(du, dd, eu, ed)
    s = np.linalg.svd(left, compute_uv=False)
    assert s[index-1] >= f-2e-11
    tail = np.ones(left.shape, dtype=bool)
    tail[:3, :3] = False
    strength = np.sum(abs(left[tail])**2)
    assert strength+1e-10 >= (index-3)*f**2
    energy = 2.5*max(eu, ed)
    rate = sum(pair_rate(mu[iu], md[jd], energy, left[iu, jd], right[iu, jd])
               for iu, jd in zip(*np.nonzero(tail)))
    kinematic = (1-4*(max(eu, ed)/energy)**2)**1.5
    assert rate+1e-10 >= kinematic*strength
    return bool(min(hard_values) < min(du, dd)/2), float(f)


def calculate():
    cfg = load("inputs.json")
    v = load(cfg["old_width_and_vev"])["v_GeV"]
    colours = load(cfg["old_loop_norm"])["Nc"]
    gate = load(cfg["old_gate"])["initial_family"]["mc"]
    families, known = cfg["family_index"], cfg["known_families"]
    y = 4*math.pi/math.sqrt(colours)
    delta = v*y/math.sqrt(2)
    control = cfg["texture_control"]
    chiral_check = coherent_texture(control["extra_uniform_size"], control["per_entry_y"],
                                    control["vector_mass_GeV"], v)
    old = load(cfg["old_fourth"])["family"]["masses_chiral_UDEN"]
    mu, md = old[:2]
    energy_list = sorted(set(cfg["readout_virtualities_GeV"] + [gate]))
    floors = {}
    for budget in cfg["weak_budget_examples"]+[y]:
        floors[str(round(budget, 8))] = {
            str(q):optimized_floor(q, v*budget/math.sqrt(2), families-known, colours)
            for q in energy_list}
    thresholds = {str(mass):dict(minimum_Y=mass*math.sqrt(2)/v,
                                 loop_strength=colours*(mass*math.sqrt(2)/v)**2/(16*math.pi**2))
                  for mass in [gate, cfg["additional_comparison_GeV"]]}
    old_rates = {str(q):colours*pair_rate(mu, md, q, 1., 0.)
                 for q in (300., 400., 500., 1000.)}
    generator = np.random.default_rng(cfg["random_seed"])
    cases = [exercise(generator, i) for i in range(cfg["random_cases"])]
    return dict(upper_on_pth_ascending_mass_GeV=delta,
                virtuality_for_first_positive_guarantee_GeV=math.sqrt(8)*delta,
                required_projected_Y_if_only_three_below=thresholds,
                optimized_extra_pair_rate_over_one_lepton=floors,
                historical_chiral_pair=dict(masses_GeV=[mu, md], rates=old_rates),
                entrywise_budget_counterexample=chiral_check,
                numerical_checks=dict(samples=len(cases),
                    common_bare_doublet_gap_below_both_half_projected_norms=sum(a for a,b in cases),
                    smallest_test_lower_singular_bound=min(b for a,b in cases)))


def serialize(value):
    if isinstance(value, dict):
        return {k:serialize(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(v) for v in value]
    if isinstance(value, (int, bool)):
        return value
    return float(f"{value:.10g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = serialize(calculate())
    if args.check:
        assert result == load("results.json")
        print("1776 OK: protected mass moments, gap-free charged-current singular values,"
              " finite timelike spectral rates and coherent counterexample")
    else:
        (HERE/"results.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")

