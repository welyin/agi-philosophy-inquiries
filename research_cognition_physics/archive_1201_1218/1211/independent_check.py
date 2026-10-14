"""Independent routing/capacity audit; no import of the round-1211 author code."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def permutation(n, d=4):
    words = list(product(range(d), repeat=n))
    index = {w: i for i, w in enumerate(words)}
    memory = n+1
    perm = []
    for w in words:
        first = next((j+1 for j, v in enumerate(w) if v == 0), 0)
        for m in range(memory):
            updated = (m+first) % memory
            v = list(w)
            if updated:
                j = updated-1
                v[0], v[j] = v[j], v[0]
            assert Counter(v) == Counter(w)
            perm.append(index[tuple(v)]*memory + updated)
    assert sorted(perm) == list(range(len(perm)))
    return words, np.array(perm, dtype=int)


def evaluate():
    p = [F(2, 5), F(3, 10), F(1, 5), F(1, 10)]
    cases = []
    for n in range(1, 5):
        words, perm = permutation(n)
        memory = n+1
        prob = [math.prod(p[i] for i in w) for w in words]
        after = F(0)
        for i, w in enumerate(words):
            out = perm[i*memory]
            v = words[out//memory]
            assert (v[0] == 0) == (0 in w)
            if v[0] == 0: after += prob[i]
        expected = 1-(1-p[0])**n
        assert after == expected < 1
        # With no pure record, iid products are scalar on occupation blocks.
        # Blockwise Ky Fan therefore leaves the target population unchanged.
        blocks = {}
        for i, w in enumerate(words):
            key = tuple(w.count(k) for k in range(4))
            blocks.setdefault(key, []).append(i)
        no_record_best = F(0)
        for ids in blocks.values():
            rank_ground = sum(words[i][0] == 0 for i in ids)
            vals = [prob[i] for i in ids]
            assert len(set(vals)) == 1
            no_record_best += sum(sorted(vals, reverse=True)[:rank_ground])
        assert no_record_best == p[0]
        cases.append(dict(n=n, dimension_full=len(perm), ground_success_exact=str(after),
                          no_record_optimum_exact=str(no_record_best),
                          residual_exact=str(1-after), all_basis_permutation=True))
    # Unknown, entangled two-material input with an inaccessible reference.
    n = 2; words, perm = permutation(n); mem = n+1; d = len(words); ref = 3
    rng = np.random.default_rng(1211007)
    pure = rng.normal(size=(d, ref))+1j*rng.normal(size=(d, ref))
    pure /= np.linalg.norm(pure)
    input_full = np.zeros((d*mem, ref), complex)
    input_full[::mem] = pure
    output = np.zeros_like(input_full); output[perm] = input_full
    recovered = output[perm]
    assert np.max(abs(recovered-input_full)) == 0
    no_ground = sum(float(np.sum(abs(pure[i])**2)) for i, w in enumerate(words) if 0 not in w)
    success = sum(float(np.sum(abs(output[i*mem:(i+1)*mem])**2)) for i, w in enumerate(words) if w[0] == 0)
    assert abs(success-(1-no_ground)) < 2e-15
    input_r = pure.conj().T @ pure
    output_r = output.conj().T @ output
    assert np.linalg.norm(input_r-output_r) < 1e-15
    # Joint correlations can supply a guaranteed ground occurrence despite
    # full-rank one-body marginals: cyclic permutations always contain 0.
    correlated = [(0, 1, 2, 3), (1, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2)]
    assert all(0 in w for w in correlated)
    assert all(sorted(w[j] for w in correlated) == list(range(4)) for j in range(4))
    # An independent pure ground stock permits exact SWAP for any unknown SR.
    sr = rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3)); sr /= np.linalg.norm(sr)
    swap_input = np.zeros((4, 4, 3), complex); swap_input[:, 0, :] = sr
    swap_output = swap_input.swapaxes(0, 1)
    assert np.sum(abs(swap_output[0])**2) > 1-1e-15
    assert np.linalg.norm(swap_output[0]-sr) == 0
    # Standard finite-spectrum entropy limit, independent of accessibility.
    spectra = [[0., 1., 2.], [0., 0., 1.]]
    entropy_limits = []
    for es in spectra:
        g = es.count(0.)
        pr = np.exp(-40*np.array(es)); pr /= pr.sum()
        entropy = -sum(x*math.log(x) for x in pr if x)
        assert abs(entropy-math.log(g)) < 1e-13
        entropy_limits.append(dict(ground_degeneracy=g, S_beta40=entropy,
                                   limit=math.log(g)))
    # A closing gap destroys a uniform limit even when every lambda>0 has g=1.
    beta = 10**6; gap = 1/beta
    pr = np.exp(-beta*np.array([0., gap])); pr /= pr.sum()
    closing_entropy = -sum(x*math.log(x) for x in pr)
    assert closing_entropy > .58
    # The actual 100-material energy-window outputs are correlated. Count only
    # their no-ground sector; never replace the pool by iid one-site marginals.
    old2 = json.loads((STAGE/'1202/results.json').read_text(encoding='utf-8'))
    old4 = json.loads((STAGE/'1204/results.json').read_text(encoding='utf-8'))
    prep = old2['preparation']; N = prep['N']; den = prep['root_denominator']
    boxes = prep['root_boxes']; pool_bounds = {}
    for name, data in [('cold', prep), ('hot', old4['new_reference'])]:
        lower, upper = map(F, data['energy_window'])
        total = 0; selected = 0
        for b in range(N+1):
            for c in range(N-b+1):
                dd = N-b-c
                lo = F(b*boxes[1][0]+c*boxes[2][0]+dd*boxes[3][0], den)
                hi = F(b*boxes[1][1]+c*boxes[2][1]+dd*boxes[3][1], den)
                inside = lower <= lo and hi <= upper
                outside = hi < lower or lo > upper
                assert inside or outside
                if inside:
                    total += math.comb(N, b)*math.comb(N-b, c)
                    selected += 1
        probability = F(total, int(data['dimension']))
        assert 0 < probability < 1
        pool_bounds[name] = dict(no_ground_words=str(total), selected_occupation_types=selected,
            pool_dimension=data['dimension'], unavoidable_failure_exact=str(probability),
            approximate_failure=float(probability),
            scope='Only the selected output pool is active; other stocks and preparation records stay passive.')
    # A faithful but non-maximally mixed record is already a purity resource.
    # Exact pair swaps in occupation blocks attain this nonzero increment.
    p0, t = F(2,5), F(3,4)
    mixed_record_gain = (2*t-1)*p0*(1-p0)
    assert mixed_record_gain == F(3,25)
    sources = ['1205/acceptance.json', '1206/acceptance.json',
               '1209/acceptance.json', '1202/results.json', '1204/results.json']
    return dict(status='PASS', all_checks_passed=True,
        source_sha256={x: sha(STAGE/x) for x in sources},
        independent_code_sha256=sha(Path(__file__)),
        complete_permutation_cases=cases,
        unknown_reference_check=dict(no_ground_probability=no_ground,
            output_ground_probability=success, inverse_recovers_full_input=True,
            passive_reference_unchanged=True),
        correlated_full_rank_marginals_but_zero_no_ground_probability=True,
        pure_same_h_stock_swap_transfers_unknown_state=True,
        finite_spectrum_entropy_examples=entropy_limits,
        closing_gap_nonuniform_entropy=closing_entropy,
        native_correlated_output_pool_bounds=pool_bounds,
        mixed_record_resource_gain_exact=str(mixed_record_gain),
        qualifications=dict(native_clock_exact_routing_proved=False,
            arbitrary_energy_battery_covered=False, new_axiom_adopted=False,
            general_third_law_derived=False, stage_complete=False))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true')
    args = parser.parse_args(); result = evaluate(); out = HERE/'independent_results.json'
    if args.write:
        assert not (HERE/'acceptance.json').exists(), 'Do not overwrite accepted results'
        out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert json.loads(out.read_text(encoding='utf-8')) == result
    print(json.dumps(dict(status='PASS', all_checks_passed=True, results_match=True,
        permutation_cases=len(result['complete_permutation_cases']))))


if __name__ == '__main__': main()
