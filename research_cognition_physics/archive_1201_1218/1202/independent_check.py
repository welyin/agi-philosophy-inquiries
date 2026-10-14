"""1202 independent candidate: exact shell counting and finite uniform unranking.
Writes only independent_results.json with --write (initial creation only).
Default recomputes and compares saved JSON, with no writes.
The conditional-SWAP controller is a declared candidate contract, not a QQ implementation.
"""
from pathlib import Path
from bisect import bisect_right
from fractions import Fraction
import argparse
import hashlib
import json
import math
import random
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OUT = HERE / "independent_results.json"


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-8, abs_tol=2e-12), (a, b)
    else:
        assert a == b, (a, b)


def run():
    N, root_den, attempts = 100, 10**18, 10
    def poly(x):
        return (32000*x**3 - 31200*x*x*root_den -
                2095*x*root_den**2 - 17*root_den**3)
    def bracket(a, b):
        a, b = a*(root_den//10**6), b*(root_den//10**6)
        sa = poly(a)
        assert sa*poly(b) < 0
        while b-a > 1:
            c = (a+b)//2
            if poly(c)*sa > 0:
                a, sa = c, poly(c)
            else:
                b = c
        return a, b
    boxes = [bracket(-54072,-54071), bracket(-9461,-9460),
             (root_den,root_den), bracket(1038532,1038533)]
    energies = np.array([(a+b)/(2*root_den) for a,b in boxes])
    factorial = [math.factorial(i) for i in range(N+1)]
    types, weights, prefix = [], [], [0]
    one, two = [0]*4, [[0]*4 for _ in range(4)]
    boundary_margin = 10**100
    for a in range(N+1):
        for b in range(N-a+1):
            for c in range(N-a-b+1):
                n = (a,b,c,N-a-b-c)
                low = sum(n[j]*boxes[j][0] for j in range(4))
                high = sum(n[j]*boxes[j][1] for j in range(4))
                inside = 2*low >= 39*root_den and 2*high <= 41*root_den
                outside = 2*high < 39*root_den or 2*low > 41*root_den
                assert inside or outside
                boundary_margin = min(boundary_margin, *[
                    abs(2*x-y*root_den) for x in (low,high) for y in (39,41)])
                if not inside:
                    continue
                w = factorial[N] // math.prod(factorial[x] for x in n)
                types.append(n)
                weights.append(w)
                prefix.append(prefix[-1]+w)
                for i in range(4):
                    one[i] += w*n[i]
                    for j in range(4):
                        two[i][j] += w*n[i]*(n[j]-(i==j))
    dimension = prefix[-1]
    assert dimension == 9399690799043690121938015947926455391852712761214920
    assert len(types) == 1728
    type_index = {n: j for j,n in enumerate(types)}

    def unrank(index):
        if not isinstance(index, int) or not 0 <= index < dimension:
            raise ValueError("Index outside accepted shell")
        block = bisect_right(prefix, index)-1
        residual = index-prefix[block]
        counts = list(types[block])
        ways, remaining, word = weights[block], N, []
        while remaining:
            for letter in range(4):
                numerator = ways*counts[letter]
                assert numerator % remaining == 0
                block_size = numerator//remaining
                if residual >= block_size:
                    residual -= block_size
                    continue
                word.append(letter)
                counts[letter] -= 1
                ways, remaining = block_size, remaining-1
                break
            else:
                raise AssertionError("No permutation block found")
        assert residual == 0 and ways == 1 and counts == [0]*4
        return tuple(word)

    def rank(word):
        if len(word) != N or any(x not in (0,1,2,3) for x in word):
            raise ValueError("Invalid level string")
        counts = [word.count(i) for i in range(4)]
        block = type_index.get(tuple(counts))
        if block is None:
            raise ValueError("String outside energy shell")
        ways, remaining, local = weights[block], N, 0
        for letter in word:
            for earlier in range(letter):
                numerator = ways*counts[earlier]
                assert numerator % remaining == 0
                local += numerator//remaining
            numerator = ways*counts[letter]
            assert numerator % remaining == 0
            ways = numerator//remaining
            counts[letter] -= 1
            remaining -= 1
        assert ways == 1
        return prefix[block]+local

    # Every type's first/last string, plus deterministic pseudorandom test indices.
    # This PRNG is only for software tests, not the proposed physical seed source.
    boundaries = {x for j in range(len(types)) for x in (prefix[j],prefix[j+1]-1)}
    rng = random.Random(1202)
    random_indices = [rng.randrange(dimension) for _ in range(256)]
    indices = sorted(boundaries | set(random_indices) | {0,1,dimension-2,dimension-1})
    max_swap_energy_residual = 0.0
    for index in indices:
        word = unrank(index)
        assert rank(word) == index
        counts = tuple(word.count(i) for i in range(4))
        assert counts in type_index
        assert 2*sum(counts[i]*boxes[i][0] for i in range(4)) >= 39*root_den
        assert 2*sum(counts[i]*boxes[i][1] for i in range(4)) <= 41*root_den
        # Per-target inventory: (e0,e1,ec,e3); B target initially |e0>.
        # The later measurement probe A, not this preparation target, uses native |c>.
        for desired in set(word):
            stock = [0,1,2,3]
            target = 0
            target, stock[desired] = stock[desired], target
            assert target == desired
            assert stock[desired] == 0
            before = energies[0]+sum(energies)
            after = energies[target]+sum(energies[j] for j in stock)
            max_swap_energy_residual = max(max_swap_energy_residual,abs(before-after))
    for bad in (-1,dimension,dimension+1):
        try:
            unrank(bad)
            raise AssertionError("Invalid rank was accepted")
        except ValueError:
            pass

    seed_bits = (dimension-1).bit_length()
    seed_space = 1 << seed_bits
    assert seed_bits == 173
    rejected = seed_space-dimension
    per_seed = Fraction(dimension,seed_space)
    failure = Fraction(rejected,seed_space)**attempts
    success = 1-failure
    assert failure < Fraction(1,4_000_000)
    # For each accepted x, first-valid sampling has this same unconditional weight.
    weight_per_x = sum((1-per_seed)**j/Fraction(seed_space,1) for j in range(attempts))
    assert dimension*weight_per_x == success
    assert weight_per_x/success == Fraction(1,dimension)

    def first_valid(seeds):
        if len(seeds) != attempts or any(not 0 <= x < seed_space for x in seeds):
            raise ValueError("Need ten valid 173-bit seed values")
        for j,x in enumerate(seeds):
            if x < dimension:
                return dict(success=True, selected_position=j, selected_index=x,
                            preserved_seed_values=list(seeds))
        return dict(success=False, selected_position=None, selected_index=None,
                    preserved_seed_values=list(seeds))
    samplers = [
        [dimension]*attempts,
        [0]+[dimension]*(attempts-1),
        [dimension]*(attempts-1)+[dimension-1],
        [dimension,dimension-1,0]+[dimension]*(attempts-3)]
    sampler_results = [first_valid(s) for s in samplers]
    assert sampler_results[0]["success"] is False
    assert sampler_results[1]["selected_index"] == 0
    assert sampler_results[2]["selected_position"] == 9
    assert sampler_results[3]["selected_position"] == 1

    r = np.array([x/(N*dimension) for x in one])
    pair = np.array([[x/(N*(N-1)*dimension) for x in row] for row in two])
    assert max(abs(pair.sum(axis=0)-r)) < 1e-14
    native = STAGE/"1201/independent_results_1201.json"
    old = json.loads(native.read_text(encoding="utf-8"))
    M = np.array(old["population_readout"]["matrix_rows_A_output_columns_B_population"])
    Mhold = np.array(old["predeclared_holdout"]["matrix"])
    joint = M@pair@Mhold.T
    assert abs(joint.sum()-1) < 1e-12
    old_coeff = old["predeclared_holdout"]["energy_coefficients_by_initial_B_population"]
    hold_energy = {k:float(np.array(v)@r) for k,v in old_coeff.items()}
    pair_h = np.diag((energies[:,None]+energies[None,:]).ravel())
    swap = np.zeros((16,16))
    for a in range(4):
        for b in range(4):
            swap[4*b+a,4*a+b] = 1
    assert np.linalg.norm(swap.T@pair_h@swap-pair_h,2) == 0
    accepted_energy = float(N*energies@r)
    initial_inventory_energy = float(N*sum(energies))
    initial_target_energy = float(N*energies[0])
    final_inventory_energy = initial_inventory_energy+initial_target_energy-accepted_energy
    table_payload = json.dumps(list(zip(types,weights)),separators=(",",":")).encode()
    examples = []
    for index in (0,1,dimension//2,dimension-2,dimension-1):
        word = unrank(index)
        examples.append(dict(index=str(index),levels="".join(map(str,word)),
                             counts=[word.count(i) for i in range(4)]))
    return dict(
        status="1202 independent candidate audit, not physical controller implementation",
        all_checks_passed=True,
        fixed_problem=dict(N=N,window=[19.5,20.5],root_denominator=root_den,
                           root_boxes=boxes,boundary_margin_lower=boundary_margin/(2*root_den)),
        exact_counting=dict(shell_dimension=str(dimension),admissible_types=len(types),
                            all_types=math.comb(N+3,3),
                            type_order="lexicographic count tuple, then lexicographic level string",
                            type_weight_table_sha256=hashlib.sha256(table_payload).hexdigest()),
        rank_unrank=dict(all_type_boundary_indices=len(boundaries),
                         deterministic_random_indices=len(random_indices),
                         total_distinct_indices_checked=len(indices),
                         invalid_indices_rejected=True,examples=examples),
        seed_sampler=dict(bits_per_seed=seed_bits,preallocated_seeds=attempts,
                          total_seed_bits=seed_bits*attempts,
                          one_seed_acceptance=float(per_seed),
                          all_failure_probability=float(failure),
                          success_probability=float(success),
                          failure_bound="less than 1/4000000, checked as exact rational",
                          conditional_uniformity_exact_fraction_check=True,
                          preserved_full_seed_records=True,
                          sampler_test_results=sampler_results),
        source_marginals=dict(one=r.tolist(),two=pair.tolist(),
                              mean_internal_energy=float(energies@r),
                              source_entropy_nats=math.log(dimension),
                              two_vs_product_TV=float(abs(pair-np.outer(r,r)).sum()/2)),
        conditional_inventory_contract=dict(
            initial_target_level=0,inventory_levels_per_target=[0,1,2,3],
            target_materials=N,inventory_materials=4*N,total_materials=5*N,
            maximum_selected_swaps=N,
            initial_target_internal_energy=initial_target_energy,
            initial_inventory_internal_energy=initial_inventory_energy,
            accepted_target_internal_energy=accepted_energy,
            accepted_inventory_internal_energy=final_inventory_energy,
            total_material_internal_energy=initial_target_energy+initial_inventory_energy,
            single_swap_bare_energy_commutator_norm=0.0,
            maximum_label_swap_energy_residual=float(max_swap_energy_residual),
            readonly_table_storage_bits_upper=173*(len(types)+1)+28*len(types),
            seed_preparation_and_controller_energy_not_derived=True,
            failed_branch_no_target_use=True),
        original_1201_connection=dict(
            source_result_sha256=hashlib.sha256(native.read_bytes()).hexdigest(),
            actual_600_output=(M@r).tolist(),actual_300_output=(Mhold@r).tolist(),
            joint_600_300=joint.tolist(),
            joint_vs_independent_same_marginal_TV=float(
                abs(joint-np.outer(M@r,Mhold@r)).sum()/2),
            holdout_300_energy=hold_energy),
        scope=dict(
            unknown_reference_assumed_gibbs=False,
            rare_4_power_N_physical_postselection_required=False,
            source_independent_of_preserved_records=False,
            same_pool_materials_independent=False,
            uniform_physical_seed_source_constructed=False,
            native_QQ_realizes_controlled_swap=False,
            internal_controller_constructed=False,
            final_1202_admission=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write:
        if OUT.exists():
            raise FileExistsError("Refusing to overwrite independent results")
        OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        compare(result,json.loads(OUT.read_text(encoding="utf-8")))
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
