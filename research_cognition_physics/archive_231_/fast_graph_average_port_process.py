"""Round 494: fast graph averaging for the retained port instrument.

Baseline 492. Same active D/G throughout; a real initial port pinch and fresh
isolated readers are extra inputs. Proof is analytic with exact finite algebra
and rational constants. Floating matrix powers are labelled diagnostics.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib, io, itertools, json, math, platform, unittest
from pathlib import Path
import numpy as np
import retained_port_instrument_mixing as prior
import branching_tree_distance_audit as trees_old

TARGET = Path(__file__).with_name("fast_graph_average_port_process_results.json")
OBS = {}
LEAVES = (0, 3, 4, 5)
# All exact numpy products below have dimension <=384 and at most four integer
# factors with entries <=32. Check BEFORE products; larger rational arithmetic
# uses Python int/Fraction. This exceeds all polynomial/projector operations.
assert 384**4 * 32**4 < 2**63

def short(x):
    return float(f"{float(x):.12g}")

@lru_cache(None)
def model():
    trees, full_h, f, adj, embedding, h, ports = prior.model()
    n = np.zeros((6, 6, 6, 6), dtype=np.int64)
    for a in range(6):
        for b in range(6):
            n[b, a] = np.diag([g[b, a] for g in adj])
    degree = np.array([1, 3, 3, 1, 1, 1], dtype=np.int64)
    h0 = np.kron(np.eye(6, dtype=np.int64), f)
    h0 += np.kron(np.diag(5-degree), np.eye(6, dtype=np.int64))
    return trees, full_h, f, adj, embedding, h, ports, n, degree, h0

def block_embed(blocks):
    ans = np.zeros((36, 36), dtype=np.result_type(*blocks))
    for a, x in enumerate(blocks):
        ans[6*a:6*a+6, 6*a:6*a+6] = x
    return ans

def blocks_of(x):
    return np.array([x[6*a:6*a+6, 6*a:6*a+6] for a in range(6)])

def b_generator(blocks):
    n, degree = model()[7:9]
    return np.array([sum((n[b,a]@blocks[a]@n[a,b]
                           for a in range(6) if a != b),
                          np.zeros_like(blocks[0]))-degree[b]*blocks[b]
                     for b in range(6)])

@lru_cache(None)
def exact_graph():
    _, _, f, adj, *_ = model()
    one = np.ones((6,6), dtype=np.int64)
    ident = np.eye(6, dtype=np.int64)
    complement = one-ident-f
    # All projectors are represented by integer numerators with denominator 24.
    numerators = [f@(f+2*ident),
                  -3*(f-4*ident)@(f+2*ident),
                  2*f@(f-4*ident)]
    return complement, numerators

def averaged_laplacian():
    l = np.diag(model()[8]).astype(np.int64)*2
    l[1,2] = l[2,1] = -2
    for a in LEAVES:
        l[a,1] = l[1,a] = l[a,2] = l[2,a] = -1
    return l  # actual Laplacian is l/2.

def classical_transition(tau):
    one = np.ones(6)
    leaf = np.array([int(a in LEAVES) for a in range(6)])
    center = np.array([0,1,-1,0,0,0])
    middle = np.array([-1,2,2,-1,-1,-1])
    p0 = np.outer(one,one)/6
    pl = np.diag(leaf)-np.outer(leaf,leaf)/4
    pm = np.outer(middle,middle)/12
    pc = np.outer(center,center)/2
    return p0+math.exp(-tau)*pl+math.exp(-3*tau)*pm+math.exp(-4*tau)*pc

def trace_graph(blocks):
    return np.einsum("aiirs->ars", blocks)

def trace_norm_cq(blocks):
    return float(sum(np.abs(np.linalg.eigvalsh((x+x.conj().T)/2)).sum()
                     for x in blocks.reshape(-1,2,2)))

@lru_cache(None)
def numeric_process():
    # One fixed diagnostic. No eigenvalue estimate supplies the theorem.
    delta = F(1,1024)
    steps = 1 << 18
    tau = steps*delta*delta
    h = model()[5]
    eig, v = np.linalg.eigh(h)
    u = (v*np.exp(-1j*float(delta)*eig))@v.conj().T
    channel = np.zeros((216,216), dtype=complex)
    for b in range(6):
        for a in range(6):
            uba = u[6*b:6*b+6,6*a:6*a+6]
            channel[36*b:36*b+36,36*a:36*a+36] = np.einsum(
                "ik,jl->ijkl",uba,uba.conj()).reshape(36,36)
    power = np.linalg.matrix_power(channel, steps)
    return delta, steps, tau, u, channel, power

def diagnostic_source():
    # Port coherence AND graph-reference entanglement. Pinch is really applied.
    psi = np.zeros((6,6,2), dtype=complex)
    psi[0,0,0] = 1
    psi[0,1,1] = 1j
    psi[3,2,0] = 1
    psi[3,4,1] = -1
    psi /= 2
    density = np.outer(psi.ravel(), psi.ravel().conj())
    blocks = np.array([density[12*a:12*a+12,12*a:12*a+12]
                       .reshape(6,2,6,2).transpose(0,2,1,3)
                       for a in range(6)])
    return psi, density, blocks

def apply_super(s, blocks):
    return (s@blocks.reshape(216,4)).reshape(6,6,6,2,2)

def rational_bounds():
    d = F(1,64)
    x = 8+6*d
    original = F(18**3,6)/(1-18*d/4)
    surrogate = x**3/(6*(1-x*d/4))+48+18*d
    assert original+surrogate < 1500
    return dict(delta_max=str(d), original_third_order_coefficient=str(original),
                surrogate_third_order_coefficient=str(surrogate),
                total_third_order_coefficient=str(original+surrogate),
                adopted_coefficient=1500, generator_A_norm_bound=8,
                generator_B_norm_bound=6, periodic_primitive_bound=48,
                averaging_bound="48*delta*(1+12*tau)",
                final_ordinary_diamond_bound="delta*(48+2076*tau)")

@lru_cache(None)
def larger_graph():
    words = sorted(set(itertools.permutations((0,0,1,1,2,2))))
    trees = sorted({trees_old.prufer_tree(8,w) for w in words},
                   key=lambda t:sorted(t))
    lookup = {t:j for j,t in enumerate(trees)}
    f = np.zeros((90,90), dtype=np.int64)
    for j,t in enumerate(trees):
        for z,count in trees_old.tree_flips(t,8).items():
            f[lookup[z],j] = count
    permutation = (1,0,2,4,3,6,5,7)
    r = np.zeros_like(f)
    for j,t in enumerate(trees):
        image = frozenset(trees_old.edge(permutation[a],permutation[b])
                          for a,b in t)
        r[lookup[image],j] = 1
    q = np.array([int((0,3) in t) for t in trees], dtype=np.int64)
    return trees, f, r, q

class Audit(unittest.TestCase):
    def test_01_original_H_and_exact_second_order_generator(self):
        _, full_h, f, _, v, h, ports, n, degree, h0 = model()
        self.assertTrue(np.array_equal(full_h@v,v@h))
        self.assertTrue(np.array_equal(h0+sum(
            (np.kron(np.eye(6,dtype=np.int64)[b:b+1].T
                     @np.eye(6,dtype=np.int64)[a:a+1], n[b,a])
             for a in range(6) for b in range(6) if a!=b),
            np.zeros((36,36),dtype=np.int64)),h))
        out_rate = np.zeros((36,36),dtype=np.int64)
        for a in range(6):
            self.assertTrue(np.array_equal(sum((n[b,a] for b in range(6)),
                    np.zeros((6,6),dtype=np.int64)),degree[a]*np.eye(6,dtype=np.int64)))
            out_rate[6*a:6*a+6,6*a:6*a+6] = degree[a]*np.eye(6,dtype=np.int64)
        checks = 0
        for a in range(6):
            for i in range(6):
                for j in range(6):
                    blocks = np.zeros((6,6,6),dtype=np.int64)
                    blocks[a,i,j] = 1
                    x = block_embed(blocks)
                    comm = h@x-x@h
                    actual_first = blocks_of(comm)
                    expected_first = np.array([f@y-y@f for y in blocks])
                    self.assertTrue(np.array_equal(actual_first,expected_first))
                    twice_actual = -blocks_of(h@comm-comm@h)
                    expected = np.array([-(f@(f@y-y@f)-(f@y-y@f)@f)
                                         for y in blocks])+2*b_generator(blocks)
                    self.assertTrue(np.array_equal(twice_actual,expected))
                    checks += 1
        OBS["exact_instrument_expansion"] = dict(
            full_H_dimension=384, active_dimension=36,
            complex_matrix_unit_checks=checks, Kraus_out_rate_diagonal=out_rate.diagonal().tolist(),
            graph_A="minus_i_commutator_F", B_is_GKSL=True,
            original_instrument_requires_external_to_H_measurement_contract=True)

    def test_02_exact_energy_blocks_and_arbitrary_graph_closure(self):
        f = model()[2]
        c, nums = exact_graph()
        ident = np.eye(6,dtype=np.int64)
        self.assertTrue(np.array_equal(c@c,ident))
        self.assertTrue(np.array_equal(f,np.ones((6,6),dtype=np.int64)-ident-c))
        self.assertTrue(np.array_equal(sum(nums),24*ident))
        for i,(r,e) in enumerate(zip(nums,(4,0,-2))):
            self.assertTrue(np.array_equal(f@r,e*r))
            self.assertTrue(np.array_equal(r@r,24*r))
            for j,s in enumerate(nums):
                if i!=j:self.assertFalse(np.any(r@s))
        n = model()[7]
        for a in LEAVES:
            q = n[1,a]
            self.assertTrue(np.array_equal(c@q@c,ident-q))
            for r in nums:
                self.assertTrue(np.array_equal(2*r@q@r,24*r))
        l2 = averaged_laplacian()
        # Adjoint closure on EVERY one of the six port population effects.
        for target in range(6):
            for a in range(6):
                operator = sum(((int(b==target)-int(a==target))*n[b,a]
                                 for b in range(6) if b!=a),
                                np.zeros((6,6),dtype=np.int64))
                avg_num = sum((r@operator@r for r in nums),
                              np.zeros((6,6),dtype=np.int64))
                self.assertTrue(np.array_equal(2*avg_num,
                           -576*l2[target,a]*ident))
        self.assertFalse(np.any(l2.sum(axis=0)))
        self.assertTrue(np.array_equal(l2,l2.T))
        # Finite role spectrum is a consequence of preselected endpoints, NOT
        # a new necessity claim and NOT evidence of a smooth 3D manifold.
        ranks = []
        for eigen in (0,2,6,8):
            ranks.append(6-prior.rank_mod(l2-eigen*ident))
        self.assertEqual(ranks,[1,3,1,1])
        OBS["exact_fast_average"] = dict(
            energies=[4,0,-2],energy_ranks=[int(np.trace(x)//24) for x in nums],
            fast_period="pi", complement_signs=[1,-1,1],
            averaged_leaf_center_effect="I/2",
            scalar_effect_closure_verified_for_all_six_ports=True,
            independent_graph_source_not_required=True,
            arbitrary_passive_reference_supported=True,
            twice_classical_laplacian=l2.tolist(),
            finite_role_eigenvalues=[0,1,3,4],finite_role_multiplicities=ranks,
            role_multiplicity_not_spatial_dimension=True)

    def test_03_rational_discrete_and_periodic_error_budget(self):
        payload = rational_bounds()
        delta = F(1,1<<22)
        n = 1<<42
        tau = n*delta*delta
        self.assertEqual(tau,F(1,4))
        single = delta*(48+2076*tau)
        readers = 2*n+1  # initial actual pinch plus both full microstep blocks.
        gamma = F(1,4000*readers)
        source = F(1,4000)
        total = 2*single+readers*gamma+source
        self.assertLess(total,F(1,1000))
        payload["finite_two_coarse_interval_contract"] = dict(
            delta=str(delta),microsteps_per_interval=n,
            tau_per_interval=str(tau),coarse_intervals=2,
            analytic_error_per_interval=str(single),analytic_error_sum=str(2*single),
            instrument_steps_and_six_state_records=readers,
            record_qubit_encoding_capacity=3*readers,
            ideal_wait_total=str(2*n*delta),
            whole_instrument_diamond_error_per_step=str(gamma),
            optional_source_trace_error=str(source),
            total_ordinary_trace_or_diamond_error=str(total),
            certified_upper_bound="1/1000",actually_run_microstep_count=0,
            active_networks=1,new_graph_sources=0,new_data_batches=0,
            detection_hardware_clocks_and_isolation_not_derived=True)
        OBS["uniform_error_certificate"] = payload

    def test_04_actual_CQ_process_with_graph_reference_and_initial_pinch(self):
        delta,n,tau,u,s,sn = numeric_process()
        psi,density,x = diagnostic_source()
        # Physical t=0 port instrument includes six records. Its nonselective
        # state is exactly x, and graph-reference traces are not reset.
        record_isometry = np.zeros((6,72),dtype=complex)
        for a in range(6):
            record_isometry[a,12*a:12*a+12] = psi[a].ravel()
        self.assertAlmostEqual(float(np.linalg.norm(record_isometry)**2),1)
        post = np.zeros_like(density)
        for a in range(6):
            post[12*a:12*a+12,12*a:12*a+12] = density[
                12*a:12*a+12,12*a:12*a+12]
        self.assertLess(np.linalg.norm(post-density),1)
        self.assertGreater(np.linalg.norm(post-density),.1)
        self.assertLess(np.linalg.norm(trace_graph(x).sum(axis=0)-np.eye(2)/2),1e-14)
        # Verify one actual U / port instrument with R left untouched.
        ur = np.kron(u,np.eye(2))
        evolved = ur@post@ur.conj().T
        direct = np.array([evolved[12*a:12*a+12,12*a:12*a+12]
                  .reshape(6,2,6,2).transpose(0,2,1,3) for a in range(6)])
        self.assertLess(np.max(np.abs(direct-apply_super(s,x))),1e-13)
        actual = trace_graph(apply_super(sn,x))
        expected = np.einsum("ba,ars->brs",classical_transition(float(tau)),trace_graph(x))
        error = trace_norm_cq(actual-expected)
        self.assertLess(error,float(delta*(48+2076*tau)))
        self.assertLess(error,.02)
        self.assertAlmostEqual(float(np.trace(actual.sum(axis=0)).real),1,places=7)
        OBS["fixed_point_numerical_diagnostic"] = dict(
            delta=str(delta),actual_microsteps_per_interval=n,tau=str(tau),
            source_has_initial_port_coherence=True,real_initial_pinch=True,
            source_graph_reference_entanglement=True,
            one_step_full_instrument_matrix_error=short(np.max(np.abs(direct-apply_super(s,x)))),
            endpoint_joint_port_R_trace_error=short(error),
            analytic_ordinary_bound=str(delta*(48+2076*tau)),
            fast_generator_reverse_used_only_in_proof=True,
            floating_power_is_diagnostic_not_general_certificate=True)

    def test_05_two_coarse_records_keep_the_same_graph(self):
        delta,n,tau,_,_,sn = numeric_process()
        _,_,x = diagnostic_source()
        first = apply_super(sn,x)
        actual = np.zeros((6,6,2,2),dtype=complex)
        for b in range(6):
            branch = np.zeros_like(first)
            branch[b] = first[b]  # retain conditional G/R, DO NOT prepare G.
            actual[b] = trace_graph(apply_super(sn,branch))
        p = classical_transition(float(tau))
        initial = trace_graph(x)
        ideal = np.zeros_like(actual)
        for b in range(6):
            reduced = sum((p[b,a]*initial[a] for a in range(6)),
                           np.zeros((2,2),dtype=complex))
            for c in range(6):ideal[b,c] = p[c,b]*reduced
        error = trace_norm_cq(actual-ideal)
        self.assertLess(error,float(2*delta*(48+2076*tau)))
        self.assertLess(error,.03)
        self.assertAlmostEqual(float(np.trace(actual.sum(axis=(0,1))).real),1,places=7)
        OBS["coarse_history_diagnostic"] = dict(
            recorded_coarse_times=[str(tau),str(2*tau)],
            joint_two_records_and_initial_R_trace_error=short(error),
            analytic_uniform_error_sum=str(2*delta*(48+2076*tau)),
            conditional_graph_retained_between_intervals=True,
            graph_traced_only_for_reported_records=True,
            full_microstep_record_history_not_claimed_Markov=True,
            full_DG_plus_records_not_claimed_classical=True,
            no_feedback_from_ignored_records=True,
            per_normalized_rare_history_error_not_claimed=True)

    def test_06_larger_recursive_graph_breaks_uniform_effect_closure(self):
        trees,f,r,q = larger_graph()
        ident = np.eye(90,dtype=np.int64)
        self.assertEqual(len(trees),90)
        self.assertTrue(np.array_equal(f,f.T))
        self.assertTrue(np.array_equal(f.sum(axis=0),np.full(90,8)))
        self.assertTrue(np.array_equal(r,r.T))
        self.assertTrue(np.array_equal(r@r,ident))
        self.assertTrue(np.array_equal(r@f,f@r))
        self.assertEqual(int(q.sum()),30)
        self.assertEqual(int(np.trace(r)),4)
        self.assertEqual(int(q@np.diag(r)),2)
        # I +/- R are twice orthogonal projections, hence positive. Both
        # normalized states commute with F, so no fast averaging changes these
        # particular edge expectations. No F spectrum or float is needed.
        for sign in (1,-1):
            a = ident+sign*r
            self.assertTrue(np.array_equal(a@a,2*a))
        plus = F(int(q.sum()+q@np.diag(r)),90+int(np.trace(r)))
        minus = F(int(q.sum()-q@np.diag(r)),90-int(np.trace(r)))
        self.assertEqual(plus,F(16,47))
        self.assertEqual(minus,F(14,43))
        self.assertEqual(plus-minus,F(30,2021))
        self.assertNotEqual(plus,F(1,3))
        OBS["recursive_size_boundary"] = dict(
            internal_vertices=[0,1,2],leaves=[3,4,5,6,7],tree_count=90,
            NNI_graph_degree=8,permutation=[1,0,2,4,3,6,5,7],
            trace_Q=30,trace_R=4,trace_QR=2,
            stationary_graph_sources=["(I+R)/94","(I-R)/86"],
            edge_03_expectations=[str(plus),str(minus)],
            difference=str(plus-minus),
            averaged_edge_effect_cannot_be_I_over_3=True,
            graph_independent_closed_port_semigroup_not_universal=True,
            no_claim_all_larger_models_or_geometries_fail=True)

def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    baseline = Path(__file__).with_name("retained_port_instrument_mixing_results.json")
    return dict(round=494,scientific_baseline_round=492,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        baseline_results_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
        scope=dict(retained_same_D_and_G=True,
            real_initial_port_pinch_and_timed_Lueders_are_extra_inputs=True,
            arbitrary_initial_graph_and_passive_reference=True,
            exact_classical_microstep_history_not_claimed=True,
            graph_averaging_is_special_six_graph_identity=True,
            no_stationary_information_erasure_claim=True,
            no_closed_H_natural_dissipation_claim=True,
            no_spatial_dimension_or_complete_position_contract_derived=True,
            full_GR_goal_completed=False,phase_closure_triggered=False),
        observations=OBS)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run",action="store_true")
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    data = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding="utf-8")) == data
    elif not args.dry_run:
        with TARGET.open("x",encoding="utf-8") as handle:
            json.dump(data,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
    print(json.dumps(data,ensure_ascii=False,indent=2))

