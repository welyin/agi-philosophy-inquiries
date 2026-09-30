"""Round 512: certified first arrival under the existing identity Hamiltonian.

All theorem certificates use integers/Fractions. Floating linear algebra is only
used for independent small instrument and payload diagnostics.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old
import uniform_identity_receipt as identity
from distributed_role_reader import exponential, opnorm

HERE = Path(__file__).resolve().parent
TARGET = HERE/'target_identity_delivery_results.json'
OBS = {}
SCALE = 2**30
ERROR_GRID = 2**50
PRIME = 1009


@lru_cache(maxsize=1)
def model():
    trees, f, _, _, base, _, _, _ = old.model(2)
    h = identity.model(trees, f, 6)
    return trees, f, base.astype(np.int64), h


def rref(matrix):
    a = [[Q(x) for x in row] for row in matrix]
    rows, cols, at = len(a), len(a[0]), 0
    pivots = []
    for col in range(cols):
        pivot = next((r for r in range(at, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[at], a[pivot] = a[pivot], a[at]
        divisor = a[at][col]
        a[at] = [x/divisor for x in a[at]]
        for row in range(rows):
            if row != at and a[row][col]:
                factor = a[row][col]
                a[row] = [x-factor*y for x,y in zip(a[row],a[at])]
        pivots.append(col)
        at += 1
        if at == rows:
            break
    return a, pivots


def integer_nullspace(matrix):
    reduced, pivots = rref(matrix)
    cols = len(matrix[0])
    out = []
    for free in range(cols):
        if free in pivots:
            continue
        v = [Q(0)]*cols
        v[free] = Q(1)
        for row, pivot in enumerate(pivots):
            v[pivot] = -reduced[row][free]
        denominator = math.lcm(*(x.denominator for x in v))
        integers = [int(x*denominator) for x in v]
        divisor = math.gcd(*integers)
        out.append([x//divisor for x in integers])
    return out


def fraction_inverse(matrix):
    d = len(matrix)
    augmented = [list(map(Q,row))+[Q(i==j) for j in range(d)] for i,row in enumerate(matrix)]
    reduced, pivots = rref(augmented)
    assert pivots == list(range(d))
    assert all(reduced[i][j] == (i==j) for i in range(d) for j in range(d))
    return [row[d:] for row in reduced]


def rank_mod(matrix):
    basis = {}
    for column in matrix.T:
        v = column.astype(np.int64).copy()%PRIME
        for pivot,b in basis.items():
            if v[pivot]:
                v = (v-v[pivot]*b)%PRIME
        nonzero = np.flatnonzero(v)
        if len(nonzero):
            pivot = int(nonzero[0])
            v = (v*pow(int(v[pivot]),-1,PRIME))%PRIME
            basis[pivot] = v
    return len(basis)


def krylov_rank(h, target, wanted):
    basis, ranks = {}, []
    block = np.eye(len(h),dtype=np.int64)[:,target*42:(target+1)*42]
    for level in range(len(h)):
        for column in block.T:
            v = column.copy()%PRIME
            for pivot,b in basis.items():
                if v[pivot]:
                    v = (v-v[pivot]*b)%PRIME
            nonzero = np.flatnonzero(v)
            if len(nonzero):
                pivot = int(nonzero[0])
                basis[pivot] = (v*pow(int(v[pivot]),-1,PRIME))%PRIME
        ranks.append(len(basis))
        if len(basis) >= wanted:
            break
        block = (h @ block)%PRIME
    assert len(basis) == wanted, ranks
    return ranks


@lru_cache(maxsize=2)
def dark_certificate(target):
    _,_,base,h = model()
    columns, energies, descriptions = [], [], []
    dark_projector = np.full((252,252),Q(0),dtype=object)
    for colour in range(6):
        local_columns = []
        allowed = [j for j in range(36) if j//6 not in {target,colour}]
        for energy in (-3,-1,1):
            kernel = integer_nullspace((base-energy*np.eye(36,dtype=np.int64))[:,allowed].tolist())
            for coefficients in kernel:
                phi = np.zeros(36,dtype=np.int64)
                phi[allowed] = coefficients
                assert np.array_equal(base @ phi,energy*phi)
                vector = np.zeros(252,dtype=np.int64)
                vector[np.arange(36)*7+colour+1] = phi
                columns.append(vector); energies.append(energy); local_columns.append(phi)
                descriptions.append(dict(colour=colour,energy=energy,base_vector=phi.tolist()))
        if local_columns:
            z = np.array(local_columns,dtype=np.int64).T
            gram = z.T @ z
            inv = np.array(fraction_inverse(gram.tolist()),dtype=object)
            project = z.astype(object) @ inv @ z.T.astype(object)
            assert np.array_equal(project @ z,z)
            assert np.array_equal(project.T,project)
            rows = np.arange(36)*7+colour+1
            dark_projector[np.ix_(rows,rows)] = project
    z = np.array(columns,dtype=np.int64).T
    assert np.array_equal(h @ z,z*np.array(energies))
    assert not np.any(z[target*42:(target+1)*42])
    assert not np.any(z[::7])
    assert rank_mod(z) == z.shape[1]
    ranks = krylov_rank(h,target,252-z.shape[1])
    bright = np.eye(252,dtype=object)-dark_projector
    return dict(target=target,dark_dimension=z.shape[1],bright_dimension=252-z.shape[1],
                prime=PRIME,krylov_ranks_by_degree=ranks,integer_basis=descriptions,
                bright=bright,z=z)


def sparse_left(h, matrix):
    out = np.empty_like(matrix,dtype=object)
    for i,row in enumerate(h):
        indices = np.flatnonzero(row)
        out[i] = sum(int(row[j])*matrix[j] for j in indices)
    return out


def nearest_ratio(numerator, denominator):
    return (2*numerator+denominator)//(2*denominator)


def upward(value):
    scaled = value*ERROR_GRID
    return Q(-(-scaled.numerator//scaled.denominator),ERROR_GRID)


@lru_cache(maxsize=1)
def fixed_unitary():
    h = model()[3]
    degree = 24
    denominator = 8**degree*math.factorial(degree)
    real = np.zeros_like(h,dtype=object)
    imag = np.zeros_like(h,dtype=object)
    power = np.eye(len(h),dtype=object)
    for k in range(degree+1):
        multiplier = denominator//(8**k*math.factorial(k))
        part = real if k%2 == 0 else imag
        part += (-1)**((k+1)//2)*multiplier*power
        if k < degree:
            power = sparse_left(h,power)
    norm_bound = max(sum(abs(int(x)) for x in row) for row in h)
    assert norm_bound == 11
    x = Q(norm_bound,8)
    tail = x**(degree+1)/(math.factorial(degree+1)*(1-x/Q(degree+2)))
    rounded = []
    for part in (real,imag):
        rounded.append(np.array([nearest_ratio(int(v)*SCALE,denominator) for v in part.flat],
                                dtype=np.int64).reshape(part.shape))
    return (*rounded,tail)


def product_bound(a,b):
    ma,mb = int(np.max(np.abs(a))),int(np.max(np.abs(b)))
    limit = 2**63-1
    assert len(a)*ma < limit and len(b)*mb < limit
    row = int(np.max(np.sum(np.abs(a),axis=1)))
    col = int(np.max(np.sum(np.abs(b),axis=0)))
    return min(row*mb,ma*col)


def fixed_product(ar,ai,br,bi):
    # Bound each dot product AND each real/imaginary addition before int64 use.
    real_bound = product_bound(ar,br)+product_bound(ai,bi)
    imag_bound = product_bound(ar,bi)+product_bound(ai,br)
    largest = max(real_bound,imag_bound)
    fast = largest < 2**63-1
    if not fast:
        ar,ai,br,bi = [a.astype(object) for a in (ar,ai,br,bi)]
    rr = ar @ br-ai @ bi
    ii = ar @ bi+ai @ br
    # Add half a unit safely before integer division; ties round toward +infinity.
    if fast and largest+SCALE//2 < 2**63:
        out = [((a+SCALE//2)//SCALE).astype(np.int64) for a in (rr,ii)]
    else:
        out = [np.array([nearest_ratio(int(v),SCALE) for v in a.flat],dtype=np.int64).reshape(a.shape)
               for a in (rr,ii)]
    return (*out,largest,not fast)


@lru_cache(maxsize=2)
def contraction_certificate(target):
    c = dark_certificate(target)
    bright = c['bright']
    br = np.array([nearest_ratio(v.numerator*SCALE,v.denominator) if isinstance(v,Q)
                   else int(v)*SCALE for v in bright.flat],dtype=np.int64).reshape(bright.shape)
    ur,ui,tail = fixed_unitary()
    ur,ui = ur.copy(),ui.copy()
    ur[target*42:(target+1)*42] = 0
    ui[target*42:(target+1)*42] = 0
    ar,ai,bound,slow = fixed_product(ur,ui,br,np.zeros_like(br))
    round_complex = Q(252,SCALE)
    eu = tail+round_complex
    eb = Q(252,2*SCALE)
    error = upward(eu*(1+eb)+eb+round_complex)
    rows = [dict(power=1,error=str(error))]
    max_product_bound, fallback = bound,int(slow)
    for step in range(16):
        ar,ai,bound,slow = fixed_product(ar,ai,ar,ai)
        max_product_bound = max(max_product_bound,bound)
        fallback += int(slow)
        error = upward(2*error+error**2+round_complex)
        rows.append(dict(power=2**(step+1),error=str(error)))
    squares = sum(int(v)**2 for v in ar.flat)+sum(int(v)**2 for v in ai.flat)
    assert 16*squares <= SCALE**2  # exact Frobenius bound <= 1/4
    assert error <= Q(1,4)
    digest = hashlib.sha256(ar.astype('<i8').tobytes()+ai.astype('<i8').tobytes()).hexdigest()
    return dict(target=target,sample_interval='1/8',block_steps=65536,
        h_norm_bound=11,taylor_degree=24,taylor_operator_tail=str(tail),
        binary_fraction_bits=30,error_grid_bits=50,initial_U_error=str(eu),
        initial_B_error=str(eb),per_product_rounding_bound=str(round_complex),
        largest_prechecked_product_bound=max_product_bound,object_integer_fallbacks=fallback,
        error_history=rows,final_integer_squared_frobenius=squares,
        final_squared_frobenius_denominator=SCALE**2,certified_frobenius_upper='1/4',
        final_operator_error=str(error),certified_contraction_upper='1/2',
        final_fixed_matrix_sha256=digest,exact_integer_certificate=True)


class Audit(unittest.TestCase):
    def test_01_existing_model_and_role_alignment(self):
        trees,f,base,h = model()
        self.assertEqual(h.shape,(252,252))
        self.assertTrue(np.array_equal(h,h.T))
        self.assertTrue(np.array_equal(h[::7,::7],base))
        self.assertTrue(all([sum(v in e for e in t) for v in range(6)] == [3,3,1,1,1,1] for t in trees))
        self.assertEqual(max(np.sum(abs(h),axis=1)),11)
        # Explicit relabeling maps round 503's saved hardware to this same model.
        ot,_,_,of = identity.old500.six_model()
        oh = identity.model(ot,of,6)
        permutation = [2,0,1,3,4,5]  # old vertex -> aligned vertex
        lookup = {t:j for j,t in enumerate(trees)}
        rows = []
        for v in range(6):
            for t in ot:
                mapped = frozenset(tuple(sorted((permutation[a],permutation[b]))) for a,b in t)
                gi = lookup[mapped]
                for colour in range(7):
                    cc = 0 if colour == 0 else permutation[colour-1]+1
                    rows.append((permutation[v]*6+gi)*7+cc)
        self.assertTrue(np.array_equal(h[np.ix_(rows,rows)],oh))
        for pair in ((0,1),(2,3),(3,4),(4,5)):
            perm = list(range(6))
            perm[pair[0]],perm[pair[1]] = perm[pair[1]],perm[pair[0]]
            rows = []
            for v in range(6):
                for t in trees:
                    mapped = frozenset(tuple(sorted((perm[a],perm[b]))) for a,b in t)
                    for colour in range(7):
                        cc = 0 if colour == 0 else perm[colour-1]+1
                        rows.append((perm[v]*6+lookup[mapped])*7+cc)
            self.assertTrue(np.array_equal(h[np.ix_(rows,rows)],h))
        OBS['existing_model'] = dict(vertices=6,graph_states=6,colours=7,dimension=252,
            internal_vertices=[0,1],leaves=[2,3,4,5],old_to_aligned_vertices=permutation,
            exact_relabel_of_503=True,all_role_permutation_generators_checked=True,
            new_propagation_terms=0)

    def test_02_exact_dark_space_and_query_sources(self):
        certificates = []
        for target,dimension in ((0,35),(2,15)):
            c = dark_certificate(target)
            self.assertEqual(c['dark_dimension'],dimension)
            self.assertEqual(c['krylov_ranks_by_degree'][-1],252-dimension)
            certificates.append({k:v for k,v in c.items() if k not in ('bright','z')})
        OBS['exact_dark_space'] = certificates

    def test_03_fixed_integer_finite_arrival_certificate(self):
        OBS['finite_arrival_certificates'] = [contraction_certificate(a) for a in (0,2)]
        for c in OBS['finite_arrival_certificates']:
            self.assertTrue(c['exact_integer_certificate'])
            self.assertLessEqual(Q(c['final_operator_error']),Q(1,4))

    def test_04_full_stopped_instrument_and_reference(self):
        h = model()[3]
        u = exponential(h/8)
        query = np.arange(36)*7
        rows = []
        for target in (0,2):
            remaining = np.eye(252,dtype=complex)[:,query]
            effect = np.zeros((36,36),complex)
            arrivals = []
            for step in range(12):
                remaining = u @ remaining
                at = remaining[target*42:(target+1)*42].copy()
                arrivals.append(at)
                effect += at.conj().T @ at
                remaining[target*42:(target+1)*42] = 0
            complete = effect+remaining.conj().T @ remaining
            error = opnorm(complete-np.eye(36))
            self.assertLess(error,2e-13)
            rows.append(dict(target=target,finite_steps=12,all_query_columns=36,
                full_instrument_completeness_error=float(f'{error:.12g}'),
                early_arrival_effect_extremes=[float(f'{v:.12g}') for v in
                    (np.linalg.eigvalsh(effect)[0],np.linalg.eigvalsh(effect)[-1])],
                timeout_and_arrival_times_retained=True))
        OBS['stopped_instrument'] = rows

    def test_05_payload_factor_of_physical_single_packet(self):
        trees,f,_,h = model()
        payload = 2
        dimension = len(h)*payload
        direct = np.zeros((dimension,dimension),dtype=np.int64)
        def index(position,graph,colour,bit):
            return ((position*6+graph)*7+colour)*payload+bit
        for v in range(6):
            for gi,tree in enumerate(trees):
                for colour in range(7):
                    for bit in range(payload):
                        col = index(v,gi,colour,bit)
                        for gj in range(6):
                            direct[index(v,gj,colour,bit),col] += f[gj,gi]
                        for a,b in tree:
                            where = b if v == a else a if v == b else v
                            direct[index(where,gi,colour,bit),col] += 1
                        if colour in (0,v+1):
                            changed = v+1 if colour == 0 else 0
                            direct[index(v,gi,changed,bit),col] += 1
        self.assertTrue(np.array_equal(direct-5*np.eye(dimension,dtype=np.int64),
                                       np.kron(h,np.eye(payload,dtype=np.int64))))
        OBS['payload_extension'] = dict(payload_dimension=payload,full_active_dimension=dimension,
            local_vacuum_plus_colour_payload_dimension=1+7*payload,
            exact_all_column_generator_factorization=True,
            shared_scalar_removed=5,capture_store_and_isolation_are_inputs=True)

    def test_06_correlated_payload_steering_and_full_instrument(self):
        h = model()[3]
        u = exponential(h/8)
        target,start = 2,3
        columns = (start*6+np.arange(6))*7
        k = u[target*42:(target+1)*42][:,columns]
        effect = k.conj().T @ k
        values,vectors = np.linalg.eigh(effect)
        low,high = vectors[:,0],vectors[:,-1]
        success = float((values[0]+values[-1])/2)
        conditional_bit_one = float(values[-1]/(values[0]+values[-1]))
        self.assertGreater(conditional_bit_one,.50001)
        # B/R correlation may be steered at an early success, despite no B coupling.
        joint = np.column_stack((low,high))/np.sqrt(2)
        good = k @ joint
        failed = u[:,columns] @ joint
        failed[target*42:(target+1)*42] = 0
        marginal = good.conj().T @ good+failed.conj().T @ failed
        self.assertLess(np.linalg.norm(marginal-np.eye(2)/2),1e-13)
        OBS['payload_information_scope'] = dict(early_success_probability=success,
            correlated_payload_conditional_one_probability=conditional_bit_one,
            unconditional_full_instrument_payload_marginal_error=float(np.linalg.norm(marginal-np.eye(2)/2)),
            every_success_preserves_correlated_payload_marginal_claimed=False,
            final_erasure_flag_channel_error_bound='2 * 4^(-m)',
            graph_and_colour_are_retained_in_environment_not_reset=True)

    def test_07_finite_resource_and_error_budget(self):
        blocks = 5
        steps = 65536*blocks
        failure = Q(1,4**blocks)
        ordinary_error = Q(1,1024*steps)
        total_channel_error = 2*failure+steps*ordinary_error
        self.assertEqual(steps,327680)
        self.assertEqual(total_channel_error,Q(3,1024))
        OBS['declared_budget'] = dict(blocks=blocks,detector_steps=steps,
            free_interval='1/8',natural_wait_upper=str(Q(steps,8)),
            ideal_timeout_upper=str(failure),ideal_delivery_erasure_diamond_upper=str(2*failure),
            per_step_ordinary_diamond_error_budget=str(ordinary_error),
            resulting_delivery_diamond_upper=str(total_channel_error),
            full_raw_binary_detection_records=steps,
            clock_capture_storage_transport_costs_are_additional=True,
            exact_ideal_expected_detector_steps_upper=str(Q(4*65536,3)),
            size_uniform_or_autonomous_budget_claimed=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(date='2026-09-28',round=512,scientific_baseline_round=511,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py','uniform_identity_receipt.py',
            'research_note_502.md','research_note_503.md','research_note_511.md')},
        scope=dict(existing_identity_hamiltonian_reused=True,
            exact_query_to_target_bright_space_certificate=True,
            exact_integer_finite_arrival_budget=True,
            arbitrary_unknown_graph_and_reference_covered=True,
            payload_identity_factor_and_erasure_channel_bound=True,
            passive_graph_free_evolution_preserved=False,
            all_conditional_payload_marginals_unchanged=False,
            all_network_sizes_proved=False,autonomous_detector_clock_capture_derived=False,
            complete_511_parallel_record_routing_implemented=False,
            dimension_three_generated=False,full_GR_goal_completed=False,
            phase_closure_triggered=False),observations=OBS)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')},ensure_ascii=False))
