"""Round 492: retain the same data and graph; mix by actual port instruments.

Scientific baseline 490. No graph/data reset, fresh network, or scalar Markov
closure is assumed. Fresh isolated reader records and timed instruments are
additional resources. The general proof is algebraic; a separate fixed-point
integer certificate gives one finite, reference-uniform mixing bound.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib, io, json, math, platform, unittest
from pathlib import Path
import numpy as np
import port_spectral_speed_ruler as prior

TARGET = Path(__file__).with_name("retained_port_instrument_mixing_results.json")
OBS = {}
DIM = 36
BLOCK = 6
SUPER = 216
SCALE = 1 << 64
LEAVES = (0, 3, 4, 5)

def short(x):
    return float(f"{float(x):.12g}")

def ceil_div(a, b):
    return -((-int(a)) // int(b))

def round_div(a, b):
    # Nearest integer, with either allowed choice at exact half-integers.
    return (int(a) + int(b)//2) // int(b)

def ceil_sqrt(a):
    r = math.isqrt(int(a))
    return r + (r*r < a)

@lru_cache(None)
def model():
    trees, full_h, f, adj, embedding, hj, hf, swap = prior.model()
    h = (hj + hf + 5*np.eye(36, dtype=np.int64)).astype(np.int64)
    ports = [np.kron(np.diag(np.eye(6, dtype=np.int64)[a]),
                     np.eye(6, dtype=np.int64)) for a in range(6)]
    return trees, full_h, f, adj, embedding, h, ports

def dephase(x):
    y = np.zeros_like(x)
    for a in range(6):
        sl = slice(6*a, 6*a+6)
        y[sl, sl] = x[sl, sl]
    return y

def step(x, u):
    return dephase(u @ x @ u.conj().T)

def rank_mod(matrix, prime=1009):
    a = np.array(matrix, dtype=np.int64) % prime
    row = 0
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, col])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        a[[row, pivot]] = a[[pivot, row]]
        a[row] = (a[row]*pow(int(a[row, col]), -1, prime)) % prime
        # Before multiplication every entry is in [0,prime); one product is
        # below prime**2, safely within int64. No integer matrix dot is used.
        assert 2*prime**2 < 2**63
        others = np.flatnonzero(a[:, col])
        others = others[others != row]
        a[others] = (a[others] - a[others, col, None]*a[row]) % prime
        row += 1
        if row == a.shape[0]:
            break
    return row

@lru_cache(None)
def hermitian_coordinates():
    # X_ij = x_re + i*x_im for i<j. Basis norms squared are 1,2,2.
    local = [(i, i, 0) for i in range(6)]
    local += [(i, j, typ) for i in range(6) for j in range(i+1, 6)
              for typ in (1, 2)]
    weights = np.array([1 if typ == 0 else 2 for _, _, typ in local]*6,
                       dtype=np.int64)
    return local, weights

def encode_blocks(blocks):
    local, _ = hermitian_coordinates()
    return np.array([float(x[i,j].real if typ != 2 else x[i,j].imag)
                     for x in blocks for i,j,typ in local])

def decode_blocks(x):
    local, _ = hermitian_coordinates()
    blocks = []
    for a in range(6):
        z = np.zeros((6,6), dtype=complex)
        for v,(i,j,typ) in zip(x[36*a:36*a+36], local):
            if typ == 0:
                z[i,i] = v
            elif typ == 1:
                z[i,j] += v; z[j,i] += v
            else:
                z[i,j] += 1j*v; z[j,i] -= 1j*v
        blocks.append(z)
    return blocks

def unitary_integer(time=F(1,16), order=24):
    """Exact polynomial -> rounded Gaussian integers / SCALE, no float input."""
    h = model()[5].astype(object)
    den = time.denominator**order * math.factorial(order)
    real = np.zeros((36,36), dtype=object)
    imag = np.zeros((36,36), dtype=object)
    hp = np.eye(36, dtype=object)
    for k in range(order+1):
        coefficient = (time.numerator**k * time.denominator**(order-k)
                       * (math.factorial(order)//math.factorial(k)))
        if k % 4 == 0: real += coefficient*hp
        if k % 4 == 1: imag -= coefficient*hp
        if k % 4 == 2: real -= coefficient*hp
        if k % 4 == 3: imag += coefficient*hp
        if k < order: hp = hp @ h  # Python arbitrary precision integers.
    vr = np.array([[round_div(SCALE*z, den) for z in row] for row in real],
                  dtype=object)
    vi = np.array([[round_div(SCALE*z, den) for z in row] for row in imag],
                  dtype=object)
    x = 9*time
    assert x < order+2
    tail = x**(order+1) / math.factorial(order+1) / (1-x/F(order+2))
    assert tail < F(1, SCALE)
    # Quantizing a 36x36 complex matrix costs <=36/SCALE in Frobenius norm.
    return vr, vi, tail, F(37, SCALE)

def complex_product(ar, ai, br, bi):
    return ar*br-ai*bi, ar*bi+ai*br

def super_integer(vr, vi):
    """Coordinate matrix of Delta Ad(V), rounded once to denominator SCALE."""
    local, weights = hermitian_coordinates()
    out = np.zeros((216,216), dtype=object)
    for b in range(6):
        for a in range(6):
            ar = vr[6*b:6*b+6, 6*a:6*a+6]
            ai = vi[6*b:6*b+6, 6*a:6*a+6]
            for row,(r,s,read) in enumerate(local):
                for col,(k,l,typ) in enumerate(local):
                    re, im = complex_product(ar[r,k],ai[r,k],ar[s,l],-ai[s,l])
                    if typ:
                        re2, im2 = complex_product(ar[r,l],ai[r,l],
                                                   ar[s,k],-ai[s,k])
                        if typ == 1: re,im = re+re2,im+im2
                        else: re,im = -(im-im2),re-re2
                    value = im if read == 2 else re
                    out[36*b+row,36*a+col] = round_div(value,SCALE)
    return out

def stationary_integer():
    p = np.zeros((216,216), dtype=object)
    diagonal = [36*a+i for a in range(6) for i in range(6)]
    value = round_div(SCALE,36)
    for i in diagonal:
        for j in diagonal: p[i,j] = value
    return p

@lru_cache(None)
def finite_certificate():
    vr,vi,tail,eu = unitary_integer()
    transition = super_integer(vr,vi)
    matrix = transition-stationary_integer()
    # Weighted HS norm: conjugate by diagonal sqrt(weights). Matrix entry
    # rounding has norm <=216/SCALE; rounded Pi costs <=36/SCALE.
    e0 = (2+eu)*eu + F(216+36,SCALE)
    error_units = ceil_div(e0.numerator*SCALE,e0.denominator)
    history = []
    for level in range(1,17):
        product = matrix @ matrix  # genuine arbitrary-precision integer dot.
        matrix = np.array([[round_div(z,SCALE) for z in row] for row in product],
                          dtype=object)
        error_units = 2*error_units + ceil_div(error_units**2,SCALE) + 216
        history.append(dict(power=2**level,error_units=int(error_units)))
    _,weights = hermitian_coordinates()
    # 2*w_i/w_j is exactly 1,2,or4. Omitting sqrt(2) below overestimates norm.
    weighted_sum = sum((2*int(weights[i])//int(weights[j]))*int(matrix[i,j])**2
                       for i in range(216) for j in range(216))
    frobenius_units = ceil_sqrt(weighted_sum)
    hs_bound = F(frobenius_units+error_units,SCALE)
    diamond = 36*hs_bound
    assert diamond < F(1,4000)
    payload = dict(delta="1/16",taylor_order=24,fixed_point_bits=64,
        taylor_operator_tail=str(tail),unitary_operator_error=str(eu),
        initial_superoperator_error=str(e0),initial_error_units=ceil_div(
            e0.numerator*SCALE,e0.denominator),
        repeated_square_levels=16,block_power=2**16,instrument_steps=2**16+1,
        error_units=int(error_units),weighted_frobenius_sum=int(weighted_sum),
        frobenius_units=int(frobenius_units),hs_operator_bound=str(hs_bound),
        ordinary_diamond_bound=str(diamond),ordinary_diamond_bound_decimal=short(diamond),
        asserted_diamond_bound="1/4000",error_recurrence=history,
        arithmetic="Python integers; Gaussian polynomial; weighted Hermitian coordinates",
        matrix_sha256=hashlib.sha256(
            json.dumps([[int(z) for z in row] for row in matrix],
                       separators=(",",":")).encode()).hexdigest())
    return payload,transition,matrix

class Audit(unittest.TestCase):
    def test_01_real_instrument_retains_data_graph_and_initial_reference(self):
        trees,full_h,f,adj,v,h,ports = model()
        self.assertTrue(np.array_equal(full_h@v,v@h))
        self.assertTrue(np.array_equal(sum(ports),np.eye(36,dtype=np.int64)))
        self.assertTrue(np.array_equal(sum(p@p for p in ports),np.eye(36,dtype=np.int64)))
        for b,p in enumerate(ports):
            raw = np.diag([int(d == 1<<(5-b)) for d in range(64)])
            self.assertTrue(np.array_equal(np.kron(raw,np.eye(6,dtype=np.int64))@v,v@p))
        delta=1/16;u=prior.unitary(h,delta)
        ks=[p@u for p in ports]
        self.assertLess(np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(36)),1e-12)
        rng=np.random.default_rng(492)
        v0=rng.normal(size=(36,2))+1j*rng.normal(size=(36,2));v0/=np.linalg.norm(v0)
        v1=rng.normal(size=(36,2))+1j*rng.normal(size=(36,2));v1/=np.linalg.norm(v1)
        histories0=[kb@ka@v0 for ka in ks for kb in ks]
        histories1=[kb@ka@v1 for ka in ks for kb in ks]
        overlap=sum(np.vdot(a,b) for a,b in zip(histories0,histories1))
        self.assertLess(abs(overlap-np.vdot(v0,v1)),1e-12)
        rho_r=sum(x.T@x.conj() for x in histories0)
        self.assertLess(np.linalg.norm(rho_r-v0.T@v0.conj()),1e-12)
        nonselective=sum(x@x.conj().T for x in histories0)
        target=step(step(v0@v0.conj().T,u),u)
        self.assertLess(np.linalg.norm(nonselective-target),1e-12)
        OBS["instrument"]=dict(dimension=36,full_dimension=384,port_outcomes=6,
            two_round_histories=36,overlap_residual=short(abs(overlap-np.vdot(v0,v1))),
            initial_reference_marginal_residual=short(np.linalg.norm(rho_r-v0.T@v0.conj())),
            data_replaced=False,graph_replaced=False,
            passive_old_records_and_new_blank_readers_required=True,
            unknown_reference_statement="all initial DG-R states supported in the one-excitation sector")

    def test_02_exact_block_algebra_and_peripheral_structure(self):
        _,_,f,adj,_,h,_=model()
        self.assertTrue(np.array_equal(f,f.T))
        self.assertTrue(np.array_equal(f@np.ones(6,dtype=np.int64),
                                       4*np.ones(6,dtype=np.int64)))
        projectors={a:np.diag([g[a,1] for g in adj]) for a in LEAVES}
        for a in LEAVES:
            self.assertTrue(np.array_equal(h[6*a:6*a+6,6:12],projectors[a]))
            self.assertTrue(np.array_equal(h[6*a:6*a+6,12:18],
                                           np.eye(6,dtype=np.int64)-projectors[a]))
        self.assertTrue(np.array_equal(h[6:12,12:18],np.eye(6,dtype=np.int64)))
        signatures=[tuple(int(projectors[a][g,g]) for a in LEAVES) for g in range(6)]
        self.assertEqual(len(set(signatures)),6)
        reached={0}
        for _ in range(6): reached |= {j for i in reached for j in range(6) if f[i,j]}
        self.assertEqual(len(reached),6)
        columns=[]
        for a in range(6):
            for i in range(6):
                for j in range(6):
                    x=np.zeros((36,36),dtype=np.int64);x[6*a+i,6*a+j]=1
                    columns.append((h@x-x@h).reshape(-1))
        rank=rank_mod(np.stack(columns,axis=1))
        self.assertEqual(rank,215)
        # h is a restriction of a sum of four-norm F and five unit SWAPs.
        # 4*9*(1/16)<6<2*pi; no floating trigonometric threshold is used.
        self.assertLess(F(4*9,16),6)
        OBS["exact_algebra"]=dict(block_operator_dimension=216,
            integer_commutator_rank_mod_1009=rank,kernel_dimension=1,
            leaf_signatures=[list(x) for x in signatures],
            graph_F_connected=True,
            sufficient_general_condition="kappa*J != 0; 0<delta<pi/(2*(4*abs(kappa)+5*abs(J)))",
            proof="HS equality, unmixed energy differences, complementary leaf edges, F connectivity",
            fixed_algebra="C I_36",unit_circle_eigenvalues=[1],
            no_floating_spectrum_used_for_theorem=True)

    def test_03_each_coupling_is_required_for_this_primitivity(self):
        _,_,f,adj,_,h,ports=model()
        hf=np.kron(np.eye(6,dtype=np.int64),f)
        hj=h-hf
        z=np.diag([g[0,1] for g in adj])
        graph_memory=np.kron(np.eye(6,dtype=np.int64),z)
        port_memory=ports[0]
        self.assertFalse(np.any(hj@graph_memory-graph_memory@hj))
        self.assertFalse(np.any(hf@port_memory-port_memory@hf))
        self.assertTrue(np.any(h@graph_memory-graph_memory@h))
        self.assertTrue(np.any(h@port_memory-port_memory@h))
        self.assertTrue(np.array_equal(dephase(graph_memory),graph_memory))
        self.assertTrue(np.array_equal(dephase(port_memory),port_memory))
        OBS["coupling_boundary"]=dict(kappa_zero_keeps_nonconstant_graph_projector=True,
            J_zero_keeps_each_port_population=True,
            no_claim_for_resonant_delta_or_full_64_dimensional_data=True)

    def test_04_population_is_not_an_exact_closed_classical_state(self):
        _,_,_,adj,_,h,_=model()
        g1=next(g for g,a in enumerate(adj) if a[0,1])
        g0=next(g for g,a in enumerate(adj) if not a[0,1])
        delta=F(1,1024)
        x=18*delta
        remainder=x**4/math.factorial(4)/(1-x/F(5))
        # Real h and basis inputs imply even probabilities, with leading
        # coefficient n_01(g). The same-source graph state is not discarded.
        self.assertLess(2*remainder,F(1,100)*delta**2)
        u=prior.unitary(h,float(delta))
        p1=float(np.sum(np.abs(u[6:12,g1])**2))
        p0=float(np.sum(np.abs(u[6:12,g0])**2))
        self.assertGreater(p1-p0,float(F(99,100)*delta**2))
        input_population=[1,0,0,0,0,0]
        OBS["population_scope"]=dict(initial_port_population=input_population,
            different_graph_basis_states=[g1,g0],
            delta=str(delta),probability_difference=short(p1-p0),
            strict_difference_lower_bound=str(F(99,100)*delta**2),
            per_probability_remainder=str(remainder),
            no_scalar_transition_on_all_CQ_inputs=True,
            generic_lumpability_is_old_not_claimed_new=True)

    def test_05_integer_finite_reference_uniform_mixing_certificate(self):
        payload,transition,power=finite_certificate()
        OBS["finite_certificate"]=payload
        # Independent implementation test: the real coordinates reproduce
        # the direct complex-block CP map on all 216 basis columns.
        local,weights=hermitian_coordinates()
        u=prior.unitary(model()[5],1/16)
        maximum=0.
        numeric=np.array(transition,dtype=float)/SCALE
        for col in range(216):
            x=np.zeros(216);x[col]=1
            blocks=decode_blocks(x)
            whole=np.zeros((36,36),dtype=complex)
            for a,z in enumerate(blocks): whole[6*a:6*a+6,6*a:6*a+6]=z
            out=step(whole,u)
            expected=encode_blocks([out[6*a:6*a+6,6*a:6*a+6] for a in range(6)])
            maximum=max(maximum,float(np.max(np.abs(expected-numeric[:,col]))))
        self.assertLess(maximum,1e-12)
        OBS["finite_certificate"]["coordinate_map_diagnostic_error"]=short(maximum)

    def test_06_reader_information_and_finite_control_ledger(self):
        payload,_,_=finite_certificate()
        m=payload["instrument_steps"]
        gamma=F(1,4000*m)
        mathematical=F(payload["ordinary_diamond_bound"])
        self.assertLess(mathematical+m*gamma,F(1,2000))
        # In the stationary state the last record equals the current port.
        # Include graph index g. Joint state vs product has exact TV 5/6.
        joint=[F(int(a==b),36) for a in range(6) for g in range(6) for b in range(6)]
        product=[F(1,216)]*216
        distance=sum(abs(x-y) for x,y in zip(joint,product))/2
        self.assertEqual(sum(joint),1)
        self.assertEqual(distance,F(5,6))
        OBS["resource_and_reference_scope"]=dict(
            finite_instrument_steps=m,total_ideal_wait_time=str(F(m,16)),
            retained_active_networks=1,data_qubits=6,graph_dimension=6,
            extra_fresh_graph_sources=0,extra_fresh_data_batches=0,
            separate_six_state_readers=m,three_qubit_record_encoding_capacity=3*m,
            maximum_whole_instrument_diamond_error_per_step=str(gamma),
            cumulative_control_error=str(m*gamma),
            actual_active_initial_reference_diamond_bound="1/2000",
            last_record_stationary_nonproduct_half_trace_distance=str(distance),
            initial_passive_R_not_all_newly_generated_records=True,
            all_old_records_and_detector_environments_isolated_and_counted=True,
            whole_history_information_not_erased=True,
            finite_control_precision_not_derived_from_original_H=True,
            no_statistical_independent_network_cloning_or_closed_H_dissipation=True)

def run():
    OBS.clear()
    stream=io.StringIO()
    result=unittest.TextTestRunner(stream=stream,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    baseline=Path(__file__).with_name("port_spectral_speed_ruler_results.json")
    return dict(round=492,scientific_baseline_round=490,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,
        baseline_results_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
        scope=dict(retained_graph_and_data_instrument_process=True,
            primitive_channel_is_one_excitation_and_extra_measurement_contract=True,
            no_fresh_network_per_step_required=True,
            arbitrary_initial_isolated_reference_supported=True,
            no_full_history_erasure_or_scalar_record_Markov_claim=True,
            no_spatial_dimension_or_endpoint_group_derived=True,
            full_GR_goal_completed=False,phase_closure_triggered=False),
        observations=OBS)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--dry-run",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding="utf-8"))==data
    elif not args.dry_run:
        with TARGET.open("x",encoding="utf-8") as handle:
            json.dump(data,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
    print(json.dumps(data,ensure_ascii=False,indent=2))

