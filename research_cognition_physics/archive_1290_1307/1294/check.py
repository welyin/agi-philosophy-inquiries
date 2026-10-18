"""1294: one author group. Exact Laurent walks and rational tails certify bounds.
Full finite electric-basis matrices cross-check Gauss, records, sources and R.
Default/--check is read-only; --write creates results.json exclusively.
"""
import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results.json'
TOL = 2e-11


def encode(x):
    if isinstance(x, F):
        return f'{x.numerator}/{x.denominator}'
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    return x


def close(a, b=0):
    e = float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    assert e < TOL, e
    return e


def bounds():
    t, beta, tmin = F(1, 512), F(1, 1024), F(1, 1024)
    # exp(x) < 1/(1-x), 0<x<1; no floating exponential enters certificates.
    exp_t, exp_b = 1/(1-4*t), 1/(1-4*beta)
    tail_over_t3 = F(4**4, 24)*t*exp_t
    assert tail_over_t3 < F(1, 40)
    lo0, lo1, hi1 = F(37, 120)**2, F(17, 120)**2, F(23, 120)**2
    assert lo0-hi1 == F(7, 120)
    assert lo1 > F(1, 50)
    thermal_tail_over_beta4 = F(16*4**5, 120)*beta*exp_b
    assert thermal_tail_over_beta4 < F(4, 15)
    assert 1-thermal_tail_over_beta4 > F(11, 15)
    return {
        'time': t, 'window': [tmin, t], 'beta': beta,
        'exp_4t_upper': exp_t, 'amplitude_tail_div_t3_upper': tail_over_t3,
        'public_tail_div_t3': F(1, 40),
        'N3_q1_probability_lower_at_t': lo1*t**6,
        'cross_group_public_lower_at_t': t**6/50,
        'cross_group_public_lower_all_window': tmin**6/50,
        'N3_q0_minus_q1_lower': F(7, 120)*t**6,
        'thermal_tail_div_beta4_upper': thermal_tail_over_beta4,
        'Zbar0_minus_Zbar1_common_lower': F(11, 15)*beta**4,
        'thermal_sector_probability_gap_lower': {
            2: F(11, 240*2)*beta**4, 3: F(11, 240*3)*beta**4},
        'physical_Z_difference_has_exp_minus_4beta': True,
        'instrument_accuracy_or_sampling_budget_certified': False,
    }


def laurent_walks():
    # Htilde entries are integer Laurent polynomials in z=exp(i phi).
    entries = {}
    def add(i, j, power, value):
        key = (i, j, power)
        entries[key] = entries.get(key, 0)+value
    for x in range(4):
        add(2*x+1, 2*x+1, 0, 1)    # pointer energy
    for x in range(3):
        for p in range(2):
            add(2*(x+1)+p, 2*x+p, 0, -1)
            add(2*x+p, 2*(x+1)+p, 0, -1)
    for p in range(2):
        add(p, 6+p, 1, -1)
        add(6+p, p, -1, -1)
    add(4, 5, 0, 1)
    add(5, 4, 0, 1)
    rows = {}
    for (i, j, k), a in entries.items():
        rows.setdefault(i, []).append((j, k, a))
    current = {(i, i, 0): 1 for i in range(8)}
    traces, vectors = [], []
    for order in range(5):
        tr = {}
        vec = {}
        for (i, j, k), a in current.items():
            if i == j:
                tr[k] = tr.get(k, 0)+a
            if j == 0 and i % 2 == 1:
                vec[(i, k)] = a
        traces.append({k: a for k, a in tr.items() if a})
        vectors.append({k: a for k, a in vec.items() if a})
        nxt = {}
        for (i, j, k), a in current.items():
            for l, power, b in rows.get(j, []):
                key = (i, l, k+power)
                nxt[key] = nxt.get(key, 0)+a*b
        current = {key: a for key, a in nxt.items() if a}
    assert all(vectors[j] == {} for j in range(3))
    assert vectors[3] == {(5, -1): 1, (5, 0): 1}
    assert all(set(tr).issubset({0}) for tr in traces[:4])
    assert traces[4][-1] == traces[4][1] == 8
    # exp(+i phi) and exp(-i phi) each contribute 8 at fourth order.
    diffs = {}
    for N in (2, 3):
        def exact_cos(k):
            r = k % N
            return F(1) if r == 0 else (F(-1) if N == 2 else F(-1, 2))
        d = sum(F(a)*(1-exact_cos(k)) for k, a in traces[4].items())
        assert d == (32 if N == 2 else 24)
        diffs[N] = d
    return {
        'record_cubic_Laurent_vector': [[i, k, a] for (i, k), a in sorted(vectors[3].items())],
        'trace_polynomials_orders_0_to_4': traces,
        'trace_H4_sector_difference': diffs,
        'thermal_beta4_difference_coefficients': {N: d/F(24) for N, d in diffs.items()},
        'all_lower_record_orders_zero': True,
        'all_lower_thermal_differences_zero': True,
    }


def ring(phi, remove_closing=False):
    K = np.zeros((4, 4), complex)
    for j in range(3):
        K[j+1, j] = K[j, j+1] = 1
    if not remove_closing:
        K[0, 3], K[3, 0] = np.exp(1j*phi), np.exp(-1j*phi)
    return K


NP = np.diag([0., 1.])
XP = np.array([[0., 1.], [1., 0.]])
P2 = np.diag([0., 0., 1., 0.])
PONE = np.kron(np.eye(4), NP)


def block(phi, g=1., ap=1., r=1., remove_closing=False):
    return -g*np.kron(ring(phi, remove_closing), np.eye(2)) + ap*PONE + r*np.kron(P2, XP)


def unitary(H, t):
    # Stable small-time Taylor matrix, diagnostic only; proof uses norm tail.
    result = np.eye(len(H), dtype=complex)
    term = result.copy()
    for n in range(1, 37):
        term = term@(-1j*t*H)/n
        result += term
    close(result.conj().T@result, np.eye(len(H)))
    return result


def full_model(N):
    basis = [(x, e, p) for x in range(4)
             for e in itertools.product(range(N), repeat=4) for p in range(2)]
    index = {v: i for i, v in enumerate(basis)}
    d = len(basis)
    forward = np.zeros((d, d), complex)
    wrong_forward = np.zeros((d, d), complex)
    pointer = np.zeros((d, d), complex)
    interaction = np.zeros((d, d), complex)
    charges = np.empty((4, d), int)
    physical_indices = []
    for col, (x, e, p) in enumerate(basis):
        ee = list(e)
        ee[x] = (ee[x]-1) % N
        forward[index[((x+1) % 4, tuple(ee), p)], col] = 1
        wrong_forward[index[((x+1) % 4, e, p)], col] = 1
        pointer[col, col] = p
        if x == 2:
            interaction[index[(x, e, 1-p)], col] = 1
        q = [(e[j]-e[(j-1) % 4]-(x == j)+(j == 0)) % N for j in range(4)]
        charges[:, col] = q
        if not any(q):
            physical_indices.append(col)
    K = forward+forward.conj().T
    H = 4*np.eye(d)-K+pointer+interaction
    J0 = np.zeros((d, 8*N), complex)
    for x in range(4):
        for ell in range(N):
            e = tuple((ell-(j < x)) % N for j in range(4))
            for p in range(2):
                J0[index[(x, e, p)], (x*N+ell)*2+p] = 1
    assert len(physical_indices) == 8*N
    assert set(np.flatnonzero(np.any(J0 != 0, axis=1))) == set(physical_indices)
    fourier = np.zeros((8*N, 8*N), complex)
    for q in range(N):
        for x in range(4):
            for ell in range(N):
                for p in range(2):
                    fourier[(x*N+ell)*2+p, q*8+x*2+p] = np.exp(2j*np.pi*q*ell/N)/math.sqrt(N)
    J = J0@fourier
    maxerr = max(close(J0.conj().T@J0, np.eye(8*N)), close(J.conj().T@J, np.eye(8*N)))
    B = np.zeros((8*N, 8*N), complex)
    for q in range(N):
        B[q*8:(q+1)*8, q*8:(q+1)*8] = 4*np.eye(8)+block(2*np.pi*q/N)
    maxerr = max(maxerr, close(H@J, J@B))
    for j in range(4):
        phase = np.exp(2j*np.pi*charges[j]/N)
        maxerr = max(maxerr, close((phase[:, None]-phase[None, :])*H))
    wrong = wrong_forward+wrong_forward.conj().T
    wrong_leak = float(np.max(np.abs(wrong@J0-J0@(J0.conj().T@wrong@J0))))
    assert wrong_leak >= .99
    # Removing the negative charge leaves no physical state in this one-particle sector.
    no_anchor_count = sum(all((e[j]-e[(j-1) % 4]-(x == j)) % N == 0 for j in range(4))
                          for x, e, p in basis)
    assert no_anchor_count == 0
    source_full = [np.eye(d), np.eye(d), -K, pointer, interaction]
    source_phys = [J.conj().T@S@J for S in source_full]
    maxerr = max(maxerr, close(B, 2*source_phys[0]+2*source_phys[1]+source_phys[2]+source_phys[3]+source_phys[4]))
    code = np.zeros((8*N, 2), complex)
    code[0, 0] = code[8, 1] = 1
    maxerr = max(maxerr, close(code.conj().T@B@code, 4*np.eye(2)),
                 close(code.conj().T@B@B@code, 18*np.eye(2)))
    P = source_phys[3]
    P0 = np.eye(8*N)-P
    maxerr = max(maxerr, close(P@B@P+P0@B@P0-B, -source_phys[4]))
    energy_parts = [4*np.eye(8*N), source_phys[2], source_phys[3], source_phys[4]]
    maxerr = max(maxerr, close(sum(1j*(B@A-A@B) for A in energy_parts)))
    # One nonmaximally entangled input diagnoses the complete R extension.
    U = unitary(B-4*np.eye(8*N), 1/512)
    inp = np.array([math.sqrt(.3), 0, 0, np.exp(.4j)*math.sqrt(.7)], complex).reshape(2, 2)
    psi = U@code@inp
    branches = [P0@psi, P@psi]
    rhoR = psi.T@psi.conj()
    maxerr = max(maxerr, close(sum(b.T@b.conj() for b in branches), inp.T@inp.conj()),
                 close(rhoR, inp.T@inp.conj()))
    completeness = sum((Q@U@code).conj().T@(Q@U@code) for Q in (P0, P))
    maxerr = max(maxerr, close(completeness, np.eye(2)))
    probs = [float(np.vdot(b, b).real) for b in branches]
    maxerr = max(maxerr, close(sum(probs), 1))
    # Full finite-sector Gibbs: no full unconstrained Fock claim.
    e, V = np.linalg.eigh(B-4*np.eye(8*N))
    Gibbs = (V*np.exp(-e/1024))@V.conj().T
    Gibbs /= np.trace(Gibbs)
    sector_weights = [float(np.trace(Gibbs[q*8:(q+1)*8, q*8:(q+1)*8]).real) for q in range(N)]
    assert sector_weights[0] > sector_weights[1]
    assert np.linalg.eigvalsh(B).min() >= 1-1e-12
    return {
        'kinematic_dimension': d, 'complete_physical_dimension': 8*N,
        'max_matrix_residual': maxerr,
        'wrong_unlinked_hopping_Gauss_leak': wrong_leak,
        'without_compensating_charge_physical_count': no_anchor_count,
        'code_compressed_energy': 4, 'code_compressed_H_squared': 18,
        'entangled_input_record_probabilities': probs,
        'thermal_sector_weights_diagnostic': sector_weights,
        'five_parameter_sources_retained': ['m_c', 'm_d', 'g', 'a_p', 'r'],
        'readout_energy_identity_verified': True,
        'R_marginal_preserved_nonselectively': True,
        'unknown_input_reversible_after_classical_read_claimed': False,
    }


def block_diagnostics():
    data, maxerr = {}, 0.
    for N in (2, 3):
        rows = []
        for q in range(N):
            phi = 2*np.pi*q/N
            H = block(phi)
            v = np.eye(8, dtype=complex)[:, 0]
            P13 = np.eye(4)[[0, 3, 2, 1]]
            D = np.diag([1, np.exp(1j*phi), np.exp(1j*phi), np.exp(1j*phi)])
            R = D@P13
            maxerr = max(maxerr, close(R@ring(phi)@R.conj().T, ring(-phi)))
            cubic = np.zeros(8, complex)
            cubic[5] = 1+np.exp(-1j*phi)
            maxerr = max(maxerr, close(PONE@H@H@H@v, cubic))
            times = []
            for t in (1/1024, 1/512):
                state = unitary(H, t)@v
                probability = float(np.vdot(PONE@state, PONE@state).real)
                if N == 2 and q == 1:
                    assert probability < 1e-28
                if N == 3 and q != 0:
                    assert probability > float(F(289, 14400)*F(t)**6)
                times.append({'t': t, 'p_record1': probability})
            eig = np.linalg.eigvalsh(H)
            rows.append({'q': q, 'records': times, 'Zbar_diagnostic': float(np.exp(-eig/1024).sum())})
        data[N] = rows
    # Exact dark subspace verified as a matrix invariant.
    Hpi = block(np.pi)
    dark = np.zeros((8, 2), complex)
    dark[0, 0] = 1
    dark[2, 1], dark[6, 1] = 1/math.sqrt(2), -1/math.sqrt(2)
    maxerr = max(maxerr, close(Hpi@dark-dark@(dark.conj().T@Hpi@dark)))
    # Delete r / g: no record; delete closing edge: all fluxes same local process.
    for kwargs in ({'r': 0}, {'g': 0}):
        state = unitary(block(2*np.pi/3, **kwargs), 1/512)[:, 0]
        maxerr = max(maxerr, close(PONE@state))
    maxerr = max(maxerr, close(block(0, remove_closing=True), block(2*np.pi/3, remove_closing=True)))
    for j in range(2):
        p0 = data[3][0]['records'][j]['p_record1']
        p1 = data[3][1]['records'][j]['p_record1']
        t = data[3][0]['records'][j]['t']
        assert p0-p1 > float(F(7, 120)*F(t)**6)
    return {'max_residual': maxerr, 'blocks': data,
            'deletion_tests_pass': True, 'small_matrices_not_certificates_of_error_bounds': True}


def bindings():
    paths = ['proof.md', 'check.py', '../_shared/round1294_admission.md',
             '../../archive_585_628/research_note_617.md',
             '../../archive_1009_1043/research_note_1021.md',
             '../../archive_554_584/research_note_557.md',
             '../../archive_554_584/research_note_559.md']
    return {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}


def calculate():
    return encode({
        'schema': 'round1294-finite-loop-local-record-common-source-v1',
        'round': 1294, 'author_check_groups': 1, 'verdict': 'PASS_AUTHOR_SCOPED',
        'source_bindings': bindings(), 'exact_bounds': bounds(),
        'exact_Laurent_walks': laurent_walks(),
        'full_electric_basis': {N: full_model(N) for N in (2, 3)},
        'block_diagnostics': block_diagnostics(),
        'scope': {
            'new_finite_graph_candidate_inputs': True,
            'noncontractible_means_selected_graph_only': True,
            'full_declared_8N_Gauss_sector': True,
            'common_q0_q1_code_and_arbitrary_R': True,
            'physical_local_pointer_formation': True,
            'declared_ideal_terminal_instrument_manufactured': False,
            'terminal_instrument_included_in_waiting_H': False,
            'terminal_readout_work_omitted': False,
            'full_five_parameter_source_menu': True,
            'four_dimensional_Hilbert_stress_certified': False,
            'thermal_reference_is_working_preparation': False,
            'full_TQFT_or_all_axiom_countermodel': False,
            'all_O_qualification_or_group_uniqueness': False,
            'new_axioms': 0,
        }})


def main():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group()
    g.add_argument('--write', action='store_true')
    g.add_argument('--check', action='store_true')
    args = p.parse_args()
    data = calculate()
    if args.write:
        with OUT.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(data, f, ensure_ascii=False, sort_keys=True, indent=2)
            f.write('\n')
        receipt = 'CREATED_EXCLUSIVELY'
    else:
        saved = json.loads(OUT.read_text(encoding='utf-8'))
        assert saved == data, 'saved result differs; no file was written'
        receipt = 'READ_ONLY_MATCH'
    print(json.dumps({'verdict': data['verdict'], 'receipt': receipt,
        'groups': 1, 'cross_group_window_lower': data['exact_bounds']['cross_group_public_lower_all_window'],
        'results_sha256': hashlib.sha256(OUT.read_bytes()).hexdigest()}, ensure_ascii=False))


if __name__ == '__main__':
    main()
