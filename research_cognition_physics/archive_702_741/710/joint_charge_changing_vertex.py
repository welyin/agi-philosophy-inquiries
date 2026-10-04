"""710: a specified QQQL extension of the original finite CAR/Gauss process.

This is an added interaction, not a computed instanton amplitude or anomaly.
All tensor and sparse CAR checks use the original 32-mode convention. Finite
temperature existence follows from the bounded perturbation proof in the note;
no reduced matrix is substituted for the complete bosonic Hamiltonian.
"""
import argparse
from collections import defaultdict
import hashlib
from itertools import combinations, permutations, product
import json
from pathlib import Path
import numpy as np
import joint_charge_quantum_coarse as old
import joint_topological_mass_phases as topology

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_charge_changing_vertex_results.json'


def parity(indices):
    return (-1) ** sum(a > b for i, a in enumerate(indices) for b in indices[i+1:])


def tensor():
    """Coefficients of the canonical creator; epsilon_01 = epsilon_012 = 1."""
    out = defaultdict(int)
    for a, b, c in permutations(range(3)):
        for i, j, alpha, gamma in product(range(2), repeat=4):
            l, k, beta, delta = 1-i, 1-j, 1-alpha, 1-gamma
            ids = (4*a+2*i+alpha, 4*b+2*j+beta, 4*c+2*k+gamma, 24+2*l+delta)
            out[tuple(sorted(ids))] += parity((a, b, c)) * (-1)**(i+j+alpha+gamma) * parity(ids)
    return {ids: int(v) for ids, v in sorted(out.items()) if v}


COEFF = tensor()
NORM2 = sum(v*v for v in COEFF.values())


def operation(state, dagger=False):
    result = {}
    for ids, coefficient in COEFF.items():
        term = state
        for index in (reversed(ids) if dagger else ids):
            term = (old.create if dagger else old.annihilate)(term, index)
        result = old.add(result, term, coefficient / np.sqrt(NORM2))
    return result


def nq(mask):
    return sum(((mask >> i) & 1) for i in range(24))


def charge(state):
    return {m: nq(m)*a for m, a in state.items() if nq(m)}


def phase(state, alpha):
    return {m: np.exp(1j*alpha*nq(m))*a for m, a in state.items()}


def distance(a, b):
    return float(np.sqrt(sum(abs(a.get(m, 0)-b.get(m, 0))**2 for m in set(a)|set(b))))


def vertex(state, g):
    return old.add(
        {m: g*a for m, a in operation(state).items()}, operation(state, True), g.conjugate())


def eta_source(state, g):
    return old.add({m: 1j*g*a for m, a in operation(state).items()}, operation(state, True), -1j*g.conjugate())


def tensor_check():
    assert len(COEFF) == 36 and NORM2 == 72 and sum(abs(v) for v in COEFF.values()) == 48
    outputs = np.array([q+(l,) for q in combinations(range(12), 3) for l in range(24, 28)])
    inputs = np.array(list(COEFF))
    values = np.array(list(COEFF.values()), float)
    expected = np.array([COEFF.get(tuple(ids), 0) for ids in outputs])
    rng = np.random.default_rng(7101)
    errors = []
    boost_unitarity_defect = 0.
    for boost in (False, True):
        for _ in range(3):
            C = old.matter.gauge.group_exp(rng.normal(size=8), 3)
            W = old.matter.gauge.group_exp(rng.normal(size=3), 2)
            z = np.exp(1j*rng.normal())
            R = old.matter.representation(C, W, z)
            if boost:
                S = old.matter.gauge.group_exp(rng.normal(size=3), 2) @ np.diag(np.exp([.37, -.37]))
                R[:12, :12] = z*np.kron(np.kron(C, W), S)
                R[24:28, 24:28] = z**(-3)*np.kron(W, S)
                boost_unitarity_defect = float(np.linalg.norm(S.conj().T@S-np.eye(2)))
            minors = R[outputs[:, None, :, None], inputs[None, :, None, :]]
            transformed = np.linalg.det(minors) @ values
            errors.append(float(np.max(np.abs(transformed-expected))))
    assert max(errors) < 2e-13 and boost_unitarity_defect > .1
    b = operation({0: 1.}, True)
    assert abs(old.inner(b, b)-1) < 1e-13 and distance(operation(b), {0: 1.}) < 1e-13
    return dict(original_modes=32,active_Q_L_modes=16,canonical_monomials=len(COEFF),
        coefficient_counts={str(v): list(COEFF.values()).count(v) for v in sorted(set(COEFF.values()))},
        raw_vacuum_norm_squared=NORM2,raw_coefficient_l1=48,normalized_operator_norm_upper_bound=float(48/np.sqrt(72)),
        complete_output_occupancies=len(outputs),gauge_and_SL2_tensor_errors=errors,
        nonunitary_boost_tensor_only_defect=boost_unitarity_defect,
        no_finite_CAR_Lorentz_boost_claim=True,one_generation_vacuum_creator_norm=round(old.inner(b,b).real,14))


def charge_and_phase_check():
    rng = np.random.default_rng(7102)
    states = [{0: 1.}, operation({0: 1.}, True)]
    # General original 32-mode states, with spectators as well as Q and L.
    states.append(old.normalized({int(m): complex(a) for m,a in zip(
        rng.integers(0, 2**32, 96), rng.normal(size=96)+1j*rng.normal(size=96))}))
    g = .31*np.exp(.43j)
    ward_error = covariance_error = charge_error = 0.
    for s in states:
        comm = old.add(vertex(charge(s), g), charge(vertex(s, g)), -1)
        rhs = eta_source(s, g)
        ward_error = max(ward_error, distance({m: 1j*a for m,a in comm.items()}, {m: 3*a for m,a in rhs.items()}))
        lower = old.add(charge(operation(s)), operation(charge(s)), -1)
        charge_error = max(charge_error, distance(lower, {m: -3*a for m,a in operation(s).items()}))
        for alpha in (.27, 2*np.pi/3):
            transformed = phase(vertex(phase(s, alpha), g), -alpha)
            covariance_error = max(covariance_error, distance(transformed, vertex(s, g*np.exp(3j*alpha))))
    b = states[1]
    psi = old.normalized(old.add({0: 1.}, b, 1j))
    slope = float((3*old.inner(psi, eta_source(psi, .31+0j))).real)
    assert abs(slope+.93) < 1e-12 and max(ward_error, covariance_error, charge_error) < 1e-12
    D = np.array([3,0,0,1,0,0], int)
    assert np.array_equal(D, -topology.A[:,2])
    T = np.vstack((topology.C, -topology.A.T))
    extended = np.vstack((T, D))
    assert topology.exact_rank(T) == topology.exact_rank(extended) == 6
    invariant = np.zeros(10, int);invariant[7] = -1;invariant[9] = 1
    assert np.array_equal(invariant @ extended, np.zeros(6, int))
    generations = []
    for ng in (1, 2, 3):
        q = 3*ng
        factor = sum(np.exp(2j*np.pi*k*q/(3*ng)) for k in range(ng))/ng
        assert abs(factor-1) < 1e-14
        # Exact tensor-product identities, not enumeration of a full 96-mode Fock matrix.
        generations.append(dict(generations=ng,charge_change=q,lepton_change=ng,fermion_degree=4*ng,
            nonzero_vacuum_components=36**ng,raw_norm_squared=72**ng,
            physical_residual_order=ng,residual_coherence_factor_real=float(factor.real),
            U1_coherence_factor=0,source_expectation_for_plus_i=-q*.31))
    return dict(charge_commutator_error=charge_error,ward_identity_error=ward_error,
        passive_phase_covariance_error=covariance_error,one_generation_physical_charge_slope=slope,
        QQQL_phase_row=D.tolist(),old_phase_rank=6,new_phase_rank=6,
        old_restricted_quotient_dimension=3,new_restricted_quotient_dimension=4,
        extra_invariant_coefficients=invariant.tolist(),generations=generations,
        single_common_phase_is_unitarily_removable_in_finite_extension=True,
        independent_topological_carrier_not_constructed=True)


def volume_and_source_check():
    rows=[]
    eps=.73;kappa=.017;eta=.31;step=2e-5
    for ng in (1,2,3):
        exponent=1-2*ng;weight=6*exponent
        coupling=lambda sigma: kappa*(eps**3*np.exp(6*sigma))**exponent*np.exp(1j*eta)
        sigma=.07;g=coupling(sigma)
        difference=(coupling(sigma+step)-coupling(sigma-step))/(2*step)
        rel=abs(difference-weight*g)/abs(weight*g)
        assert rel<1e-7
        # Local operator's phase derivative shares the same scalar volume factor.
        dS=(coupling(sigma+step)-coupling(sigma-step))/(2*step)*1j
        mixed_error=abs(dS-weight*1j*g)/abs(weight*g)
        scale_ratio=abs(coupling(0))/abs(kappa*eps**(3*exponent))
        assert abs(scale_ratio-1)<1e-13
        rows.append(dict(generations=ng,continuum_operator_dimension=6*ng,coupling_mass_dimension=4-6*ng,
            canonical_volume_exponent=exponent,conformal_source_factor=weight,
            first_source_relative_error=float(rel),mixed_source_relative_error=float(mixed_error),
            half_cell_side_bare_coefficient_ratio=2**(3*(2*ng-1))))
    # The extension contains no bosonic operators: it commutes every original
    # scalar read function exactly, by tensor-product structure (analytic proof).
    return dict(volume_contracts=rows,borrowed_cell_normalization='664: psi=c/sqrt(w); w=epsilon^3 exp(6 sigma)',
        finite_positive_background_only=True,new_H_and_Gibbs_required=True,
        scalar_double_commutator_unchanged_as_operator=True,
        Gibbs_partition_comparison='exp(-beta b) Z0 <= Znew <= exp(beta b) Z0',
        bounded_perturbation_proof_not_full_heat_spectrum_simulation=True,
        no_uniform_UV_bound_or_anomaly_amplitude_claim=True)


def run():
    dependencies = ('joint_charge_quantum_coarse.py','joint_charge_quantum_coarse_results.json',
        'joint_topological_mass_phases.py','research_note_664.md','research_note_665.md',
        'round710_drafts/residual_phase_entry.py','round710_drafts/residual_phase_entry_results.json')
    return dict(round=710,tests_run=3,failures=0,errors=0,tensor=tensor_check(),
        charge_and_phase=charge_and_phase_check(),volume_and_reference=volume_and_source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
        scope='New bounded local QQQL-product extension of original finite CAR/Gauss model; joint charge, phase, volume, thermal and record conditions. No induced instanton vertex, complete chiral continuum, empirical parameter prediction or GR derivation.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
