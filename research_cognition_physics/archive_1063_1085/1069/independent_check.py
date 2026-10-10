"""Independent numerical witnesses for round 1069. No author imports or writes.
All universal claims rely on the report proof; these finite witnesses test formulas.
The overlap a=3/4 differs from the author witness a=9/25.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
pauli = [X, Y, Z]
kets = [np.array([1, 0], complex), np.array([0, 1], complex)]
a = 3.0 / 4.0
qkets = [np.array([np.sqrt(a), np.sqrt(1-a)], complex),
         np.array([-np.sqrt(1-a), np.sqrt(a)], complex)]
outer = lambda x: np.outer(x, x.conj())
P = [outer(x) for x in kets]
Q = [outer(x) for x in qkets]
W = [Y, X]
V = [X, Z]
MP = [np.kron(P[r], W[r]) for r in range(2)]
MQ = [np.kron(Q[r], V[r]) for r in range(2)]
checks = {}

def check(group, left, right):
    residual = float(np.linalg.norm(left-right))
    checks.setdefault(group, []).append(residual)
    if not residual < 2e-12:
        raise AssertionError((group, residual))

# Sharp effects and repeatability do not require private-memory neutrality.
for projections, operators in [(P, MP), (Q, MQ)]:
    for z, op in zip(projections, operators):
        project = np.kron(z, I)
        check('effects', op.conj().T @ op, project)
        check('repeatability', project @ op, op)
        for observable in [I, X, Y, Z]:
            check('branch_closure', op.conj().T @ np.kron(observable, I) @ op,
                  np.kron(z @ observable @ z, I))

# A known prior P outcome s is available. This is not a recovery with unknown s.
# The correction is a full unitary, hence CPTP on the entire physical carrier.
transition = []
for s in range(2):
    row = []
    for r in range(2):
        dunitary = np.column_stack([kets[s], kets[1-s]]) @ np.column_stack([qkets[r], qkets[1-r]]).conj().T
        correction = np.kron(dunitary, V[r].conj().T)
        check('recovery_unitarity', correction.conj().T @ correction, np.eye(4))
        amplitude = np.vdot(qkets[r], kets[s])
        row.append(float(abs(amplitude)**2))
        check('recovery_on_code', correction @ MQ[r] @ np.kron(P[s], I),
              amplitude * np.kron(P[s], I))
        check('recovery_after_actual_P_branch', correction @ MQ[r] @ MP[s], amplitude * MP[s])
    transition.append(row)
    check('recovery_probability_sum', np.array([sum(row)]), np.array([1.0]))

# Keep a nontrivial K-reference entangled state; recorded correction preserves it.
phi = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
rhoKR = outer(phi)
for s in range(2):
    rho = np.kron(P[s], rhoKR)
    recovered = np.zeros((8, 8), complex)
    for r in range(2):
        dunitary = np.column_stack([kets[s], kets[1-s]]) @ np.column_stack([qkets[r], qkets[1-r]]).conj().T
        corrected_branch = np.kron(np.kron(dunitary, V[r].conj().T) @ MQ[r], I)
        recovered += corrected_branch @ rho @ corrected_branch.conj().T
    check('reference_recovery', recovered, rho)

# A mixed separable task resource has the required opposite-role coherence.
w = 1.0 / 4.0
singlet = np.array([0, 1, -1, 0], complex)/np.sqrt(2)
omega = (1-w)*np.eye(4)/4 + w*outer(singlet)
decomposition = 0.25*np.kron(I/2, I/2)
for sigma in pauli:
    for sign in [-1, 1]:
        decomposition += 0.125*np.kron((I+sign*sigma)/2, (I-sign*sigma)/2)
check('werner_separable_decomposition', decomposition, omega)
check('werner_normalization', np.array([np.trace(omega)]), np.array([1.0]))
check('werner_purity', np.array([np.trace(omega@omega)]), np.array([19/64]))

def phase(z, t):
    return I + (np.exp(1j*t)-1)*z

cross_norms = []
for z in [P[0], Q[0]]:
    cross = np.kron(z, I-z) @ omega @ np.kron(I-z, z)
    cross_norms.append(float(np.linalg.norm(cross)))
    check('werner_cross_block', np.array([np.linalg.norm(cross)]), np.array([w/2]))
    for t in [0.17, 0.71, 1.2, -0.43]:
        u = phase(z, t)
        joint = np.kron(u, u)
        check('werner_same_revision', joint @ omega @ joint.conj().T, omega)

# Preserve an arbitrary entangled private K_A K_B relation in physical A/B order.
# Initial ordering D_A,D_B,K_A,K_B is changed to D_A,K_A,D_B,K_B.
full = np.kron(omega, rhoKR).reshape((2,)*8).transpose(0,2,1,3,4,6,5,7).reshape(16,16)
for z in [P[0], Q[0]]:
    for t in [0.37, 0.91]:
        local = np.kron(phase(z,t), I)
        joint = np.kron(local, local)
        check('private_joint_resource', joint @ full @ joint.conj().T, full)

# A strict contraction of coherence violates preservation of the same resource.
lam = 3.0/5.0
out = omega.copy()
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                out[2*i+j, 2*k+l] *= (1 if i==k else lam)*(1 if j==l else lam)
check('dephasing_cross_formula', np.array([out[1,2]]), np.array([-9/200]))
expected_defect = np.sqrt(2)*w*(1-lam**2)/2
check('dephasing_defect_formula', np.array([np.linalg.norm(out-omega)]), np.array([expected_defect]))
if not np.linalg.norm(out-omega) > 0.1:
    raise AssertionError('The nonunitary revision must fail the same-resource test')

pt = omega.reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4)
result = {
    'round': 1069,
    'independent_algorithm': 'standalone_numpy_operator_identities_and_separable_mixture',
    'author_code_imported': False,
    'overlap_a': a,
    'memory_dimension': 2,
    'known_prior_P_outcome_required': True,
    'group_check_counts': {k: len(v) for k,v in checks.items()},
    'matrix_checks': sum(len(v) for v in checks.values()),
    'group_max_residuals': {k: max(v) for k,v in checks.items()},
    'max_residual': max(max(v) for v in checks.values()),
    'transition_probabilities_P_to_Q': transition,
    'werner_w': w,
    'werner_eigenvalues': np.linalg.eigvalsh(omega).tolist(),
    'werner_partial_transpose_eigenvalues': np.linalg.eigvalsh(pt).tolist(),
    'werner_cross_block_norms': cross_norms,
    'full_carrier_dimension_per_side': 4,
    'full_joint_resource_rank': int(np.linalg.matrix_rank(full)),
    'dephasing_lambda': lam,
    'cross_coherence_before': float(omega[1,2].real),
    'cross_coherence_after': float(out[1,2].real),
    'dephasing_defect_frobenius': float(np.linalg.norm(out-omega)),
    'dephasing_expected_defect': float(expected_defect),
    'passed': True,
    'spatial_dimension_derived': False,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
