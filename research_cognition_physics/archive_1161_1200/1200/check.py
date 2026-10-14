"""Exact parity audit of the frozen round-1199 finite L-type object.

Finite Weyl algebra and rational soft-pole coefficient are certificates.
No finite-box graviton spectrum, scattering matrix, instrument or GR is inferred.
The local matrix calculation below is a diagnostic, not the exact proof.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / '1199'


def frozen_module():
    path = OLD / 'check.py'
    accepted = json.loads((OLD / 'acceptance.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(path.read_bytes()).hexdigest() == accepted['files_sha256']['check.py']
    spec = importlib.util.spec_from_file_location('round1199_frozen', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def paired_rows(rows, weights, prime):
    coeff = defaultdict(F)
    for row, weight in zip(rows, weights):
        for sign in (-1, 1):
            coeff[tuple(int(x) for x in sign * row % prime)] += weight
    mismatches = sum(coeff[tuple(-x % prime for x in key)] != val
                     for key, val in list(coeff.items()))
    assert mismatches == 0
    return {'original_rows': len(rows), 'collected_weyl_terms': len(coeff),
            'inversion_coefficient_mismatches': mismatches}


def minkowski(a, b):
    return a[0]*b[0] - sum(a[i]*b[i] for i in range(1, 4))


def check():
    old = frozen_module()
    (s, v, r, t), rw, tw = old.lattice(4)
    prime, n = 5, s.shape[1]
    previous = json.loads((OLD / 'results.json').read_text(encoding='utf-8'))
    hashes = {name: hashlib.sha256(arr.astype('<i8').tobytes()).hexdigest()
              for name, arr in zip(('S', 'V', 'R', 'T'), (s, v, r, t))}
    assert hashes == previous['stencils_sha256']
    # Inspect each independent real coupling separately. The constant in HU
    # is central, so only its nonconstant Weyl coefficients enter this test.
    pairing = {
        'HU_scalar_per_U2': paired_rows(s, [F(-prime, 2)]*len(s), prime),
        'HU_vector_per_U1': paired_rows(v, [F(-prime, 2)]*len(v), prime),
        'Hg_per_g': paired_rows(r, [prime*w for w in rw], prime),
        'HJ_per_J': paired_rows(t, [prime*w for w in tw], prime),
    }
    assert not np.any(s @ v.T)
    sx, sp = old.rref(s, prime)
    vz, vp = old.rref(v, prime)
    k = n-len(sx)-len(vz)
    assert k == 132
    # CSS basis labels are ker(V)/row(S). Inversion fixes exactly the zero
    # coset: 2[x]=0 implies [x]=0 over F5. Thus Tr_code(C)=1 exactly.
    assert pow(2, -1, prime) == 3
    dimension = prime**k
    even, odd = (dimension+1)//2, (dimension-1)//2
    assert even+odd == dimension and even-odd == 1 and odd > 0
    nonzero_odd = {}
    for name, rows, basis, pivots in (('R', r, sx, sp), ('T', t, vz, vp)):
        reduced = old.quotient(rows, basis, pivots, prime)
        i = int(np.flatnonzero(np.any(reduced, axis=1))[0])
        assert np.any(2*reduced[i] % prime)
        # O=(W-W^dagger)/(2i); tau(W^2)=tau(W^-2)=0 gives tau(O^2)=1/2.
        nonzero_odd[name] = {'row': i, 'W_and_W_squared_nonidentity': True,
                             'normalized_code_trace_O_squared': '1/2'}

    # Exact finite action on computational basis, including Z phase exponents.
    for j in range(prime):
        assert -(-j % prime) % prime == j
        assert (-((-j % prime)+1)) % prime == (j-1) % prime
        assert (-j) % prime == (prime-1)*j % prime

    # Small numerical diagnostic only; the full tensor C proof uses the above
    # exact local identities and does not enumerate 5**384 states.
    c = np.zeros((prime, prime))
    for j in range(prime):
        c[-j % prime, j] = 1
    x = np.roll(np.eye(prime), 1, axis=0)
    z = np.diag(np.exp(2j*np.pi*np.arange(prime)/prime))
    diagnostics = {'unitary_residual': float(np.linalg.norm(c.T@c-np.eye(prime))),
                   'X_inverse_residual': float(np.linalg.norm(c@x@c-x.T)),
                   'Z_inverse_residual': float(np.linalg.norm(c@z@c-z.conj().T))}
    assert max(diagnostics.values()) < 1e-13

    # Independent of the finite lattice: test the leading-soft coefficient
    # required of an Einstein scattering branch (round 326's theorem).
    hard = [(1, 0, 0, 1), (1, 0, 0, -1), (1, 1, 0, 0), (1, -1, 0, 0)]
    signs = [-1, -1, 1, 1]
    qhat = (1, 0, 1, 0)
    pol = (0, 1, 0, -1)  # spacetime diagonal, spatial norm squared=2
    assert minkowski(qhat, qhat) == 0
    assert all(minkowski(p, p) == 0 for p in hard)
    assert all(sum(eta*p[a] for eta, p in zip(signs, hard)) == 0 for a in range(4))
    assert all(pol[a]*qhat[a] == 0 for a in range(4))
    assert pol[0]-sum(pol[1:]) == 0
    terms = [F(eta*sum(pol[a]*p[a]**2 for a in range(4)), minkowski(p, qhat))
             for eta, p in zip(signs, hard)]
    soft = sum(terms)
    assert soft == 4
    mandelstam = [minkowski(tuple(hard[0][i]+hard[1][i] for i in range(4)),
                            tuple(hard[0][i]+hard[1][i] for i in range(4)))]
    for other in (hard[2], hard[3]):
        delta = tuple(hard[0][i]-other[i] for i in range(4))
        mandelstam.append(minkowski(delta, delta))
    assert mandelstam == [4, -2, -2]
    return {
        'schema': 'finite_L_parity_and_conditional_selfcoupling_v1',
        'status': 'PASS_EXACT_PARITY_AND_CONDITIONAL_SOFT_WITNESS',
        'check_groups': 2, 'box': [4, 4, 4], 'prime': prime, 'qudits': n,
        'frozen_stencils_sha256': hashes, 'full_H_pairing': pairing,
        'local_inversion': 'unitary involution; X->X^-1 and Z->Z^-1',
        'code': {'scalar_rank': len(sx), 'vector_rank': len(vz),
                 'dimension': f'5^{k}', 'C_trace': 1,
                 'even_dimension': f'(5^{k}+1)/2', 'odd_dimension': f'(5^{k}-1)/2',
                 'fixed_cosets': 1},
        'nonzero_C_odd_operators': nonzero_odd,
        'finite_exact_result': 'P_minus exp(-it H) P_plus = 0 for all real t, '
                              'including arbitrary passive reference; odd multitime '
                              'correlators vanish for any C-invariant state.',
        'perturbation_bound': {'hypothesis': 'norm of C-odd part of perturbed full H <= delta',
                               'opposite_parity_probability_upper': 'min(1,(abs(t)*delta)^2)',
                               'physical_readout_certified': False},
        'soft_witness': {'hard_momenta_at_soft_limit': hard, 'in_out_signs': signs,
                         'q_over_omega': qhat, 'polarization_diagonal': pol,
                         'polarization_norm_squared': 2, 's_t_u': mandelstam,
                         'omega_times_each_geometric_term': list(map(str, terms)),
                         'omega_times_soft_sum': str(soft),
                         'omega_M5_over_kappa_M4_limit': str(soft/2),
                         'hypotheses': 'unbroken C, C-odd one-graviton branch, '
                                       'appropriate scattering/soft limit, nonzero kappa and hard M4; '
                                       'not certified for the finite lattice'},
        'diagnostic_only': diagnostics,
        'scope': 'Excludes the unchanged unbroken C-odd branch from an Einstein soft '
                 'selfcoupling identification, conditionally on scattering qualifications. '
                 'Does not exclude symmetry-broken or C-even collective encodings, '
                 'new rules, finite-error effective alternatives, GR or dark sectors.',
        'new_adopted_axioms': 0, 'gr_hard_deliverables_completed': 0,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = check()
    raw = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    if args.write:
        path = HERE/'results.json'
        if path.exists() and json.loads(path.read_text(encoding='utf-8')) != result:
            raise RuntimeError('Refusing to replace a different result')
        if not path.exists():
            path.write_text(raw, encoding='utf-8')
    print(raw, end='')
