"""Read-only finite check for task_selection_after1063.md, not a new round.

--save-exclusive creates the adjacent result once; default only prints.
No physical instrument preparation or resource implementation is claimed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np


def swap(d):
    value = np.zeros((d*d, d*d), complex)
    for i in range(d):
        for j in range(d):
            value[j*d+i, i*d+j] = 1
    return value


def trace_qa(value, d):
    return np.einsum('ijkijl->kl', value.reshape((d,)*6))


def check_dimension(d):
    eye = np.eye(d)
    sw = swap(d)
    minus = (np.eye(d*d)-sw)/2
    plus = (np.eye(d*d)+sw)/2
    resource = 2*minus/(d*(d-1))
    projects = [np.kron(minus, eye), np.kron(plus, eye)]
    # Generalized Weyl realization of the covariant CP endpoint. Its individual
    # Kraus operators need not be invariant; only the whole channel is covariant.
    omega = np.exp(2j*np.pi/d)
    ks = []
    for a in range(d):
        for b in range(d):
            if a == b == 0:
                continue
            w = np.zeros((d, d), complex)
            for j in range(d):
                w[(j+a) % d, j] = omega**(b*j)
            ks.append(w/np.sqrt(d*d-1))

    def endpoint(x):
        return sum((k@x@k.conj().T for k in ks), np.zeros_like(x))

    norm_residual = np.max(np.abs(sum(k.conj().T@k for k in ks)-eye))
    lam = -1/(d*d-1)
    eta = d/(2*(d-1)**2*(d+1))
    choi = np.zeros((d*d, d*d), complex)
    branch_residual = endpoint_residual = full_residual = 0.0
    for i in range(d):
        for j in range(d):
            unit = np.zeros((d,d), complex)
            unit[i,j] = 1
            joint = np.kron(unit, resource)
            ns = [trace_qa(p@joint@p, d) for p in projects]
            expected = [((d-2)*np.trace(unit)*eye+unit)/(2*d*(d-1)),
                        (d*np.trace(unit)*eye-unit)/(2*d*(d-1))]
            branch_residual = max(branch_residual, *(float(np.max(np.abs(x-y)))
                                                    for x,y in zip(ns, expected)))
            endpoint_residual = max(endpoint_residual, float(np.max(np.abs(
                endpoint(unit)-(lam*unit+(1-lam)*np.trace(unit)*eye/d)))))
            output = ns[0]+endpoint(ns[1])
            full_residual = max(full_residual, float(np.max(np.abs(
                output-(eta*unit+(1-eta)*np.trace(unit)*eye/d)))))
            choi += np.kron(unit, output)/d
    phi = np.eye(d).reshape(-1)/np.sqrt(d)
    fe = float(np.vdot(phi, choi@phi).real)
    expected_fe = (1+(d*d-1)*eta)/(d*d)
    # Finite vertices of the algebraic optimization: sharp P± with endpoints
    # lambda_- and lambda_+. General POVMs follow analytically from a,b>=0.
    vertex_eta = [(x-y)/(2*d*(d-1)) for x in [lam,1] for y in [lam,1]]
    assert abs(max(vertex_eta)-eta) < 1e-13
    assert abs(fe-expected_fe) < 1e-13
    assert np.linalg.eigvalsh(choi).min() > -1e-13
    trace_b = np.einsum('ijkj->ik', choi.reshape(d,d,d,d))
    tp_residual = float(np.max(np.abs(trace_b-eye/d)))
    residual = max(branch_residual, endpoint_residual, full_residual,
                   float(norm_residual), tp_residual)
    assert residual < 1e-12
    if d == 2:
        # Heralded success is distinct from the deterministic average channel.
        assert np.linalg.matrix_rank(minus, tol=1e-12) == 1
        assert abs(fe-0.5) < 1e-13
    return dict(d=d, eta_max=eta, entanglement_fidelity=fe,
                pure_input_fidelity=(1+(d-1)*eta)/d,
                branch_map_residual=branch_residual,
                correction_channel_residual=endpoint_residual,
                full_channel_residual=full_residual,
                cp_tp_residual=max(float(norm_residual),tp_residual),
                choi_minimum_eigenvalue=float(np.linalg.eigvalsh(choi).min()),
                fixed_contract_only=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--save-exclusive', action='store_true')
    args = parser.parse_args()
    result = dict(passed=True, formal_round=None, new_scientific_groups=0,
                  new_accepted_cognitive_axioms=0,
                  scope='One-way scalar messages, invariant QA effects, all covariant B channels, fixed antisymmetric resource; deterministic B output only.',
                  dimensions=[check_dimension(d) for d in range(2,6)])
    output = json.dumps(result, ensure_ascii=False, indent=2)
    if args.save_exclusive:
        with Path(__file__).with_name('scalar_handoff_results.json').open('x', encoding='utf-8') as f:
            f.write(output+'\n')
    print(output)


if __name__ == '__main__':
    main()
