"""Optional floating candidate discovery for four independent decoder angles.

Requires SciPy; none of its optimizer outputs alone is a rigorous bound.
The saved rational certificates are verified by separate default scripts.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

from joint_history_channel import I2, X, Y, Z
from soft_history_decoder import branch_weights, schur_root


MODIFIED = (0, 2, 3, 5)
P = np.array(branch_weights(), dtype=float)
SCALE_VECTOR = np.r_[np.sqrt(P), 1.]


def unitaries(angles):
    a, b, c, d = angles
    return np.array((np.kron(Z, np.kron(np.cos(a)*Z+np.sin(a)*X, I2)),
        np.kron(Z, np.eye(4)), np.kron(Z, np.kron(np.cos(b)*X+np.sin(b)*Z, I2)),
        np.kron(I2, np.kron(I2, np.cos(c)*Z+np.sin(c)*X)), np.eye(8),
        np.kron(I2, np.kron(I2, np.cos(d)*X+np.sin(d)*Z))))


TRUE = unitaries(np.zeros(4))
PAULI_PRODUCTS = np.array([[np.kron(a, b.conj()).real/4
                          for b in (I2, X, Y, Z)] for a in (I2, X, Y, Z)])


def kernel(angles):
    repair = unitaries(angles)
    ws = np.array([repair[g].T @ TRUE[r] for r in range(6) for g in range(6)]+[np.eye(8)])
    products = np.einsum('aki,bkj->abij', ws, ws).reshape(37, 37, 2, 4, 2, 4)
    partial = np.trace(products, axis1=2, axis2=4)
    bare = np.eye(4)[None, None]/2+(partial+partial.transpose(0, 1, 3, 2))/8
    return bare*SCALE_VECTOR[:, None, None, None]*SCALE_VECTOR[None, :, None, None], ws


def angle_gradient(angles, ws, sigma, coefficients, denominator):
    a = coefficients*np.sqrt(P)
    residual = np.eye(8)-np.einsum('j,jab->ab', a, ws[:-1])
    gradient = []
    for index, g in enumerate(MODIFIED):
        c, s = np.cos(angles[index]), np.sin(angles[index])
        local = -s*Z+c*X if g % 3 == 0 else -s*X+c*Z
        derivative = np.kron(Z, np.kron(local, I2)) if g < 3 else np.kron(I2, np.kron(I2, local))
        dr = -sum(a[6*r+g]*derivative.T @ TRUE[r] for r in range(6))
        dn = residual.T @ dr+dr.T @ residual
        partial = np.trace(dn.reshape(2, 4, 2, 4), axis1=0, axis2=2)
        gradient.append(np.trace(sigma @ partial)/(4*denominator))
    return np.array(gradient)


def relaxed_candidate(angles, start=None):
    from scipy.optimize import minimize
    basis, ws = kernel(angles)
    def objective(vector):
        norm = vector @ vector
        v = vector/np.sqrt(norm)
        sigma = np.outer(v, v)
        gram = np.einsum('abij,ji->ab', basis, sigma)
        value, c, d = schur_root(gram)
        z = np.r_[-c, 1.]
        h = np.einsum('a,b,abij->ij', z, z, basis)/d
        return -value, -2*(h-value*np.eye(4)) @ v/np.sqrt(norm)
    if start is None:
        start = np.array([.1577, .4126, .4887, -.7524])
    result = minimize(objective, start, method='BFGS', jac=True, options={'gtol':2e-11, 'maxiter':100})
    v = result.x/np.linalg.norm(result.x); sigma = np.outer(v, v)
    gram = np.einsum('abij,ji->ab', basis, sigma)
    value, c, d = schur_root(gram)
    z = np.r_[-c, 1.]
    h = np.einsum('a,b,abij->ij', z, z, basis)/d
    return {'value':value, 'upper_diagnostic':np.linalg.eigvalsh(h)[-1],
            'vector':v, 'coefficients':c*np.sqrt(P),
            'angle_gradient':angle_gradient(angles, ws, sigma, c, d)}


def product_candidate(angles, start=None):
    from scipy.optimize import minimize
    basis, ws = kernel(angles)
    gram_basis = np.einsum('abij,uvji->uvab', basis, PAULI_PRODUCTS)
    def jet(point):
        u, v = np.r_[1., point[:3]], np.r_[1., point[3:]]
        gram = np.einsum('i,j,ijab->ab', u, v, gram_basis)
        value, c, d = schur_root(gram)
        z = np.r_[-c, 1.]
        first = np.concatenate((np.einsum('j,ijab->iab', v, gram_basis[1:]),
                                np.einsum('i,ijab->jab', u, gram_basis[:, 1:])))
        gradient = np.einsum('a,iab,b->i', z, first, z)/d
        return value, gradient, c, d
    def objective(point):
        value, gradient, _, _ = jet(point)
        return -value, -gradient
    if start is None:
        from soft_decoder_certificate import witness_states
        start = np.array(witness_states(), dtype=float).flatten()
    constraints = [{'type':'ineq', 'fun':lambda p:1-p[:3] @ p[:3],
                    'jac':lambda p:np.r_[-2*p[:3], np.zeros(3)]},
                   {'type':'ineq', 'fun':lambda p:1-p[3:] @ p[3:],
                    'jac':lambda p:np.r_[np.zeros(3), -2*p[3:]]}]
    result = minimize(objective, start, method='SLSQP', jac=True,
                      constraints=constraints, options={'ftol':2e-14,'maxiter':150})
    point = result.x.copy()
    point[:3] /= max(1., np.linalg.norm(point[:3])); point[3:] /= max(1., np.linalg.norm(point[3:]))
    value, _, c, d = jet(point)
    sigma = np.einsum('i,j,ijab->ab', np.r_[1., point[:3]], np.r_[1., point[3:]], PAULI_PRODUCTS)
    return {'value':value, 'point':point, 'coefficients':c*np.sqrt(P),
            'angle_gradient':angle_gradient(angles, ws, sigma, c, d)}


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'.research_runtime'))
    from scipy.optimize import minimize
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--product', action='store_true')
    args = parser.parse_args()
    def objective(angles):
        result = finder(angles)
        return result['value'], result['angle_gradient']
    finder = product_candidate if args.product else relaxed_candidate
    start = np.full(4, 2*np.arctan(3/40))
    print('start', finder(start), flush=True)
    result = minimize(objective, start, jac=True, method='L-BFGS-B',
                      bounds=[(0, np.pi/2)]*4, options={'gtol':1e-10, 'ftol':1e-15, 'maxiter':100})
    candidate = finder(result.x)
    report = {'angles':result.x.tolist(), 'tan_half_angles':np.tan(result.x/2).tolist(),
              'error_diagnostic':float(candidate['value']),
              'input_domain':'product' if args.product else 'all real relaxed',
              'vector_or_point':candidate.get('vector',candidate.get('point')).tolist(),
              'coefficients':candidate['coefficients'].tolist(),
              'gradient':candidate['angle_gradient'].tolist(),
              'optimizer_message':str(result.message), 'certified':False}
    if not args.product:
        report['relaxed_upper_diagnostic'] = float(candidate['upper_diagnostic'])
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
