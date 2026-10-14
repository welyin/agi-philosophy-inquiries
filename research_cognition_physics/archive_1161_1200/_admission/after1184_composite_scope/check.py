"""Finite algebraic admission check; not a GR realization or a research round.

The four-dimensional mediator replaces the proposed unbounded oscillator.
No truncated Fock simulation, source instrument, or geometry is inferred.
"""
from pathlib import Path
import argparse
import json
import numpy as np

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
A = (I + X) / 2
n = (I - Z) / 2
rot = np.array([[np.sqrt(2/3), np.sqrt(1/3)],
                [-np.sqrt(1/3), np.sqrt(2/3)]])
W = np.kron(I-A, I) + np.kron(A, rot)
E0 = np.kron(I, np.array([[1.], [0.]]))
E = W @ E0
Hm = np.eye(4) + np.kron(I, n) + np.sqrt(2)*np.kron(A, X)
Hc = I-A
expected = np.kron(Hc, I) + np.kron(I+2*A, n)
norm = lambda a: float(np.linalg.norm(a, ord=2))
checks = {
    'unitary_residual': norm(W.T@W-np.eye(4)),
    'full_diagonalization_residual': norm(W.T@Hm@W-expected),
    'encoded_generator_residual': norm(Hm@E-E@Hc),
    'transported_port_invariance_residual': norm((W@np.kron(Z,I)@W.T)@E-E@Z),
    'bare_port_compression_residual': norm(E.T@np.kron(Z,I)@E-np.sqrt(2/3)*Z),
    'bare_port_leakage': norm((np.eye(4)-E@E.T)@np.kron(Z,I)@E),
    'bare_port_leakage_exact': float(1/np.sqrt(3)),
}
rows = []
for alpha,beta in [(1.,1.), (0.,0.), (0.,.75), (5.,2.)]:
    J = np.eye(4)+alpha*np.kron(I,n)+np.sqrt(2)*beta*np.kron(A,X)
    expected_J = I+(alpha-4*beta)*A/3
    leak = (np.eye(4)-E@E.T)@J@E
    rows.append({'alpha':alpha, 'beta':beta,
                 'compressed_source_residual':norm(E.T@J@E-expected_J),
                 'source_leakage_identity_residual':norm(leak.T@leak-2*(alpha-beta)**2*A/9),
                 'source_leakage_norm':norm(leak),
                 'source_minus_energy_norm':norm(expected_J-Hc),
                 'soft_Y_initial_derivative_A0':.3,
                 'soft_Y_initial_derivative_A1':.3*(1+(alpha-4*beta)/3)})
assert max(v for k,v in checks.items() if k.endswith('residual')) < 1e-13
assert abs(checks['bare_port_leakage']-checks['bare_port_leakage_exact']) < 1e-13
assert max(row['compressed_source_residual'] for row in rows) < 1e-13
assert max(row['source_leakage_identity_residual'] for row in rows) < 1e-13
out = {
    'status':'ADMISSION_ALGEBRA_ONLY_NO_NEW_ROUND',
    'full_dimension':4, 'encoded_dimension':2,
    'full_spectrum':np.linalg.eigvalsh(Hm).tolist(),
    'checks':checks, 'source_rows':rows,
    'analytic_source_condition':'4 beta - alpha = 3 (encoded initial response only)',
    'limits':['No spatial or metric dictionary for a newly introduced soft spin',
              'Complete transport into an old parent is a further input',
              'A transported exact embedding alone selects no new effective dynamics',
              'Floating checks support the exact identities in README, not a GR certificate']
}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--write', action='store_true', help='explicitly refresh this mutable preview result')
args = parser.parse_args()
if args.write:
    Path(__file__).with_name('results.json').write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,ensure_ascii=False,indent=2))
