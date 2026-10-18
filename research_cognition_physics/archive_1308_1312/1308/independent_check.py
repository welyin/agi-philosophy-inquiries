"""Independent exact arithmetic/on-shell and numerical flavor cross-check of 1308.

Default is read-only. --write exclusively creates independent_results.json.
No author routines are imported; the analytic proof carries general quantifiers.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'independent_results.json'


def beta(nc, ng):
    # Gauge contribution + each Weyl doublet + one complex doublet.
    return -F(22, 3) + ng * (nc + 1) * F(1, 3) + F(1, 6)


def on_shell(nc):
    k = (nc - 1) // 2
    # There are k+1 legs along each of the two independent null directions.
    # The first of each group is incoming, the other k are outgoing.
    signed = [-k] + [1] * k
    momenta = [(e, 0, 0, s * e) for s in (-1, 1) for e in signed]
    assert all(p[0] ** 2 - sum(x*x for x in p[1:]) == 0 for p in momenta)
    assert tuple(sum(p[i] for p in momenta) for i in range(4)) == (0, 0, 0, 0)
    assert sum(p[0] < 0 for p in momenta) == 2
    # After matching the distinct colors to the external legs, the exterior
    # sign cancels epsilon_color's sign. The two antisymmetric contractions
    # in each pair have the same sign; the symmetric polynomial is
    # [2(x*w-y*z)]^k (x*Lw-y*Lz-z*Ly+w*Lx).
    # Its x^(k+1) w^k Lw component is positive, including assignments.
    component = (2 ** k) * factorial(k + 1) * factorial(k)
    # Square of its magnitude includes the two incoming energies k, with
    # spinors in independent directions; common normalization is omitted.
    magnitude_square = component ** 2 * k ** 2
    assert magnitude_square > 0
    return {'Nc': nc, 'incoming': 2, 'outgoing': 2*k,
            'signed_energy_per_null_direction': signed,
            'momentum_residual': [0, 0, 0, 0],
            'external_component': component,
            'unnormalized_squared_magnitude': magnitude_square}


def flavor(n):
    v = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    m = np.diag(np.arange(1, n+1)).astype(complex)
    m[:3, :3] = v @ np.diag([1, 2, 3]) @ v.conj().T
    # New spectator links differ from the author's rational rotations.
    for j in range(2, n-1):
        m[j, j+1] = m[j+1, j] = 0.0001
    eig = np.linalg.eigvalsh(m)
    assert eig[0] > 0.99 and np.min(np.diff(eig)) > 0.99
    h = m @ m
    phase = float(np.imag(h[0,1] * h[1,2] * h[2,0]))
    assert abs(phase) > 1
    visited = {0}
    while True:
        new = visited | {j for i in visited for j in range(n) if abs(h[i,j]) > 1e-12}
        if new == visited: break
        visited = new
    assert len(visited) == n
    return {'ng': n, 'positive_simple_mass': True, 'full_connected_flavor_graph': True,
            'abs_H_triangle_imag_greater_than_one': True,
            'beta_weak': str(beta(3, n))}


def run():
    for nc in range(3, 62):
        # Independently list (Y,B,L) in Q,u^c,d^c,L,e^c,N,H order.
        table = [(F(1,2*nc), F(1,nc), 0),
                 (F(-nc-1,2*nc), F(-1,nc), 0),
                 (F(nc-1,2*nc), F(-1,nc), 0),
                 (-F(1,2), 0, 1), (1, 0, -1), (0, 0, -1), (F(1,2), 0, 0)]
        assert [(nc+2)*b-l-4*y for y,b,l in table] == [1,1,-3,1,-3,1,-2]
    # Analytic P+W imply Nc<=6; R+E gives Nc=3 mod4. Hence finite exact
    # exhaustion below is already a complete integer check after that bound.
    candidates = [(nc,n) for nc in range(3,7) for n in range(3,6)
                  if (nc+1)%4 == 0 and beta(nc,n) <= 0]
    assert candidates == [(3,3),(3,4),(3,5)]
    phases = {'removed_E_5_3': (3*(5+1))%4,
              'retained_3_3': (3*(3+1))%4,
              'retained_3_4': (4*(3+1))%4,
              'retained_3_5': (5*(3+1))%4}
    assert phases == {'removed_E_5_3':2,'retained_3_3':0,'retained_3_4':0,'retained_3_5':0}
    files = [HERE/'proof.md', HERE/'check.py', HERE/'results.json',
             HERE.parent/'_admission/round1308_admission.md', Path(__file__).resolve()]
    return {'schema':'round1308_independent_v1', 'status':'PASS_SCOPED_CONDITIONAL_CHAIN',
            'survivors':candidates, 'on_shell_witnesses':[on_shell(n) for n in (3,5,7)],
            'flavor_checks':[flavor(n) for n in (3,4,5)],
            'SU2_instanton_r_exponent_mod4':phases,
            'E_independence_proved':False, 'new_axioms':0,
            'same_interacting_instrument_claim':False,
            'source_sha256':{p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in files}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert json.loads(OUT.read_text(encoding='utf-8')) == json.loads(json.dumps(result))
    print(json.dumps({'status':result['status'],'survivors':result['survivors'],
                      'on_shell_checks':3,'source_bindings':len(result['source_sha256'])}))
