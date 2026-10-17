"""Independent 1246 checks: exact quadratic-field Taylor identities and
constraint-reduced Decimal integration in the scale factor a.

The Decimal trajectories are diagnostics, not interval certificates.  The
finite enclosure is the author's separately reviewed Picard/Taylor proof.
Default execution is read-only; --write exclusively creates the result.
"""
from __future__ import annotations

import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULT = HERE / 'independent_results.json'
RAD = F(219, 200)
ORDER = 10


class Q:
    """Exact a+b*sqrt(219/200), independently implemented arithmetic."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        if isinstance(a, Q):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, other):
        other = Q(other)
        return Q(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-Q(other))

    def __rsub__(self, other):
        return Q(other) - self

    def __mul__(self, other):
        other = Q(other)
        return Q(self.a*other.a + RAD*self.b*other.b,
                 self.a*other.b + self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = F(other)
        return Q(self.a/other, self.b/other)

    def __eq__(self, other):
        other = Q(other)
        return self.a == other.a and self.b == other.b

    def pack(self):
        return {'rational': str(self.a), 'sqrt_219_200_coefficient': str(self.b)}


class P:
    """Truncated formal polynomial; no numerical interval library imported."""
    __slots__ = ('c',)

    def __init__(self, data=0):
        if isinstance(data, P):
            self.c = list(data.c)
        elif isinstance(data, (tuple, list)):
            self.c = [Q(a) for a in data[:ORDER+1]]
            self.c += [Q()]*(ORDER+1-len(self.c))
        else:
            self.c = [Q(data)] + [Q()]*ORDER

    def __add__(self, other):
        other = P(other)
        return P([a+b for a, b in zip(self.c, other.c)])

    __radd__ = __add__

    def __neg__(self):
        return P([-a for a in self.c])

    def __sub__(self, other):
        return self + (-P(other))

    def __rsub__(self, other):
        return P(other) - self

    def __mul__(self, other):
        other = P(other)
        out = [Q()]*(ORDER+1)
        for n in range(ORDER+1):
            out[n] = sum((self.c[k]*other.c[n-k] for k in range(n+1)), Q())
        return P(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return P([a/other for a in self.c])

    def __pow__(self, n):
        ans = P(1)
        for _ in range(n):
            ans = ans*self
        return ans

    def derivative(self):
        return P([self.c[k]*k for k in range(1, ORDER+1)])


FNAMES = ('u', 'H', 'p', 'v', 'J', 'q', 'r', 'z', 'L', 'psi')
TIME = P([0, 1])
W = 64*TIME**3 - 192*TIME**4 + 192*TIME**5 - 64*TIME**6
WP = W.derivative()


def formal_energy(y, alpha):
    u, H, p, v, J, q, r, z, light, psi = y
    mass = P(F(7, 10)) + J + (q*q+r*r)/2 + alpha*W*q*p/4
    dust = (mass+F(1, 50)+p/100)*u**3
    rods = F(3, 20)*(1+p/5)*u*u
    rho = (v*v+p*p)/2 + rods + dust
    rho_plus_pressure = v*v + rods*F(2, 3) + dust
    return mass, rho, rho_plus_pressure


def formal_rhs(y, alpha):
    u, H, p, v, J, q, r, z, light, psi = y
    mass, rho, rho_plus_pressure = formal_energy(y, alpha)
    return [-H*u, -F(3, 2)*rho_plus_pressure, v,
            -3*H*v-p-F(3, 100)*u*u-(F(1, 100)+alpha*W*q/4)*u**3,
            -alpha*WP*q*p/4, r, -q-alpha*W*p/4,
            p/10, u, W*q*p/4]


def formal_branch(alpha):
    y = [P(a) for a in (1, Q(0, 1), 0, 0, F(1, 10), F(1, 2), 0, 0, 0, 0)]
    for n in range(ORDER):
        rhs = formal_rhs(y, alpha)
        for state, rate in zip(y, rhs):
            state.c[n+1] = rate.c[n]/(n+1)
    mass, rho, rho_plus_pressure = formal_energy(y, alpha)
    by_name = dict(zip(FNAMES, y))
    mass_balance = mass.derivative()-alpha*W*by_name['q']*by_name['v']/4
    continuity = rho.derivative()+3*by_name['H']*rho_plus_pressure
    constraint = by_name['H']**2-rho
    assert all(c == 0 for c in mass_balance.c[:ORDER])
    assert all(c == 0 for c in continuity.c[:ORDER])
    assert all(c == 0 for c in constraint.c)
    return by_name | {'M_c': mass, 'rho': rho,
                      'Brod': by_name['u']**2/10}


def formal_group():
    off, on = formal_branch(0), formal_branch(1)
    differences = {key: on[key]-off[key] for key in off}
    expected = {'p': (5, Q(-F(2, 5))),
                'M_c': (5, Q(-F(8, 125))),
                'rho': (6, Q(0, -F(21, 500))),
                'H': (6, Q(-F(21, 1000))),
                'u': (7, Q(F(3, 1000))),
                'Brod': (7, Q(F(3, 5000))),
                'L': (8, Q(F(3, 8000))),
                'z': (6, Q(-F(1, 150)))}
    leading = {}
    for name, (power, value) in expected.items():
        assert all(c == 0 for c in differences[name].c[:power]), name
        assert differences[name].c[power] == value, name
        leading[name] = dict(power=power, coefficient=value.pack())
    assert off['rho'].c[0] == Q(RAD) == on['rho'].c[0]
    assert off['M_c'].c[0] == Q(F(37, 40)) == on['M_c'].c[0]
    assert differences['rho'].c[5] == 0
    assert on['psi'].c[6] != 0 and off['psi'].c[6] != 0
    # A pointwise force deletion identity, not a changed-model trajectory.
    t, q, p, v, u = F(1, 4), F(2, 5), -F(1, 50), -F(1, 40), F(4, 5)
    w = 64*t**3*(1-t)**3
    force_deleted_continuity_defect = w*q*v*u**3/4
    assert force_deleted_continuity_defect == -F(27, 50000)
    return dict(status='PASS', method='exact Q(sqrt(219/200)) formal series; order 10',
                initial_mass='37/40', initial_rho='219/200',
                complete_rho_difference_coefficient_order5='0',
                first_nonzero_coefficients=leading,
                complete_mass_exchange_residual_coefficients_through_order9='0',
                complete_continuity_residual_coefficients_through_order9='0',
                constraint_residual_coefficients_through_order10='0',
                program_phase_preserved_in_both_branches=True,
                missing_pointer_force_defect=str(force_deleted_continuity_defect),
                claim='Exact local coefficients and conserved-source identities, not a finite-time enclosure.')


DNAMES = ('a', 'p', 'v', 'J', 'q', 'r', 'z', 'L', 'psi')


def decimal_energy(y, t, alpha):
    a, p, v, J, q, r, z, light, psi = y
    w = 64*t**3*(1-t)**3
    mass = D('.7')+J+(q*q+r*r)/2+D(alpha)*w*q*p/4
    rods = D('.15')*(1+p/5)/a**2
    dust = (mass+D('.02')+p/100)/a**3
    rho = (v*v+p*p)/2+rods+dust
    H = rho.sqrt()
    return mass, rods, dust, rho, H


def decimal_rhs(t, y, alpha):
    a, p, v, J, q, r, z, light, psi = y
    w = 64*t**3*(1-t)**3
    wp = 192*t*t*(1-t)**2*(1-2*t)
    mass, rods, dust, rho, H = decimal_energy(y, t, alpha)
    # H is reconstructed from the full constraint, not integrated via H'.
    return [a*H, v, -3*H*v-p-D('.03')/a**2-(D('.01')+D(alpha)*w*q/4)/a**3,
            -D(alpha)*wp*q*p/4, r, -q-D(alpha)*w*p/4,
            p/10, 1/a, w*q*p/4]


def decimal_branch(alpha, steps):
    y = list(map(D, ('1', '0', '0', '.1', '.5', '0', '0', '0', '0')))
    dt = D(1)/steps
    minima = {'M_c': D(10), 'J': D(10), 'F': D(10), 'M_r': D(10)}
    derivative_residual = D(0)
    for n in range(steps):
        t = D(n)/steps
        k1 = decimal_rhs(t, y, alpha)
        k2 = decimal_rhs(t+dt/2, [a+dt*b/2 for a, b in zip(y, k1)], alpha)
        k3 = decimal_rhs(t+dt/2, [a+dt*b/2 for a, b in zip(y, k2)], alpha)
        k4 = decimal_rhs(t+dt, [a+dt*b for a, b in zip(y, k3)], alpha)
        # An independent chain rule for H=sqrt(rho) checks Raychaudhuri.
        a, p, v, J, q, r, z, light, psi = y
        ap, pp, vp, Jp, qp, rp, zp, Lp, psip = k1
        mass, rods, dust, rho, H = decimal_energy(y, t, alpha)
        w = 64*t**3*(1-t)**3
        wp = 192*t*t*(1-t)**2*(1-2*t)
        mdot = Jp+q*qp+r*rp+D(alpha)*(wp*q*p+w*qp*p+w*q*pp)/4
        rhodot = v*vp+p*pp+D('.03')*pp/a**2-2*H*rods
        rhodot += (mdot+pp/100)/a**3-3*H*dust
        residual = rhodot/(2*H)+D('1.5')*(v*v+D(2)*rods/3+dust)
        derivative_residual = max(derivative_residual, abs(residual))
        for name, val in {'M_c': mass, 'J': J, 'F': 1+p/5, 'M_r': D('.2')+p/10}.items():
            minima[name] = min(minima[name], val)
        y = [a+dt*(b+2*c+2*d+e)/6 for a, b, c, d, e in zip(y, k1, k2, k3, k4)]
    mass, rods, dust, rho, H = decimal_energy(y, D(1), alpha)
    assert min(minima.values()) > 0
    assert derivative_residual < D('1e-50')
    result = dict(zip(DNAMES, y)) | {'u': 1/y[0], 'H': H, 'Brod': D('.1')/y[0]**2, 'M_c': mass}
    return result, minima, derivative_residual


def decimal_group():
    rows = []
    results = []
    with localcontext() as ctx:
        ctx.prec = 60
        for steps in (512, 1024, 2048):
            off, min0, res0 = decimal_branch(0, steps)
            on, min1, res1 = decimal_branch(1, steps)
            delta = {key: on[key]-off[key] for key in ('Brod', 'L', 'z', 'M_c', 'u', 'H')}
            rows.append(dict(steps=steps, off={k: str(v) for k, v in off.items()},
                             on={k: str(v) for k, v in on.items()},
                             difference={k: str(v) for k, v in delta.items()},
                             sampled_minima_off={k: str(v) for k, v in min0.items()},
                             sampled_minima_on={k: str(v) for k, v in min1.items()},
                             chain_rule_residual_max=str(max(res0, res1))))
            results.append((off, on, delta))
        fine = results[-1][2]
        convergence = {}
        for key in fine:
            e1 = abs(results[1][2][key]-results[0][2][key])
            e2 = abs(results[2][2][key]-results[1][2][key])
            assert e2 < e1/8 and e2 < D('1e-11'), key
            assert abs(fine[key]) > 100*e2, key
            convergence[key] = {'coarse_difference': str(e1), 'fine_difference': str(e2),
                                'ratio': str(e1/e2)}
        assert fine['Brod'] > 0 and fine['L'] > 0 and fine['z'] < 0
    return dict(status='PASS', method='60-digit Decimal RK4 in a with H=sqrt(full rho)',
                T='1', rows=rows, step_halving_diagnostics=convergence,
                distinguishes_geometry_from_direct_pointer_channel=True,
                no_enclosure_claim=True,
                sampled_positivity_is_not_an_all_time_certificate=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(diagnose=False):
    result = dict(round=1246, status='PASS', scientific_groups=2,
                  formal_group=formal_group(), decimal_group=decimal_group(),
                  actual_quantum_readout_certified=False,
                  CT_adopted=False, whole_GR_proved=False)
    if not diagnose:
        files = ('admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md')
        result['author_sha256'] = {name: sha(HERE/name) for name in files}
        result['independent_code_sha256'] = sha(Path(__file__))
        author = json.loads((HERE/'results.json').read_text('utf8'))
        assert author['status'] == 'PASS'
        comparisons = {}
        with localcontext() as ctx:
            ctx.prec = 60
            fine = result['decimal_group']['rows'][-1]['difference']
            mesh = result['decimal_group']['step_halving_diagnostics']
            for name, author_key in (('Brod', 'delta_Brod'), ('L', 'delta_L'), ('z', 'delta_z')):
                lo, hi = map(D, author['group_finite_geometry'][author_key])
                value = D(fine[name])
                distance = max(lo-value, value-hi, D(0))
                last_mesh_change = D(mesh[name]['fine_difference'])
                # A consistency diagnostic only: the mesh difference is not a
                # certified bound on the Decimal integration's truncation error.
                assert distance < 2*last_mesh_change
                comparisons[name] = dict(author_enclosure=[str(lo), str(hi)],
                    Decimal_distance_from_enclosure=str(distance),
                    last_mesh_change=str(last_mesh_change),
                    diagnostic_tolerance=str(2*last_mesh_change))
        result['author_certificate_comparison'] = {
            'author_status': author['status'],
            'author_interval_method_reviewed_separately': True,
            'independent_Decimal_convergence_not_rigorous_error_bound': True,
            'finite_signal_consistency_diagnostics': comparisons}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--diagnose', action='store_true')
    args = parser.parse_args()
    assert not (args.write and args.diagnose)
    value = run(args.diagnose)
    if args.write:
        with RESULT.open('x', encoding='utf8') as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        mode = 'INITIAL_SAVE'
    elif args.diagnose:
        mode = 'DIAGNOSTIC_ONLY'
    else:
        assert json.loads(RESULT.read_text('utf8')) == value, 'read-only result mismatch'
        mode = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status=value['status'], mode=mode, scientific_groups=2,
                          formal_leading=value['formal_group']['first_nonzero_coefficients'],
                          finite_difference=value['decimal_group']['rows'][-1]['difference']),
                     ensure_ascii=False))
