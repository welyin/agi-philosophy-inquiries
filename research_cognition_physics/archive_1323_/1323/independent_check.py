"""Independent finite-Fock checks for round 1323. No author code is imported.

Default/--check recomputes and compares every saved field; missing results fail.
--write-new exclusively creates independent_results.json after all checks pass.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / 'independent_results.json'
ADMISSION = HERE.parent / '_admission/round1323_admission.md'
TOL = 3e-10


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def occupations(n):
    return [(k, n-k) for k in range(n+1)] if n >= 0 else []


def fermion(mask, mode, create=False):
    occupied = bool(mask & (1 << mode))
    if occupied == create:
        return None
    phase = -1 if (mask & ((1 << mode)-1)).bit_count() % 2 else 1
    return mask ^ (1 << mode), phase


def exponential(h, t):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j*t*vals)) @ vecs.conj().T


def spin(n):
    basis = occupations(n)
    size = n+1
    raising = np.zeros((size, size), complex)
    for col, (a, b) in enumerate(basis):
        if b:
            raising[col+1, col] = np.sqrt((a+1)*b)
    return [(raising+raising.conj().T)/2,
            (raising-raising.conj().T)/(2j),
            np.diag([(a-b)/2 for a, b in basis]).astype(complex)]


def extract_direction(n):
    """Analytic symmetric-tensor map, constructed separately from parent CAR."""
    out = [(a, b) for a in range(2) for b in occupations(n-1)]
    lookup = {state: i for i, state in enumerate(out)}
    w = np.zeros((2*n, n+1), complex)
    for col, occ in enumerate(occupations(n)):
        for a in range(2):
            if occ[a]:
                lowered = list(occ)
                lowered[a] -= 1
                w[lookup[a, tuple(lowered)], col] = np.sqrt(occ[a]/n)
    return w


def parent(nref, flavors, matrix):
    # Single actual fermion: right modes 0..f-1, left modes f..3f-1.
    f = flavors
    basis = [(1 << mode, *occ) for mode in range(3*f)
             for occ in occupations(nref-(mode >= f))]
    index = {state: k for k, state in enumerate(basis)}
    right = [k for k, (mask, _, _) in enumerate(basis)
             if (mask.bit_length()-1) < f]
    left = [k for k in range(len(basis)) if k not in right]
    t = np.zeros((len(basis), len(basis)), complex)
    for col, (mask, b0, b1) in enumerate(basis):
        for j in range(f):
            first = fermion(mask, j)
            if first is None:
                continue
            for a in range(2):
                occ = [b0, b1]
                if not occ[a]:
                    continue
                factor = np.sqrt(occ[a])
                occ[a] -= 1
                for i in range(f):
                    second = fermion(first[0], f+a*f+i, create=True)
                    if second is None:
                        continue
                    target = (second[0], *occ)
                    if target not in index:
                        raise AssertionError('A reachable state was omitted')
                    t[index[target], col] += first[1]*second[1]*factor*matrix[i, j]
    v = t+t.conj().T
    bare = v.copy()
    if nref:
        v /= np.sqrt(nref)
    pr = np.eye(len(basis), dtype=complex)[:, right]
    pl = np.eye(len(basis), dtype=complex)[:, left]
    nr = pr @ pr.conj().T
    nl = pl @ pl.conj().T
    nb = np.diag([b0+b1 for _, b0, b1 in basis]).astype(complex)
    # Independent direction/flavor intertwiner with the actual occupancy order.
    w = np.zeros((len(left), len(right)), complex)
    li = {basis[k]: row for row, k in enumerate(left)}
    if nref:
        for col, k in enumerate(right):
            mask, b0, b1 = basis[k]
            j = mask.bit_length()-1
            for a, count in enumerate((b0, b1)):
                if count:
                    occ = [b0, b1]
                    occ[a] -= 1
                    w[li[(1 << (f+a*f+j), *occ)], col] = np.sqrt(count/nref)
    return dict(basis=basis, v=v, bare=bare, pr=pr, pl=pl, w=w,
                nr=nr, nl=nl, nb=nb)


def two_parent(nref, active):
    # Each of two actual participants has one right and two left CAR modes.
    modes = [(a, b) for a in range(3) for b in range(3, 6)]
    basis = [(1 << a | 1 << b, *occ) for a, b in modes
             for occ in occupations(nref-(a != 0)-(b != 3))]
    index = {state: k for k, state in enumerate(basis)}
    t = np.zeros((len(basis), len(basis)), complex)
    rmode, lmodes = (0, [1, 2]) if active == 0 else (3, [4, 5])
    for col, (mask, b0, b1) in enumerate(basis):
        first = fermion(mask, rmode)
        if first is None:
            continue
        # This participant is right initially, hence N_b+n_L(active)=N_b.
        total = b0+b1
        for a, lmode in enumerate(lmodes):
            occ = [b0, b1]
            if not occ[a]:
                continue
            factor = np.sqrt(occ[a]/total)
            occ[a] -= 1
            second = fermion(first[0], lmode, create=True)
            target = (second[0], *occ)
            if target not in index:
                raise AssertionError('Sequential reachable state omitted')
            t[index[target], col] += first[1]*second[1]*factor
    return basis, t+t.conj().T


class Checks:
    def __init__(self):
        self.count = 0
        self.residuals = {}

    def near(self, name, a, b=0, tol=TOL):
        err = float(np.max(np.abs(np.asarray(a)-np.asarray(b)), initial=0))
        if not np.isfinite(err) or err > tol:
            raise AssertionError((name, err, tol))
        self.count += 1
        self.residuals[name] = err

    def truth(self, name, assertion):
        if not assertion:
            raise AssertionError(name)
        self.count += 1


def compute():
    c = Checks()
    sx = np.array([[0, 1], [1, 0]], complex)/2
    sy = np.array([[0, -1j], [1j, 0]], complex)/2
    sz = np.diag([.5, -.5]).astype(complex)
    direction_data = []
    for n in range(1, 5):
        w = extract_direction(n)
        js = spin(n)
        jo = [np.kron(s, np.eye(n))+np.kron(np.eye(2), j)
              for s, j in zip([sx, sy, sz], spin(n-1))]
        c.near(f'N{n}.W_isometry', w.conj().T @ w, np.eye(n+1))
        constraints = []
        for a, (left, right) in enumerate(zip(jo, js)):
            c.near(f'N{n}.Lie{a}', left @ w-w @ right)
            constraints.append(np.kron(np.eye(n+1), left)
                               -np.kron(right.T, np.eye(2*n)))
        singular = np.linalg.svd(np.vstack(constraints), compute_uv=False)
        nullity = int(np.sum(singular < 1e-10))
        c.truth(f'N{n}.Hom_dimension', nullity == 1)
        axis = np.array([.31, -.42, .57])
        ui = exponential(sum(x*j for x, j in zip(axis, js)), .73)
        uo = exponential(sum(x*j for x, j in zip(axis, jo)), .73)
        c.near(f'N{n}.finite_SU2', uo @ w, w @ ui)
        bound = 1/(n+1)
        c.truth(f'N{n}.exact_capacity', 2*(n+1) > 2*n)
        c.near(f'N{n}.capacity_bound', 1-2*n/(2*(n+1)), bound)
        c.truth(f'N{n}.center_obstruction', (-1)**(n+1) != (-1)**n)
        direction_data.append(dict(N=n, direction_input=n+1, direction_output=2*n,
                                   intertwiner_dimension=nullity,
                                   right_doublet_error_lower_bound=bound,
                                   next_singular_value=float(sorted(singular)[1])))

    gamma = .7
    clock = np.pi/(2*gamma)
    single = []
    for f in (3, 4, 5):
        for n in range(5):
            p = parent(n, f, gamma*np.eye(f))
            v, pr, pl, w = [p[k] for k in ('v', 'pr', 'pl', 'w')]
            d = len(p['basis'])
            c.truth(f'f{f}.N{n}.dimension', d == (3*n+1)*f)
            c.near(f'f{f}.N{n}.Hermitian', v, v.conj().T)
            c.near(f'f{f}.N{n}.particle_total', p['nr']+p['nl'], np.eye(d))
            c.near(f'f{f}.N{n}.conserved_N', p['nb']+p['nl'], n*np.eye(d))
            er, omega = 1.3, .4
            h0 = er*p['nr']+(er+omega)*p['nl']+omega*p['nb']
            c.near(f'f{f}.N{n}.free_energy_shell', h0, (er+omega*n)*np.eye(d))
            c.near(f'f{f}.N{n}.energy_commutator', h0 @ v-v @ h0)
            if n == 0:
                c.near(f'f{f}.N0.no_extract', v)
                c.near(f'f{f}.N0.no_left', p['nl'])
                continue
            c.near(f'f{f}.N{n}.compressed_creation', pl.conj().T @ v @ pr, gamma*w)
            c.near(f'f{f}.N{n}.code_isometry', w.conj().T @ w, np.eye(len(pr.T)))
            dark = p['nl']-pl @ w @ w.conj().T @ pl.conj().T
            c.near(f'f{f}.N{n}.dark_projector', dark @ dark, dark)
            c.near(f'f{f}.N{n}.dark_rank', np.trace(dark), f*(n-1))
            c.near(f'f{f}.N{n}.dark_no_coupling', v @ dark)
            un = exponential(v, clock)
            ub = exponential(p['bare'], clock)
            c.near(f'f{f}.N{n}.perfect_full_input', un @ pr, -1j*pl @ w)
            pn = float(np.linalg.norm(pl.conj().T @ un @ pr, 'fro')**2/pr.shape[1])
            pb = float(np.linalg.norm(pl.conj().T @ ub @ pr, 'fro')**2/pr.shape[1])
            c.near(f'f{f}.N{n}.bare_probability', pb, np.sin(gamma*np.sqrt(n)*clock)**2)
            if n == 4:
                c.near(f'f{f}.N4.bare_failure', pb)
                c.near(f'f{f}.N4.normalized_success', pn, 1)
            # Maximally entangled with a partner of full input dimension.
            choi = un @ pr/np.sqrt(pr.shape[1])
            ideal = -1j*pl @ w/np.sqrt(pr.shape[1])
            c.near(f'f{f}.N{n}.full_partner_vector', choi, ideal)
            # All short times use the same compressed conversion table.
            for k, time in enumerate((.017, .13, .39)):
                evolved = exponential(v, time) @ pr
                expected = np.cos(gamma*time)*pr-1j*np.sin(gamma*time)*pl @ w
                c.near(f'f{f}.N{n}.short{k}', evolved, expected)
            c.near(f'f{f}.N{n}.right_short_generator', pl.conj().T @ v @ pr/gamma, w)
            c.near(f'f{f}.N{n}.bare_short_generator', pl.conj().T @ p['bare'] @ pr/gamma,
                   np.sqrt(n)*w)
            # Same endpoint including phase does not fix a short-time generator.
            alias = 4*(n-1)+1
            c.near(f'f{f}.N{n}.endpoint_alias_same_phase',
                   exponential(alias*v, clock) @ pr, -1j*pl @ w)
            if n > 1:
                c.truth(f'f{f}.N{n}.endpoint_alias_distinct_generator',
                        np.linalg.norm((alias-1)*pl.conj().T @ v @ pr) > 1)
            # Same complete parent and arbitrary unknown input/partner preparation.
            rng = np.random.default_rng(1000+10*f+n)
            z = rng.normal(size=(pr.shape[1], 3))+1j*rng.normal(size=(pr.shape[1], 3))
            z /= np.linalg.norm(z)
            t = .47
            psi0, psit = pr @ z, exponential(h0+v, t) @ pr @ z
            energy = lambda a: float(np.vdot(a, (h0+v) @ a).real)
            c.near(f'f{f}.N{n}.complete_energy', energy(psit), energy(psi0))
            pop = float(np.vdot(psit, p['nl'] @ psit).real)
            c.near(f'f{f}.N{n}.expected_source_transfer', pop, np.sin(gamma*t)**2)
            c.near(f'f{f}.N{n}.reference_count', np.vdot(psit, p['nb'] @ psit), n-pop)
            c.near(f'f{f}.N{n}.right_count', np.vdot(psit, p['nr'] @ psit), 1-pop)
            # Free energies remain independent sources before imposing resonance.
            de, dw = .21, -.17
            source_delta = de*p['nl']+dw*p['nb']
            new_h0 = er*p['nr']+(er+omega+de)*p['nl']+(omega+dw)*p['nb']
            c.near(f'f{f}.N{n}.independent_source_identity', new_h0-h0, source_delta)
            single.append(dict(N=n, flavors=f, parent_dimension=d,
                               geometric_dark_dimension=f*(n-1),
                               normalized_success=pn, bare_success=pb))

    general = []
    for f in (3, 4, 5):
        j = np.arange(f)
        fourier = np.exp(2j*np.pi*np.outer(j, j)/f)/np.sqrt(f)
        matrix = fourier @ np.diag(np.linspace(.24, 1.13, f)) @ np.diag(np.exp(.17j*j))
        um, sm, vh = np.linalg.svd(matrix)
        for n in (1, 2, 4):
            p = parent(n, f, matrix)
            v, pr, pl, w = [p[k] for k in ('v', 'pr', 'pl', 'w')]
            tm = .83
            a = vh.conj().T @ np.diag(np.cos(tm*sm)) @ vh
            b = -1j*um @ np.diag(np.sin(tm*sm)) @ vh
            aa, bb = np.kron(a, np.eye(n+1)), np.kron(b, np.eye(n+1))
            c.near(f'general.f{f}.N{n}.all_input_instrument',
                   exponential(v, tm) @ pr, pr @ aa+pl @ w @ bb)
            c.near(f'general.f{f}.N{n}.all_failure_retained', aa.conj().T @ aa+bb.conj().T @ bb,
                   np.eye((n+1)*f))
            min_success = float(np.min(np.linalg.eigvalsh(bb.conj().T @ bb)))
            c.near(f'general.f{f}.N{n}.worst_success', min_success, np.min(np.sin(tm*sm)**2))
            c.truth(f'general.f{f}.N{n}.not_automatically_perfect', min_success < .5)
            # Full partner-dimensional Choi vector includes BOTH locations.
            got = exponential(v, tm) @ pr/np.sqrt(pr.shape[1])
            expected = (pr @ aa+pl @ w @ bb)/np.sqrt(pr.shape[1])
            c.near(f'general.f{f}.N{n}.Choi_complete', got, expected)
            for kind in ('real', 'imaginary'):
                direction = np.zeros((f, f), complex)
                direction[0, 1] = 1 if kind == 'real' else 1j
                shift = .017
                changed = parent(n, f, matrix+shift*direction)['v']
                expected_source = parent(n, f, direction)['v']
                c.near(f'general.f{f}.N{n}.{kind}_M_source',
                       (changed-v)/shift, expected_source)
            general.append(dict(N=n, flavors=f, minimum_success=min_success,
                                singular_values=[float(x) for x in sm]))

    sequential = []
    for n in range(1, 5):
        basis, ha = two_parent(n, 0)
        b2, hb = two_parent(n, 1)
        c.truth(f'sequential.N{n}.same_parent', basis == b2)
        idx = {s: k for k, s in enumerate(basis)}
        source = np.zeros((len(basis), n+1), complex)
        for k, occ in enumerate(occupations(n)):
            source[idx[(1 << 0 | 1 << 3, *occ)], k] = 1
        got = exponential(hb, np.pi/2) @ exponential(ha, np.pi/2) @ source
        if n == 1:
            both = [k for k, (mask, _, _) in enumerate(basis)
                    if not (mask & 1) and not (mask & (1 << 3))]
            c.near('sequential.N1.no_second_quantum', got[both])
            c.near('sequential.N1.full_norm', got.conj().T @ got, np.eye(n+1))
            continue
        target = np.zeros_like(source)
        for k, occ in enumerate(occupations(n)):
            for a in range(2):
                first = list(occ)
                if not first[a]:
                    continue
                fac = np.sqrt(first[a])
                first[a] -= 1
                for b in range(2):
                    last = list(first)
                    if not last[b]:
                        continue
                    factor = fac*np.sqrt(last[b])/np.sqrt(n*(n-1))
                    last[b] -= 1
                    target[idx[(1 << (a+1) | 1 << (b+4), *last)], k] = -factor
        c.near(f'sequential.N{n}.actual_two_even_CAR_steps', got, target)
        c.near(f'sequential.N{n}.full_partner_isometry', target.conj().T @ target, np.eye(n+1))
        combined = np.kron(np.eye(2), extract_direction(n-1)) @ extract_direction(n)
        c.near(f'sequential.N{n}.analytic_composition_isometry', combined.conj().T @ combined,
               np.eye(n+1))
        permutation = np.zeros((4*(n-1), 4*(n-1)))
        for a in range(2):
            for b in range(2):
                for k in range(n-1):
                    permutation[(2*b+a)*(n-1)+k, (2*a+b)*(n-1)+k] = 1
        c.near(f'sequential.N{n}.no_extra_order_coefficient', permutation @ combined, combined)
        count = np.diag([b0+b1 for _, b0, b1 in basis])
        c.near(f'sequential.N{n}.actual_two_consumed', target.conj().T @ count @ target,
               (n-2)*np.eye(n+1))
        sequential.append(dict(N=n, complete_parent_dimension=len(basis),
                               coefficient=1/np.sqrt(n*(n-1))))

    return dict(
        schema='round1323_independent_math_v1', status='PASS_INDEPENDENT_CHECKS_ONLY',
        scientific_acceptance=False, author_code_imported=False,
        checks=c.count, max_absolute_residual=max(c.residuals.values()),
        numerical_tolerance=TOL, gamma=gamma, common_clock=clock,
        direction=direction_data, single_parent=single,
        general_flavor_matrix=general, sequential=sequential,
        residuals=c.residuals,
        scope={
            'reference_window': [1, 2, 3, 4], 'flavors': [3, 4, 5],
            'zero_reference_sector_kept': True,
            'full_unknown_reference_not_known_singlet': True,
            'capacity_lemma_reused_not_new': True,
            'TW_extension_is_extra_task': True,
            'all_SM_mass_columns_or_dynamic_gauge_instrument_certified': False,
            'two_step_CAR_numeric_flavors': 1,
            'two_step_arbitrary_flavors_by_identity_tensor_extension': True,
        },
        source_sha256={
            'independent_check.py': sha(Path(__file__)),
            '../_admission/round1323_admission.md': sha(ADMISSION),
        })


def compare(current, saved, path='result'):
    if isinstance(current, dict):
        if not isinstance(saved, dict) or set(current) != set(saved):
            raise AssertionError((path, 'dictionary fields changed'))
        for key in current:
            compare(current[key], saved[key], path+'.'+key)
    elif isinstance(current, list):
        if not isinstance(saved, list) or len(current) != len(saved):
            raise AssertionError((path, 'sequence shape changed'))
        for i, (x, y) in enumerate(zip(current, saved)):
            compare(x, y, f'{path}[{i}]')
    elif isinstance(current, float):
        if isinstance(saved, bool) or not isinstance(saved, (int, float)):
            raise AssertionError((path, 'numeric type changed'))
        if not np.isfinite(saved) or abs(current-saved) > TOL*(1+abs(current)):
            raise AssertionError((path, current, saved))
    elif type(current) is not type(saved) or current != saved:
        raise AssertionError((path, current, saved))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--write-new', action='store_true')
    args = parser.parse_args()
    if not args.write_new and not RESULT.is_file():
        raise FileNotFoundError('Saved independent_results.json is required for read-only --check')
    if args.write_new and RESULT.exists():
        raise FileExistsError('--write-new never replaces existing results')
    result = compute()
    if args.write_new:
        with RESULT.open('x', encoding='utf8', newline='\n') as out:
            json.dump(result, out, ensure_ascii=False, indent=2, allow_nan=False)
            out.write('\n')
    else:
        compare(result, json.loads(RESULT.read_text(encoding='utf8')))
    print(json.dumps(dict(status='PASS', mode='write-new' if args.write_new else 'read-only-check',
                          checks=result['checks'], max_residual=result['max_absolute_residual'],
                          result_path=str(RESULT)), ensure_ascii=False))


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    main()
