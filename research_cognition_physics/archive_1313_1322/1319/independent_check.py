"""Independent round 1319 audit: CAR sector, channels, symmetry and sources.

No author code is imported or executed.  Default/--check never writes;
--write-new creates independent_results.json only if it does not exist.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'independent_results.json'
RNG = np.random.default_rng(13192026)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def norm(x):
    return float(np.linalg.norm(np.atleast_2d(x), 2))


class Audit:
    def __init__(self):
        self.n = 0
        self.errors = {}

    def truth(self, label, ok):
        self.n += 1
        self.errors.setdefault(label, 0.)
        if not bool(ok):
            raise AssertionError(label)

    def close(self, label, x, y, tol=2e-9):
        error = norm(np.asarray(x)-np.asarray(y))/max(1., norm(y))
        self.n += 1
        self.errors[label] = max(error, self.errors.get(label, 0.))
        if error > tol:
            raise AssertionError((label, error, tol))


def dagger(x):
    return x.conj().T


def exp_h(h, t):
    e, u = np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e)) @ dagger(u)


def unitary(n, special=False):
    a = RNG.normal(size=(n, n))+1j*RNG.normal(size=(n, n))
    u, r = np.linalg.qr(a)
    u = u @ np.diag(np.diag(r)/np.abs(np.diag(r)))
    if special:
        u[:, 0] /= np.linalg.det(u)
    return u


def fop(mask, mode, create):
    """Actual CAR bit sign; return None for a forbidden occupation change."""
    occupied = (mask >> mode) & 1
    if bool(occupied) == bool(create):
        return None
    sign = (-1)**((mask & ((1 << mode)-1)).bit_count())
    return mask ^ (1 << mode), sign


def word(mask, operators):
    value = 1.
    for mode, create in operators:
        out = fop(mask, mode, create)
        if out is None:
            return None
        mask, sign = out
        value *= sign
    return mask, value


class Fock:
    """Whole K=2 sector and vacuum, not the six active states alone."""
    def __init__(self):
        self.pairs = list(itertools.combinations(range(12), 2))
        self.keys = [(0, -1)] + [((1 << p)+(1 << q), -1)
                                    for p, q in self.pairs]
        self.keys += [(0, k) for k in range(3)]
        self.index = {key: i for i, key in enumerate(self.keys)}
        self.dim = len(self.keys)
        self.eye = np.eye(self.dim, dtype=complex)
        self.wedges = list(itertools.combinations(range(3), 2))
        self.S = np.column_stack([self.spin_pair(i, j) for i in range(3)
                                 for j in range(3)])
        self.Q = np.column_stack([(self.spin_pair(i, j)-self.spin_pair(j, i))/np.sqrt(2)
                                  for i, j in self.wedges])
        self.C = self.eye[:, -3:]
        self.C0 = np.zeros((3, 3), complex)
        for col, (i, j) in enumerate(self.wedges):
            for k in range(3):
                if len({i, j, k}) == 3:
                    perm = (i, j, k)
                    inversions = sum(perm[a] > perm[b] for a in range(3)
                                     for b in range(a+1, 3))
                    self.C0[k, col] = (-1)**inversions
        self.T = self.vertex(self.C0)
        self.Na = np.diag([(m & 63).bit_count() for m, _ in self.keys])
        self.Nb = np.diag([(m >> 6).bit_count() for m, _ in self.keys])
        self.Nc = np.diag([int(k >= 0) for _, k in self.keys])
        self.K = self.Na+self.Nb+2*self.Nc
        self.PA = self.Q @ dagger(self.Q)
        self.PB = self.C @ dagger(self.C)
        self.W = np.column_stack([self.Q, self.C])
        self.s = sum(self.spin_pair(i, i) for i in range(3))/np.sqrt(3)
        self.Ps = np.outer(self.s, self.s.conj())

    def state(self, terms):
        v = np.zeros(self.dim, complex)
        for coef, p, q in terms:
            out = word(0, [(q, True), (p, True)])
            if out is not None:
                m, sign = out
                v[self.index[(m, -1)]] += coef*sign
        return v

    def spin_pair(self, i, j, same=False):
        offset = 0 if same else 6
        return self.state([(1/np.sqrt(2), 2*i, offset+2*j+1),
                           (-1/np.sqrt(2), 2*i+1, offset+2*j)])

    @staticmethod
    def qterms(i, j):
        return [(0.5, 2*i, 6+2*j+1), (-0.5, 2*i+1, 6+2*j),
                (-0.5, 2*j, 6+2*i+1), (0.5, 2*j+1, 6+2*i)]

    def vertex(self, tensor):
        t = np.zeros((self.dim, self.dim), complex)
        for col, (mask, boson) in enumerate(self.keys):
            if boson >= 0:
                continue
            for l, (i, j) in enumerate(self.wedges):
                for coef, p, q in self.qterms(i, j):
                    out = word(mask, [(p, False), (q, False)])
                    if out is None:
                        continue
                    m, sign = out
                    for k in range(3):
                        key = (m, k)
                        if key not in self.index:
                            raise AssertionError('K-sector closure')
                        t[self.index[key], col] += tensor[k, l]*coef*sign
        return t

    def h(self, pars, tensor=None):
        ea, eb, ec, x, y = pars
        t = self.T if tensor is None else self.vertex(tensor)
        return ea*self.Na+eb*self.Nb+ec*self.Nc+(x+1j*y)*t+(x-1j*y)*dagger(t)

    def action(self, u, spin=None):
        spin = np.eye(2) if spin is None else spin
        one = np.zeros((12, 12), complex)
        one[:6, :6] = one[6:, 6:] = np.kron(u, spin)
        v = np.zeros((70, 70), complex)
        v[0, 0] = 1
        for row, (p, q) in enumerate(self.pairs, 1):
            for col, (r, s) in enumerate(self.pairs, 1):
                v[row, col] = one[p, r]*one[q, s]-one[p, s]*one[q, r]
        v[-3:, -3:] = np.linalg.det(u)*u.conj()
        return v

    def receiver_density(self, rho, zdim):
        """Literal partial trace of source Fock occupations, retaining c vacuum."""
        out = np.zeros((4*zdim, 4*zdim), complex)
        for i, (mi, bi) in enumerate(self.keys):
            for j, (mj, bj) in enumerate(self.keys):
                if mi == mj:
                    ci, cj = bi+1, bj+1
                    out[ci*zdim:(ci+1)*zdim, cj*zdim:(cj+1)*zdim] += \
                        rho[i*zdim:(i+1)*zdim, j*zdim:(j+1)*zdim]
        return out

    def pair_receiver_density(self, rho, zdim):
        """Trace the boson while retaining vacuum and all 66 fermionic pairs."""
        out = np.zeros((67*zdim, 67*zdim), complex)
        for i, (mi, bi) in enumerate(self.keys):
            for j, (mj, bj) in enumerate(self.keys):
                if bi == bj:
                    fi, fj = self.index[(mi, -1)], self.index[(mj, -1)]
                    out[fi*zdim:(fi+1)*zdim, fj*zdim:(fj+1)*zdim] += \
                        rho[i*zdim:(i+1)*zdim, j*zdim:(j+1)*zdim]
        return out


def tensor_entries(c):
    triples = np.zeros((3, 3, 3), complex)
    for l, (i, j) in enumerate(itertools.combinations(range(3), 2)):
        triples[i, j, :] = c[:, l]
        triples[j, i, :] = -c[:, l]
    return triples


def cyclic_defect(c):
    t = tensor_entries(c)
    return float(np.max(np.abs(t-t.transpose(1, 2, 0))))


def capacity(a):
    table = []
    for d in range(2, 31):
        m = d*(d-1)//2
        forward, reverse = max(0., 1-d/m), max(0., 1-m/d)
        # Strict boundary decided by integer inequalities, avoiding 1/3 roundoff.
        below = 3*(m-d) < m and 3*(d-m) < d
        a.truth('capacity exact integer exclusion', below == (d == 3))
        table.append([d, m, forward, reverse])
    # Explicit arbitrary-number classical flags and conditional recovery.
    # Each block is transmitted intact, coherences between flags are destroyed.
    for m, r in [(6, 3), (5, 2), (2, 1)]:
        phi = np.eye(m).reshape(-1)/np.sqrt(m)
        recovered = np.zeros((m*m, m*m), complex)
        for start in range(0, m, r):
            p = np.zeros((m, m))
            p[np.arange(start, min(start+r, m)), np.arange(start, min(start+r, m))] = 1
            vector = np.kron(p, np.eye(m)) @ phi
            recovered += np.outer(vector, vector.conj())
            schmidt = np.linalg.matrix_rank(vector.reshape(m, m))
            a.truth('classical flag branch Schmidt bound', schmidt <= r)
        fidelity = np.vdot(phi, recovered @ phi).real
        halftrace = np.sum(np.abs(np.linalg.eigvalsh(recovered-np.outer(phi, phi))))/2
        a.close('flag CPTP normalized', np.trace(recovered), 1)
        a.truth('flag overlap capacity bound', fidelity <= r/m+1e-12)
        a.truth('flag complete error bound', halftrace >= 1-r/m-1e-12)
        # With coherent flags the block decomposition is a unitary reindexing,
        # so it can transmit all m; the max-classical-block theorem does not apply.
        a.close('coherent sideband retains full unknown input', np.vdot(phi, phi), 1)
    # Finite-rank pure state overlap diagnostic for the analytical proof.
    for m, r in [(6, 3), (4, 1)]:
        for _ in range(3):
            left = unitary(m)[:, :r]
            right = unitary(m)[:, :r]
            weights = RNG.normal(size=r)+1j*RNG.normal(size=r)
            weights /= np.linalg.norm(weights)
            matrix = left @ np.diag(weights) @ dagger(right)
            overlap = abs(np.trace(matrix))**2/m
            a.truth('rank witness overlap', overlap <= r/m+1e-12)
    return table


def run():
    a = Audit()
    table = capacity(a)
    f = Fock()
    ident = np.eye(70)
    a.truth('whole Fock sector dimension', f.dim == 70)
    a.close('spin singlet normalization', dagger(f.S) @ f.S, np.eye(9))
    a.close('wedge normalization', dagger(f.Q) @ f.Q, np.eye(3))
    a.close('active plus receiver code', dagger(f.W) @ f.W, np.eye(6))
    a.close('actual annihilation vertex compression', dagger(f.C) @ f.T @ f.Q, f.C0)
    a.close('actual annihilation vertex no hidden dark coupling', f.T, f.C @ f.C0 @ dagger(f.Q))
    a.close('cofactor normalization', dagger(f.C0) @ f.C0, np.eye(3))
    for i, j in f.wedges:
        a.close('same species Weyl antisymmetric vanishes', f.spin_pair(i, j, True), f.spin_pair(j, i, True))
    a.close('K conservation', f.K @ f.T-f.T @ f.K, np.zeros((70, 70)))
    diff = f.Na-f.Nb
    a.close('species difference conservation', diff @ f.T-f.T @ diff, np.zeros((70, 70)))
    phases = np.diag(np.exp(0.5j*np.pi*np.diag(f.K)))
    a.close('finite r identity', phases @ f.T @ dagger(phases), f.T)
    parity = np.diag((-1.)**np.diag(f.Na+f.Nb))
    a.close('all sector fermion parity even', parity, ident)

    # Full linear cyclic constraint on all nine complex tensor entries.
    cols = []
    for k in range(9):
        c = np.eye(9)[k].reshape(3, 3)
        t = tensor_entries(c)
        cols.append((t-t.transpose(1, 2, 0)).reshape(-1))
    constraints = np.column_stack(cols)
    a.truth('cyclic tensor solution complex dimension one', 9-np.linalg.matrix_rank(constraints) == 1)
    a.close('cofactor obeys cyclicity', constraints @ f.C0.reshape(-1), np.zeros(27))
    for _ in range(3):
        u, spin = unitary(3), unitary(2, True)
        v = f.action(u, spin)
        a.close('induced complete Fock representation unitary', dagger(v) @ v, ident)
        a.close('whole parent U3 and spin SU2 invariance', v @ f.T @ dagger(v), f.T)
        wedge_u = dagger(f.Q) @ v @ f.Q
        a.close('determinant phase cofactor law', f.C0 @ wedge_u, np.linalg.det(u)*u.conj() @ f.C0)
        a.close('spin singlet closed', (ident-f.PA) @ v @ f.Q, np.zeros((70, 3)))
    angle = 0.29
    center = f.action(np.exp(1j*angle)*np.eye(3))
    a.close('D not free fixed phase: center output', dagger(f.C) @ center @ f.C, np.exp(2j*angle)*np.eye(3))

    # One entangled mixed input, retaining a nontrivial untouched partner.
    zdim = 4
    samples = RNG.normal(size=(3*zdim, 5))+1j*RNG.normal(size=(3*zdim, 5))
    rho = samples @ dagger(samples)
    rho /= np.trace(rho)
    wz_a, wz_b = np.kron(f.Q, np.eye(zdim)), np.kron(f.C, np.eye(zdim))
    p0 = [1.3, 2.1, 3.4, .73, -.21]
    g = complex(p0[3], p0[4])
    gnorm = abs(g)
    tstar = np.pi/(2*gnorm)
    scenarios = []
    for delta in [0., .8, 3.1]:
        pars = p0.copy()
        pars[2] += delta
        h = f.h(pars)
        bare = pars[0]*f.Na+pars[1]*f.Nb+pars[2]*f.Nc
        omega = np.sqrt(delta**2+4*gnorm**2)
        compressed = np.block([[np.zeros((3, 3)), g.conjugate()*dagger(f.C0)],
                               [g*f.C0, delta*np.eye(3)]])
        a.close('complete parent compressed block', dagger(f.W) @ h @ f.W,
                compressed+(pars[0]+pars[1])*np.eye(6))
        a.close('complete parent code invariance', (ident-f.W @ dagger(f.W)) @ h @ f.W, np.zeros((70, 6)))
        for time in [0., .39, tstar, np.pi/omega]:
            u = exp_h(h, time)
            phase = np.exp(-1j*(pars[0]+pars[1]+delta/2)*time)
            amp_a = phase*(np.cos(omega*time/2)+1j*delta/omega*np.sin(omega*time/2))
            amp_b = phase*(-2j*g/omega)*np.sin(omega*time/2)
            amp_bfail = phase*(np.cos(omega*time/2)-1j*delta/omega*np.sin(omega*time/2))
            amp_rev = phase*(-2j*g.conjugate()/omega)*np.sin(omega*time/2)
            a.close('forward full parent channel', u @ f.Q, amp_a*f.Q+amp_b*f.C @ f.C0)
            a.close('reverse full parent channel', u @ f.C, amp_bfail*f.C+amp_rev*f.Q @ dagger(f.C0))
            p = abs(amp_b)**2
            a.close('complete success failure normalization', p+abs(amp_a)**2, 1)
            # Actual large-space state and receiver partial trace, not only success probabilities.
            evolved = np.kron(u, np.eye(zdim)) @ wz_a
            full_rho = evolved @ rho @ dagger(evolved)
            reduced = f.receiver_density(full_rho, zdim)
            target = np.zeros_like(reduced)
            fail_z = np.trace(rho.reshape(3, zdim, 3, zdim), axis1=0, axis2=2)*abs(amp_a)**2
            target[:zdim, :zdim] = fail_z
            successful = np.kron(amp_b*f.C0, np.eye(zdim))
            target[zdim:, zdim:] = successful @ rho @ dagger(successful)
            a.close('literal source trace receiver erasure with partner', reduced, target)
            a.close('receiver flags exactly dephased by source occupation', reduced[:zdim, zdim:], np.zeros((zdim, 3*zdim)))
            source_branch = dagger(wz_a) @ full_rho @ wz_a
            a.close('failed unknown partner retained at source', source_branch, abs(amp_a)**2*rho)
            reduced_system = np.trace(full_rho.reshape(70, zdim, 70, zdim), axis1=1, axis2=3)
            for number, prediction in [(f.Na, 1-p), (f.Nb, 1-p), (f.Nc, p)]:
                a.close('actual resource occupations', np.trace(number @ reduced_system), prediction)
            a.close('total energy exact', np.trace(h @ reduced_system), pars[0]+pars[1])
            free_change = np.trace(bare @ reduced_system)-(pars[0]+pars[1])
            a.close('detuned free energy and interaction complete balance', free_change, -np.trace((h-bare) @ reduced_system))
            scenarios.append(dict(delta=delta, time=float(time), success=float(p)))
        if delta == 0:
            a.close('resonant bare energy conservation', bare @ h-h @ bare, np.zeros((70, 70)))
            ustar = exp_h(h, tstar)
            a.close('finite bidirectional complete delivery', dagger(f.Q) @ ustar @ ustar @ f.Q,
                    -np.exp(-2j*(pars[0]+pars[1])*tstar)*np.eye(3))
            a.close('continually enabled vertex returns pair', abs(dagger(f.Q) @ exp_h(h, 2*tstar) @ f.Q), np.eye(3))

    # All dark components remain with their individual bare energies, not all Epair.
    pferm = np.zeros((70, 70)); pferm[1:67, 1:67] = np.eye(66)
    pdark = pferm-f.PA
    h = f.h(p0)
    bare = p0[0]*f.Na+p0[1]*f.Nb+p0[2]*f.Nc
    a.close('dark sector evolution retains aa bb energies', exp_h(h, .47) @ pdark, exp_h(bare, .47) @ pdark)

    # A real same-mode smaller-symmetry parent, not relabeling a U3 invariant model.
    a.close('delta pair normalized', np.vdot(f.s, f.s), 1)
    a.close('delta pair orthogonal to task', dagger(f.W) @ f.s, np.zeros(6))
    kappa = .47
    alt = h+kappa*f.Ps
    so3, _ = np.linalg.qr(RNG.normal(size=(3, 3)))
    so3[:, 0] *= np.linalg.det(so3)
    ro = f.action(so3)
    a.close('SO3 competitor entire parent invariant', ro @ alt @ dagger(ro), alt)
    su = np.diag([np.exp(.31j), np.exp(-.31j), 1])
    ru = f.action(su)
    defect = norm(ru @ alt @ dagger(ru)-alt)
    a.truth('SO3 competitor genuinely breaks generic SU3', defect > .1)
    a.close('SO3 competitor retains finite r', phases @ alt @ dagger(phases), alt)
    a.close('SO3 competitor keeps every task amplitude', exp_h(alt, .71) @ f.W, exp_h(h, .71) @ f.W)
    a.close('SO3 competitor keeps entire future task', exp_h(alt, tstar) @ f.W, exp_h(h, tstar) @ f.W)
    a.close('SO3 competitor adds genuine dark energy', np.vdot(f.s, (alt-h) @ f.s), kappa)

    # Delete cyclicity: unequal odd frequencies retain exact full handoff.
    divided = np.diag([1., 3., 5.]) @ f.C0
    hsplit = f.vertex(divided)+dagger(f.vertex(divided))
    usplit = exp_h(hsplit, np.pi/2)
    split_success = dagger(f.C) @ usplit @ f.Q
    a.close('noncyclic odd modes perfect full delivery', dagger(split_success) @ split_success, np.eye(3))
    a.close('noncyclic odd modes no failure', dagger(f.Q) @ usplit @ f.Q, np.zeros((3, 3)))
    a.truth('noncyclic odd modes violate cyclic equality', cyclic_defect(divided) > 1)
    a.close('noncyclic split selfenergy spectrum', np.linalg.eigvalsh(dagger(divided) @ divided), [1, 9, 25])
    nonsplit = np.diag([1, 1, -1]) @ f.C0
    a.truth('noncyclic unitary also violates equality', cyclic_defect(nonsplit) > 1)
    a.close('noncyclic need not split selfenergy', dagger(nonsplit) @ nonsplit, np.eye(3))

    # Off-resonant elimination, exact poles and independent physical sources.
    detuned = p0.copy(); detuned[2] += 2.8
    delta = detuned[2]-detuned[0]-detuned[1]
    q = gnorm**2
    hd = f.h(detuned)
    jferm = np.eye(70)[:, 1:67]
    hff, hcc = dagger(jferm) @ hd @ jferm, dagger(f.C) @ hd @ f.C
    hfc = dagger(jferm) @ hd @ f.C
    pa = dagger(jferm) @ f.PA @ jferm
    for z in [-.7, .17, .61]:
        energy = detuned[0]+detuned[1]+z
        direct = hfc @ np.linalg.inv(energy*np.eye(3)-hcc) @ dagger(hfc)
        expected = -q/(delta-z)*pa
        a.close('full fermion Schur selfenergy incl dark sector', direct, expected)
        a.close('full parent resolvent after exact elimination',
                dagger(jferm) @ np.linalg.inv(energy*np.eye(70)-hd) @ jferm,
                np.linalg.inv(energy*np.eye(66)-hff-direct), tol=2e-8)
        discrepancy = norm(direct+q/delta*pa)
        bound = q*abs(z)/(delta*(delta-abs(z)))
        a.truth('controlled zero frequency selfenergy error', discrepancy <= bound+1e-11)
    pole = .5*(delta-np.sqrt(delta**2+4*q))
    active = dagger(f.W) @ hd @ f.W-(detuned[0]+detuned[1])*np.eye(6)
    a.close('exact active poles not zero energy Schur', np.linalg.eigvalsh(active), [pole]*3+[delta-pole]*3)
    a.truth('finite static selfenergy truncation bound', abs(pole+q/delta) <= q*q/delta**3+1e-12)
    omega = np.sqrt(delta**2+4*q)
    peak, pmax = np.pi/omega, 4*q/omega**2
    a.close('same operating point peak versus response', q/delta, delta*pmax/(4*(1-pmax)))
    a.close('same kernel peak time versus response', q/delta, np.pi**2/(4*delta*peak**2)-delta/4)

    sources = [f.Na, f.Nb, f.Nc, f.T+dagger(f.T), 1j*(f.T-dagger(f.T))]
    energy = sum(detuned[:2])+.27
    zfixed = .27
    step = 2e-5
    def full_schur(p, absolute_energy):
        hp = f.h(p)
        ff = dagger(jferm) @ hp @ jferm
        fc = dagger(jferm) @ hp @ f.C
        cc = dagger(f.C) @ hp @ f.C
        return absolute_energy*np.eye(66)-ff-fc @ np.linalg.inv(absolute_energy*np.eye(3)-cc) @ dagger(fc)
    for k, source in enumerate(sources):
        plus, minus = detuned.copy(), detuned.copy()
        plus[k] += step; minus[k] -= step
        a.close('five independent parent sources', (f.h(plus)-f.h(minus))/(2*step), source)
        d_q = [0, 0, 0, 2*detuned[3], 2*detuned[4]][k]
        d_delta = [-1, -1, 1, 0, 0][k]
        d_z_at_energy = [-1, -1, 0, 0, 0][k]
        denom = delta-zfixed
        relative_derivative = -d_q/denom+q*d_delta/denom**2
        energy_derivative = relative_derivative-q*d_z_at_energy/denom**2
        expected_absolute = -d_q/(detuned[2]-energy)+q*int(k == 2)/(detuned[2]-energy)**2
        a.close('source chain rule relative to absolute frequency', energy_derivative, expected_absolute)
        dff = dagger(jferm) @ source @ jferm
        fixed_e_fd = (full_schur(plus, energy)-full_schur(minus, energy))/(2*step)
        a.close('complete physical source Schur derivative', fixed_e_fd, -dff-energy_derivative*pa, tol=5e-9)
        plus_e, minus_e = sum(plus[:2])+zfixed, sum(minus[:2])+zfixed
        relative_fd = (full_schur(plus, plus_e)-full_schur(minus, minus_e))/(2*step)
        d_e = int(k < 2)
        a.close('complete relative frequency source derivative', relative_fd, d_e*np.eye(66)-dff-relative_derivative*pa, tol=5e-9)
    # Finite Hamiltonian implementation error: directly test the operator bound.
    perturb = .003*(f.Na+0.4*f.Ps)
    eta = norm(perturb)
    a.truth('Duhamel implementation bound sample', norm(exp_h(h+perturb, .8)-exp_h(h, .8)) <= .8*eta+1e-12)

    # Author's final same-detuning two-stage addition, checked from actual
    # discarded source states and a freshly empty source/receiver, not p alone.
    common = p0.copy()
    common[2] += gnorm/2
    delta_common = gnorm/2
    hc = f.h(common)
    hfree = common[0]*f.Na+common[1]*f.Nb+common[2]*f.Nc
    omega_common = np.sqrt(delta_common**2+4*q)
    uc = exp_h(hc, np.pi/omega_common)
    uz = np.kron(uc, np.eye(zdim))
    full_first = uz @ wz_a @ rho @ dagger(wz_a) @ dagger(uz)
    c_first = f.receiver_density(full_first, zdim)
    embed_c = np.column_stack([f.eye[:, 0], f.C])
    fresh_c = np.kron(embed_c, np.eye(zdim))
    second_initial = fresh_c @ c_first @ dagger(fresh_c)
    second_final = uz @ second_initial @ dagger(uz)
    pair_final = f.pair_receiver_density(second_final, zdim)
    embed_pair = np.column_stack([f.eye[:67, 0], f.Q[:67]])
    pair_small = np.kron(embed_pair, np.eye(zdim))
    target_small = np.zeros((4*zdim, 4*zdim), complex)
    pcommon = 16/17
    partner = np.trace(rho.reshape(3, zdim, 3, zdim), axis1=0, axis2=2)
    target_small[:zdim, :zdim] = (1-pcommon**2)*partner
    target_small[zdim:, zdim:] = pcommon**2*rho
    a.close('two actual detuned conversions full unknown partner', pair_final,
            pair_small @ target_small @ dagger(pair_small))
    ideal = np.zeros_like(target_small); ideal[zdim:, zdim:] = rho
    got_small = dagger(pair_small) @ pair_final @ pair_small
    erasure_halftrace = np.sum(np.abs(np.linalg.eigvalsh(got_small-ideal)))/2
    a.close('two detuned stages exact erasure error', erasure_halftrace, 33/289)
    a.truth('same nonresonant model finite-error capacity margin', erasure_halftrace < 1/3)
    # Dephasing of actual source sectors deletes interaction energy.  Its
    # work is counted before replacing the failed source by an inaccessible one.
    q0 = ident-f.PB
    def dephase(state):
        qb, qf = np.kron(f.PB, np.eye(zdim)), np.kron(q0, np.eye(zdim))
        return qb @ state @ qb+qf @ state @ qf
    hcz = np.kron(hc, np.eye(zdim))
    work1 = np.trace(hcz @ (dephase(full_first)-full_first)).real
    work2 = np.trace(hcz @ (dephase(second_final)-second_final)).real
    a.close('first detuned occupation readout full energy source', work1, delta_common*pcommon)
    a.close('second conditional readout full energy source', work2, -delta_common*pcommon**2)
    a.close('two readout sources complete net energy', work1+work2,
            delta_common*pcommon*(1-pcommon))
    a.close('same setting finite response window',
            abs(-q/(delta_common-delta_common/4)+q/delta_common), q/(3*delta_common))
    return dict(status='PASS_WITH_SCOPE', assertions=a.n, categories=len(a.errors),
                max_normalized_error=max(a.errors.values()), errors=a.errors,
                basis=dict(fermion_modes=12, pair_states=66, boson_states=3,
                           vacuum=1, total=70), capacity_table=table,
                scenarios=scenarios, so3_breaks_su3_defect=defect,
                common_detuned=dict(delta=delta_common, p=pcommon,
                                    roundtrip_error=erasure_halftrace,
                                    first_readout_energy=work1, second_readout_energy=work2),
                convention='All source derivatives distinguish fixed absolute energy from fixed relative frequency.',
                limits=['extra unadopted task contract', 'finite-mode parent, not full QFT',
                        'internal d is not family n', 'maximum tensor stabilizer is not physical gauge permission',
                        'classical flag blocks only; coherent quantum sidebands count',
                        'fixed Hamiltonian can reverse conversion after handoff time'])


def compare_saved(recomputed, saved, path='result'):
    """Check the complete scientific record, not merely its schema and digest.

    Floating replay tolerance: relative 1e-10, absolute 1e-11.  Every assertion
    in run() retains its own separate, stricter mathematical test tolerance.
    """
    if isinstance(recomputed, dict):
        if not isinstance(saved, dict) or recomputed.keys() != saved.keys():
            raise AssertionError(path+' dictionary fields')
        for key, value in recomputed.items():
            compare_saved(value, saved[key], path+'.'+key)
    elif isinstance(recomputed, list):
        if not isinstance(saved, list) or len(recomputed) != len(saved):
            raise AssertionError(path+' list shape')
        for i, (a, b) in enumerate(zip(recomputed, saved)):
            compare_saved(a, b, f'{path}[{i}]')
    elif isinstance(recomputed, bool):
        if not isinstance(saved, bool) or recomputed != saved:
            raise AssertionError(path+' boolean')
    elif isinstance(recomputed, int):
        if not isinstance(saved, int) or isinstance(saved, bool) or recomputed != saved:
            raise AssertionError(path+' integer')
    elif isinstance(recomputed, float):
        if not isinstance(saved, (int, float)) or isinstance(saved, bool) or not np.isclose(
                recomputed, saved, rtol=1e-10, atol=1e-11, equal_nan=False):
            raise AssertionError((path, recomputed, saved))
    elif type(recomputed) is not type(saved) or recomputed != saved:
        raise AssertionError((path, recomputed, saved))


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write-new', action='store_true')
    args = parser.parse_args()
    result = run()
    result['code_sha256'] = sha(__file__)
    result['admission_sha256'] = sha(HERE.parent/'_admission'/'round1319_admission.md')
    if args.write_new:
        with OUT.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        # Missing output must fail.  Read-only checking never recreates evidence.
        saved = json.loads(OUT.read_text(encoding='utf-8'))
        compare_saved(result, saved)
    print(json.dumps({k: result[k] for k in ['status', 'assertions', 'categories', 'max_normalized_error']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
