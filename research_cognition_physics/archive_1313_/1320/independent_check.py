"""Independent 1320 audit.  Actual CAR sector, open roles and matching dictionary.

No author code imports.  Default/--check is read-only and requires saved evidence;
--write-new exclusively creates the result.  Replay compares every stored field.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE/'independent_results.json'
RNG = np.random.default_rng(13202026)


def adj(x): return x.conj().T
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fnorm(x): return float(np.linalg.norm(np.asarray(x)))


class Audit:
    def __init__(self): self.n, self.errors = 0, {}
    def truth(self, key, test):
        self.n += 1; self.errors.setdefault(key, 0.)
        if not bool(test): raise AssertionError(key)
    def close(self, key, x, y, tol=2e-9):
        err = fnorm(np.asarray(x)-np.asarray(y))/max(1., fnorm(y))
        self.n += 1; self.errors[key] = max(err, self.errors.get(key, 0.))
        if not np.isfinite(err) or err > tol: raise AssertionError((key, err, tol))


def evolve(h, t):
    e, u = np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e)) @ adj(u)


def unitary(n, special=False):
    u, r = np.linalg.qr(RNG.normal(size=(n, n))+1j*RNG.normal(size=(n, n)))
    u = u @ np.diag(np.diag(r)/np.abs(np.diag(r)))
    if special: u[:, 0] /= np.linalg.det(u)
    return u


def action(mask, mode, create):
    occupied = bool(mask & (1 << mode))
    if occupied == create: return None
    return mask ^ (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())


def pair_vector(terms):
    """Sparse two-particle state from ordered creation words."""
    ans = {}
    for coef, p, q in terms:
        one = action(0, q, True)
        two = action(one[0], p, True)
        if two is not None:
            ans[two[0]] = ans.get(two[0], 0.)+coef*one[1]*two[1]
    return {m: v for m, v in ans.items() if abs(v) > 1e-15}


def eps3(i, j, k):
    if len({i, j, k}) != 3: return 0
    t = (i, j, k)
    return (-1)**sum(t[a] > t[b] for a in range(3) for b in range(a+1, 3))


def qq_terms(k, flavor1=0, flavor2=0, full=False):
    """full=True is the unsymmetrized flavor polynomial, useful for CAR tests."""
    ans = []
    colors = itertools.permutations(range(3), 2) if full else itertools.combinations(range(3), 2)
    for i, j in colors:
        color = eps3(i, j, k)
        if not color: continue
        for weak, spin in itertools.product(range(2), repeat=2):
            p = 12*flavor1+4*i+2*weak+spin
            q = 12*flavor2+4*j+2*(1-weak)+(1-spin)
            ans.append((color*(-1)**(weak+spin)/2, p, q))
    return ans


class Fock:
    def __init__(self):
        self.pairs = list(itertools.combinations(range(28), 2))
        self.keys = [(0, -1)]+[((1 << i)+(1 << j), -1) for i, j in self.pairs]+[(0, k) for k in range(3)]
        self.index = {x: i for i, x in enumerate(self.keys)}
        self.dim = len(self.keys)
        self.eye = np.eye(self.dim, dtype=complex)
        self.termsA = [qq_terms(k) for k in range(3)]
        self.termsC = [[((-1)**(w+s)/2, 12+4*k+2*w+s, 24+2*(1-w)+(1-s))
                        for w, s in itertools.product(range(2), repeat=2)] for k in range(3)]
        self.A = np.column_stack([self.state(t) for t in self.termsA])
        self.C = np.column_stack([self.state(t) for t in self.termsC])
        self.B = self.eye[:, -3:]
        self.W = np.column_stack([self.A, self.B, self.C])
        self.TA, self.TC = self.vertex(self.termsA), self.vertex(self.termsC)
        self.N = [np.diag([((m >> start) & ((1 << length)-1)).bit_count()
                           for m, _ in self.keys]) for start, length in [(0, 12), (12, 12), (24, 4)]]
        self.N.append(np.diag([int(s >= 0) for _, s in self.keys]))
        self.K = sum(self.N[:3])+2*self.N[3]

    def state(self, terms):
        v = np.zeros(self.dim, complex)
        for m, coef in pair_vector(terms).items(): v[self.index[(m, -1)]] = coef
        return v

    def vertex(self, terms):
        v = np.zeros((self.dim, self.dim), complex)
        for col, (m, s) in enumerate(self.keys):
            if s >= 0: continue
            for k in range(3):
                for coef, p, q in terms[k]:
                    first = action(m, p, False)
                    if first is None: continue
                    second = action(first[0], q, False)
                    if second is None: continue
                    key = (second[0], k)
                    if key not in self.index: raise AssertionError('outside full K sector')
                    v[self.index[key], col] += np.conj(coef)*first[1]*second[1]
        return v

    def h(self, pars):
        e_q, e_aq, e_al, e_s, ar, ai, cr, ci = pars
        a, c = ar+1j*ai, cr+1j*ci
        return sum(e*n for e, n in zip(pars[:4], self.N))+a*self.TA+a.conjugate()*adj(self.TA)+c*self.TC+c.conjugate()*adj(self.TC)

    def bare(self, pars): return sum(e*n for e, n in zip(pars[:4], self.N))

    def evolution(self, pars, t):
        """Exact complete evolution; the verified dark subspace remains free."""
        h = self.h(pars)
        bare = self.bare(pars)
        ua = evolve(adj(self.W) @ h @ self.W, t)
        u0 = np.diag(np.exp(-1j*np.diag(bare)*t))
        return u0+self.W @ (ua-adj(self.W) @ u0 @ self.W) @ adj(self.W)

    def transformed_code(self, code, one):
        """Induced wedge action without assuming a code representation."""
        out = np.zeros_like(code)
        for col in range(code.shape[1]):
            for j in np.flatnonzero(abs(code[:, col]) > 1e-14):
                m, s = self.keys[j]
                if s >= 0 or not m: continue
                r, t = [i for i in range(28) if m & (1 << i)]
                for i, (p, q) in enumerate(self.pairs, 1):
                    out[i, col] += code[j, col]*(one[p, r]*one[q, t]-one[p, t]*one[q, r])
        return out

    def receiver_trace(self, psi, receive_q):
        """True partial trace on a pure system-partner state.

        All reachable receiving occupations are retained, including single
        particles in unused dark states. Bosonic occupation is part of source.
        """
        mask_q = (1 << 12)-1
        def split(m): return (m & mask_q, m >> 12) if receive_q else (m >> 12, m & mask_q)
        recvkeys = sorted({split(m)[0] for m, _ in self.keys})
        index = {m: j for j, m in enumerate(recvkeys)}
        zdim = psi.shape[1]
        by_source = {}
        for j, (m, b) in enumerate(self.keys):
            receive, source = split(m)
            key = (source, b)
            if key not in by_source: by_source[key] = np.zeros((len(recvkeys), zdim), complex)
            by_source[key][index[receive]] += psi[j]
        rho = sum(np.outer(v.reshape(-1), v.reshape(-1).conj()) for v in by_source.values())
        code = self.A if receive_q else self.C
        rec_code = np.zeros((len(recvkeys), 3), complex)
        for j, (m, b) in enumerate(self.keys):
            if b == -1 and split(m)[1] == 0: rec_code[index[split(m)[0]]] += code[j]
        empty = np.zeros(len(recvkeys)); empty[index[0]] = 1
        return rho, rec_code, empty


def roles(a):
    for w in range(1, 20):
        # 3w input, 3-dimensional receiver; strict < 1/2 requires w < 2.
        a.truth('open W capacity excludes all larger factors', (2*(w-1) < w) == (w == 1))
    omega = np.exp(2j*np.pi/3)
    rep = np.diag([1, omega, omega**2])
    twist = omega
    intertwiner = np.zeros((3, 3), complex)
    for j in range(3): intertwiner[(j+1) % 3, j] = 1
    a.close('disconnected fullrank character self twist', intertwiner @ (twist*rep), rep @ intertwiner)
    a.close('disconnected twist cube one', twist**3, 1)
    a.truth('disconnected twist not trivial', abs(twist-1) > 1)
    # U(1), V* weights -1,0,1 and D=det(V)=1; W* weight +1.
    partial = np.zeros((3, 3)); partial[1, 0] = partial[2, 1] = 1
    for theta in [.19, .61]:
        out = np.diag(np.exp(1j*theta*np.array([-1, 0, 1])))
        inp = np.exp(1j*theta)*out
        a.close('connected partial intertwiner loophole', partial @ inp, out @ partial)
    a.truth('connected partial map not full rank', np.linalg.matrix_rank(partial) == 2)
    # Lie derivative of determinant condition forces any continuous weight to zero.
    a.truth('connected character infinitesimal cubic identity', all(3*k != 0 for k in [-3, -2, -1, 1, 2, 3]))


def run():
    a = Audit(); roles(a); f = Fock()
    a.truth('whole parent dimension', f.dim == 382 and len(f.pairs) == 378)
    a.close('three actual code normalizations', adj(f.W) @ f.W, np.eye(9))
    a.close('QQ annihilation actual vertex', f.TA, f.B @ adj(f.A))
    a.close('antiQL annihilation actual vertex', f.TC, f.B @ adj(f.C))
    for k in range(3):
        one, two = pair_vector(qq_terms(k, 0, 1, True)), pair_vector(qq_terms(k, 1, 0, True))
        a.truth('same Q symmetric flavor coefficient from CAR', one == two and bool(one))
    for vertex in [f.TA, f.TC]:
        a.close('complete K conservation', f.K @ vertex-vertex @ f.K, np.zeros_like(vertex))
    qcharge = f.N[0]-f.N[1]+3*f.N[2]+2*f.N[3]
    rdiag = np.exp(.5j*np.pi*np.diag(f.N[0]-f.N[1]-f.N[2]+2*f.N[3]))
    for vertex in [f.TA, f.TC]:
        a.close('old hypercharge relation at actual vertex', qcharge @ vertex-vertex @ qcharge, np.zeros_like(vertex))
        a.close('finite old r compatibility', rdiag[:, None]*vertex*rdiag.conj()[None, :], vertex)
    # No different flavor role or smaller physical group is eliminated here.
    for _ in range(2):
        u, weak, spin = unitary(3), unitary(2, True), unitary(2, True)
        one = np.zeros((28, 28), complex)
        one[:12, :12] = np.kron(u, np.kron(weak, spin))
        one[12:24, 12:24] = np.kron(u.conj(), np.kron(weak.conj(), spin.conj()))
        one[24:, 24:] = np.linalg.det(u)*np.kron(weak.conj(), spin.conj())
        rb = np.linalg.det(u)*u.conj()
        a.close('QQ full color weak spin covariance', f.transformed_code(f.A, one), f.A @ rb)
        a.close('antiQL character and weak spin covariance', f.transformed_code(f.C, one), f.C @ rb)

    pars = [1.1, .9, 1.3, 2.2, .7, .2, -.3, .6]
    ga, gc = complex(*pars[4:6]), complex(*pars[6:8])
    h, bare = f.h(pars), f.bare(pars)
    a.close('complete active invariant subspace', h @ f.W, f.W @ (adj(f.W) @ h @ f.W))
    dark = np.eye(382)-f.W @ adj(f.W)
    a.close('all dark sectors retained with bare energy', (h-bare) @ dark, np.zeros_like(h))
    a.close('factorized evolution equals full parent exponential', f.evolution(pars, .31), evolve(h, .31), tol=8e-9)
    strengths, gaps = [], []
    phi = np.eye(3, dtype=complex)/np.sqrt(3)  # maximal entangled Choi witness
    for ratio, detuning in [(1., 0.), (1.3, 0.), (2.2, 0.), (1.3, .6)]:
        p = pars.copy(); p[3] = 2*p[0]+detuning
        p[6], p[7] = ratio*pars[6], ratio*pars[7]
        aa, cc = complex(*p[4:6]), complex(*p[6:8])
        g2 = abs(aa)**2+abs(cc)**2
        gamma = (abs(aa)**2-abs(cc)**2)/g2
        pbound = 1-gamma**2
        times = [np.pi/np.sqrt(g2)] if detuning == 0 else [.21, 1.1, 2.7]
        for t in times:
            u = f.evolution(p, t)
            ac, ca = adj(f.C) @ u @ f.A, adj(f.A) @ u @ f.C
            ptransfer = float(np.trace(adj(ac) @ ac).real/3)
            a.truth('all time endpoint transfer bound', ptransfer <= pbound+2e-11)
            a.close('equal forward and reverse probabilities', adj(ca) @ ca, ptransfer*np.eye(3))
            a.close('endpoint transfer is scalar on internal data', ac, np.trace(ac)/3*np.eye(3))
            a.close('all actual failure branches retained', adj(f.A) @ adj(u) @ u @ f.A, np.eye(3))
            if detuning == 0:
                a.close('resonant endpoint maximum', ptransfer, pbound)
                a.close('mediator returns empty at transfer', adj(f.B) @ u @ f.A, np.zeros((3, 3)))
                predicted = -np.exp(-2j*p[0]*t)*2*aa*cc.conjugate()/g2
                a.close('orientation and conjugation of actual transfer', ac, predicted*np.eye(3))
                # Actual source trace yields Choi state of complete erasure.
                for input_code, receive_q in [(f.A, False), (f.C, True)]:
                    psi = u @ input_code @ phi
                    reduced, recvcode, empty = f.receiver_trace(psi, receive_q)
                    success = (recvcode @ phi).reshape(-1)
                    target = ptransfer*np.outer(success, success.conj())+(1-ptransfer)*np.kron(np.outer(empty, empty), np.eye(3)/3)
                    a.close('full receiver trace with unknown partner and failure', reduced, target, tol=8e-9)
                error_tr, error_rt = 1-ptransfer, 1-ptransfer**2
                a.close('direct error equals strength imbalance squared', error_tr, gamma**2)
                a.close('roundtrip bound distinct from direct error', gamma**2, 1-np.sqrt(1-error_rt))
                strengths.append(dict(ratio=ratio, delta=detuning, transfer=ptransfer,
                                      imbalance=gamma, direct_error=error_tr, roundtrip_error=error_rt))
            # Endpoint energies equal. Free energy in a detuned mediator is
            # compensated by interaction energy, before any actual readout.
            psi = u @ f.A @ np.array([1, 1j, -1])/np.sqrt(3)
            current_h, current_bare = f.h(p), f.bare(p)
            initial_e = 2*p[0]
            a.close('full parent energy conserved', np.vdot(psi, current_h @ psi), initial_e)
            a.close('detuned mediator complete energy balance', np.vdot(psi, current_bare @ psi)-initial_e,
                    -np.vdot(psi, (current_h-current_bare) @ psi))
        # A dark combination survives; the mediator never carries all 6
        # endpoint unknown amplitudes. It is not part of the one-endpoint task.
        darkcode = (cc*f.A-aa*f.C)/np.sqrt(g2)
        a.close('simultaneous endpoint dark code is real', f.h(p) @ darkcode, 2*p[0]*darkcode)
        a.close('dark code normalized', adj(darkcode) @ darkcode, np.eye(3))
    # Precisely balanced complex couplings provide the nonempty witness.
    balanced = pars.copy(); balanced[6:8] = [-pars[5], pars[4]]
    tbal = np.pi/np.sqrt(2*(pars[4]**2+pars[5]**2))
    ubal = f.evolution(balanced, tbal)
    a.close('balanced complete endpoint exchange', adj(f.A) @ adj(ubal) @ f.C @ adj(f.C) @ ubal @ f.A, np.eye(3))

    # Actual two-pass recovery from the first receiver's state only.  There is
    # no access to the departed Q source, and failures are never postselected.
    p = pars.copy(); p[6] *= 1.3; p[7] *= 1.3
    gsum = sum(abs(complex(*p[j:j+2]))**2 for j in [4, 6])
    u = f.evolution(p, np.pi/np.sqrt(gsum))
    received, receiving_code, empty = f.receiver_trace(u @ f.A @ phi, False)
    small_code = np.column_stack([empty, receiving_code])
    compress = np.kron(small_code, np.eye(3))
    small = adj(compress) @ received @ compress
    eigenvalues, eigenvectors = np.linalg.eigh(small)
    fresh = np.column_stack([f.eye[:, 0], f.C])
    recovered = None
    for j, value in enumerate(eigenvalues):
        if value < 1e-13: continue
        wave = u @ fresh @ (np.sqrt(value)*eigenvectors[:, j].reshape(4, 3))
        out, target_code, target_empty = f.receiver_trace(wave, True)
        recovered = out if recovered is None else recovered+out
    ptrans = float(np.trace(adj(f.A) @ adj(u) @ f.C @ adj(f.C) @ u @ f.A).real/3)
    success_wave = (target_code @ phi).reshape(-1)
    target = ptrans**2*np.outer(success_wave, success_wave.conj())+(1-ptrans**2)*np.kron(np.outer(target_empty, target_empty), np.eye(3)/3)
    a.close('two actual receiver-only stages complete partner channel', recovered, target, tol=8e-9)

    # Deleting simultaneous constant portals removes the balance restriction.
    # The extra independent controls are public; they are not silently available.
    first, second = pars.copy(), pars.copy()
    first[6:8] = [0., 0.]
    second[4:6] = [0., 0.]; second[6] *= 3; second[7] *= 3
    t_first = np.pi/(2*abs(complex(*first[4:6])))
    t_second = np.pi/(2*abs(complex(*second[6:8])))
    sequential = f.evolution(second, t_second) @ f.evolution(first, t_first)
    transfer = adj(f.C) @ sequential @ f.A
    a.close('independent portal scheduling bypasses fixed-H balance', adj(transfer) @ transfer, np.eye(3))
    a.truth('scheduled counter uses genuinely unequal strengths', abs(abs(complex(*first[4:6]))-abs(complex(*second[6:8]))) > 1)

    # Complete finite energy sources, retaining their independent mode identity.
    sources = f.N+[f.TA+adj(f.TA), 1j*(f.TA-adj(f.TA)), f.TC+adj(f.TC), 1j*(f.TC-adj(f.TC))]
    step = 2e-5
    for k, source in enumerate(sources):
        pp, pm = pars.copy(), pars.copy(); pp[k] += step; pm[k] -= step
        a.close('eight full mode sources', (f.h(pp)-f.h(pm))/(2*step), source)
    # Exact finite energy Schur kernel is NOT the relativistic mass-squared kernel.
    ends = np.column_stack([f.A, f.C])
    hs = adj(f.B) @ h @ f.B
    coupling = adj(f.B) @ h @ ends
    for energy in [.37, 1.41, 3.1]:
        sigma = adj(coupling) @ np.linalg.inv(energy*np.eye(3)-hs) @ coupling
        expected = np.kron(np.outer(np.conj([ga, gc]), [ga, gc])/(energy-pars[3]), np.eye(3))
        a.close('finite energy shared Schur kernel', sigma, expected)

    # Fixed, declared mode overlap dictionary to dimensionless current couplings.
    lam = np.array([.4+.2j, -.3+.7j])
    kappa = np.array([1.7-.2j, .8+.3j])  # energy units; includes statistics factors
    mass2 = 11.0
    def gram(ll, mm, s=0.): return np.outer(ll.conj(), ll)/(mm-s)
    kernel = gram(lam, mass2)
    a.close('raw current kernel Hermitian', kernel, adj(kernel))
    a.truth('raw current kernel positive', np.linalg.eigvalsh(kernel).min() >= -1e-13)
    a.close('single shared mediator saturation', abs(kernel[0, 1])**2, kernel[0, 0]*kernel[1, 1])
    gmode = kappa*lam
    weighted = abs(kappa)**2*np.diag(kernel).real
    imbalance = (weighted[0]-weighted[1])/sum(weighted)
    direct = imbalance**2
    field_ratio = weighted[0]/weighted[1]
    root = np.sqrt(direct)
    a.truth('task to fixed overlap Gram ratio bound', (1-root)/(1+root)-1e-12 <= field_ratio <= (1+root)/(1-root)+1e-12)
    a.close('fixed mode dictionary same imbalance', imbalance,
            (abs(gmode[0])**2-abs(gmode[1])**2)/sum(abs(gmode)**2))
    for invariant in [-2., .7, 2.5]:
        k = gram(lam, mass2, invariant)
        bound = abs(invariant)/(mass2-abs(invariant))*fnorm(kernel)
        a.truth('mass squared finite momentum window', fnorm(k-kernel) <= bound+1e-12)
        a.close('all current coefficients share same denominator', k, mass2/(mass2-invariant)*kernel)
    # Independent sources in the field current kernel: mass-squared and both complex couplings.
    base = np.array([mass2, lam[0].real, lam[0].imag, lam[1].real, lam[1].imag])
    def from_parameters(pp): return gram(np.array([complex(pp[1], pp[2]), complex(pp[3], pp[4])]), pp[0], .8)
    for j in range(5):
        pp, pm = base.copy(), base.copy(); pp[j] += step; pm[j] -= step
        if j == 0: derivative = -gram(lam, mass2, .8)/(mass2-.8)
        else:
            dl = np.zeros(2, complex); dl[(j-1)//2] = 1 if j % 2 else 1j
            derivative = (np.outer(dl.conj(), lam)+np.outer(lam.conj(), dl))/(mass2-.8)
        a.close('complete relativistic current kernel sources', (from_parameters(pp)-from_parameters(pm))/(2*step), derivative)
    # A physical variation changes scalar frequency, overlap normalization and
    # both field couplings together.  Do not freeze these chained sources.
    dlam = np.array([.12-.07j, -.08+.04j])
    alpha = np.array([.09+.03j, -.05+.02j])
    scalar_momentum2 = .6
    def physical_curve(x):
        mm = mass2+x
        ll = lam+x*dlam
        kk = kappa*np.exp(alpha*x)*(mass2/mm)**.25
        gg = kk*ll
        pp = pars.copy()
        pp[3] = np.sqrt(mm+scalar_momentum2)
        pp[4:8] = [gg[0].real, gg[0].imag, gg[1].real, gg[1].imag]
        return pp
    dg = kappa*dlam+lam*kappa*(alpha-1/(4*mass2))
    source_curve = f.N[3]/(2*np.sqrt(mass2+scalar_momentum2))
    source_curve = source_curve+dg[0]*f.TA+dg[0].conjugate()*adj(f.TA)+dg[1]*f.TC+dg[1].conjugate()*adj(f.TC)
    a.close('field to finite mode complete chained source',
            (f.h(physical_curve(step))-f.h(physical_curve(-step)))/(2*step), source_curve)
    j0 = np.array([[[1, .2j], [.3, -.1]], [[.2, -.4], [.1j, .7]]], complex)
    dj = np.array([[[0, .1], [.2j, 0]], [[.03, 0], [0, -.06j]]], complex)
    def observable_curve(x):
        kk = gram(lam+x*dlam, mass2+x, .8)
        jj = j0+x*dj
        return sum(kk[i, j]*adj(jj[i]) @ jj[j] for i in range(2) for j in range(2))
    dkernel = (np.outer(dlam.conj(), lam)+np.outer(lam.conj(), dlam))/(mass2-.8)-np.outer(lam.conj(), lam)/(mass2-.8)**2
    current_kernel = gram(lam, mass2, .8)
    complete_d = sum(dkernel[i, j]*adj(j0[i]) @ j0[j]+current_kernel[i, j]*(adj(dj[i]) @ j0[j]+adj(j0[i]) @ dj[j])
                     for i in range(2) for j in range(2))
    a.close('field response keeps both current variation terms',
            (observable_curve(step)-observable_curve(-step))/(2*step), complete_d)
    # Removal of fixed overlaps allows arbitrary lambda magnitude ratio despite
    # balanced mode transfers; a second independent mediator spoils rank one.
    for ratio in [.03, 20.]:
        free_lam = np.array([ratio, 1.])
        adjustable_overlap = np.array([1/ratio, 1.])
        a.close('free overlaps do not select field coupling ratio', adjustable_overlap*free_lam, [1, 1])
    other = gram(np.array([.6, -.2j]), 7.)
    a.truth('extra mediator breaks one-kernel saturation', np.linalg.det(kernel+other).real > 1e-6)
    contact = kernel+np.diag([.1, 0])
    a.truth('independent contact breaks one-kernel saturation', np.linalg.det(contact).real > 1e-6)
    return dict(status='PASS_WITH_SCOPE', round=1320, assertions=a.n, categories=len(a.errors),
                max_normalized_frobenius_error=max(a.errors.values()), errors=a.errors,
                basis=dict(fermion_modes=28, pair_states=378, scalar_states=3, vacuum=1, total=382),
                task_samples=strengths, fixed_overlap_gram=dict(direct_error=direct, weighted_ratio=field_ratio),
                scope=['open W selected only under full-code capacity contract',
                       'connected full-rank intertwiner required for exact character identity',
                       'single declared flavor code, not all n squared flavor information',
                       'joint family-number qualifications are not checked in this round',
                       'finite energy kernel differs from relativistic mass squared kernel',
                       'fixed overlap bridge conditional; no full QFT instrument certification',
                       'single positive shared kernel and no independent contact at the same order'])


def compare(x, y, path='result'):
    if isinstance(x, dict):
        if not isinstance(y, dict) or x.keys() != y.keys(): raise AssertionError(path+' fields')
        for k in x: compare(x[k], y[k], path+'.'+k)
    elif isinstance(x, list):
        if not isinstance(y, list) or len(x) != len(y): raise AssertionError(path+' sequence')
        for i, (a, b) in enumerate(zip(x, y)): compare(a, b, f'{path}[{i}]')
    elif isinstance(x, bool):
        if type(y) is not bool or x != y: raise AssertionError(path)
    elif isinstance(x, int):
        if type(y) is not int or x != y: raise AssertionError(path)
    elif isinstance(x, float):
        if not isinstance(y, (int, float)) or isinstance(y, bool) or not np.isclose(x, y, rtol=1e-10, atol=1e-11, equal_nan=False): raise AssertionError((path, x, y))
    elif type(x) is not type(y) or x != y: raise AssertionError((path, x, y))


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true'); group.add_argument('--write-new', action='store_true')
    args = parser.parse_args()
    result = run()
    result['code_sha256'] = sha(__file__)
    result['admission_sha256'] = sha(HERE.parent/'_admission'/'round1320_admission.md')
    if args.write_new:
        with OUT.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write('\n')
    else:
        compare(result, json.loads(OUT.read_text(encoding='utf-8')))
    print(json.dumps({k: result[k] for k in ['status', 'assertions', 'categories', 'max_normalized_frobenius_error']}, ensure_ascii=False))


if __name__ == '__main__': main()
