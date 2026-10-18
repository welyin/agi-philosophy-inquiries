"""Independent finite Fock-sector audit; does not import the author package.

Default/--check are read-only. --write-new exclusively creates this reviewer's
independent_results.json. All matrices are built from the mode monomials before
the logical codes are used. Finite diagnostics do not replace all-d/all-partner
proofs or grant preparation/measurement permissions.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(13130071)
TOL = 4e-11


def unitary(d: int, special: bool = False) -> np.ndarray:
    z = RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d))
    q, r = np.linalg.qr(z)
    q = q @ np.diag(np.exp(-1j * np.angle(np.diag(r))))
    if special:
        q[:, 0] /= np.linalg.det(q)
    return q


def epsilon(a: int, b: int, c: int) -> int:
    if len({a, b, c}) != 3:
        return 0
    t = (a, b, c)
    return (-1) ** sum(t[i] > t[j] for i in range(3) for j in range(i + 1, 3))


def occ(d: int, index: int | None = None) -> tuple[int, ...]:
    return tuple(int(k == index) for k in range(d))


def sector(kind: str, d: int, n: int):
    """State = (fermion side, activity index/-1, flavor, boson occupation)."""
    basis = []
    if kind in ('A', 'B'):
        for a, b, j in itertools.product(range(d), range(d), range(n)):
            basis.append(('L', a, j, occ(d, b)))
        if kind == 'A':
            basis += [('R', -1, i, occ(d)) for i in range(n)]
        else:
            for a, b, i in itertools.product(range(d), range(d), range(n)):
                basis.append(('R', a, i, occ(d, b)))
    else:
        assert kind == 'C' and d == 3
        for a, b, c, j in itertools.product(range(d), range(d), range(d), range(n)):
            basis.append(('L', a, j, occ(d, b) + occ(d, c)))
        basis += [('R', -1, i, (0,) * (2 * d)) for i in range(n)]
    return basis


def build(kind: str, d: int, m: np.ndarray):
    n = m.shape[0]
    assert m.shape == (n, n)
    basis = sector(kind, d, n)
    index = {s: k for k, s in enumerate(basis)}
    h = np.zeros((len(basis), len(basis)), complex)
    transitions = 0
    # Each monomial is r^dagger (boson annihilators) l and its adjoint.
    terms = []
    for i, j in itertools.product(range(n), range(n)):
        if kind == 'A':
            for a in range(d):
                terms.append((a, -1, i, j, (a,), m[i, j] / np.sqrt(d)))
        elif kind == 'B':
            for a in range(d):
                terms.append((a, a, i, j, (), m[i, j]))
        else:
            for a, b, c in itertools.permutations(range(3)):
                terms.append((a, -1, i, j, (b, d + c),
                              epsilon(a, b, c) * m[i, j] / np.sqrt(6)))
    for col, state in enumerate(basis):
        side, active, flavor, bos = state
        for left_active, right_active, i, j, annih, coefficient in terms:
            if side == 'L' and active == left_active and flavor == j:
                counts = list(bos)
                amp = coefficient
                for b in annih:
                    amp *= np.sqrt(counts[b])
                    counts[b] -= 1
                if amp != 0:
                    target = ('R', right_active, i, tuple(counts))
                    assert target in index, (kind, state, target)
                    h[index[target], col] += amp
                    transitions += 1
            if side == 'R' and active == right_active and flavor == i:
                counts = list(bos)
                amp = coefficient.conjugate()
                for b in reversed(annih):
                    counts[b] += 1
                    amp *= np.sqrt(counts[b])
                if amp != 0:
                    target = ('L', left_active, j, tuple(counts))
                    assert target in index, (kind, state, target)
                    h[index[target], col] += amp
                    transitions += 1
    jl = np.zeros((len(basis), n), complex)
    jr = np.zeros_like(jl)
    for row, (side, active, flavor, bos) in enumerate(basis):
        if side == 'L':
            if kind in ('A', 'B'):
                if bos[active] == 1:
                    jl[row, flavor] = 1 / np.sqrt(d)
            else:
                b, c = bos[:d].index(1), bos[d:].index(1)
                jl[row, flavor] = epsilon(active, b, c) / np.sqrt(6)
        elif kind in ('A', 'C'):
            jr[row, flavor] = 1
        elif bos[active] == 1:
            jr[row, flavor] = 1 / np.sqrt(d)
    nr = np.diag([int(s[0] == 'R') for s in basis]).astype(complex)
    nl = np.eye(len(basis)) - nr
    nb1 = np.diag([sum(s[3][:d]) for s in basis]).astype(complex)
    nb2 = None if kind != 'C' else np.diag([sum(s[3][d:]) for s in basis]).astype(complex)
    return basis, h, jl, jr, nl, nr, nb1, nb2, transitions


def block_diag(a, b):
    out = np.zeros((a.shape[0] + b.shape[0], a.shape[1] + b.shape[1]), complex)
    out[:a.shape[0], :a.shape[1]] = a
    out[a.shape[0]:, a.shape[1]:] = b
    return out


def norm(x):
    return float(np.linalg.norm(x))


def logical_blocks(m, t):
    u, s, vh = np.linalg.svd(m)
    a = (vh.conj().T * np.cos(t * s)) @ vh
    b = -1j * (u * np.sin(t * s)) @ vh
    return a, b, s


def evolve(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def report():
    residuals = {}
    diagnostics = []
    def record(label, value):
        value = float(value)
        residuals[label] = max(value, residuals.get(label, 0.0))
        assert np.isfinite(value) and value < TOL, (label, value)
    for n in (1, 2, 3):
        left, right = unitary(n), unitary(n)
        s = np.linspace(0.6, 1.4, n)
        m = left @ np.diag(s) @ right.conj().T
        for kind, d in [('A', 2), ('A', 3), ('A', 4), ('B', 2), ('B', 3), ('B', 4), ('C', 3)]:
            basis, h, jl, jr, nl, nr, nb1, nb2, transitions = build(kind, d, m)
            dim = len(basis)
            expected_dim = (d*d + 1)*n if kind == 'A' else (2*d*d*n if kind == 'B' else 28*n)
            assert dim == expected_dim
            ident = np.eye(dim)
            record('parent_hermiticity', norm(h-h.conj().T))
            record('input_isometry', norm(jl.conj().T@jl-np.eye(n)))
            record('output_isometry', norm(jr.conj().T@jr-np.eye(n)))
            record('code_orthogonality', norm(jl.conj().T@jr))
            record('actual_forward_intertwining', norm(h@jl-jr@m))
            record('actual_reverse_intertwining', norm(h@jr-jl@m.conj().T))
            p = jl@jl.conj().T + jr@jr.conj().T
            record('code_no_leakage', norm((ident-p)@h@p))
            if kind in ('A', 'C'):
                record('dark_complement', norm(h@(ident-p)))
            record('fermion_number_conserved', norm(h@(nl+nr)-(nl+nr)@h))
            charge1 = nb1 + nr if kind != 'B' else nb1
            record('first_boson_charge_conserved', norm(h@charge1-charge1@h))
            if nb2 is not None:
                record('second_boson_charge_conserved', norm(h@(nb2+nr)-(nb2+nr)@h))
            q = unitary(d, special=True)
            representation = np.kron(np.kron(q, q.conj()), np.eye(n))
            if kind == 'C':
                representation = np.kron(np.kron(np.kron(q, q), q), np.eye(n))
            right_rep = representation if kind == 'B' else np.eye(n)
            full_rep = block_diag(representation, right_rep)
            record('full_group_invariance', norm(full_rep@h-h@full_rep))
            record('input_group_singlet', norm(full_rep@jl-jl))
            record('output_group_singlet', norm(full_rep@jr-jr))
            e_l, e_b1, e_b2 = 2.0, 0.7, 0.4
            e_r = e_l + (e_b1 if kind != 'B' else 0) + (e_b2 if kind == 'C' else 0)
            h0 = e_l*nl + e_r*nr + e_b1*nb1
            if nb2 is not None:
                h0 += e_b2*nb2
            record('free_energy_resonance', norm(h0@h-h@h0))
            # A deliberate detuning must be caught by the energy commutator.
            assert norm((h0+0.17*nr)@h-h@(h0+0.17*nr)) > 0.05
            z = n + 1
            v = RNG.normal(size=n*z) + 1j*RNG.normal(size=n*z)
            v /= np.linalg.norm(v)
            # A Choi input tests the full linear map; a distinct random partner
            # witnesses that the physical branch output is handled directly.
            phi = np.eye(n, dtype=complex).reshape(-1)/np.sqrt(n)
            for t in (0.0, 0.31, 0.87, 1.47):
                ut = evolve(h, t)
                a, b, spec = logical_blocks(m, t)
                target = jl@a + jr@b
                actual = ut@jl
                record('full_time_isometry', norm(actual-target))
                record('instrument_completeness', norm(a.conj().T@a+b.conj().T@b-np.eye(n)))
                record('success_whole_branch', norm(nr@actual-jr@b))
                record('failure_whole_branch', norm(nl@actual-jl@a))
                record('arbitrary_partner_witness', norm(np.kron(actual,np.eye(z))@v-np.kron(target,np.eye(z))@v))
                record('choi_map_identity', norm(np.kron(actual,np.eye(n))@phi-np.kron(target,np.eye(n))@phi))
                for projector, j, op in ((nr,jr,b),(nl,jl,a)):
                    physical = np.kron(projector@actual, np.eye(z))@v
                    formula = np.kron(j@op, np.eye(z))@v
                    record('branch_density_partner', norm(np.outer(physical,physical.conj())-np.outer(formula,formula.conj())))
                success = nr@actual
                expected_ref = 1.0 if kind == 'B' else 0.0
                record('reference_success_number', norm(nb1@success-expected_ref*success))
                if nb2 is not None:
                    record('second_reference_consumed', norm(nb2@success))
                record('failure_reference_retained', norm(nb1@(nl@actual)-nl@actual))
                if kind == 'B':
                    state = np.kron(success,np.eye(z))@v
                    rightstate = state[d*d*n*z:].reshape(d,d,n,z)
                    # Trace the b index; right activity remains maximally mixed.
                    reduced = np.einsum('abin,cbjm->aincjm',rightstate,rightstate.conj()).reshape(d*n*z,d*n*z)
                    logical = np.kron(b,np.eye(z))@v
                    expected = np.kron(np.eye(d)/d,np.outer(logical,logical.conj()))
                    record('charged_receiver_flavor_access', norm(reduced-expected))
                initial = np.kron(jl,np.eye(z))@v
                evolved = np.kron(actual,np.eye(z))@v
                total_energy = np.kron(h0+h,np.eye(z))
                record('complete_energy_expectation', abs(np.vdot(evolved,total_energy@evolved)-np.vdot(initial,total_energy@initial)))
            diagnostics.append({'branch':kind,'d':d,'flavors':n,'full_sector_dimension':dim,
                                'mode_transitions':transitions,
                                'complement_coupling_norm':norm(h@(ident-p))})
    # Strong task distinctions: isometric transfer, generic filtering, rank loss.
    task_checks=[]
    for n in (2,3):
        u=unitary(n)
        m=0.8*u
        t=np.pi/(2*0.8)
        for kind,d in [('A',3),('B',3),('C',3)]:
            _,h,jl,jr,_,_,_,_,_=build(kind,d,m)
            record('complete_transfer_all_roles', norm(evolve(h,t)@jl+1j*jr@u))
        a,b,_=logical_blocks(m,t)
        record('complete_transfer_zero_failure', norm(a))
        record('complete_transfer_unity_success', norm(b.conj().T@b-np.eye(n)))
    a,b,_=logical_blocks(np.diag([1.0,2.0]),np.pi/2)
    record('full_rank_but_bad_time_probability', norm(b.conj().T@b-np.diag([1.0,0.0])))
    task_checks.append('full rank does not imply complete transfer at a fixed time')
    a,b,_=logical_blocks(np.diag([1.0,0.0]),0.9)
    record('rank_deficient_dark_flavor', norm(b[:,1]))
    task_checks.append('Yukawa rank loss leaves dark flavor in this parent model')
    a,b,_=logical_blocks(np.diag([0.5,1.2]),0.8)
    success_gram=b.conj().T@b
    assert abs(success_gram[0,0]-success_gram[1,1])>0.3
    task_checks.append('unequal success singular values filter unknown states')
    # Receiver erasure distance, not the coherent full-output distance.
    for n in (2, 3):
        z=n+1
        u=unitary(n)
        m=0.8*u
        v=RNG.normal(size=n*z)+1j*RNG.normal(size=n*z)
        v/=np.linalg.norm(v)
        rho=np.outer(v,v.conj())
        rho_z=np.trace(rho.reshape(n,z,n,z),axis1=0,axis2=2)
        ideal_v=np.kron(u,np.eye(z))@v
        ideal=block_diag(np.outer(ideal_v,ideal_v.conj()),np.zeros((z,z)))
        for displacement in (-0.31,0.0,0.24):
            t=(np.pi/2+displacement)/0.8
            a,b,_=logical_blocks(m,t)
            success_v=np.kron(b,np.eye(z))@v
            failure_v=np.kron(a,np.eye(z))@v
            failure_rho=np.outer(failure_v,failure_v.conj())
            failure_z=np.trace(failure_rho.reshape(n,z,n,z),axis1=0,axis2=2)
            out=block_diag(np.outer(success_v,success_v.conj()),failure_z)
            distance=0.5*np.sum(np.abs(np.linalg.eigvalsh(out-ideal)))
            q=np.cos(0.8*t)**2
            record('receiver_erasure_distance_witness',abs(distance-q))
            record('failure_partner_state',norm(failure_z-q*rho_z))
    # Independent unstructured finite-matrix perturbations test the stated
    # Duhamel estimate even when the perturbation mixes code and complement.
    _,h,jl,_,_,_,_,_,_=build('A',3,np.diag([0.8,1.1]))
    z=RNG.normal(size=h.shape)+1j*RNG.normal(size=h.shape)
    delta=z+z.conj().T
    eta=0.03
    delta*=eta/np.linalg.norm(delta,2)
    for t in (0.2,1.3,2.7):
        difference=evolve(h+delta,t)-evolve(h,t)
        record('duhamel_operator_excess',max(0.0,np.linalg.norm(difference,2)-t*eta))
        record('duhamel_code_excess',max(0.0,np.linalg.norm(difference@jl,2)-t*eta))
    task_checks.append('receiver erasure and bounded parent perturbations retain declared error bounds')
    # Perturb one group contraction coefficient: the audit must reject it.
    _,h,jl,jr,_,_,_,_,_=build('A',3,np.eye(2))
    h_bad=h.copy(); h_bad[18,0]+=0.11; h_bad[0,18]+=0.11
    q=unitary(3,special=True)
    rep=block_diag(np.kron(np.kron(q,q.conj()),np.eye(2)),np.eye(2))
    assert norm(rep@h_bad-h_bad@rep)>0.02
    task_checks.append('broken contraction is detected by group invariance test')
    sources={}
    for p in (HERE/'independent_check.py',HERE.parent/'_admission/round1313_admission.md'):
        sources[str(p.relative_to(HERE.parent))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return {'status':'PASS','tolerance':TOL,'max_residual':max(residuals.values()),
            'residuals':residuals,'cases':diagnostics,'task_distinctions':task_checks,
            'source_sha256':sources,
            'scope':['Finite global invariant mode model with declared relational preparation.',
                     'No arbitrary bare activity-state decoding or free preparation instrument.',
                     'No Standard Model, anomaly, local Gauss or field-continuum certification.',
                     'A/C dark complement differs from B; only code instruments are compared.',
                     'Preparation, timing, final effects and post-handoff isolation remain explicit contracts.']}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    assert not (args.check and args.write_new)
    out=report()
    destination=HERE/'independent_results.json'
    if args.write_new:
        with destination.open('x',encoding='utf-8') as f:
            json.dump(out,f,ensure_ascii=False,indent=2); f.write('\n')
    elif destination.exists():
        saved=json.loads(destination.read_text(encoding='utf-8'))
        assert saved==out, 'Independent result drift; do not overwrite the saved review.'
    print(json.dumps({'status':out['status'],'cases':len(out['cases']),
                      'checks':len(out['residuals']),'max_residual':out['max_residual'],
                      'written':args.write_new},ensure_ascii=False))


if __name__=='__main__':
    main()
