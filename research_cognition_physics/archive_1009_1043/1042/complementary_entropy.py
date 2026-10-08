"""Finite complementary recovery, central entropy and residual encoding freedom.

No geometric interpretation is used by this calculation. Default: reproduce and
compare; --write exclusively creates the result. Python + existing NumPy only.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'complementary_entropy_results.json'
INPUTS = [
    'archive_370_428/research_note_391.md',
    'archive_301_341/research_note_305.md',
    'archive_301_341/research_note_312.md',
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_990_1008/research_note_1008.md',
    'archive_1009_/1041/dependency_delta_1041.json',
]


def entropy(rho):
    ev = np.linalg.eigvalsh((rho + rho.conj().T)/2)
    assert ev.min() > -1e-12
    ev = ev[ev > 1e-14]
    return float(-np.dot(ev, np.log(ev)))


def log_matrix(rho):
    ev, u = np.linalg.eigh((rho + rho.conj().T)/2)
    assert ev.min() > 1e-12, ev.min()
    return (u * np.log(ev)) @ u.conj().T


def relative_entropy(rho, sigma):
    return float(-entropy(rho) - np.trace(rho @ log_matrix(sigma)).real)


def binary_entropy(p):
    return entropy(np.diag([p, 1-p]))


def encoding(ps):
    v = np.zeros((64, 8), complex)
    for alpha in range(2):
        for a in range(2):
            for b in range(2):
                col = 4*alpha + 2*a + b
                for k in range(2):
                    ia = 4*alpha + 2*a + k
                    ib = 4*alpha + 2*b + k
                    v[8*ia + ib, col] = np.sqrt(ps[alpha] if k == 0 else 1-ps[alpha])
    return v


def boundary(rho, side):
    t = rho.reshape(8, 8, 8, 8)
    return np.einsum('abcb->ac', t) if side == 'A' else np.einsum('abad->bd', t)


def logical_reduction(rho, side):
    """The canonical direct-sum algebra state, normalized in dimension four."""
    out = np.zeros((4, 4), complex)
    for alpha in range(2):
        block = rho[4*alpha:4*alpha+4, 4*alpha:4*alpha+4].reshape(2,2,2,2)
        reduced = np.einsum('abcb->ac', block) if side == 'A' else np.einsum('abad->bd', block)
        out[2*alpha:2*alpha+2, 2*alpha:2*alpha+2] = reduced
    return out


def decode(rho):
    out = np.zeros((4,4), complex)
    for alpha in range(2):
        block = rho[4*alpha:4*alpha+4, 4*alpha:4*alpha+4].reshape(2,2,2,2)
        out[2*alpha:2*alpha+2, 2*alpha:2*alpha+2] = np.einsum('akbk->ab', block)
    return out


def algebra_modular(sigma, side):
    red = logical_reduction(sigma, side)
    out = np.zeros((8,8), complex)
    for alpha in range(2):
        km = -log_matrix(red[2*alpha:2*alpha+2, 2*alpha:2*alpha+2])
        out[4*alpha:4*alpha+4, 4*alpha:4*alpha+4] = (
            np.kron(km, np.eye(2)) if side == 'A' else np.kron(np.eye(2), km))
    return out


def kraus_products(v, side):
    t = v.reshape(8,8,8)
    ks = [t[:,j,:] for j in range(8)] if side == 'A' else [t[j,:,:] for j in range(8)]
    flat = np.stack([(a.conj().T @ b).reshape(-1) for a in ks for b in ks])
    _, s, vh = np.linalg.svd(flat, full_matrices=False)
    basis = vh[s > 1e-10].reshape(-1,8,8)
    # Column-major vectorization: vec(XG-GX)=(G^T x I-I x G)vec(X).
    constraints = np.concatenate([np.kron(g.T,np.eye(8))-np.kron(np.eye(8),g) for g in basis])
    singular = np.linalg.svd(constraints, compute_uv=False)
    return int(np.sum(s > 1e-10)), int(64-np.sum(singular > 1e-10))


def run():
    rng = np.random.default_rng(1042)
    states = []
    for _ in range(14):
        x = rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
        rho = x@x.conj().T + np.eye(8)
        states.append(rho/np.trace(rho))
    psi = np.zeros(8,complex); psi[0]=1/np.sqrt(2); psi[7]=1j/np.sqrt(2)
    states += [np.outer(psi,psi.conj()), np.diag([1.,0,0,0,0,0,0,0])]
    sigma = states[0]
    families = [(0.5,0.25),(0.25,0.5),(0.125,0.375)]
    errors = dict(isometry=0., decoded_channel=0., entropy_identity=0.,
                  relative_entropy=0., modular_identity=0., matrix_unit_reconstruction=0.)
    records = []
    for ps in families:
        v = encoding(ps)
        errors['isometry'] = max(errors['isometry'], float(np.linalg.norm(v.conj().T@v-np.eye(8))))
        ls = [binary_entropy(p) for p in ps]
        area = np.diag(np.repeat(ls,4))
        ranks = {}
        for side in ('A','B'):
            ranks[side] = list(kraus_products(v,side))
            assert ranks[side] == [8,8]
            sig_out = boundary(v@sigma@v.conj().T,side)
            kphys = -log_matrix(sig_out)
            lifted = np.kron(kphys,np.eye(8)) if side=='A' else np.kron(np.eye(8),kphys)
            errors['modular_identity'] = max(errors['modular_identity'], float(np.linalg.norm(
                v.conj().T@lifted@v-area-algebra_modular(sigma,side))))
            for alpha in range(2):
                for i in range(2):
                    for j in range(2):
                        o = np.zeros((2,2)); o[i,j]=1
                        logical = np.zeros((8,8)); physical = np.zeros((8,8))
                        logical[4*alpha:4*alpha+4,4*alpha:4*alpha+4] = (
                            np.kron(o,np.eye(2)) if side=='A' else np.kron(np.eye(2),o))
                        physical[4*alpha:4*alpha+4,4*alpha:4*alpha+4] = np.kron(o,np.eye(2))
                        lifted_o = np.kron(physical,np.eye(8)) if side=='A' else np.kron(np.eye(8),physical)
                        errors['matrix_unit_reconstruction'] = max(errors['matrix_unit_reconstruction'],
                            float(np.linalg.norm(lifted_o@v-v@logical)))
            for rho in states:
                red = logical_reduction(rho,side)
                phys = boundary(v@rho@v.conj().T,side)
                errors['decoded_channel'] = max(errors['decoded_channel'],float(np.linalg.norm(decode(phys)-red)))
                errors['entropy_identity'] = max(errors['entropy_identity'],abs(
                    entropy(phys)-entropy(red)-np.trace(rho@area).real))
                errors['relative_entropy'] = max(errors['relative_entropy'],abs(
                    relative_entropy(phys,sig_out)-relative_entropy(red,logical_reduction(sigma,side))))
        records.append(dict(p=list(ps), central_entropy=ls, kraus_product_rank_and_correctable_dimension=ranks))
    assert max(errors.values()) < 2e-12, errors

    # Same representations/decoders, different physical edge spectra.
    rho = states[-1]
    vp,vq = encoding(families[0]),encoding(families[1])
    ap,aq = boundary(vp@rho@vp.conj().T,'A'), boundary(vq@rho@vq.conj().T,'A')
    distance = float(np.sum(abs(np.linalg.eigvalsh(ap-aq)))/2)
    assert abs(distance-0.25)<1e-14
    entropy_gap = entropy(ap)-entropy(aq)
    exact_gap = 0.75*np.log(3)-np.log(2)
    assert abs(entropy_gap-exact_gap)<1e-14 and exact_gap>0
    # An independently specified edge energy; its mean changes with the resource.
    edge_h = np.diag(np.tile([0.,1.],4))
    energy_gap = 2*np.trace((aq-ap)@edge_h).real
    assert abs(energy_gap-0.5)<1e-14

    # An actual noisy encoding with a known distance to this fixed isometry.
    noise_checks, worst_ratio = 0,0.
    for eps in (0.001,0.01,0.05):
        bound = eps*np.log(7)+binary_entropy(eps)
        for rho in states:
            full = vp@rho@vp.conj().T
            noisy = (1-eps)*full+eps*np.eye(64)/64
            for side in ('A','B'):
                ideal, out = boundary(full,side), boundary(noisy,side)
                err = abs(entropy(out)-entropy(ideal))
                assert err <= bound+1e-12
                worst_ratio = max(worst_ratio,err/bound)
                noise_checks += 1
    return dict(round=1042, date='2026-10-08', hypothesis='complementary_recovery_selects_form_not_central_entropy_value',
                code_dimension=8, boundary_dimensions=[8,8], families=records,
                states_per_family=len(states), errors=errors,
                physical_edge_trace_distance=distance, central_entropy_gap=entropy_gap,
                exact_gap='3*log(3)/4-log(2)', independently_declared_total_edge_energy_gap=energy_gap,
                noisy_entropy_checks=noise_checks, largest_error_to_bound=worst_ratio,
                new_scientific_calibration_groups=1, new_adopted_cognitive_axioms=0,
                full_common_parent_model_established=False, geometry_certified=False, goal_completed=False,
                historical_sha256={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in INPUTS})


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False); f.write('\n')
    else:
        assert result == json.loads(OUT.read_text('utf8')), 'Scientific result changed'
    print(json.dumps(dict(round=1042,passed=True,maximum_residual=max(result['errors'].values()),
                         physical_distance=result['physical_edge_trace_distance'],
                         entropy_gap=result['central_entropy_gap'], noisy_checks=result['noisy_entropy_checks'])))
