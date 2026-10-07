"""772: same physical covariance, paired BRST representation and mean dictionary.

Matrices calibrate the full inherited Cauchy symbol, the algebraic quartet,
and an original-potential vacuum jet. Infinite-dimensional claims are proved
in the note; no interacting charge, physical preparation or PDE is simulated.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_microlocal_gauge_projection as old
import joint_brst_relative_source as jets
import joint_scalar_propagation_matching as scalar

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_brst_quartet_matching_results.json'


def mx(a):
    return float(np.max(np.abs(a)))


def herm(a):
    return (a+a.conj().T)/2


def cauchy_pairing():
    K, c, c0, q, q0 = old.original_principal_setup()
    _, _, H, H0 = old.make_projection(K, c, c0)
    T = np.linalg.solve(H, q)
    T0 = np.linalg.solve(H0, q0)
    sharp = np.linalg.solve(H0, K.conj().T@H)
    M = sharp@T@T@K
    L0 = T@K@np.linalg.solve(M, T0)
    L = L0-K@np.linalg.solve(q0, L0.conj().T@q@L0)/2
    U = np.linalg.solve(q0, L.conj().T@q)
    V = np.linalg.solve(q0, K.conj().T@q)
    E = np.eye(len(q))-K@U-L@V
    errors = dict(dual_pair=mx(K.conj().T@q@L-q0), dual_isotropic=mx(L.conj().T@q@L),
                  slice_idempotent=mx(E@E-E), slice_charge_adjoint=mx(E.conj().T@q-q@E),
                  left_inverse=mx(U@K-np.eye(len(q0))), constraint_of_slice=mx(V@E),
                  split_identity=mx(K@U+L@V+E-np.eye(len(q))), frequency_compatibility=mx(c@L-L@c0))
    physical = E.conj().T@(q@c+.1*H)@E
    gauge = U.conj().T@q0@c0@V+V.conj().T@q0@c0@U
    paired = herm(physical+gauge)
    rng = np.random.default_rng(772)
    A = herm(.03*rng.normal(size=(32,32)))
    B = .015*rng.normal(size=(32,126))
    off = herm(V.conj().T@A@V+V.conj().T@B@E+E.conj().T@B.T@V)
    prior = paired+off
    delta_form = paired-prior
    qi = np.linalg.inv(q)
    delta_raw = qi@delta_form@qi
    C = U@delta_raw-.5*U@delta_raw@U.conj().T@K.conj().T
    errors.update(unchanged_physical_covariance=mx(E.conj().T@delta_form@E),
                  unchanged_gauge_intertwining=mx(delta_form@K),
                  exact_new_intertwining=mx(np.linalg.solve(q,paired)@K-K@c0),
                  raw_constraint=mx(V@delta_raw),
                  gauge_factorization=mx(delta_raw-K@C-C.conj().T@K.conj().T))
    assert max(errors.values()) < 3e-12
    assert mx(delta_form) > .001
    return dict(errors=errors, changed_unphysical_covariance_norm=mx(delta_form),
                physical_rank=int(np.linalg.matrix_rank(E, tol=1e-10)),
                scope='Inherited full 126x32 Cauchy principal calibration plus a finite smooth-block analogue; it does not compute the original background covariance numerically.')


def quartet_fock():
    # Three bosonic modes (u,v,p) and two odd modes (c,b); p is physical.
    # Finite TOTAL degree is invariant under Q, so there is no cutoff defect.
    cutoff = 4
    basis = [(u,v,p,c,b) for u,v,p in itertools.product(range(cutoff+1), repeat=3)
             for c,b in itertools.product(range(2), repeat=2) if u+v+p+c+b <= cutoff]
    lookup = {x:i for i,x in enumerate(basis)}
    Q = np.zeros((len(basis),len(basis)))
    J = np.zeros_like(Q)
    for j,(u,v,p,c,b) in enumerate(basis):
        if c:
            Q[lookup[(u+1,v,p,0,b)],j] += np.sqrt(u+1)
        if v and not b:
            Q[lookup[(u,v-1,p,c,1)],j] += (-1)**c*np.sqrt(v)
        J[lookup[(v,u,p,b,c)],j] = (-1)**(c*b)
    N = np.diag([u+v+c+b for u,v,p,c,b in basis])
    physical = np.diag([float(u+v+c+b == 0) for u,v,p,c,b in basis])
    h = Q.T@np.diag([1/n if n else 0 for n in np.diag(N)])
    errors = dict(nilpotent=mx(Q@Q), Krein_self_adjoint=mx(Q.T@J-J@Q),
                  Hodge_number=mx(Q@Q.T+Q.T@Q-N), contracting_homotopy=mx(Q@h+h@Q-np.eye(len(basis))+physical),
                  fundamental_symmetry=mx(J@J-np.eye(len(basis))))
    _, s, vh = np.linalg.svd(Q)
    rank = int(np.sum(s > 1e-10))
    kernel = vh[rank:].T
    gram = herm(kernel.T@J@kernel)
    eig = np.linalg.eigvalsh(gram)
    positive = int(np.sum(eig>1e-10))
    null = int(np.sum(np.abs(eig)<1e-10))
    assert max(errors.values()) < 2e-14 and eig[0] > -2e-14
    assert positive == cutoff+1 and null == rank
    return dict(errors=errors, total_degree_cutoff=cutoff, vector_count=len(basis),
                image_rank=rank, kernel_dimension=kernel.shape[1], physical_positive_dimension=positive,
                kernel_null_dimension=null, minimum_kernel_norm_eigenvalue=float(eig[0]),
                scope='An exact finite-degree invariant sector of one quartet and one physical mode. The note proves the arbitrary-mode finite-particle identity; this is not a continuum interacting Hilbert-space completion.')


def original_source_dictionary():
    _, u, _ = scalar.parameters()
    phi = np.array([np.sqrt(u[0]),0.,0.,0.,np.sqrt(u[1])])
    t = np.zeros((5,5)); t[1,0]=.5; t[0,1]=-.5
    k = t@phi
    mean = t@k
    e, hk, half_third = jets.potential_gradient_jet(phi,k)
    _, hm, _ = jets.potential_gradient_jet(phi,mean)
    source_change = 2*half_third  # delta W=2 k tensor k, so delta J=S'''(k,k).
    F = 2-phi@phi/6
    h2 = phi[:4]@phi[:4]; k2=k@k
    record_mean = 2*(phi[:4]@mean[:4])/F+h2*(phi@mean)/(3*F**2)
    record_covariance = 2*k2/F+h2*k2/(3*F**2)
    errors = dict(stationary_original_potential=mx(e), gauge_is_a_linear_solution_of_potential=mx(hk),
                  source_plus_mean=mx(source_change+hm), relational_record_compensation=abs(record_mean+record_covariance))
    assert max(errors.values()) < 2e-14
    assert mx(source_change)>.001 and abs(record_covariance)>.01 and F>0
    return dict(errors=errors, original_vacuum=phi.tolist(), F=float(F),
                source_change=source_change.tolist(), mean_dictionary_change=mean.tolist(),
                invariant_record='(sum of four real Higgs components squared)/F',
                record_mean_change=float(record_mean), record_covariance_change=float(record_covariance),
                scope='Original U=V/F^2 at its stationary potential vacuum, with a pure-gauge covariance jet. This is a normalization/source-dictionary calibration, not the dynamical 753 background or a new actual instrument.')


def run():
    result = dict(round=772, tests_run=3, failures=0, errors=0,
                  complete_Cauchy_pairing=cauchy_pairing(), algebraic_quartet=quartet_fock(),
                  original_source_dictionary=original_source_dictionary())
    deps = ('research_note_753.md','research_note_765.md','research_note_766.md','research_note_767.md',
            'research_note_768.md','research_note_771.md','joint_microlocal_gauge_projection.py',
            'joint_brst_relative_source.py','round772_drafts/research_note_772_working.md')
    result['dependency_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps}
    result['scope']=('An adapted unphysical representative preserves the inherited physical Hadamard state and admits algebraic free BRST quartet positivity. Its smooth source change is compensated by the second-order mean dictionary for invariant observables. It is not a proof of interacting BRST/QME, a full completed charge domain, an actual instrument, finite-strength self-consistency, or original Q-to-continuum equivalence.')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--write-results',action='store_true')
    args=p.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
