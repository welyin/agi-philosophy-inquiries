"""767: calibrations for the physical slice and common smooth-noise repair.

The infinite-dimensional proof is in note767. Matrices check algebra, not a
numerical solution of the original coupled pseudodifferential operators.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_microlocal_gauge_projection as old

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_physical_hadamard_positivity_results.json'


def mx(a):
    return float(np.max(np.abs(a)))


def herm(a):
    return (a+a.conj().T)/2


def root(a):
    v, u = np.linalg.eigh(herm(a))
    assert v[0] > 0
    return (u*np.sqrt(v))@u.conj().T


def physical_slice():
    K, c, c0, q, q0 = old.original_principal_setup()
    B, Pi, H, H0 = old.make_projection(K, c, c0)
    T = np.linalg.solve(H, q)
    sharp = np.linalg.solve(H0, K.conj().T@H)
    second_gram = sharp@T@T@K
    R = Pi-T@K@np.linalg.solve(second_gram, sharp@T)
    errors = dict(projection=mx(R@R-R), H_self_adjoint=mx(R.conj().T@H-H@R),
                  gauge_removed=mx(R@K), constraint=mx(sharp@T@R),
                  frequency_commutation=mx(R@c-c@R), T_commutation=mx(T@c-c@T))
    U = root(H); Ui = np.linalg.inv(U)
    rh = herm(U@R@Ui)
    vv, basis = np.linalg.eigh(rh)
    E = basis[:, vv > .5]
    ch = U@c@Ui
    ce = herm(E.conj().T@ch@E)
    cv, cu = np.linalg.eigh(ce)
    up, um = E@cu[:, cv>.5], E@cu[:, cv<.5]
    lp = herm(Ui.conj().T@q@c@Ui)
    lm = herm(-Ui.conj().T@q@(np.eye(len(q))-c)@Ui)
    minp = float(np.linalg.eigvalsh(herm(up.conj().T@lp@up))[0])
    minm = float(np.linalg.eigvalsh(herm(um.conj().T@lm@um))[0])
    # lp need not vanish as a full-space vector on physical minus vectors in
    # arbitrary representatives; use the physical compression for the test.
    errors['opposite_branch_plus'] = mx(E.conj().T@lp@um)
    errors['opposite_branch_minus'] = mx(E.conj().T@lm@up)
    assert max(errors.values()) < 3e-13
    assert minp > .01 and minm > .01
    assert E.shape[1] == 62 and up.shape[1] == um.shape[1] == 31
    return dict(errors=errors, physical_representative_rank=E.shape[1],
                positive_frequency_rank=up.shape[1], negative_frequency_rank=um.shape[1],
                positive_branch_minima=[minp, minm],
                scope='One inherited principal calibration of the new exact physical-slice formula. Counts are checks, not a dimension derivation.')


def row_majorant():
    # The corresponding infinite kernel is exp(-0.4*(i+j)) times the bounded
    # Hermitian phase below. Finite truncation checks the row-sum inequality.
    n = 24
    i, j = np.meshgrid(np.arange(1, n+1), np.arange(1, n+1), indexing='ij')
    envelope = np.exp(-.4*(i+j))
    r1 = envelope*(np.cos(.7*(i-j))+1j*np.sin(.31*(i-j)))
    r2 = -.8*envelope*(np.cos(.23*(i-j))+1j*np.sin(.61*(i-j)))
    assert mx(r1-r1.conj().T)<1e-14 and mx(r2-r2.conj().T)<1e-14
    d = np.diag(np.sum(np.abs(r1)+np.abs(r2), axis=1))
    minima = [float(np.linalg.eigvalsh(d+sign*r)[0])
              for r in (r1, r2) for sign in (-1, 1)]
    assert min(minima) > 0
    p = np.diag(np.r_[np.ones(n//2), np.zeros(n//2)])
    aplus = p+r1; aminus = np.eye(n)-p+r2
    before = [float(np.linalg.eigvalsh(a)[0]) for a in (aplus, aminus)]
    after = [float(np.linalg.eigvalsh(a+d)[0]) for a in (aplus, aminus)]
    assert min(before) < -1e-5 and min(after)>0
    error = mx((aplus+d)-(aminus+d)-(aplus-aminus))
    assert error < 1e-14
    return dict(majorant_minima=minima, variance_minima_before=before,
                variance_minima_after=after, unchanged_difference_error=error,
                scope='Finite check of a proven infinite-rank smoothing domination lemma, not an eigenvalue bound for the full original PDE.')


def exact_gauge_reality_repair():
    q = np.array([[0, 1, 0, 0], [1, 0, 0, 0],
                  [0, 0, 1, 0], [0, 0, 0, -1.]], complex)
    K = np.eye(4, dtype=complex)[:, :1]
    c = np.diag([1., 1., 1., 0.]).astype(complex)
    c0 = np.ones((1, 1), complex)
    B, Pi, *_ = old.make_projection(K, c, c0)
    v = np.array([1., -1., 0., 0.])/np.sqrt(2)
    w = np.array([0., 0., 0., 1.]); theta=.37
    U = (np.eye(4)+(np.cos(theta)-1)*(np.outer(v,v)+np.outer(w,w))
         +np.sin(theta)*(np.outer(w,v)-np.outer(v,w)))
    ct, _ = old.repair(U@c@U.T, c0, K, B, Pi, q, np.ones((1,1)))
    lp, lm = herm(q@ct), herm(-q@(np.eye(4)-ct))
    E = np.eye(4)[:, 2:]
    # Use a common dominating variance; no exact-Hamiltonian operation is
    # asserted by this covariance choice.
    s = float(np.sin(theta)**2)
    S = np.diag([0., 0., s, s]).astype(complex)
    b = old.block
    qp = b(q, -q.conj()); kp=b(K,K.conj())
    plus=b(lp,lm.conj()); minus=b(lm,lp.conj()); noise=b(S,S.conj())
    eye=np.eye(4); C=np.block([[np.zeros((4,4)),eye],[eye,np.zeros((4,4))]])
    ep=b(E,E)
    repaired_plus, repaired_minus = plus+noise, minus+noise
    errors=dict(full_CCR=mx(repaired_plus-repaired_minus-qp),
                noise_annihilates_gauge=mx(noise@kp),
                real_structure=mx(C@repaired_plus.conj()@C-repaired_minus),
                same_gauge_identity=mx((repaired_plus-plus)@kp),
                physical_gauge_pairing=mx(ep.conj().T@repaired_plus@kp))
    before=[float(np.linalg.eigvalsh(herm(ep.conj().T@a@ep))[0]) for a in (plus,minus)]
    after=[float(np.linalg.eigvalsh(herm(ep.conj().T@a@ep))[0]) for a in (repaired_plus,repaired_minus)]
    assert max(errors.values())<1e-14 and min(before)<-.1 and min(after)>-1e-14
    return dict(errors=errors, physical_variance_minima_before=before,
                physical_variance_minima_after=after,
                scope='Doubled real-structure completion of the 766 finite failure block. It checks simultaneous repair, not a numerical Hadamard existence proof.')


def run():
    result=dict(round=767, tests_run=3, failures=0, errors=0,
                physical_slice=physical_slice(), smoothing_majorant=row_majorant(),
                exact_gauge_reality_repair=exact_gauge_reality_repair())
    deps=('research_note_753.md','research_note_765.md','research_note_766.md',
          'joint_microlocal_gauge_projection.py','round767_drafts/STATUS.md')
    result['dependency_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps}
    result['scope']=('Analytic existence of a positive physical Hadamard quasifree state for the '
                     'inherited complete linear bosonic theory on the fixed compact irreducible '
                     'background, by common smoothing covariance addition. Not purity, a unique '
                     'state, a resource-bounded preparation, nonlinear quantum gravity, all-sector '
                     'renormalized Ward identities, or equivalence with the original Q process. '
                     'Numerics calibrate algebraic lemmas only.')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
