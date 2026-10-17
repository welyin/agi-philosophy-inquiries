"""Round 1281: finite positive comparison and independently fixed calibration.

Default is a read-only replay. --write creates results.json only if absent.
All universal statements are proved in proof.md; finite matrices audit them.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]


def norm(x):
    return float(np.linalg.norm(x, 2))


def root_positive(x, inverse=False):
    w, v = np.linalg.eigh(x)
    assert w.min() > 0
    return (v * (1 / np.sqrt(w) if inverse else np.sqrt(w))) @ v.conj().T


def encoding(r):
    d = len(r)
    q = np.eye(d) + r @ r
    return np.vstack((1j * r, np.eye(d))) @ root_positive(q, True)


def block(a, c, k):
    return np.block([[a, -1j*c], [1j*c.conj().T, k]])


def exp_unitary(h, t):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j*t*w)) @ v.conj().T


def energy(h, b, rho):
    return float(np.trace(rho @ b.conj().T @ h @ b).real)


def assert_close(x, y, tol=2e-10):
    assert norm(np.atleast_2d(x-y)) <= tol


def calculate():
    eye = np.eye(2, dtype=complex)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    z = np.diag([1, -1]).astype(complex)
    r, a, rho = x, eye, eye/2
    b = encoding(r)
    h = block(a, r, eye)
    assert_close(b.conj().T @ b, eye)
    assert_close(h @ b, np.zeros((4, 2)))
    assert np.linalg.eigvalsh(h).min() >= -1e-12
    assert_close(exp_unitary(h, 0.7) @ b, b)
    # Explicit passive partner: no tracing out the logical encoding.
    bell = np.eye(2).reshape(-1)/np.sqrt(2)
    encoded = np.kron(b, eye) @ bell
    evolved = np.kron(exp_unitary(h, 0.7), eye) @ encoded
    assert np.linalg.norm(evolved - encoded) < 1e-12

    bad_c = x + z/4
    bad_k = 9*eye/8
    bad_h = block(a, bad_c, bad_k)
    assert np.linalg.eigvalsh(bad_h).min() > 0
    assert abs(energy(bad_h, b, rho) - 1/16) < 1e-12
    assert_close(bad_k-bad_c.conj().T@bad_c, eye/16)
    evolved_bad = exp_unitary(bad_h, 0.7) @ b
    for vec in (np.array([1, 0]), np.array([1, 1])/np.sqrt(2)):
        p = b @ vec
        p1 = exp_unitary(bad_h, 0.7) @ p
        assert abs(np.vdot(p, bad_h @ p) - np.vdot(p1, bad_h @ p1)) < 1e-12

    scaling=[]
    for t in (1/2, 1/4, 1/8, 1/16):
        c = x + t*z
        k = (1+t*t)*eye
        ht = block(eye, c, k)
        nt = b.conj().T @ ht @ b
        assert_close(nt, (t*t/2)*eye)
        assert abs(norm(c-r) - np.sqrt(2*energy(ht,b,rho))) < 1e-12
        scaling.append(dict(t=t, epsilon=t*t/2, coupling_error=t))

    # Removing positivity permits cancellation at exactly zero mean.
    indefinite = block(eye, np.zeros((2, 2)), -eye)
    assert_close(b.conj().T@indefinite@b, np.zeros((2, 2)))
    assert norm(indefinite@b) > 0.9
    assert np.linalg.eigvalsh(indefinite).min() == -1
    # A rank-one input misses a genuine second logical direction.
    rr = np.diag([1.,2.]); cc = np.diag([1.,3.])
    br = encoding(rr); hr = block(eye, cc, cc@cc)
    nr = br.conj().T@hr@br
    assert_close(nr, np.diag([0.,1/5]))
    assert abs(energy(hr,br,np.diag([1.,0.]))) < 1e-12
    # If R is allowed to follow C, incompatible fixed references both pass.
    for chosen in (x,z):
        own_b=encoding(chosen)
        own_h=block(eye,chosen,eye)
        assert_close(own_h@own_b,np.zeros((4,2)))
    assert energy(block(eye,z,eye),b,rho) > 0.9

    rng = np.random.default_rng(1281)
    def herm(d):
        m=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
        return (m+m.conj().T)/2
    diagnostics=[]
    for d in (2,3,4):
        for sample in range(8):
            rr=herm(d); raw=herm(d)
            aa=raw@raw+np.eye(d)/3
            cc=aa@rr+(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))/7
            mm=herm(d); mm=mm@mm/50
            kk=cc.conj().T@np.linalg.solve(aa,cc)+mm
            hh=block(aa,cc,kk); bb=encoding(rr)
            pp=herm(d); pp=pp@pp+np.eye(d)/2; pp/=np.trace(pp)
            mu=float(np.linalg.eigvalsh(pp).min())
            eps=energy(hh,bb,pp)
            tau=eps*norm(np.eye(d)+rr@rr)/mu
            delta=np.sqrt(norm(aa)*tau)
            e=cc-aa@rr
            yy=np.vstack((1j*rr,np.eye(d)))
            dd=yy.conj().T@hh@yy
            assert_close(dd,mm+e.conj().T@np.linalg.solve(aa,e))
            assert norm(e) <= delta+1e-9
            assert norm(kk-rr@aa@rr) <= tau+2*norm(rr)*delta+1e-9
            assert norm(mm) <= tau+1e-9
            assert np.linalg.eigvalsh(hh).min() >= -1e-10
            diagnostics.append((norm(e)/delta,norm(mm)/tau))
    # Exact theorem also covers singular A, hence inverse A is not a premise.
    singular_a=np.diag([1.,0.])
    singular_h=block(singular_a,singular_a@x,x@singular_a@x)
    assert_close(singular_h@b,np.zeros((4,2)))
    assert np.linalg.eigvalsh(singular_h).min() >= -1e-12

    sources=[HERE/'proof.md',Path(__file__).resolve(),
        RESEARCH/'archive_923_934/research_note_933.md',
        RESEARCH/'archive_956_989/research_note_959.md',
        RESEARCH/'archive_956_989/research_note_962.md',
        RESEARCH/'archive_1009_1043/research_note_1016.md']
    return dict(round=1281,status='PASS_WITH_SCOPE',new_axioms=0,
        theorem='positive comparison + independently fixed full-rank coherent zero calibration => C=AR,K=RAR',
        exact_witness=dict(d=2,reference='Pauli X',comparison_eigenvalues=[0,0,2,2],
                           preserves_arbitrary_partner=True,comparison_cost_zero=True),
        positive_mismatch=dict(c='X+Z/4',k='9I/8',schur_remainder='I/16',
                               exact_calibration_energy='1/16',
                               finite_evolution_code_change=round(norm(evolved_bad-b),12)),
        necessity=dict(without_positivity_min_eigenvalue=-1,
                       rank_deficient_missed_energy='1/5',
                       reference_cannot_be_fitted_to_candidate=True),
        square_root_order=scaling,
        numerical_diagnostics=dict(noncommuting_general_cases=len(diagnostics),
             largest_c_error_bound_ratio=round(max(q[0] for q in diagnostics),12),
             largest_schur_bound_ratio=round(max(q[1] for q in diagnostics),12),singular_a_test='PASS'),
        scope=dict(finite_comparison_model=True,all_local_field_variations_assumed=False,
                   full_field_theory_generated=False,gauge_dynamics_generated=False,
                   full_cognitive_axiom_model_certified=False,preparation_and_records_free=False),
        source_sha256={str(p.relative_to(RESEARCH)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=calculate();dest=HERE/'results.json'
    if args.write:
        with dest.open('x',encoding='utf-8') as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED'
    else:
        assert json.loads(dest.read_text('utf-8'))==data,'Saved result/source mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status=data['status'],mode=mode,round=1281,
        calibration_energy=data['positive_mismatch']['exact_calibration_energy'],
        general_cases=data['numerical_diagnostics']['noncommuting_general_cases'],new_axioms=0),ensure_ascii=False))


if __name__=='__main__':main()
