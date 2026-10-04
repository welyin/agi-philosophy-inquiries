"""608: locality and independent composition of the declared 607 completion.

Finite qutrits test the *completion rule*, not a replacement standard model.
All-sector counterexamples are distinguished from admissible-sector dynamics.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "joint_projection_local_composition_results.json"
I = np.eye(3, dtype=complex)
P0 = np.diag([1, 1, 0]).astype(complex)
Z = np.diag([1, 0, -1]).astype(complex)
D = np.diag([1, -1, 0]).astype(complex)
X = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], complex)
A = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], complex)


def tensor(xs):
    out = np.ones((1, 1), complex)
    for x in xs:
        out = np.kron(out, x)
    return out


def at(x, i, n):
    return tensor([x if j == i else I for j in range(n)])


def norm(x):
    return float(np.linalg.norm(x, 2))


def comm(a, b):
    return a @ b - b @ a


def projections(n):
    return [at(P0, i, n) for i in range(n)]


def pinch(H, P):
    Q = np.eye(len(H)) - P
    return P @ H @ P + Q @ H @ Q


def local_completion(H, ps):
    for p in ps:
        H = pinch(H, p)
    return H


def unitary(H, t):
    e, v = np.linalg.eigh(H)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def embedding(n):
    v = np.eye(3, dtype=complex)[:, :2]
    return tensor([v] * n)


def receiver(psi, n):
    m = psi.reshape(3, 3 ** (n - 1))
    return m @ m.conj().T


def counterexample():
    rows = []
    t = .275
    for n in (2, 3, 4):
        ps = projections(n)
        P = tensor([P0] * n)
        H = n * np.eye(3**n) - at(X, 0, n)
        Hg = pinch(H, P)
        Hl = local_completion(H, ps)
        R = tensor([I] + [P0] * (n - 1))
        expected = n * np.eye(3**n) - at(X, 0, n) @ (np.eye(3**n) - R)
        assert norm(Hg - expected) < 1e-12
        original = norm(comm(comm(H, at(Z, 0, n)), at(X, n-1, n)))
        global_defect = norm(comm(comm(Hg, at(Z, 0, n)), at(X, n-1, n)))
        local_defect = norm(comm(comm(Hl, at(Z, 0, n)), at(X, n-1, n)))
        assert original < 1e-12 and local_defect < 1e-12
        assert abs(global_defect - 2) < 1e-12
        psi = np.zeros(3**n, complex)
        psi[0] = 1
        remote_flip = at(X + np.diag([0, 1, 0]), n-1, n)
        assert norm(remote_flip.conj().T @ remote_flip - np.eye(3**n)) < 1e-12
        U = unitary(Hg, t)
        diff = receiver(U @ psi, n) - receiver(U @ remote_flip @ psi, n)
        td = float(np.sum(abs(np.linalg.eigvalsh(diff))) / 2)
        assert abs(td - abs(np.sin(t))) < 1e-12
        Q = np.eye(3**n) - P
        rows.append(dict(sites=n,code_dimension=2**n,separation=n-1,
            old_double_commutator=original,global_double_commutator=global_defect,
            local_double_commutator=local_defect,receiver_trace_distance=td,
            global_Q_block_retained_error=norm(Q @ (Hg-H) @ Q),
            local_Q_block_change=norm(Q @ (Hl-H) @ Q),
            common_P_block_error=norm(P @ (Hg-Hl) @ P)))
    return dict(rows=rows,waiting_time=t,exact_trace_distance=float(abs(np.sin(t))),
        encoding_leaves_admissible_sector=True,
        not_a_signal_inside_restricted_physical_sector=True)


def interacting(n):
    dim = 3**n
    terms = []
    for i in range(n):
        terms.append(((i,), 3*np.eye(dim)+.2*at(D,i,n)+.13*at(A,i,n)-.31*at(X,i,n)))
    for i in range(n-1):
        terms.append(((i,i+1),.12*at(A,i,n)@at(A,i+1,n)+.09*at(X,i,n)@at(X,i+1,n)))
    H = sum(h for _,h in terms)
    return H, terms


def composition_check():
    rows = []
    for n in (2,3,4):
        H, terms = interacting(n)
        ps=projections(n);J=embedding(n);P=J@J.conj().T
        Hg=pinch(H,P);Hl=local_completion(H,ps);h=J.conj().T@H@J
        sumlocal=np.zeros_like(H);support_error=0.;growth=0.
        for sites,term in terms:
            local=local_completion(term,[ps[i] for i in sites])
            allsites=local_completion(term,ps)
            assert norm(local-allsites)<1e-12
            sumlocal+=local
            growth=max(growth,norm(local)-norm(term))
            for j in range(n):
                if j not in sites:
                    support_error=max(support_error,norm(comm(local,at(X,j,n))),
                                      norm(comm(local,at(A,j,n))))
        assert norm(sumlocal-Hl)<1e-11 and support_error<1e-12 and growth<1e-12
        mapping=max(norm(Hg@J-J@h),norm(Hl@J-J@h))
        assert mapping<1e-12
        minimum=float(np.linalg.eigvalsh(Hl)[0]);assert minimum>0
        # Completion of independent systems agrees with the tensor sum.
        H1,_=interacting(1);H2,_=interacting(n-1)
        independent=np.kron(H1,np.eye(3**(n-1)))+np.kron(np.eye(3),H2)
        local_ind=local_completion(independent,ps)
        factors=np.kron(local_completion(H1,projections(1)),np.eye(3**(n-1)))+np.kron(np.eye(3),local_completion(H2,projections(n-1)))
        composition=norm(local_ind-factors);assert composition<1e-12
        beta=.3
        Zg=float(np.sum(np.exp(-beta*np.linalg.eigvalsh(Hg))))
        Zl=float(np.sum(np.exp(-beta*np.linalg.eigvalsh(Hl))))
        assert abs(Zg-Zl)>1e-5
        rows.append(dict(sites=n,intertwining_error=mapping,support_error=support_error,
            independent_composition_error=composition,maximum_term_norm_increase=growth,
            minimum_local_energy=minimum,ambient_partition_global=Zg,
            ambient_partition_local=Zl,ambient_completions_differ=norm(Hg-Hl)))
    return dict(rows=rows,sector_Gibbs_same_analytically=True,ambient_Gibbs_not_same=True,
        locality_theorem_requires_actual_local_projectors=True)


def record_and_source_check():
    n=3;gamma=.12;kappa=np.exp(-2*gamma)
    H,_=interacting(n);ps=projections(n);J=embedding(n);P=J@J.conj().T
    Hg=kappa*pinch(H,P);Hl=kappa*local_completion(H,ps)
    h=kappa*(J.conj().T@H@J)
    instruments=[]
    for site in (0,2):
        diag=np.real(np.diag(.5*np.eye(3**n)+.2*at(D,site,n)))
        instruments.append([np.diag(np.sqrt(diag)),np.diag(np.sqrt(1-diag))])
    v=np.zeros((2**n,2),complex)
    v[0,0]=1/np.sqrt(2);v[-1,1]=np.exp(.37j)/np.sqrt(2)
    max_error=0.;state_error=0.;energy_error=0.;total=0.
    for a,b in itertools.product(range(2),repeat=2):
        L1=instruments[0][a];L2=instruments[a][b]
        def branch(K,l1,l2):
            return unitary(K,.06)@l2@unitary(K,.09+.04*a)@l1@unitary(K,.14)
        Kg=branch(Hg,L1,L2);Kl=branch(Hl,L1,L2)
        small=branch(h,J.conj().T@L1@J,J.conj().T@L2@J)
        max_error=max(max_error,norm(Kg@J-J@small),norm(Kl@J-J@small))
        vg=Kg@J@v;vl=Kl@J@v
        diff=np.outer(vg.ravel(),vg.ravel().conj())-np.outer(vl.ravel(),vl.ravel().conj())
        state_error+=float(np.sum(abs(np.linalg.eigvalsh(diff))))
        energy_error=max(energy_error,abs(float(np.trace(vg.conj().T@Hg@vg).real-np.trace(vl.conj().T@Hl@vl).real)))
        total+=float(np.linalg.norm(vl)**2)
    assert max_error<1e-11 and state_error<1e-11 and energy_error<1e-10 and abs(total-1)<1e-11
    source=-2*Hl;eps=1e-6
    fd=(np.exp(-2*(gamma+eps))-np.exp(-2*(gamma-eps)))/(2*eps)*local_completion(H,ps)
    fd_error=norm(fd-source);assert fd_error<2e-8
    L=instruments[0]
    injection=sum(x@Hl@x for x in L)-Hl
    original_injection=sum(x@(kappa*H)@x for x in L)-kappa*H
    injerr=norm(injection-local_completion(original_injection,ps))
    source_injection=sum(x@source@x for x in L)-source
    assert injerr<1e-11 and norm(source_injection+2*injection)<1e-11
    assert norm(injection)>1e-3
    return dict(branches=4,reference_dimension=2,operator_history_error=max_error,
        complete_record_state_trace_norm=state_error,total_probability=total,
        energy_error=energy_error,source_finite_difference_error=fd_error,
        injection_commuting_square_error=injerr,injection_norm=norm(injection),
        arbitrary_reference_guarantee_is_analytic=True)


def run():
    data=dict(counterexample=counterexample(),composition=composition_check(),
              history_and_sources=record_and_source_check())
    deps=('research_note_360.md','research_note_363.md','research_note_466.md',
          'research_note_606.md','research_note_607.md',
          'joint_projected_process_completion.py','cognitive_foundation_bridge_605.md')
    return dict(round=608,tests_run=3,failures=0,errors=0,**data,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(global_completion_not_monoidal=True,
            bounded_local_factor_completion_is_conditional=True,
            original_sector_histories_sources_preserved=True,
            not_actual_GW_factorization_or_unbounded_gauge_LR=True,
            not_derived_dimension_chiral_measure_or_GR=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=608,tests=3,all_passed=True,
        signal=result['counterexample']['exact_trace_distance'],
        histories=result['history_and_sources'])))
