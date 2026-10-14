"""Finite internal-budget lift: exact permutations, rational bounds, NumPy diagnostics.

No large witness matrix is materialized. --write-results is explicit; default verifies.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import factorial, log, sqrt
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entropy(matrix):
    eigen = np.linalg.eigvalsh(matrix)
    assert eigen.min() > -1e-12
    eigen = eigen[eigen > 1e-14]
    return float(-np.sum(eigen * np.log(eigen)))


def exp_bounds(x, terms=80):
    """Exact positive Taylor sum plus a geometric majorant of its tail."""
    x = F(x)
    assert x >= 0 and x < terms + 2
    s = sum((x**k / factorial(k) for k in range(terms + 1)), F(0))
    first = x**(terms + 1) / factorial(terms + 1)
    return s, s + first / (1 - x / (terms + 2))


def small_case(blocks, length, with_record, beta=None):
    energies = (0, 100, 201)  # units 1/100; x=02, y=11 differ by 1
    words = list(product(range(3), repeat=2*blocks))
    lookup = {w:i for i,w in enumerate(words)}
    ns, nb, nk = len(words), length+2, 2 if with_record else 1
    na = nb*nk
    target, category = [], []
    for w in words:
        v = list(w); label = 0
        for j in range(blocks):
            pair = tuple(w[2*j:2*j+2])
            if pair in [(0,2), (1,1)]:
                label = 1 if pair == (0,2) else -1
                v[2*j:2*j+2] = (1,1) if label == 1 else (0,2)
                break
        category.append(label); target.append(lookup[tuple(v)])
    energy = [sum(energies[j] for j in w) for w in words]
    # Complete global permutation, including unused ladder endpoints and both bits.
    perm = []
    for s in range(ns):
        for k in range(nb):
            for bit in range(nk):
                ss, kk = s, k
                shift = category[s]
                if shift and 0 <= k+shift < nb:
                    ss, kk = target[s], k+shift
                outbit = bit ^ int(category[ss] == 1) if with_record else 0
                perm.append((ss*nb+kk)*nk+outbit)
                assert energy[s]+k == energy[ss]+kk
    assert sorted(perm) == list(range(ns*na))
    # Isometry from arbitrary material input; pointer starts at zero.
    V = np.zeros((ns*na, ns))
    J = np.zeros_like(V)
    for s in range(ns):
        for k in range(1,length+1):
            V[perm[(s*nb+k)*nk],s] = 1/sqrt(length)
            bit = int(category[target[s]] == 1) if with_record else 0
            J[(target[s]*nb+k)*nk+bit,s] = 1/sqrt(length)
    expected = np.diag([F(2,length) if cat else 0 for cat in category]).astype(float)
    assert np.max(np.abs((V-J).T@(V-J)-expected)) < 1e-12
    assert np.max(np.abs(V.T@V-np.eye(ns))) < 1e-12
    if beta is None:
        p0 = [F(1,2),F(1,4),F(1,4)]
        probs = [np.prod([p0[j] for j in w]) for w in words]
    else:
        raw = np.exp(-beta*np.array(energies)/100)
        p0 = raw/raw.sum()
        probs = [float(np.prod([p0[j] for j in w])) for w in words]
    p = np.array(probs,dtype=float)
    output = np.zeros(ns)
    aux = np.zeros((na,na))
    for s in range(ns):
        z = V[:,s].reshape(ns,na)
        output += p[s]*np.sum(z*z,axis=1)
        aux += p[s]*(z.T@z)
    desired = np.zeros(ns)
    for s in range(ns): desired[target[s]] = p[s]
    assert np.max(np.abs(output-desired)) < 1e-12
    initial = np.zeros(na)
    for k in range(1,length+1): initial[k*nk] = 1/sqrt(length)
    initialrho = np.outer(initial,initial)
    eig = np.linalg.eigvalsh(aux)
    assert np.count_nonzero(eig>1e-12) <= 3
    se = entropy(aux)
    da = float(np.abs(np.linalg.eigvalsh(aux-initialrho)).sum()/2)
    minority = sum(p[s] for s in range(ns) if category[s] == -1)
    e = 2/length + (minority if with_record else 0)
    assert float(initial@aux@initial) >= 1-e-1e-12
    sb = log(3) if e >= 2/3 else -e*log(e)-(1-e)*log(1-e)+e*log(2)
    assert se <= sb+1e-12 and da <= sqrt(e)+1e-12
    ea = np.array([k/100 for k in range(nb) for bit in range(nk)])
    deltaa = float(np.diag(aux)@ea-initialrho.diagonal()@ea)
    deltas = float((output-p)@np.array(energy)/100)
    assert abs(deltaa+deltas)<1e-12 and abs(deltaa)<=.01+1e-12
    drift = float(np.abs(output-p).sum()/2)
    assert abs(float(-np.sum(output*np.log(output)))+float(np.sum(p*np.log(p)))) < 1e-12
    if with_record:
        read = sum(aux[j,j] for j in range(na) if j%nk==1)
        assert abs(read-minority)<1e-12
    if beta is None:
        a,b = F(1,8),F(1,16)
        exact = (1-(1-a-b)**blocks)*(a-b)/(a+b)
        assert abs(drift-float(exact)) < 1e-12
    else:
        exact = None
        rel = float(np.sum(output*np.log(output/p)))
        assert abs(rel+beta*deltaa) < 1e-11
        assert drift <= sqrt(max(0,(abs(beta)*abs(deltaa)+se)/2))+1e-10
    return dict(blocks=blocks,L=length,record_included=with_record,beta=beta,
                complete_basis_size=ns*na,exact_energy_and_bijection=True,
                exact_material_distance=str(exact) if exact is not None else None,
                material_distance=round(drift,12),auxiliary_entropy=round(se,12),
                auxiliary_trace_distance=round(da,12),auxiliary_energy_change=round(deltaa,12),
                auxiliary_rank=int(np.count_nonzero(eig>1e-12)),
                unknown_input_isometry_squared_norm_exact=str(F(2,length)))


def native_certificate():
    d=json.loads((STAGE/'1202/results.json').read_text(encoding='utf-8'))
    r=list(map(F,d['one_marginal_rational']))
    x=[0,120,0,80]; y=[112,0,88,0]
    a=np.prod([z**k for z,k in zip(r,x)])
    b=np.prod([z**k for z,k in zip(r,y)])
    assert isinstance(a,F) and a>b>0
    q=a+b; M=(q.denominator+q.numerator-1)//q.numerator
    assert (M-1)*q<1<=M*q
    L=10**6; e=F(2,L)+b/q
    assert b/q<F('.000447')
    assert e<F('0.000449') and e<F('0.0212')**2
    old=json.loads((STAGE/'1208/results.json').read_text(encoding='utf-8'))['native_parity_witness']
    elo,ehi=[4*F(v) for v in old['energy_interval']]
    assert 0<elo<ehi<F('.00339')
    # e^-1 < 60/163 follows from exp(1)>163/60 (strict beyond 5th term).
    assert exp_bounds(F(1))[0]>F(163,60)
    dl=F(103,163)*(a-b)/q
    assert dl>F('.6313')
    assert exp_bounds(F('7.71'))[0]>1/F('.000449')
    assert exp_bounds(F('.694'))[0]>2
    sup=F('.000449')*(F('7.71')+1+F('.694'))
    assert sup<F('.00423')
    meanhi=ehi*F(L+1,2); varhi=ehi**2*F(L*L-1,12)
    assert meanhi<1695 and varhi<960000
    return dict(reference='1202 cold one-material marginal, independently prepared copies',
                x_occupation=x,y_occupation=y,block_length=200,blocks=str(M),
                total_materials=str(200*M),ladder_L=L,ladder_dimension=L+2,record_dimension=2,
                energy_gap_interval=[str(elo),str(ehi)],
                energy_change_absolute_upper='.00339',auxiliary_infidelity_upper='.000449',
                auxiliary_trace_distance_upper='.0212',auxiliary_entropy_upper='.00423',
                entropy_upper_exact=str(sup),material_distance_strict_lower='.6313',
                initial_ladder_mean_energy_upper='1695',initial_ladder_variance_upper='960000',
                first_x_output_probability_upper='.000447',
                original_QQ_realization_certified=False,physical_preparation_certified=False,
                all_matrices_materialized=False,finite_integer_certificate=True)


def evaluate():
    cases=[small_case(m,L,rec) for m in [1,2] for L in [4,16] for rec in [False,True]]
    gibbs=[small_case(1,8,True,beta) for beta in [-1.,0.,1.]]
    r=np.array([.4,.2,.4]);tau=np.ones(3)/3
    entropy_delta=float(-np.sum(r*np.log(r))-log(3))
    assert entropy_delta<0 and abs((r-tau)@np.arange(3))<1e-15
    assert abs(np.abs(r-tau).sum()/2-2/15)<1e-15
    sources=['1202/results.json','1207/proof.md','1207/results.json','1208/results.json','1208/acceptance.json']
    return dict(schema='internal_budget_reference_1209_v1',status='AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE',
                source_sha256={s:sha(STAGE/s) for s in sources},
                exact_small_models=cases,gibbs_budget_diagnostics=gibbs,
                one_sided_entropy_deletion=dict(energy_change=0,entropy_change=round(entropy_delta,12),
                    material_and_auxiliary_distance_exact='2/15',same_auxiliary_return=False),
                native_record_inclusive_certificate=native_certificate(),
                general_theorem_requires_analytic_proof=True,all_checks_passed=True,
                new_adopted_axioms=0,actual_native_operation_permission_certified=False)


if __name__=='__main__':
    result=evaluate(); path=HERE/'results.json'
    if '--write-results' in sys.argv:
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert result==json.loads(path.read_text(encoding='utf-8')),'Saved results mismatch'
    print(json.dumps(dict(all_checks_passed=True,results_match=True,
        exact_small_models=len(result['exact_small_models']),gibbs_cases=3,
        native_material_distance_lower='.6313',record_included=True),ensure_ascii=False))
