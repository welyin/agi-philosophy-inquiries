"""Finite local-CPTP normalization tests for bipartite qubit processes.

Convention: legs AI, AO, BI, BO; unnormalized process W=4*w; Tr(w)=1.
Default is read-only. This does not simulate an autonomous physical apparatus.
"""
from fractions import Fraction
from functools import reduce
from itertools import product
from pathlib import Path
import argparse
import json
import math
import numpy as np

PAULI = [np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex), np.diag([1, -1]).astype(complex)]
R_TYPES = [(a, 0, b, 0) for a, b in product(range(4), repeat=2)]
P_TYPES = [(a, o, b, 0) for a, o, b in product(range(4), range(1, 4), range(1, 4))]
Q_TYPES = [(a, 0, b, o) for a, b, o in product(range(1, 4), range(4), range(1, 4))]


def kron_all(matrices):
    return reduce(np.kron, matrices)


def pauli_word(word):
    return kron_all([PAULI[k] for k in word])


def pairing(x, y):
    z = np.einsum('ij,ji->', x, y)
    assert abs(z.imag) < 2e-11
    return float(z.real)


def reorder(matrix, order):
    n = len(order)
    return matrix.reshape([2]*(2*n)).transpose(list(order)+[n+i for i in order]).reshape(2**n, 2**n)


def twirl(matrix, axes):
    """Replace selected qubits by normalized identity, in original leg order."""
    n = int(round(math.log2(matrix.shape[0])))
    kept = [i for i in range(n) if i not in axes]
    order = kept + list(axes)
    d, e = 2**len(kept), 2**len(axes)
    reordered = reorder(matrix, order).reshape(d, e, d, e)
    partial = np.einsum('iaja->ij', reordered)
    return reorder(np.kron(partial, np.eye(e)/e), np.argsort(order))


def trace_distance(x, y):
    return float(np.sum(np.abs(np.linalg.eigvalsh(x-y)))/2)


def channel(u, v, sign):
    """Two-qubit input -> two-qubit output, with all measurement branches."""
    branches = []
    for outcome in (-1, 1):
        effect = (np.eye(4) + outcome*u.T)/2
        state = (np.eye(4) + sign*outcome*v)/4
        assert np.linalg.eigvalsh(effect).min() > -1e-12
        assert np.linalg.eigvalsh(state).min() > -1e-12
        assert abs(np.trace(state)-1) < 1e-12
        branches.append(np.kron(effect.T, state))
    choi = sum(branches)
    assert np.max(np.abs(choi-(np.eye(16)+sign*np.kron(u, v))/4)) < 1e-12
    assert np.max(np.abs(np.einsum('iaja->ij', choi.reshape(4,4,4,4))-np.eye(4))) < 1e-12
    return choi, branches


def random_channel(rng):
    iso = rng.normal(size=(4,2)) + 1j*rng.normal(size=(4,2))
    iso, _ = np.linalg.qr(iso)
    # Two Kraus operators, output row then input column; CJ input first.
    result = np.zeros((4,4), complex)
    for k in (iso[:2], iso[2:]):
        vec = k.T.reshape(-1)
        result += np.outer(vec, vec.conj())
    assert np.max(np.abs(np.einsum('iaja->ij',result.reshape(2,2,2,2))-np.eye(2))) < 1e-12
    return result


def calculate():
    assert len(R_TYPES)==16 and len(P_TYPES)==len(Q_TYPES)==36
    words = R_TYPES + P_TYPES + Q_TYPES
    assert len(set(words))==88
    basis = np.array([pauli_word(w) for w in words])
    gram = np.einsum('aij,bji->ab', basis, basis)
    assert np.max(np.abs(gram-16*np.eye(88))) < 1e-12
    rng = np.random.default_rng(1057)
    max_norm_error = max_branch_sum = max_projection_error = max_hs_error = 0.0
    max_single_error = max_twirl_error = 0.0
    min_branch = 1.0
    diagnostic_rows = []
    direct_settings = 0
    for case in range(12):
        weights = rng.normal(size=87)
        perturbation = np.einsum('a,aij->ij', weights, basis[1:])
        perturbation *= 0.75 / np.max(np.abs(np.linalg.eigvalsh(perturbation)))
        w = (np.eye(16)+perturbation)/16
        assert np.linalg.eigvalsh(w).min() > 0
        coefficients = np.einsum('ij,aji->a', w, basis).real
        cp, cq = coefficients[16:52], coefficients[52:]
        x = w-twirl(w,(1,))
        y = w-twirl(w,(3,))
        cmax, dmax = float(np.max(np.abs(cp))), float(np.max(np.abs(cq)))
        menu_defect = float(np.max(np.abs(np.outer(cp,cq))))
        assert abs(menu_defect-cmax*dmax) < 1e-14
        ordered_p = twirl(w,(3,))
        ordered_q = twirl(w,(1,))
        assert min(np.linalg.eigvalsh(ordered_p).min(),np.linalg.eigvalsh(ordered_q).min()) > -1e-12
        distance = min(trace_distance(w,ordered_p),trace_distance(w,ordered_q))
        assert distance <= 3*min(cmax,dmax)+1e-12
        assert distance <= 3*math.sqrt(menu_defect)+1e-12
        for _ in range(8):
            ma, mb = random_channel(rng), random_channel(rng)
            test = np.kron(ma,mb)
            max_single_error = max(max_single_error, abs(pairing(4*w,test)-1))
            max_twirl_error = max(max_twirl_error,abs(pairing(4*ordered_p,test)-1),abs(pairing(4*ordered_q,test)-1))
        ww = np.kron(w,w)
        forbidden = ww-twirl(ww,(1,5))-twirl(ww,(3,7))+twirl(ww,(1,3,5,7))
        explicit = np.kron(x,y)+np.kron(y,x)
        max_projection_error = max(max_projection_error,float(np.linalg.norm(forbidden-explicit)))
        max_hs_error = max(max_hs_error,abs(float(np.linalg.norm(forbidden))-math.sqrt(2)*np.linalg.norm(x)*np.linalg.norm(y)))
        grouped = 16*reorder(ww,[0,4,1,5,2,6,3,7])
        # Dense Born contractions are independent of the coefficient formula.
        pairs = [(int(np.argmax(abs(cp))),int(np.argmax(abs(cq))))]
        pairs += [(int(rng.integers(36)),int(rng.integers(36))) for _ in range(5)]
        for i,j in pairs:
            p,q = P_TYPES[i],Q_TYPES[j]
            ua,va = np.kron(PAULI[p[0]],PAULI[q[0]]),np.kron(PAULI[p[1]],PAULI[q[1]])
            ub,vb = np.kron(PAULI[p[2]],PAULI[q[2]]),np.kron(PAULI[p[3]],PAULI[q[3]])
            for s,t in product((-1,1),repeat=2):
                ma,ba = channel(ua,va,s)
                mb,bb = channel(ub,vb,t)
                z = pairing(grouped,np.kron(ma,mb))
                expected = 1+s*t*cp[i]*cq[j]
                branch = [pairing(grouped,np.kron(a,b)) for a,b in product(ba,bb)]
                min_branch = min(min_branch,min(branch))
                max_norm_error = max(max_norm_error,abs(z-expected))
                max_branch_sum = max(max_branch_sum,abs(sum(branch)-z))
                direct_settings += 1
        diagnostic_rows.append({'case':case,'finite_menu_defect':menu_defect,'chosen_twirl_trace_distance':distance,
                                'certified_formula_bound':3*math.sqrt(menu_defect)})
    sharp=[]
    p=pauli_word((0,3,3,0)); q=pauli_word((3,0,0,3))
    for t in (Fraction(1,16),Fraction(1,8),Fraction(1,4),Fraction(1,2)):
        w=(np.eye(16)+float(t)*(p+q))/16
        assert np.linalg.eigvalsh(w).min() > -1e-12
        assert abs(pairing(w,p)-float(t)) < 1e-12 and abs(pairing(w,q)-float(t)) < 1e-12
        dist=min(trace_distance(w,twirl(w,(1,))),trace_distance(w,twirl(w,(3,))))
        assert abs(dist-float(t/2)) < 1e-12
        rows=[]
        for k in (1,2,4,8,16):
            defect=(1+t*t)**k-1
            rows.append({'blocks':k,'independent_copies':2*k,'exact_normalization_defect':str(defect)})
        sharp.append({'t':str(t),'exact_nearest_order_trace_distance':str(t/2),
                      'exact_two_copy_menu_defect':str(t*t),'amplification':rows})
    k=1
    while (Fraction(257,256))**k-1 <= Fraction(1,10):
        k+=1
    assert k==25
    assert Fraction(257,256)**24-1 <= Fraction(1,10) < Fraction(257,256)**25-1
    for value in (max_norm_error,max_branch_sum,max_projection_error,max_hs_error,max_single_error,max_twirl_error):
        assert value < 2e-11
    assert min_branch >= -2e-11
    return {'round':1057,'all_checks_passed':True,
            'scope':'finite four-qubit-leg process; declared independent copies and cross-copy laboratory permissions',
            'exact_counts':{'allowed_types_including_identity':88,'each_signalling_family':36,
                            'fixed_signed_two_copy_settings':5184},
            'normalization_tolerance_is_not_a_probability_metric':True,
            'state_trace_distance_is_not_a_process_strategy_norm':True,
            'matrix_diagnostics':{'processes':12,'direct_two_copy_settings':direct_settings,
                'max_dense_born_formula_error':max_norm_error,'max_instrument_branch_sum_error':max_branch_sum,
                'min_raw_branch_weight':min_branch,'max_forbidden_projection_error':max_projection_error,
                'max_hilbert_schmidt_identity_error':max_hs_error,'max_single_copy_normalization_error':max_single_error,
                'max_ordered_twirl_normalization_error':max_twirl_error},
            'random_process_diagnostics':diagnostic_rows,'sharp_family':sharp,
            'finite_amplification_certificate':{'t':'1/16','total_error_threshold':'1/10','first_blocks':25,
                                                'independent_copies':50,
                                                'previous_defect':str(Fraction(257,256)**24-1),
                                                'first_exceeding_defect':str(Fraction(257,256)**25-1)},
            'new_cognitive_axioms':0,'physical_instrument_or_source_certified':False,
            'experimental_data':False,'independent_review_included':False}


def compare(saved, actual):
    if isinstance(actual,float):
        assert math.isclose(saved,actual,rel_tol=2e-10,abs_tol=2e-12),(saved,actual)
    elif isinstance(actual,dict):
        assert saved.keys()==actual.keys()
        for k in actual: compare(saved[k],actual[k])
    elif isinstance(actual,list):
        assert len(saved)==len(actual)
        for a,b in zip(saved,actual): compare(a,b)
    else:
        assert saved==actual,(saved,actual)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=Path(__file__).with_name('results.json')
    if args.save:
        with path.open('x',encoding='utf-8') as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write('\n')
    else:
        compare(json.loads(path.read_text(encoding='utf-8')),result)
    print(json.dumps({'round':1057,'all_checks_passed':True,'mode':'save_initial' if args.save else 'read_only',
                      'exact_counts':result['exact_counts'],'matrix_diagnostics':result['matrix_diagnostics']},ensure_ascii=False))


if __name__=='__main__':
    main()
