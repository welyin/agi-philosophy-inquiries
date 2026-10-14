"""1195 internal, non-blind independent review. No author module imports.

Read-only by default. --write-results uses exclusive creation of review_checks.json.
Grades are computed by Pascal DP and explicit subset counting. The independent
Haar polynomial (mu=0) is contracted with exact rational Grassmann signs.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
AUTHOR_FILES = ('baseline.json', 'proof.md', 'check.py', 'results.json',
                'admission.md', 'sources.md', 'NEXT.md', '../research_note_1195.md')
BASIS = (0, 3, 5, 10, 12, 15)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def mm(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))), F(0))
             for j in range(len(B[0]))] for i in range(len(A))]

def tr(A):
    return sum((A[i][i] for i in range(len(A))), F(0))

def transpose(A):
    return [list(row) for row in zip(*A)]

def gproduct(left, right, bits=8):
    if left & right:
        return 0, 0
    L = [i for i in range(bits) if left & (1 << i)]
    R = [i for i in range(bits) if right & (1 << i)]
    inversions = sum(x > y for x in L for y in R)
    return left | right, (-1 if inversions % 2 else 1)

def haar_coefficients(l1, l2):
    # Obtained independently by expanding original A.1/A.2, before review.
    # The two-site variables are p1,p2,b1,b2 at each endpoint; mu=0.
    q = 2*(l1*l1 + l2*l2)
    return {0: F(1), 0x0f: -l1, 0xf0: -l1, 0xff: q,
            0x3c: l2, 0x5a: l2, 0xa5: l2, 0xc3: l2}

def contraction(l1, l2):
    T = [[F(0) for _ in BASIS] for _ in BASIS]
    for col, state in enumerate(BASIS):
        for kernel_mask, value in haar_coefficients(l1, l2).items():
            mask, sign = gproduct(kernel_mask, state)
            if not sign or mask & 15 != 15:
                continue
            out = mask >> 4
            assert out in BASIS
            T[BASIS.index(out)][col] -= sign*value  # A.4 orientation
    return T

def template(l1, l2):
    T = [[F(0) for _ in BASIS] for _ in BASIS]
    T[0][0] = T[5][5] = l1
    T[0][5] = -F(1)
    T[5][0] = -2*(l1*l1+l2*l2)
    for j, sign in ((1,-1), (2,1), (3,1), (4,-1)):
        T[j][j] = sign*l2
    return T

def grading_check():
    counts = [1]
    for _ in range(16):
        updated = [0]*(len(counts)+1)
        for degree, n in enumerate(counts):
            updated[degree] += n
            updated[degree+1] += n
        counts = updated
    direct = [0]*17
    for mask in range(1 << 16):
        direct[mask.bit_count()] += 1
    assert counts == direct
    by_residue = [sum(counts[j] for j in range(r,17,12)) for r in range(12)]
    pairs = [(a,b) for a in range(0,17,2) for b in range(0,17,2)
             if (96+16-a+b) % 12 == 0]
    allowed = {r: set() for r in range(0,12,2)}
    for a,b in pairs:
        allowed[a % 12].add(b % 12)
    assert all(allowed[r] == {(r+8) % 12} for r in allowed)
    capacities = [min(by_residue[r], by_residue[next(iter(allowed[r]))])
                  for r in range(0,12,2)]
    even_dimension = sum(counts[::2])
    bound = sum(capacities)
    grade8_output = sum(counts[b] for a,b in pairs if a == 8)
    assert even_dimension == 32768 and bound == 13951
    assert even_dimension-bound == 18817
    assert counts[8]-grade8_output == 11049
    full_grade_only = sum(min(by_residue[r],by_residue[(r+8)%12]) for r in range(12))
    assert full_grade_only == 24991
    return {'method':'Pascal DP, independently checked by all 65536 subsets',
            'degree_counts':counts, 'residue_dimensions':by_residue,
            'even_block_rank_capacities':capacities,
            'even_dimension':even_dimension, 'rank_upper_even':bound,
            'nullity_lower_even':even_dimension-bound,
            'grade8_nullity_lower':counts[8]-grade8_output,
            'full_grade_only_loose_bound':full_grade_only,
            'scope':'single-cell unmodified amplitude; actual rank not computed; odd inputs additionally vanish'}

def kernel_and_block_check():
    samples = ((F(1),F(1,2)), (F(0),F(1)), (F(-2),F(1,3)), (F(3,5),F(-2,7)))
    for l1,l2 in samples:
        T = contraction(l1,l2)
        assert T == template(l1,l2)
        q = 2*(l1*l1+l2*l2)
        W = [[F(0) for _ in BASIS] for _ in BASIS]
        for i,w in enumerate((q,F(1),F(1),F(1),F(1),F(1))): W[i][i]=w
        assert mm(W,T) == mm(transpose(T),W)
        assert q > l1*l1 and l2 != 0
    T = contraction(F(1),F(1,2))
    T2 = mm(T,T); T3 = mm(T2,T); T4 = mm(T2,T2)
    traces = [tr(T2),tr(T3),tr(T4)]
    assert traces == [F(8),F(17),F(179,4)]
    Q = [[F(0) for _ in BASIS] for _ in BASIS]
    for i,x in enumerate(BASIS):
        for j,y in enumerate(BASIS):
            mask,sg = gproduct(x,y,4)
            if mask == 15: Q[i][j]=F(sg)
    C = [[F(0) for _ in BASIS] for _ in BASIS]
    for mask,v in haar_coefficients(F(1),F(1,2)).items():
        C[BASIS.index(mask&15)][BASIS.index(mask>>4)]=v
    CQ = mm(transpose(C),Q)
    assert T == [[-x for x in row] for row in CQ]
    arr = np.array(T,dtype=float)
    metric = np.array([2.5,1,1,1,1,1],float)
    A = np.sqrt(metric)[:,None]*arr/np.sqrt(metric)[None,:]
    assert np.max(abs(A-A.T)) < 1e-14
    values, V = np.linalg.eigh(A)
    exact_spectrum = np.sort([1-np.sqrt(2.5),1+np.sqrt(2.5),-.5,-.5,.5,.5])
    assert np.max(abs(values-exact_spectrum)) < 1e-14
    assert sum(values<0) == 3
    assert max(abs(values))-min(abs(values)) > 2
    signs = np.sign(values)
    energies = -0.5*np.log(values**2)
    H = (V*energies)@V.T
    P = (V*signs)@V.T
    positive_step = (V*np.exp(-energies))@V.T
    reconstructed = P@positive_step
    block_reconstruction = (V*np.exp(-2*energies))@V.T
    errors = {'signed_one_step':float(np.max(abs(reconstructed-A))),
              'positive_two_step':float(np.max(abs(block_reconstruction-A@A))),
              'involution':float(np.max(abs(P@P-np.eye(6)))),
              'commutator_PH':float(np.max(abs(P@H-H@P)))}
    assert max(errors.values()) < 1e-13
    assert np.linalg.eigvalsh(A@A).min() > .24999999999999
    qspec = np.linalg.eigvalsh(np.array(Q,dtype=float))
    assert np.max(abs(qspec-np.array([-1,-1,-1,1,1,1]))) < 1e-14
    # Removing the sign changes even the single-step kernel and odd trace.
    assert np.linalg.norm(positive_step-A,ord=2) > 1
    return {'haar_coefficients_at_fixed_point':{hex(k):str(v) for k,v in haar_coefficients(F(1),F(1,2)).items()},
            'exact_matrix':[[str(v) for v in row] for row in T],
            'exact_traces_2_3_4':[str(x) for x in traces],
            'eigenvalues':values.tolist(),
            'exact_positive_metric_diagonal':['5/2','1','1','1','1','1'],
            'berezin_pairing_inertia':[3,3,0],
            'two_step_eigenvalues':np.linalg.eigvalsh(A@A).tolist(),
            'reconstruction_errors':errors,
            'parameter_samples_checked':len(samples),
            'scope':'Berezin contraction is not an identified CP superoperator; new W is not original reflection positivity'}

def run():
    grade = grading_check()
    kernel = kernel_and_block_check()
    author = json.loads((HERE/'results.json').read_text(encoding='utf-8-sig'))
    assert author['wetterich_cell']['rank_upper_even'] == grade['rank_upper_even']
    assert author['wetterich_cell']['nullity_lower_even'] == grade['nullity_lower_even']
    assert author['wetterich_cell']['grade_dimensions'] == grade['residue_dimensions']
    assert author['dv_original_kernel']['exact_traces'] == kernel['exact_traces_2_3_4']
    assert np.array_equal(np.array(author['dv_original_kernel']['matrix']),np.array(contraction(F(1),F(1,2)),float))
    for name in AUTHOR_FILES:
        assert (HERE/name).is_file(), 'Missing final author file: '+name
    return {'status':'PASS_WITH_SCOPE', 'review_type':'internal non-blind independent mathematical/code review',
            'author_check_groups':3, 'independent_review_groups':3,
            'grading':grade, 'berezin_and_blocking':kernel,
            'author_files_sha256':{name:sha(HERE/name) for name in AUTHOR_FILES},
            'reviewer_files_sha256':{name:sha(HERE/name) for name in ('reviewer_check.py','review.md')},
            'not_certified':['full four-dimensional transfer or physical quotient',
                             'original reflection positivity or physical observable adjoint',
                             'original real-time arbitrary-reference quantum instrument',
                             'GR, dark matter, dark energy, or exclusion of all spinor candidates']}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    obj=run()
    result=HERE/'review_checks.json'
    if args.write_results:
        with result.open('x',encoding='utf-8',newline='\n') as handle:
            json.dump(obj,handle,ensure_ascii=False,indent=2)
            handle.write('\n')
    elif result.exists():
        saved=json.loads(result.read_text(encoding='utf-8'))
        assert obj==saved, 'Read-only rerun or frozen bindings differ from saved review'
    print(json.dumps(obj,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
