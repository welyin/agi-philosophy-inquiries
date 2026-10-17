"""Independent finite-query reciprocity check (one mathematical group).

No author-code imports. Exact rational frame algebra and direct probability
queries supplement, but do not replace, the all-unitary proof. Default is
read-only comparison; --write-once creates the result file exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations, product
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 2e-11


def qpoly_add(a, b):
    c = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a): c[i] += x
    for i, x in enumerate(b): c[i] += x
    return c


def qpoly_mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] += x*y
    return c


def det_charpoly(matrix):
    n = len(matrix)
    ans = [F(0)] * (n+1)
    for p in permutations(range(n)):
        inv = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        term = [F((-1) ** inv)]
        for i, j in enumerate(p):
            term = qpoly_mul(term, [matrix[i][j], F(-1)] if i == j else [matrix[i][j]])
        ans = qpoly_add(ans, term)
    return ans


def pure_reduction(vector, dims, keep):
    rest = tuple(i for i in range(len(dims)) if i not in keep)
    tensor = vector.reshape(dims).transpose(tuple(keep) + rest)
    rows = math.prod(dims[i] for i in keep)
    mat = tensor.reshape(rows, -1)
    return mat @ mat.conj().T


def hsnorm(a):
    return float(np.linalg.norm(a))


def choi_test(U, g, m):
    d = g*m
    assert hsnorm(U.conj().T @ U - np.eye(d)) < TOL
    # Matrix entries U[(a,b),(c,d)] give normalized pure Choi components.
    v = U.reshape(-1) / math.sqrt(d)
    dims = (g, m, g, m)
    X = pure_reduction(v, dims, (0, 2))
    Y = pure_reduction(v, dims, (1, 3))
    rb = pure_reduction(v, dims, (0, 2, 3)) - np.kron(X, np.eye(m)/m)
    rf = pure_reduction(v, dims, (1, 3, 2)) - np.kron(Y, np.eye(g)/g)
    p = float(np.trace(X @ X).real)
    pb = float(np.trace(Y @ Y).real)
    rhs = 1-p
    out_g = pure_reduction(v, dims, (0,))
    out_m = pure_reduction(v, dims, (1,))
    errors = [
        abs(m*hsnorm(rb)**2-rhs), abs(g*hsnorm(rf)**2-rhs),
        abs(p-pb), hsnorm(out_g-np.eye(g)/g), hsnorm(out_m-np.eye(m)/m)
    ]
    assert max(errors) < TOL
    return rb, rf, p, max(errors)


def prob(U, Fout, rho, sigma):
    out = U @ np.kron(rho, sigma) @ U.conj().T
    z = np.trace(Fout @ out)
    assert abs(z.imag) < TOL
    return float(z.real)


def haar_diagnostic(rng, d):
    a = rng.normal(size=(d,d)) + 1j*rng.normal(size=(d,d))
    q, r = np.linalg.qr(a)
    phases = np.diag(r)
    return q @ np.diag(phases / np.abs(phases))


def calculate():
    I = np.eye(2, dtype=complex)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Y = np.array([[0,-1j],[1j,0]], dtype=complex)
    Z = np.diag([1,-1]).astype(complex)
    pauli = [X,Y,Z]
    effects = [(I+q)/2 for q in pauli]
    preps = [I/2] + effects
    differences = pauli
    basis = [I/math.sqrt(2)] + [q/math.sqrt(2) for q in pauli]
    coords = lambda a: np.array([np.trace(b @ a).real for b in basis])
    af = np.array([coords(e-np.trace(e)*I/2)[1:] for e in effects])
    ap = np.array([coords(r) for r in preps])
    ad = np.array([coords(s)[1:] for s in differences])
    sf = float(np.linalg.svd(af, compute_uv=False)[-1])
    sp = float(np.linalg.svd(ap, compute_uv=False)[-1])
    sd = float(np.linalg.svd(ad, compute_uv=False)[-1])
    gram = [
        [F(2),F(1,2),F(1,2),F(1,2)],
        [F(1,2),F(1,2),F(0),F(0)],
        [F(1,2),F(0),F(1,2),F(0)],
        [F(1,2),F(0),F(0),F(1,2)]
    ]
    direct = det_charpoly(gram)
    factored = qpoly_mul(qpoly_mul([F(-1,2),F(1)], [F(-1,2),F(1)]),
                         [F(1,4),F(-5,2),F(1)])
    assert direct == factored
    assert hsnorm(ap.T @ ap - np.array(gram, dtype=float)) < TOL
    sp2 = (5-math.sqrt(21))/4
    assert abs(sp*sp-sp2) < TOL
    assert abs(sf*sf-F(1,2)) < TOL and abs(sd*sd-2) < TOL

    U = np.zeros((4,4),dtype=complex)
    for a,b in product(range(2), repeat=2):
        U[2*a+(a^b), 2*a+b] = 1
    rb,rf,p,identity_error = choi_test(U,2,2)
    rho0, rho1 = effects[2], I-effects[2]
    sigma = effects[2]
    EM = effects[2]
    df = prob(U,np.kron(I,EM),rho0,sigma)-prob(U,np.kron(I,EM),rho1,sigma)
    Af = hsnorm(EM)*hsnorm(rho0-rho1)*hsnorm(sigma)
    assert abs(df-1) < TOL and abs(Af-math.sqrt(2)) < TOL

    reverse = []
    reconstruction_errors = []
    for f,r,s in product(effects,preps,differences):
        # s is the difference of the two actual normal qubit states (I +/- s)/2.
        delta = prob(U,np.kron(f,I),r,(I+s)/2)-prob(U,np.kron(f,I),r,(I-s)/2)
        rhs = 4*np.trace(np.kron(np.kron(f,r.T),s.T) @ rb)
        reconstruction_errors.append(abs(delta-rhs))
        reverse.append(delta)
    N = len(reverse)
    assert N == 36
    assert max(reconstruction_errors) < TOL
    reverse_max = max(map(abs,reverse))
    assert abs(reverse_max-1) < TOL
    # The forward contraction follows the (Mout,Min,Gin) order chosen above.
    forward_choi = 4*np.trace(np.kron(np.kron(EM,sigma.T),(rho0-rho1).T) @ rf)
    assert abs(forward_choi-df) < TOL
    assert np.linalg.norm(reverse)+TOL >= 4*sf*sp*sd*hsnorm(rb)
    c = sf*sp*sd/(math.sqrt(N)*Af)
    assert abs(c*c-(5-math.sqrt(21))/288) < TOL
    # Exact proof of c > 1/30: sqrt(21) < 117/25 = 5 - 288/900.
    assert F(21) < F(117,25)**2
    assert F(5)-F(288,900) == F(117,25)
    assert c > 1/30
    e = F(1,1000)
    robust_lower = (1-2*e)/30-2*e
    assert robust_lower == F(469,15000) > 0

    # Deleting access invalidates the positive frame premise, even for this U.
    diag_preps = [I/2,effects[2],I-effects[2]]
    deleted_input = [
        prob(U,np.kron(f,I),r,(I+s)/2)-prob(U,np.kron(f,I),r,(I-s)/2)
        for f,r,s in product(effects,diag_preps,differences)
    ]
    z_readout = [
        prob(U,np.kron(effects[2],I),r,(I+s)/2)-prob(U,np.kron(effects[2],I),r,(I-s)/2)
        for r,s in product(preps,differences)
    ]
    assert max(map(abs,deleted_input)) < TOL
    assert max(map(abs,z_readout)) < TOL
    assert np.linalg.matrix_rank(np.array([coords(r) for r in diag_preps]),tol=TOL) == 2

    rng = np.random.default_rng(12710037)
    diagnostics = []
    for g,m in [(2,3),(3,2),(2,4)]:
        for repetition in range(3):
            V = haar_diagnostic(rng,g*m)
            _,_,purity,err = choi_test(V,g,m)
            diagnostics.append({"g":g,"m":m,"sample":repetition,
                                "operator_schmidt_purity":purity,"identity_error":err})
    # Product unitaries provide the zero-signalling endpoint independently.
    product_U = np.kron(haar_diagnostic(rng,2),haar_diagnostic(rng,3))
    rb0,rf0,pp,_ = choi_test(product_U,2,3)
    assert hsnorm(rb0) < TOL and hsnorm(rf0) < TOL and abs(pp-1)<TOL

    return {
        "verdict":"PASS_FINITE_QUERY_RECIPROCITY_MATHEMATICAL_ONLY",
        "independent_groups":1,
        "independent_code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "normalized_choi_identity":"m*||Rb||_2^2 = g*||Rf||_2^2 = 1-Tr(X^2)",
        "frame_characteristic_polynomial_ascending":[str(x) for x in direct],
        "frame_squares":{"F":"1/2","P":"(5-sqrt(21))/4","D":"2"},
        "cnot":{"operator_schmidt_purity":p,"identity_error":identity_error,
                "forward_difference":df,"forward_norm_factor":Af,
                "cartesian_difference_queries":N,"underlying_probability_queries":2*N,
                "reverse_max_difference":reverse_max,
                "max_direct_vs_choi_error":float(max(reconstruction_errors)),
                "c_display":c,"c_squared_exact":"(5-sqrt(21))/288",
                "c_rational_lower":"1/30",
                "per_probability_error":"1/1000","robust_reverse_lower":str(robust_lower)},
        "deletion_checks":{"only_diagonal_G_preparation_reverse_max":float(max(map(abs,deleted_input))),
                           "only_Z_readout_reverse_max":float(max(map(abs,z_readout))),
                           "diagonal_preparation_linear_rank":2,
                           "required_preparation_rank":4},
        "unequal_dimensions_diagnostics_only":diagnostics,
        "scope":{
            "closed_full_unitary_required":True,
            "complete_actual_input_and_output_frames_certified":False,
            "fixed_environment_may_be_discarded":False,
            "choi_is_a_proof_device_not_a_free_preparation":True,
            "new_physical_geometry_port":False,
            "curvature_or_GR_claim":False,
            "round_1271_flat_model_excluded":False,
            "new_numbered_round":False
        }
    }


def same(a,b):
    if isinstance(a,float) and isinstance(b,(int,float)):
        return abs(a-b) <= TOL
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return a==b


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-once",action="store_true")
    args=parser.parse_args()
    result=calculate()
    out=HERE/"independent_results.json"
    if args.write_once:
        with out.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
        action="WRITE_ONCE_PASS"
    else:
        saved=json.loads(out.read_text(encoding="utf-8-sig"))
        assert same(result,saved),"saved independent result mismatch"
        action="READ_ONLY_MATCH"
    print(json.dumps({"verdict":result["verdict"],"action":action,
                      "independent_groups":1,
                      "robust_reverse_lower":result["cnot"]["robust_reverse_lower"]}))


if __name__=="__main__":
    main()