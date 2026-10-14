"""Independent review of 1192. Does not import or execute author check.py.
Exact inverse formula, generic Christoffel/Ricci tensor computation, beta moments,
and rational finite budget. --write creates review_checks.json once; default verifies.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ZERO = F(0)
HALF = F(1, 2)

def enc(x):
    return {"exact": str(x), "float": float(x)}

def basis_vector(n, j):
    return [F(i == j) for i in range(n)]

def ten_velocities():
    vectors = [[F(1), ZERO, ZERO, ZERO]]
    for axis in (1, 2, 3):
        for sign in (1, -1):
            v = [F(5, 4), ZERO, ZERO, ZERO]
            v[axis] = sign * F(3, 4)
            vectors.append(v)
    for i, j in ((1, 2), (1, 3), (2, 3)):
        v = [F(3, 2), ZERO, ZERO, ZERO]
        v[i], v[j] = F(1), HALF
        vectors.append(v)
    return vectors

def explicit_reconstruction(y):
    """Construct the symmetric matrix directly from ten quadratic contractions."""
    c = [[ZERO for j in range(4)] for i in range(4)]
    c[0][0] = y[0]
    for i in (1, 2, 3):
        yp, ym = y[2*i-1], y[2*i]
        c[0][i] = c[i][0] = F(4, 15) * (yp-ym)
        c[i][i] = F(8, 9)*(yp+ym)-F(25, 9)*y[0]
    for slot, (i, j) in enumerate(((1, 2), (1, 3), (2, 3)), start=7):
        c[i][j] = c[j][i] = (y[slot]-F(9, 4)*c[0][0]-3*c[0][i]
                            -F(3, 2)*c[0][j]-c[i][i]-F(1, 4)*c[j][j])
    return c

def contraction(c, u):
    return sum((c[i][j]*u[i]*u[j] for i in range(4) for j in range(4)), ZERO)

def review_design():
    us = ten_velocities()
    assert all(u[0]**2-sum((z*z for z in u[1:]), ZERO) == 1 for u in us)
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    columns = [explicit_reconstruction(basis_vector(10, k)) for k in range(10)]
    for k, c in enumerate(columns):
        assert [contraction(c, u) for u in us] == basis_vector(10, k)
    row_norms = [sum((abs(c[i][j]) for c in columns), ZERO) for i, j in pairs]
    assert max(row_norms) == F(24, 5)
    full = [[F((i+2)*(j+3), 31) for j in range(4)] for i in range(4)]
    full = [[full[min(i,j)][max(i,j)] for j in range(4)] for i in range(4)]
    assert explicit_reconstruction([contraction(full, u) for u in us]) == full
    missing = columns[-1]
    assert [contraction(missing, u) for u in us] == [ZERO]*9+[F(1)]
    assert missing[2][3] == 1 and missing[3][2] == 1
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))
    assert list(map(str, row_norms)) == author["clock_design"]["inverse_row_norms"]
    return {"inverse_infinity_norm": "24/5", "exact_inverse_rows_verified": 10,
            "row_norms": list(map(str, row_norms)),
            "nine_direction_null_verified": True, "signature": "+---"}

def zeros(*shape):
    if not shape:
        return ZERO
    return [zeros(*shape[1:]) for _ in range(shape[0])]

def ricci_from_metric_jet(x, eps):
    """Generic index contractions; no closed Ricci formula enters the computation."""
    a, ap, app = 1+eps*x*x, 2*eps*x, 2*eps
    g = zeros(4,4)
    gi = zeros(4,4)
    g[0][0], gi[0][0] = a*a, 1/(a*a)
    for i in (1,2,3):
        g[i][i] = gi[i][i] = F(-1)
    dg, ddg = zeros(4,4,4), zeros(4,4,4,4)
    dg[1][0][0] = 2*a*ap
    ddg[1][1][0][0] = 2*(ap*ap+a*app)
    dgi = zeros(4,4,4)
    for r in range(4):
        for i in range(4):
            for j in range(4):
                dgi[r][i][j] = -sum((gi[i][m]*dg[r][m][n]*gi[n][j]
                                     for m in range(4) for n in range(4)), ZERO)
    gamma, dgamma = zeros(4,4,4), zeros(4,4,4,4)
    for up in range(4):
        for m in range(4):
            for n in range(4):
                for j in range(4):
                    bracket = dg[m][j][n]+dg[n][j][m]-dg[j][m][n]
                    gamma[up][m][n] += HALF*gi[up][j]*bracket
                    for r in range(4):
                        dbracket = ddg[r][m][j][n]+ddg[r][n][j][m]-ddg[r][j][m][n]
                        dgamma[r][up][m][n] += HALF*(dgi[r][up][j]*bracket+gi[up][j]*dbracket)
    ricci = zeros(4,4)
    for m in range(4):
        for n in range(4):
            ricci[m][n] = sum((dgamma[r][r][m][n]-dgamma[n][r][m][r]
                              for r in range(4)), ZERO)
            ricci[m][n] += sum((gamma[r][r][s]*gamma[s][m][n]
                               -gamma[r][n][s]*gamma[s][m][r]
                               for r in range(4) for s in range(4)), ZERO)
    K = [[F(m == 0 and n == 0)-g[0][0]*gi[m][n]/4 for n in range(4)] for m in range(4)]
    assert sum((g[m][n]*K[m][n] for m in range(4) for n in range(4)), ZERO) == 0
    expected = zeros(4,4)
    expected[0][0], expected[1][1] = a*app, -app/a
    assert ricci == expected
    density = a*sum((K[m][n]*ricci[m][n] for m in range(4) for n in range(4)), ZERO)
    assert density == a*a*app/2
    B, L = 2*eps*(1+eps), 2*eps+6*eps*eps
    assert max(abs(gamma[r][m][n]) for r in range(4) for m in range(4) for n in range(4)) <= B
    assert max(abs(dgamma[s][r][m][n]) for s in range(4) for r in range(4)
               for m in range(4) for n in range(4)) <= L
    assert abs(g[0][0]-1) < HALF
    return density

def beta_moment(power_half):
    """Integral x^(2m)(1-x^2)^4 dx over [-1,1], using B(m+1/2,5)."""
    a = F(2*power_half+1, 2)
    den = F(1)
    for j in range(5):
        den *= a+j
    return F(24)/den

def review_curvature():
    eps = F(1,1000)
    vals = []
    for x in (F(-1), F(-1,2), ZERO, F(1,2), F(1)):
        vals.append({"x": str(x), "density": str(ricci_from_metric_jet(x, eps))})
    q0, q2, q4 = [beta_moment(m) for m in (0,1,2)]
    assert q0 == F(256,315)
    direct = (eps*q0+2*eps*eps*q2+eps**3*q4)*q0**2
    # Independent integration by parts: integrate q' x^r by -r * moment(q).
    derivative_term = eps*q0+6*eps*eps*q2+5*eps**3*q4
    quadratic_term = -4*eps*eps*q2-4*eps**3*q4
    by_parts = (derivative_term+quadratic_term)*q0**2
    assert direct == by_parts > 0
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))
    assert str(direct) == author["curved_diagnostic"]["exact_Ricci_pairing"]["exact"]
    d, w, dw = 4, 6, 72
    c1, c2 = 2*d*d*w, 2*d**3*dw+2*d**4*w
    assert (c1,c2) == (192,12288)
    # All benchmark trajectories have conserved a(x)^2 v^0.
    energy_upper = F(3,2)*(1+eps*F(1,16))
    assert energy_upper < 2
    assert energy_upper**2-1 < F(5,4)**2
    assert F(1,4)+F(1,5)*F(5,4) < 1
    return {"metric_jet_points": vals, "generic_Ricci_components_exact": True,
            "beta_moments": list(map(str,(q0,q2,q4))),
            "weak_pairing": enc(direct), "integration_by_parts_pairing": enc(by_parts),
            "IBP_constants": [c1,c2],
            "benchmark_velocity_buffer_via_Killing_energy": True,
            "scope": "specified smooth metric diagnostic, not a coupled physical-history construction"}

def review_budget():
    B0, L, V, tau = F(1,100), F(1,100), F(2), F(1,1000)
    dx, du, hb, h = F(1,10**12), F(1,10**12), F(1,10000), F(1,20)
    connection_change = 16*L*V**3
    velocity_change = 32*B0*V*(16*B0*V**2)
    J = connection_change+velocity_change
    b = 2*(dx+tau*du)/tau**2+32*B0*V*du+J*tau/3
    B = min(B0, F(24,5)*(b+16*V**2*L*hb)+L*h)
    Q = beta_moment(0)
    I0, I1 = Q**3, 8*Q**2
    weak = lambda z: 192*z*I1+12288*z*z*I0
    source_lower = F(3,8)*25*F(1,2)
    es, eE, kappa = F(1,8), F(1,8), F(1)
    margin = kappa*(source_lower-es)-weak(B)-eE
    threshold = (weak(B)+eE)/(source_lower-es)
    assert margin > F(3,4)
    assert kappa*(source_lower-es)-weak(B0)-eE < 0
    assert F(1,10)*(source_lower-es) < weak(B)+eE
    values = {"J":J,"b":b,"Gamma_updated":B,"I0":I0,"I1":I1,
              "weak_Ricci_updated_bound":weak(B),"weak_Ricci_prior_bound":weak(B0),
              "source_reference_lower":source_lower,"strict_margin":margin,
              "kappa_exclusion_threshold":threshold}
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))
    for key,value in values.items():
        assert str(value) == author["finite_certificate"][key]["exact"], key
    return {key:enc(value) for key,value in values.items()}

def review_fluid():
    """Direct full tensor contraction, including a nontrivial lapse and boost."""
    for lapse in (F(1), F(6,5)):
        g = [[ZERO for _ in range(4)] for _ in range(4)]
        gi = [[ZERO for _ in range(4)] for _ in range(4)]
        g[0][0], gi[0][0] = lapse*lapse, 1/(lapse*lapse)
        for j in (1,2,3):
            g[j][j] = gi[j][j] = F(-1)
        K = [[F(m == 0 and n == 0)-g[0][0]*gi[m][n]/4
              for n in range(4)] for m in range(4)]
        for ut, ux in ((F(1),ZERO),(F(5,4),F(3,4))):
            u = [ut/lapse, ux, ZERO, ZERO]
            lowered = [sum((g[m][n]*u[n] for n in range(4)),ZERO) for m in range(4)]
            assert sum((lowered[m]*u[m] for m in range(4)),ZERO) == 1
            kuu = sum((K[m][n]*lowered[m]*lowered[n] for m in range(4) for n in range(4)),ZERO)
            assert kuu == lapse*lapse*(ut*ut-F(1,4))
            for rho, pressure in ((F(1),ZERO),(F(1),F(-1)),(F(1),F(-2,3))):
                T = [[(rho+pressure)*lowered[m]*lowered[n]-pressure*g[m][n]
                      for n in range(4)] for m in range(4)]
                measured = sum((K[m][n]*T[m][n] for m in range(4) for n in range(4)),ZERO)
                assert measured == (rho+pressure)*kuu
                if ux == 0:
                    assert measured == F(3,4)*g[0][0]*(rho+pressure)
    cases = {"dust":str(F(3,4)), "vacuum":str(ZERO), "negative_pressure":str(F(1,4))}
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))
    assert cases == author["fluid_tracefree_pairs"]
    return {"rest_eta_pairs":cases, "boosted_tensor_contraction_verified":True,
            "scope":"ideal-fluid algebra only; no dark-sector identity or acceleration diagnosis"}

def bindings():
    required = [HERE/name for name in
                ("admission.md","proof.md","sources.md","check.py",
                 "results.json","baseline.json","NEXT.md")]
    required.append(HERE.parent/"research_note_1192.md")
    result = {}
    for p in required:
        assert p.is_file(), str(p)
        key = str(p.relative_to(HERE.parent)).replace("\\","/")
        result[key] = hashlib.sha256(p.read_bytes()).hexdigest()
    return result

def run():
    return {"status":"PASS_WITH_SCOPE",
            "independent_method":"explicit inverse; generic metric-jet Ricci; beta integrals; exact rational budget",
            "design":review_design(), "curvature":review_curvature(), "finite_budget":review_budget(),
            "fluid":review_fluid(), "author_files_sha256":bindings(),
            "scope":{"actual_clock_accuracy_certified":False,
                     "same_history_resources_realized":False,
                     "smooth_complete_source_matching_certified":False,
                     "all_permitted_effective_dictionaries_excluded":False,
                     "gr_goal_complete":False,
                     "mathematical_conditional_implication_verified":True}}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    data = run()
    target = HERE/"review_checks.json"
    if args.write:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        assert json.loads(target.read_text(encoding="utf-8-sig")) == data
    print(json.dumps({"status":data["status"],"strict_margin":data["finite_budget"]["strict_margin"],
                      "bound_author_files":len(data["author_files_sha256"])},ensure_ascii=False))
