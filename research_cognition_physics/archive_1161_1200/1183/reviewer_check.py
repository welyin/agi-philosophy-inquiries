"""Independent exact review of 1183; never imports the author implementation.

Four checks review the same four scientific groups; no added scientific count.
Fraction polynomials in z=sin(x) replace the author's nonlinear differences.
--math-only does not bind files. --write creates the result once; default
recomputes every exact assertion and compares without writing.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


class P:
    def __init__(self, coeff=0):
        a = list(coeff) if isinstance(coeff, (tuple, list)) else [coeff]
        self.a = [F(x) for x in a]
        while len(self.a) > 1 and self.a[-1] == 0:
            self.a.pop()

    def __add__(self, other):
        b = other if isinstance(other, P) else P(other)
        return P([(self.a[i] if i < len(self.a) else 0)
                  + (b.a[i] if i < len(b.a) else 0)
                  for i in range(max(len(self.a), len(b.a)))])

    __radd__ = __add__

    def __neg__(self):
        return P([-a for a in self.a])

    def __sub__(self, other):
        return self + -(other if isinstance(other, P) else P(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        b = other if isinstance(other, P) else P(other)
        out = [F(0)] * (len(self.a) + len(b.a) - 1)
        for i, a in enumerate(self.a):
            for j, c in enumerate(b.a):
                out[i+j] += a*c
        return P(out)

    __rmul__ = __mul__

    def __truediv__(self, x):
        return self * (1/F(x))

    def __pow__(self, n):
        out = P(1)
        for _ in range(n):
            out = out*self
        return out

    def __eq__(self, other):
        return self.a == (other if isinstance(other, P) else P(other)).a

    def dz(self):
        return P([i*self.a[i] for i in range(1, len(self.a))] or [0])

    def dxx(self):
        return (1-Z*Z)*self.dz().dz() - Z*self.dz()

    def avg(self):
        # Exact normalized circle integrals of sin(x)^n.
        return sum((c*F(math.comb(n, n//2), 2**n)
                    for n, c in enumerate(self.a) if n % 2 == 0), F(0))

    def at(self, z):
        return sum((c*F(z)**i for i,c in enumerate(self.a)), F(0))

    def record(self):
        return [str(x) for x in self.a]


Z = P([0,1])


def calculate(lam, amplitude=1):
    pi = [P(2), 2+amplitude*Z, P(2)]
    trace = sum(pi)
    h = sum(x*x for x in pi) - lam*trace*trace
    p2 = -2*h
    assert h+p2/2 == 0

    # Independent derivatives of (sum q_i^2 pi_i^2-lam*(sum q_i*pi_i)^2)
    # /sqrt(prod q_i) and p^2/(2 sqrt(prod q_i)), evaluated at q_i=1.
    gq = [2*x*x-2*lam*trace*x-h/2 for x in pi]
    mq = [-p2/4]*3
    qdot = [2*(x-lam*trace) for x in pi]
    pidot = [-(a+b) for a,b in zip(gq,mq)]
    assert all(r == -2*x*x+2*lam*trace*x for r,x in zip(pidot,pi))

    # Direct first-order jet of the full nonlinear algebraic density along
    # q_i(t)=1+t*qdot_i, pi_i(t)=pi_i+t*pidot_i. This is not an evolution.
    trace_dot = sum(dq*x+r for dq,x,r in zip(qdot,pi,pidot))
    numerator_dot = sum(2*dq*x*x+2*x*r
                        for dq,x,r in zip(qdot,pi,pidot)) - 2*lam*trace*trace_dot
    gravity_dot = numerator_dot - h*sum(qdot)/2
    matter_dot = -p2*sum(qdot)/4
    assert gravity_dot+matter_dot == 0
    rdot = qdot[0].dxx() - sum(qdot).dxx()
    hdot = gravity_dot+matter_dot-rdot
    assert rdot == 2*(1-2*lam)*amplitude*Z
    assert hdot == -rdot

    # Full derivative of D_1: pi^{ij} is a contravariant tensor density
    # of weight one, so div(pi)^i=partial_j pi^{ij}+Gamma^i_jk pi^{jk}.
    # Divide all terms by cos(x); all remaining coefficients are polynomials.
    connection_contraction = (
        qdot[0].dz()*pi[0]-qdot[1].dz()*pi[1]-qdot[2].dz()*pi[2])/2
    d_g_dot_over_cos = -2*(pidot[0].dz()+connection_contraction)
    d_m_dot_over_cos = p2.dz()/2  # p*partial_x(phi_dot)=p*partial_x p
    assert d_g_dot_over_cos+d_m_dot_over_cos == 0

    # H[1], H[M] bracket directly from the two functional gradients.
    n, m = P(1), 2+Z
    curvature_q_n = [P(0)]*3
    curvature_q_m = [P(0), -Z, -Z]
    bracket_density = sum(
        (n*(g+mq_i)+cn)*(m*dq) -
        (n*dq)*(m*(g+mq_i)+cm)
        for g,mq_i,cn,cm,dq in zip(gq,mq,curvature_q_n,curvature_q_m,qdot))
    bracket = bracket_density.avg()
    assert bracket == amplitude*(1-2*lam)
    assert (m*hdot).avg() == -bracket
    cross_gm = -sum(n*dq*m*a for dq,a in zip(qdot,mq))
    cross_mg = sum(n*a*m*dq for dq,a in zip(qdot,mq))
    assert cross_gm+cross_mg == 0

    # Negative control: remove matter q-feedback, retaining its energy.
    bad_pidot = [-g for g in gq]
    bad_trace_dot = sum(dq*x+r for dq,x,r in zip(qdot,pi,bad_pidot))
    bad_num_dot = sum(2*dq*x*x+2*x*r
                      for dq,x,r in zip(qdot,pi,bad_pidot))-2*lam*trace*bad_trace_dot
    omitted_feedback_residual = bad_num_dot-(h+p2/2)*sum(qdot)/2
    assert omitted_feedback_residual == h*sum(qdot)/2
    assert omitted_feedback_residual != 0

    return {
        "lambda": str(lam), "amplitude": amplitude,
        "h_coefficients": h.record(), "p_squared_coefficients": p2.record(),
        "gravity_q_gradients": [p.record() for p in gq],
        "matter_q_gradients": [p.record() for p in mq],
        "qdot": [p.record() for p in qdot],
        "pidot": [p.record() for p in pidot],
        "gravity_algebraic_density_dot": gravity_dot.record(),
        "matter_density_dot": matter_dot.record(),
        "Rdot": rdot.record(), "Hdot": hdot.record(),
        "Dg_dot_over_cos": d_g_dot_over_cos.record(),
        "Dm_dot_over_cos": d_m_dot_over_cos.record(),
        "full_Ddot": "0",
        "cross_gm_plus_mg": "0",
        "bracket_H1_HM": str(bracket),
        "smeared_HM_dot": str((m*hdot).avg()),
        "omitted_matter_q_feedback_residual": omitted_feedback_residual.record(),
    }


def exact_math():
    cases = [calculate(lam) for lam in (F(2,5), F(1,2), F(3,5))]
    h_low_lambda = P([F(-12,5), F(-4,5), F(3,5)])
    h_high_lambda = P([F(-48,5), F(-16,5), F(2,5)])
    assert -1-h_low_lambda == (1+Z)*(7-3*Z)/5
    assert h_high_lambda+F(62,5) == (1-Z)*(14-2*Z)/5
    assert h_low_lambda.at(-1) == -1
    assert h_high_lambda.at(1) == F(-62,5)
    # h_lambda decreases with lambda since -(6+z)^2 is strictly negative.
    assert P([-36,-12,-1]) == -(6+Z)**2

    lapse_cases = []
    for a in (F(0), F(1,4), F(1,2)):
        n = 1+a*Z
        obstruction = (n*n*(1-Z*Z)).avg()
        div_n2_grad_pi = (n*n).dz()*(1-Z*Z)-n*n*Z
        ibp = ((6+Z)*div_n2_grad_pi).avg()
        assert ibp == -obstruction
        assert obstruction == F(1,2)+a*a/8 > 0
        assert obstruction >= (1-a)**2/2
        lapse_cases.append({"a":str(a), "integral":str(obstruction),
                            "integration_by_parts":str(ibp),
                            "minimum_lapse_lower_bound":str(1-a)})
    cmc = [calculate(lam,0) for lam in (F(2,5),F(1,2),F(3,5))]
    assert all(c["Rdot"] == ["0"] and c["Hdot"] == ["0"]
               and c["bracket_H1_HM"] == "0" for c in cmc)
    return {
        "group1_positive_source_initial_constraints": {
            "h_bounds":["-62/5","-1"], "rho_bounds":["1","62/5"],
            "p_squared_bounds":["2","124/5"],
            "upper_slack_factorization":"(1+z)(7-3z)/5",
            "lower_slack_factorization":"(1-z)(14-2z)/5",
            "initial_H":"0", "initial_D":["0","0","0"]},
        "group2_complete_tangent": cases,
        "group3_positive_lapse": {
            "arbitrary_smooth_positive_lapse_argument":
                "0=-integral N^2 |grad pi|^2 <= -n0^2/2 < 0",
            "exact_test_integrals":lapse_cases},
        "group4_cmc_deletion": cmc,
    }


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_author(math_result):
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8-sig"))
    assert author["round"] == 1183 and author["groups"] == 4 and author["status"] == "PASS"
    assert author["exact"]["gravity_density_bounds"] == ["-62/5","-1"]
    expected = {r["lambda"]:r["smeared_HM_dot"]
                for r in math_result["group2_complete_tangent"]}
    assert author["exact"]["smeared_unit_lapse_derivative"] == expected
    for row in author["nonlinear_tangent_checks"]:
        expected_value = 2*row["lam"]-1
        assert row["constraint_density_derivative_error"] < 2e-8
        assert abs(row["smeared_derivative"]-expected_value) < 2e-8
        assert row["momentum_derivative_error"] < 1e-10
        assert row["curvature_derivative_error"] < 1e-8
    for row, exact in zip(author["positive_lapse_checks"],
                          math_result["group3_positive_lapse"]["exact_test_integrals"]):
        assert abs(row["integrated_positive_obstruction"]-float(F(exact["integral"]))) < 1e-12
    for row in author["cmc_deletion_checks"]:
        assert row["first_tangent_residual"] < 2e-7
    for relative,digest in author["source_sha256"].items():
        assert sha(STAGE/relative) == digest, relative
    return {"author_groups":4, "source_sha256":author["source_sha256"],
            "comparison":"Exact independent derivatives agree with saved author diagnostics",
            "diagnostic_not_finite_time_certificate":True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--math-only", action="store_true")
    args = parser.parse_args()
    result = {"round":1183, "status":"PASS_WITH_SCOPE", "groups":4,
              "new_scientific_groups":0, "blockers":0,
              "method":"Independent exact Fraction polynomials, no author imports",
              "exact":exact_math()}
    if args.math_only:
        print("MATH_ONLY_PASS 4 groups; exact initial constraints, full tangent, lapse, CMC")
        return
    result["author_comparison"] = verify_author(result["exact"])
    result["author_core_sha256"] = {
        name:sha(HERE/name) for name in (
            "../research_note_1183.md","admission.md","proof.md",
            "check.py","results.json","NEXT.md")}
    result["reviewer_code_sha256"] = sha(Path(__file__))
    destination = HERE/"review_checks.json"
    if args.write:
        # Exclusive creation: never replace a previous scientific result.
        with destination.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
        print("WRITE_ONCE_PASS 4 groups; blockers 0")
    else:
        saved = json.loads(destination.read_text(encoding="utf-8-sig"))
        assert saved == result, "saved independent result or binding differs"
        print("READ_ONLY_MATCH 4 groups; blockers 0")


if __name__ == "__main__":
    main()
