"""1292 exact finite checks. Default/--check are read-only; --write creates once."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "results.json"
KS = (-3, -1, 0, 1, 2)


def frac(x):
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def encode(x):
    if isinstance(x, F):
        return frac(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def mod1(x):
    return F(x) % 1


def add(x, y):
    return (mod1(x[0] + y[0]), (x[1] + y[1]) % 4)


def multiple(n, x):
    return (mod1(n * x[0]), (n * x[1]) % 4)


def order(x, maximum):
    for n in range(1, maximum + 1):
        if multiple(n, x) == (F(0), 0):
            return n
    raise AssertionError("order bound insufficient")


def cover_fields(k):
    # (name, SU(3) triality, weak-diagonal weight, integer q, r/2pi phase)
    # W=(-1)^a diag(exp(3i theta),exp(-3i theta)); weak central sign
    # is included separately in cover_phase for the doublet entries.
    return [
        ("Q_up", 1, 1, 1 + 6*k, F(1,4)),
        ("Q_down", 1, -1, 1 + 6*k, F(1,4)),
        ("uc", -1, 0, -4 - 6*k, F(1,4)),
        ("dc", -1, 0, 2 - 6*k, F(1,4)),
        ("L_up", 0, 1, -3 - 18*k, F(1,4)),
        ("L_down", 0, -1, -3 - 18*k, F(1,4)),
        ("ec", 0, 0, 6 + 18*k, F(1,4)),
        ("N", 0, 0, 18*k, F(1,4)),
        ("H_up", 0, 1, 3, F(1,2)),
        ("H_down", 0, -1, 3, F(1,2)),
        ("S", 0, 0, -36*k, F(1,2)),
    ]


def cover_phase(field, color, theta, a):
    _, triality, weak, q, rphase = field
    return mod1(F(color*triality,3) + (q+3*weak)*theta
                + a*rphase + (F(a,2) if weak else 0))


def group_checks(k):
    fields = cover_fields(k)
    expected_kernel = {(n % 3, mod1(F(n,6)), 0) for n in range(6)}
    # For k=0 a finite root grid checks the fixed central kernel, not U(1).
    scale = max(1, abs(k))
    # Original covering theta has period 2pi; this grid contains all
    # possible centers in the H,S stabilizer for nonzero k.
    grid = 72*scale
    stabilizer = []
    kernel = set()
    Hdown, S = fields[9], fields[10]
    for c in range(3):
        for j in range(grid):
            theta = F(j, grid)
            for a in range(4):
                if cover_phase(Hdown,c,theta,a) == 0 and cover_phase(S,c,theta,a) == 0:
                    stabilizer.append((c,theta,a))
                    if all(cover_phase(f,c,theta,a) == 0 for f in fields):
                        kernel.add((c,theta,a))
    assert kernel == expected_kernel
    # Four Yukawa charge equations plus SNN, checked using the original table.
    qQ, qu, qd, qL, qe, qN, qH, qS = (
        1+6*k,-4-6*k,2-6*k,-3-18*k,6+18*k,18*k,3,-36*k)
    assert qQ+qH+qu == qQ-qH+qd == qL-qH+qe == qL+qH+qN == qS+2*qN == 0
    assert mod1(F(1,2)+2*F(1,4)) == 0
    result = {"k":k, "original_cover_kernel_order":len(kernel),
              "Yukawa_and_SNN_charges_zero":True}
    if k == 0:
        assert all(a % 2 == 0 for _,_,a in stabilizer)
        for t in (F(0),F(1,7),F(2,5),F(11,13)):
            for a in (0,2):
                assert mod1(-18*k*t+F(a,2)) == 0
        result.update({"continuous_U1_remains":True, "K0":"U(1) x Z2^F",
                       "complete_stabilizer_connected":True,
                       "N_pure_internal_character":"1"})
        return result
    n = abs(k)
    # K after the theta~theta+pi kernel is removed: t=phi/(2pi).
    K = {(F(j,36*n),a) for j in range(36*n) for a in range(4)
         if mod1(-18*k*F(j,36*n)+F(a,2)) == 0}
    x, z = (mod1(F(1,36*k)),1), (F(0),2)
    represented = {add(multiple(j,x),multiple(e,z))
                   for j in range(36*n) for e in range(2)}
    assert represented == K and len(K) == 72*n
    assert order(x,36*n) == 36*n and order(z,2) == 2
    assert z not in {multiple(j,x) for j in range(36*n)}
    color = multiple(12*k,x)
    assert color == (F(1,3),0)
    Hc = {multiple(j,color) for j in range(3)}
    Hcz = {add(c,multiple(e,z)) for c in Hc for e in range(2)}
    assert len(Hc) == 3 and len(Hcz) == 6
    def quotient_count(subgroup):
        remaining = set(K)
        count = 0
        while remaining:
            g = next(iter(remaining))
            orbit = {add(g,h) for h in subgroup}
            assert orbit <= K
            remaining.difference_update(orbit)
            count += 1
        return count
    before_spin, after_spin = quotient_count(Hc), quotient_count(Hcz)
    assert before_spin == 24*n and after_spin == 12*n
    pure = {F(j,18*n) for j in range(18*n)}
    assert all(mod1(-18*k*t) == 0 for t in pure)
    # Internal components additionally quotient the color-center subgroup.
    assert len(pure)//3 == 6*n
    N_character = mod1(9*k*F(1,18*n))
    assert N_character == F(1,2)
    result.update({
        "K_order":len(K), "x_order":36*n, "z_order":2,
        "internal_component_order":6*n,
        "with_r_component_order_before_Spin":before_spin,
        "complete_component_order_after_Spin":after_spin,
        "N_pure_internal_character":"-1",
        "negative_k_supported":k<0,
        "delete_color_quotient_wrong_factor":3,
        "delete_Spin_identification_wrong_factor":2,
    })
    return result


def dot(v,w):
    return sum((a*b for a,b in zip(v,w)),F(0))


def quadratic_parent(k,g,gp,v,w,A):
    qh, qs = (-g/2,gp/2),(F(0),-6*k*gp)
    # Direct |D H|^2+|D S|^2 at the background; E=1/2 A^T M^2 A.
    return v*v*dot(qh,A)**2/2 + w*w*dot(qs,A)**2/2


def matrix_from_parent(k,g,gp,v,w):
    zero=(F(0),F(0)); basis=((F(1),F(0)),(F(0),F(1)))
    e=lambda a:quadratic_parent(k,g,gp,v,w,a)
    a=e(basis[0])+e((-F(1),F(0)))-2*e(zero)
    d=e(basis[1])+e((F(0),-F(1)))-2*e(zero)
    b=e((F(1),F(1)))-e(basis[0])-e(basis[1])+e(zero)
    return [[a,b],[b,d]]


def det(M):
    return M[0][0]*M[1][1]-M[0][1]*M[1][0]


def tr(M):
    return M[0][0]+M[1][1]


def matmul(A,B):
    return [[sum((A[i][h]*B[h][j] for h in range(2)),F(0))
             for j in range(2)] for i in range(2)]


def trans(A):
    return [[A[j][i] for j in range(2)] for i in range(2)]


def inverse(M):
    d=det(M)
    assert d
    return [[M[1][1]/d,-M[0][1]/d],[-M[1][0]/d,M[0][0]/d]]


def shifted(M,s):
    return [[M[0][0]+s,M[0][1]],[M[1][0],M[1][1]+s]]


def response(M,s):
    return s*inverse(shifted(M,s))[1][1]


def matrix_checks():
    rows=[]
    parameter_sets=[(F(1),F(1),F(1),F(1)),
                    (F(2),F(3,2),F(4,3),F(2,5))]
    for g,gp,v,w in parameter_sets:
        qh,qs0=(-g/2,gp/2),(F(0),F(0))
        a=g*g/(g*g+gp*gp)
        mZ=v*v*(g*g+gp*gp)/4
        Mzero=matrix_from_parent(0,g,gp,v,w)
        delta_star=det(matrix_from_parent(1,g,gp,v,w))/tr(matrix_from_parent(1,g,gp,v,w))
        for k in KS:
            M=matrix_from_parent(k,g,gp,v,w)
            C=[[v*v*qh[0],v*v*qh[1]],[F(0),-6*k*gp*w*w]]
            assert M[0][0]==v*v*g*g/4 and M[0][1]==-v*v*g*gp/4
            assert M[1][1]==v*v*gp*gp/4+36*k*k*gp*gp*w*w
            assert det(M)==9*k*k*g*g*gp*gp*v*v*w*w
            assert det(C)==3*k*g*gp*v*v*w*w
            rank=1 if not k else 2
            assert (det(C)==0)==(k==0)
            for s in (a*delta_star/8,a*delta_star/6,a*delta_star/4):
                G=inverse(shifted(M,s))
                assert matmul(shifted(M,s),G)==[[F(1),F(0)],[F(0),F(1)]]
                current=matmul(matmul(C,G),trans(C))
                assert current[0][0]>0 and det(current)>=0
                if k:
                    assert det(current)>0
                    assert det(M)/tr(M)>=delta_star
                    assert response(M,s)<=s/(s+delta_star)<=a/4
                    assert response(Mzero,s)-response(M,s)>=3*a/4
                else:
                    assert response(M,s)==a+(1-a)*s/(s+mZ)
            rows.append({"parameters":[g,gp,v,w],"k":k,"mass_matrix":M,
                         "current_matrix":C,"mass_determinant":det(M),
                         "current_rank":rank})
        assert matrix_from_parent(-1,g,gp,v,w)==matrix_from_parent(1,g,gp,v,w)
    return rows


def polynomial_mul(p,q):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            e=tuple(i+j for i,j in zip(a,b))
            out[e]=out.get(e,F(0))+x*y
    return {e:c for e,c in out.items() if c}


def derivative(p,axis):
    out={}
    for e,c in p.items():
        if e[axis]:
            f=list(e); f[axis]-=1
            out[tuple(f)]=c*e[axis]
    return out


def evaluate(p,x):
    return sum((c*prod(t**n for t,n in zip(x,e)) for e,c in p.items()),F(0))


def prod(xs):
    out=F(1)
    for x in xs: out*=x
    return out


def potential_checks():
    # Lower Higgs radial/phase, complex S radial/phase. The other Higgs
    # angular directions are gauge-related and are not counted here.
    v,w,lH,lS=F(1),F(1),F(1),F(1)
    fH={(1,0,0,0):2*v,(2,0,0,0):F(1),(0,2,0,0):F(1)}
    fS={(0,0,1,0):2*w,(0,0,2,0):F(1),(0,0,0,2):F(1)}
    V={}
    for polynomial,weight in [(fH,lH/4),(fS,lS/4)]:
        for e,c in polynomial_mul(polynomial,polynomial).items():
            V[e]=V.get(e,F(0))+weight*c
    origin=(F(0),)*4
    assert evaluate(V,origin)==0
    assert all(evaluate(derivative(V,i),origin)==0 for i in range(4))
    hessian=[[evaluate(derivative(derivative(V,i),j),origin) for j in range(4)] for i in range(4)]
    assert hessian==[[F(2),F(0),F(0),F(0)],[F(0)]*4,
                     [F(0),F(0),F(2),F(0)],[F(0)]*4]
    # Positivity is analytically a sum of squares; these are implementation checks.
    for x in [(F(1,2),F(1,3),F(-1,4),F(2,3)),(F(-2),F(0),F(1),F(0))]:
        assert evaluate(V,x)>=0
    return {"Hessian_radial_phase":hessian,"bounded_below_sum_of_squares":True,
            "k0_complex_S_phase_not_deleted":True}


def source_checks():
    # A homogeneous timelike covariant-derivative test of the scalar source.
    # N and a are varied independently BEFORE any equations of motion.
    N,a=F(3,2),F(5,4)
    Vpot=F(2,5)
    A=(F(1,3),F(1,2))
    rows=[]
    for k in KS:
        W=quadratic_parent(k,F(1),F(1),F(1),F(1),A)
        # reduced L=N a^3 [W/N^2 - V], a Laurent polynomial.
        L={(-1,3):W,(1,3):-Vpot}
        rho=-evaluate(derivative(L,0),(N,a))/(a**3)
        pressure=evaluate(derivative(L,1),(N,a))/(3*N*a*a)
        assert rho==W/(N*N)+Vpot
        assert pressure==W/(N*N)-Vpot
        Wdeleted=quadratic_parent(0,F(1),F(1),F(1),F(1),A)
        omitted=(W-Wdeleted)/(N*N)
        assert omitted==18*k*k*A[1]*A[1]/(N*N)
        assert (omitted>0)==(k!=0)
        # The same source quadratic form gives the gauge-current variation.
        M=matrix_from_parent(k,F(1),F(1),F(1),F(1))
        q=lambda b:quadratic_parent(k,F(1),F(1),F(1),F(1),(A[0],b))
        derivative_B=(q(A[1]+1)-q(A[1]-1))/2
        assert derivative_B==M[1][0]*A[0]+M[1][1]*A[1]
        rows.append({"k":k,"rho_scalar":rho,"P_scalar":pressure,
                     "omitted_S_source_if_dropped":omitted,
                     "gauge_current_derivative_B":derivative_B})
    return {"test":"off-shell homogeneous scalar source variation only",
            "rows":rows,"full_quantum_source_budget_certified":False}


def finite_window():
    one=F(1)
    M0=matrix_from_parent(0,one,one,one,one)
    M1=matrix_from_parent(1,one,one,one,one)
    delta=det(M1)/tr(M1)
    a=F(1,2); lo=F(9,584); hi=F(9,292)
    assert delta==F(18,73) and hi==a*delta/4
    assert response(M0,hi)==F(82,155)
    eps0=epsk=F(1,64)
    ideal=3*a/4
    surviving=ideal-eps0-epsk
    assert surviving==F(11,32)
    # All k bound is analytic. Samples verify the response calculated
    # by inversion rather than substituting the bound itself.
    rows=[]
    for k in KS:
        M=matrix_from_parent(k,one,one,one,one)
        for s in (lo,(lo+hi)/2,hi):
            R=response(M,s)
            rows.append({"k":k,"s":s,"R_B":R,
                         "R0_minus_Rk":response(M0,s)-R})
            if k: assert response(M0,s)-R>=ideal
    # Delete w>0: the two responses coincide at w=0. No common window
    # is claimed as the symmetry-breaking scale is sent to zero.
    Mdeleted=matrix_from_parent(1,one,one,one,F(0))
    assert Mdeleted==M0
    # Analytic monotonicity of A x/(B+C x): cross multiply for sample
    # positive x; the proof derivative is A B/(B+C x)^2 > 0.
    A=F(9); B=F(1,2); C=F(36)
    for x,y in ((F(1),F(4)),(F(1),F(9)),(F(4),F(9))):
        dx=A*x/(B+C*x);dy=A*y/(B+C*y)
        assert dy-dx==A*B*(y-x)/((B+C*x)*(B+C*y))>0
    return {"g":1,"gprime":1,"v":1,"w":1,"lambda_H":1,"lambda_S":1,
            "s_window":[lo,hi],"delta_star":delta,"a":a,
            "uniform_ideal_gap_all_nonzero_integer_k":ideal,
            "per_branch_assumed_total_response_error":[eps0,epsk],
            "conditional_surviving_gap":surviving,
            "root_upper_endpoint_R0":response(M0,hi),
            "root_upper_endpoint_R1":response(M1,hi),
            "samples":rows,
            "delete_S_condensate_removes_gap":True,
            "Euclidean_window_is_not_compact_spacetime_instrument":True}


def source_bindings():
    paths=["proof.md","check.py","../_shared/round1292_admission.md",
           "../_admission/neutrino_phase_selection_screen.md",
           "../../archive_1009_1043/1035/drafts/joint_selection_derivation.md",
           "../../archive_1046_1062/1060/proof.md",
           "../../archive_1288_/1289/proof.md"]
    return {p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}


def calculate():
    return encode({
        "schema":"round1292-snn-residual-response-v1",
        "round":1292,"author_check_groups":1,"verdict":"PASS_AUTHOR_SCOPED",
        "source_bindings":source_bindings(),
        "cover_and_residual_checks":[group_checks(k) for k in KS],
        "kinetic_variation_and_current_checks":matrix_checks(),
        "potential_and_retained_phase":potential_checks(),
        "source_variation_and_deletion":source_checks(),
        "finite_invariant_response":finite_window(),
        "scope":{
            "all_integer_k_analytic_proof":True,
            "negative_k_finite_checks":True,
            "same_complex_S_field_for_k0":True,
            "full_old_fermion_anomaly_inherited_not_reproved":True,
            "gauge_background_is_not_Elitzur_observable":True,
            "leading_translation_invariant_vacuum_Euclidean_response":True,
            "finite_response_errors_are_assumed_contract":True,
            "actual_instrument_or_unknown_quantum_code_certified":False,
            "exact_all_order_quantum_poles_certified":False,
            "all_bundles_have_defect_free_VEV":False,
            "k_to_minus_k_full_model_equivalence_claimed":False,
            "all_cognitive_axioms_select_k0":False,
            "new_axioms":0,
        },
    })


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=calculate()
    if args.write:
        raw=(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+"\n").encode("utf-8")
        with OUT.open("xb") as stream: stream.write(raw)
        receipt="CREATED_EXCLUSIVELY"
    else:
        saved=json.loads(OUT.read_text(encoding="utf-8"))
        if saved!=data:
            raise AssertionError("saved result differs; no file was written")
        receipt="READ_ONLY_MATCH"
    print(json.dumps({"verdict":data["verdict"],"receipt":receipt,
                      "uniform_ideal_gap":"3/8","conditional_gap":"11/32",
                      "results_sha256":hashlib.sha256(OUT.read_bytes()).hexdigest()},
                     ensure_ascii=False))

if __name__=="__main__":
    main()
