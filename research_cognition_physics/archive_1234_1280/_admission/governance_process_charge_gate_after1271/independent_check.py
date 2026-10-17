"""Independent finite diagnostics; no Type-II or physical implementation simulation.
Default: recalculate and compare frozen JSON. --write: exclusive first save.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
from math import factorial

HERE = Path(__file__).resolve().parent
ZERO = F(0)
def mat(rows): return [[F(x) for x in row] for row in rows]
def ident(n): return [[F(i == j) for j in range(n)] for i in range(n)]
def diag(values): return [[F(x) if i == j else ZERO for j in range(len(values))] for i, x in enumerate(values)]
def tr(a): return sum(a[i][i] for i in range(len(a)))
def adj(a): return list(map(list, zip(*a)))
def add(a,b): return [[x+y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]
def sub(a,b): return [[x-y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]
def mul(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def scale(c,a): return [[c*x for x in row] for row in a]
def sum_mats(xs): 
    xs=list(xs)
    z=diag([0]*len(xs[0]))
    for x in xs: z=add(z,x)
    return z
def channel(ks,rho): return sum_mats(mul(mul(k,rho),adj(k)) for k in ks)
def heis(ks,b): return sum_mats(mul(mul(adj(k),b),k) for k in ks)
def comm(a,b): return sub(mul(a,b),mul(b,a))
def is_zero(a): return all(x==0 for row in a for x in row)
def kron(a,b): return [[a[i][j]*b[k][l] for j in range(len(a[0])) for l in range(len(b[0]))] for i in range(len(a)) for k in range(len(b))]
def frac(x): return str(x)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def matrix_group():
    C=diag([0,0,1,1]); I=ident(4)
    K0=diag([1,F(3,5),1,F(3,5)])
    K1=mat([[0,F(4,5),0,0],[0,0,0,0],[0,0,0,F(4,5)],[0,0,0,0]])
    ks=[K0,K1]
    assert heis(ks,I)==I
    assert all(is_zero(comm(k,C)) for k in ks)
    for B in [C, sub(I,C), add(scale(F(7,3),sub(I,C)),scale(F(-5,7),C))]:
        assert heis(ks,B)==B
    defect=sub(channel(ks,I),I)
    assert defect==diag([F(16,25),F(-16,25),F(16,25),F(-16,25)])
    rho=diag([F(1,8),F(1,8),F(1,4),F(1,2)])
    branches=[mul(mul(k,rho),adj(k)) for k in ks]
    probs=[tr(r) for r in branches]
    cond=[tr(mul(C,r))/p for r,p in zip(branches,probs)]
    assert probs==[F(3,5),F(2,5)]
    assert cond==[F(43,60),F(4,5)]
    assert sum(p*c for p,c in zip(probs,cond))==tr(mul(C,rho))==F(3,4)
    # Distinct conditional feedback, retaining both records.
    swap=mat([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
    fks=[K0,mul(swap,K1)]
    assert heis(fks,I)==I and heis(fks,C)==C
    # Exact operator identity, hence arbitrary states on system plus reference.
    R=mat([[1,F(2,3)],[F(2,3),-1]])
    extended=[kron(k,ident(2)) for k in fks]
    assert heis(extended,kron(C,R))==kron(C,R)
    # q rises by one using another same-region energy carrier.
    Q=diag([0,0,1,1]); A=diag([0,1,0,1]); total=add(Q,A)
    U=mat([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    r0=diag([0,1,0,0]); r1=channel([U],r0)
    assert is_zero(comm(U,total))
    dq=tr(mul(Q,sub(r1,r0))); da=tr(mul(A,sub(r1,r0)))
    assert (dq,da)==(1,-1)
    # Same restricted A=B(C^2) x I state, different commutant charge extension.
    q0=kron(diag([F(1,2),F(1,2)]),diag([1,0]))
    q1=kron(diag([F(1,2),F(1,2)]),diag([0,1]))
    for i in range(2):
        for j in range(2):
            e=diag([0,0]); e[i][j]=F(1)
            assert tr(mul(kron(e,ident(2)),q0))==tr(mul(kron(e,ident(2)),q1))
    Ccomm=kron(ident(2),diag([0,1]))
    assert tr(mul(Ccomm,sub(q1,q0)))==1
    return {
      "status":"PASS_EXACT_FRACTION",
      "density_channel_unital":False,
      "unital_defect_diagonal":[frac(defect[i][i]) for i in range(4)],
      "branch_probabilities":list(map(frac,probs)),
      "branch_charge_means":list(map(frac,cond)),
      "unconditional_charge_mean":"3/4",
      "all_charge_spectral_functions_fixed":True,
      "feedback_and_arbitrary_reference_identity":True,
      "q_change":frac(dq),"aux_energy_change":frac(da),
      "same_visible_state_commutant_charge_difference":"1"
    }


def entropy_group():
    # Different elementary certificates from author's upper-tail expansion:
    # e^x < (1-x/n)^(-n), from -log(1-y)>y.
    exp_14_9_upper=(F(144,137))**32
    exp_1_upper=F(6,5)**6
    exp_1_lower=sum(F(1) / factorial(j) for j in range(4))
    exp_3_lower=sum(F(3)**j / factorial(j) for j in range(9))
    exp_07_lower=sum(F(7,10)**j / factorial(j) for j in range(5))
    assert exp_14_9_upper < 5
    assert exp_1_upper < 3 and exp_1_lower == F(8,3)
    assert exp_3_lower > 20 and exp_3_lower**3 > 3001
    assert exp_07_lower > 2
    # arctan(u)<u, pi>3 imply xstar<3/2+1/18.
    xstar_upper=F(3,2)+F(1,18)
    assert xstar_upper==F(14,9)
    support_upper=F(3,8)**2*(1-F(3,8))
    probability_upper=10*support_upper
    gap_lower=1-probability_upper
    assert support_upper==F(45,512)
    assert gap_lower==F(31,256)
    # delta<=1/100, alpha=1/100. log3001<9; log(1+alpha)<alpha.
    eta_upper=F(1,100)*9+F(1,100)
    cost_lower=1-eta_upper
    assert eta_upper==F(1,10) and cost_lower==F(9,10)>F(7,10)
    return {
      "status":"PASS_EXACT_FRACTION",
      "xstar_strict_upper":frac(xstar_upper),
      "exp_14_9_strict_upper":frac(exp_14_9_upper),
      "exp_1_strict_upper":frac(exp_1_upper),
      "exp_1_strict_lower":frac(exp_1_lower),
      "exp_9_lower_via_cube":frac(exp_3_lower**3),
      "exp_07_strict_lower":frac(exp_07_lower),
      "support_trace_strict_upper":frac(support_upper),
      "forward_half_L1_strict_lower":frac(gap_lower),
      "entropy_eta_strict_upper_at_delta_001":frac(eta_upper),
      "aux_entropy_cost_strict_lower_at_delta_001":frac(cost_lower),
      "two_dimensional_auxiliary_entropy_upper":"7/10",
      "finite_corner_product_invariant":"(t_B*C_0)*(tau(e_1)/t_B)=C_0*tau(e_1)"
    }

def source_bindings():
    names=["README.md","proof.md","check.py","results.json","source_baseline.json"]
    hashes={name:sha(HERE/name) for name in names}
    baseline=json.loads((HERE/"source_baseline.json").read_text(encoding="utf-8"))
    archive=HERE.parent.parent
    for relative,expected in baseline["sources_sha256"].items():
        assert sha(archive/relative)==expected,relative
    return hashes

def calculate():
    out={"schema":"governance_process_charge_independent_v1",
         "status":"PASS_WITH_SCOPE","new_scientific_rounds":0,
         "matrix_diagnostic":matrix_group(),
         "entropy_bounds":entropy_group(),
         "author_files_sha256":source_bindings(),
         "independent_diagnostic_groups":2,
         "new_formal_check_groups":0,
         "scope":"Finite rational matrices check algebra only; no Type-II realization, instruments, thermodynamic heat cost or new geometry certificate.",
         "self_sha256":sha(Path(__file__))}
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    args=ap.parse_args()
    result=calculate()
    dest=HERE/"independent_results.json"
    if args.write:
        with dest.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
        print("WRITE_ONCE_PASS")
    else:
        saved=json.loads(dest.read_text(encoding="utf-8"))
        assert saved==result,"independent result mismatch"
        print("READ_ONLY_MATCH")
    print(json.dumps({"groups":2,"status":result["status"]},ensure_ascii=False))
if __name__=="__main__": main()
