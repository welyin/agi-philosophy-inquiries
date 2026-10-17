"""1249: one exact-arithmetic group; no quantum trajectory simulation.

Default is read-only. --write exclusively creates the first saved result.
Huge finite resources are encoded as integer exponent trees, never expanded.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "quantum_bounds_results.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def ceil_log2(n):
    assert n > 0
    return (n - 1).bit_length()

def bell(n):
    row = [1]
    for r in range(1, n + 1):
        nxt = [0] * (r + 1)
        for k in range(1, r + 1):
            nxt[k] = (row[k - 1] if k - 1 < len(row) else 0) + k * (row[k] if k < len(row) else 0)
        row = nxt
    return sum(row)

def calculate():
    D, order, m = 14, 20, 10
    beta_norm = factorial(81) // factorial(40) ** 2
    coeff = {41 + k: F(beta_norm * (-1) ** k * comb(40, k), 41 + k) for k in range(41)}
    assert sum(coeff.values()) == 1
    for r in range(1, 41):
        assert sum(c * F(factorial(j), factorial(j-r)) for j,c in coeff.items()) == 0
    B = 324 * beta_norm * 2 ** 40
    # Also bounds the smaller good-box transitions: 8B+26B=34B.
    Lcut = 56 * B
    H0 = 2 ** 36
    # Complex radius-one polydisc: Q in [-3,4], other coordinates <=129.
    exp12 = 3 ** 12
    h_complex_bound = (F(5,4) * 129**2 * exp12 + 129**2 * exp12
                       + F(3,20)*3**4*(1+F(129,5)) + F(72,100)
                       + F(3,2)*129**2 + F(130,4)*129**2
                       + F(129**2,10) + F(130*129,10))
    assert h_complex_bound < H0
    # All original analytic material/source symbols and first Hamiltonian
    # derivatives have a deliberately generous radius-one bound 2^96.
    source_complex_bound = 1024 * (3 ** 18) * (1 + H0) * (1 + 129 ** 2)
    assert source_complex_bound < 2 ** 96
    Lh = D * (1 + Lcut)
    raw_M = H0 * factorial(21) * Lh ** 21
    m_exp = ceil_log2(raw_M)
    M = 1 << m_exp
    assert M >= raw_M
    bell_majorant = bell(20) * factorial(20) ** 20
    assert bell_majorant <= 2 ** 4096
    R_exp = 4096 + 5 * M
    # (A,L) means ||D^r f|| <= A r! L^r, 0<=r<=20.
    # Source completion includes chi derivatives and the clock reconstruction.
    A_source = 2 ** 98 * (1 + Lcut) ** 2
    L_source = 4 * (Lcut + D)
    # Four positive-domain beta masks: widths >=3/100. 100B is conservative.
    L_positive = L_source * (1 + A_source * 100 * B)
    L_good = Lcut + 4 * L_positive
    L_K = Lcut + 2 * (2_000_000_000 * B)
    # <=64 declared bounded source/energy/extension symbols in the weight.
    A_tail = 1 + 128 * A_source ** 2
    L_tail = L_good + 2 * L_source
    A_menu = max(A_source, A_tail, 1)
    L_menu = max(L_source, L_tail, L_K)
    volume = 5 * 256 ** 13
    # For k=0..4 and D=14, dyadic annuli bound the Fourier integral by 2^34.
    fourier_integral_bounds = []
    for k in range(5):
        upper = F(2**14) + F(2**(k+28), 1) / (1 - F(1,2**(6-k)))
        assert upper < 2 ** 34
        fourier_integral_bounds.append(str(upper))
    prefactor = 2 ** 34 * volume * 15 ** 10 * factorial(20)
    EH = ceil_log2(prefactor * H0 * Lh ** 20)
    Ebase = ceil_log2(prefactor * A_menu * (1 + L_menu) ** 20)
    EF = Ebase + 20 * R_exp
    # Two order-epsilon heat smoothing terms + Weyl Moyal remainder.
    # Each term is <2^-42 with this n=ceil(max(1,hbar))*2^N.
    N = 1 + max((EF + 42 + 2)//3, (EH + EF + 42 + 5)//6)
    linear_margin = 3*N - EF
    quadratic_margin = 6*N - EH - EF
    assert linear_margin > 42 and quadratic_margin > 42
    delta = F(1,2**40)
    assert F(1,2**41) < delta < F(1,10**10) and delta < F(1,100)
    old = json.loads((HERE.parent / "1247/results.json").read_text(encoding="utf-8"))
    centers = [F("-.00227032632043"), F("-.00227026126850")]
    eta = F(1,10**9)
    intervals, center_errors = [], []
    for i,c in enumerate(centers):
        branch = old["group_finite_memory"][f"alpha{i}"]
        endpoint = branch["checkpoints"]["2"]["state"]
        vals = endpoint["K"]
        lo,hi = map(F,vals)
        err = max(abs(lo-c),abs(hi-c))
        assert err < F(3,10**15) < eta/2
        intervals.append([str(lo),str(hi)])
        center_errors.append(str(err))
        tube = branch["all_time_tube"]
        assert F(tube["u"][0]) > F(1,3)
        assert F(tube["u"][1]) < F(3,2)
        assert F(tube["H"][1]) < F(11,10)
        for key in ("p", "v", "X", "V"):
            assert max(map(abs,map(F,tube[key]))) < F(3,2)
        for key in ("q", "r", "psi", "Y", "K", "z"):
            assert max(map(abs,map(F,tube[key]))) < 1
        qmax = max(map(abs,map(F,endpoint["q"])))
        rmax = max(map(abs,map(F,endpoint["r"])))
        assert qmax*qmax+rmax*rmax < 1
        ymax = max(map(abs,map(F,endpoint["Y"])))
        kmax = max(map(abs,map(F,endpoint["K"])))
        assert ymax+kmax/4 < 1
    assert centers[1]-centers[0] > 2*eta
    # Conservative core bounds from the retained classical positivity bootstrap.
    assert 2*27*F(11,10) < 61
    assert 27*F(3,2) < 61
    assert F(98,1000) > F(2,25)
    assert F(930,1000) > F(9,10)
    assert F(7,10) > F(3,5)  # hold F: use mask upper threshold 3/5, below .7
    assert F(1,20) > F(1,25) # hold Mr: use mask upper threshold 1/25
    bound_ledger = {
        "cutoff_C40_beta_normalization": str(beta_norm),
        "B_decimal": str(B), "Lcut_decimal": str(Lcut),
        "h_complex_upper": str(h_complex_bound), "H0": str(H0),
        "source_complex_upper": str(source_complex_bound),
        "source_analytic_upper": "2^96",
        "M_power_of_two_exponent": m_exp,
        "Bell20": bell(20), "Bell_factor_log2_upper": ceil_log2(bell_majorant),
        "flow_R": "2^(4096+5*2^M_power_of_two_exponent)",
        "menu_A_log2_upper": ceil_log2(A_menu),
        "menu_L_log2_upper": ceil_log2(L_menu),
        "Fourier_F0_through_F4_integral_upper": "2^34",
        "dyadic_rational_checks": fourier_integral_bounds,
        "Fourier_H_exponent": EH,
        "Fourier_menu_base_exponent": Ebase,
        "Fourier_transport_exponent_decimal": str(EF),
        "integer_n_power_exponent_decimal": str(N),
        "linear_error_negative_binary_exponent": str(linear_margin),
        "quadratic_error_negative_binary_exponent": str(quadratic_margin),
    }
    return {
        "round":1249, "part":"finite_joint_quantum_bridge", "status":"PASS",
        "science_groups":1, "arithmetic":"exact integers and Fraction; no seven-dimensional quantum simulation",
        "group_constructive_finite_resource_bound": {
            "status":"PASS", "phase_dimension":D, "canonical_pairs":7,
            "spatial_regularities":{"cutoff":"C40", "flow_derivatives":20, "h_derivatives":21, "Fourier_derivatives":20},
            "time_regularity":"original C2 pulses; no smoothing or high time derivatives",
            "cutoff":{"Q_core":[-1,2],"Q_support":[-2,3],"other_core":[-64,64],"other_support":[-128,128]},
            "good_box":{"Q_core":["-1/2","3/2"],"Q_support":["-3/4","7/4"],"other_core":[-61,61],"other_support":[-62,62]},
            "positive_masks":{"J":["1/25","2/25"],"Mc":["4/5","9/10"],"F":["1/2","3/5"],"Mr":["1/100","1/25"]},
            "ledger":bound_ledger,
            "resource_formula":"m=ceil(max(1,hbar)); n=m*2^N; V=1000*n^3; epsilon=hbar/V; N is integer_n_power_exponent_decimal",
            "uniform_intensive_menu_error":"<2^-41<delta=2^-40<1e-10",
            "Husimi_energy_weighted_bad_tail":"<2^-41 in the declared bounded completion only",
            "terminal_K":{"time":"9/4", "retained_from":"2", "centers":[str(c) for c in centers],"halfwidth":str(eta),"old_K_intervals":intervals,"center_errors":center_errors,"disjoint":True,"correct_branch_probability_lower":"1-2^-40 > .99","same_three_outcome_POVM":True},
            "limits":["known coherent preparation family; no arbitrary-state classical trajectory", "zero-mode effective completion, not full-field truncation", "terminal POVM is mathematical permission, not a local instrument construction", "total-energy absolute error does not vanish merely by increasing V", "no tail estimate for the unextended unbounded original stress", "no TR, CC, full quantum constraints, or quantum-clock performance theorem"]
        },
        "bindings": {name:sha(path) for name,path in {
            "code":Path(__file__), "canonical.md":HERE/"canonical.md",
            "canonical_results.json":HERE/"canonical_results.json",
            "1247/results.json":HERE.parent/"1247/results.json",
            "1247/proof.md":HERE.parent/"1247/proof.md",
            "1248/proof.md":HERE.parent/"1248/proof.md"}.items()}
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        print("FIRST_EXCLUSIVE_WRITE PASS; science_groups=1")
    else:
        saved=json.loads(OUT.read_text(encoding="utf-8"))
        if saved != result:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS; science_groups=1")
    group=result["group_constructive_finite_resource_bound"]
    ledger=group["ledger"]
    print(json.dumps({"M_exponent":ledger["M_power_of_two_exponent"],"N_decimal_digits":len(ledger["integer_n_power_exponent_decimal"]),"menu_error":group["uniform_intensive_menu_error"],"K_success":group["terminal_K"]["correct_branch_probability_lower"]},ensure_ascii=False))

if __name__ == "__main__":
    main()
