"""1278 one fixed diagnostic group. Default is read-only; --write creates once."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ORDER = 64
T = F(1, 8)
R0 = F(3, 20)
EPS = F(3, 20)
NU = F(7, 20)
V = 65**3
ELL = F(16, 65)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def convolution(a, b, n):
    return sum((a[j]*b[n-j] for j in range(n+1)), F(0))

def series(c):
    # ux, uy, uz, gx, gy, gz; r=(3/20)u.
    a = [[F(x)] for x in (0, 0, -1, 1, 0, 0)]
    for n in range(ORDER):
        ux, uy, uz, gx, gy, gz = a
        rhs = (
            -2*convolution(gz, uy, n),
            2*convolution(gz, ux, n)-2*c*uz[n],
            2*c*uy[n],
            -2*R0*convolution(uz, gy, n),
            2*R0*convolution(uz, gx, n)+2*gz[n],
            -2*gy[n],
        )
        for j in range(6):
            a[j].append(rhs[j]/(n+1))
    assert all(abs(row[n]) <= 4**n for row in a for n in range(ORDER+1))
    for n in range(ORDER+1):
        norm_u = sum((convolution(a[j], a[j], n) for j in range(3)), F())
        norm_g = sum((convolution(a[j], a[j], n) for j in range(3, 6)), F())
        assert norm_u == norm_g == (1 if n == 0 else 0)
        energy = c*R0*a[0][n]-a[3][n]+R0*convolution(a[2], a[5], n)
        assert energy == (-1 if n == 0 else 0)
    return a

def evaluate(a):
    return sum((v*T**n for n, v in enumerate(a)), F())

def integral(a):
    return sum((v*T**(n+1)/(n+1) for n, v in enumerate(a)), F())

def bounds(value, radius):
    return (value-radius, value+radius)

def outward(interval, digits=24):
    scale = 10**digits
    lo = interval[0]*scale
    hi = interval[1]*scale
    l = lo.numerator//lo.denominator
    h = -((-hi.numerator)//hi.denominator)
    return {"lower": str(F(l, scale)), "upper": str(F(h, scale)),
            "decimal_midpoint": float((interval[0]+interval[1])/2)}

def atan_interval(x, n=50):
    s = sum(((-1)**j*x**(2*j+1)/(2*j+1) for j in range(n)), F())
    nex = (-1)**n*x**(2*n+1)/(2*n+1)
    return min(s, s+nex), max(s, s+nex)

def sin_interval(x, n=25):
    assert 0 <= x <= 1
    s = sum(((-1)**j*x**(2*j+1)/math.factorial(2*j+1) for j in range(n)), F())
    nex = (-1)**n*x**(2*n+1)/math.factorial(2*n+1)
    return min(s, s+nex), max(s, s+nex)

def positive_sin_interval(interval):
    # Directed rational reboxing prevents nested transcendental approximants
    # from creating huge denominators. This is an outward enclosure.
    scale = 10**30
    lower = interval[0]*scale
    upper = interval[1]*scale
    lo = F(lower.numerator//lower.denominator, scale)
    hi = F(-((-upper.numerator)//upper.denominator), scale)
    assert 0 <= lo <= interval[0] <= interval[1] <= hi <= 1
    return sin_interval(lo)[0], sin_interval(hi)[1]

def tensor(*xs):
    out = np.array([[1.0]], dtype=complex)
    for x in xs:
        out = np.kron(out, x)
    return out

def matrices():
    i2 = np.eye(2)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    z = np.diag([1, -1])
    n = np.diag([0, 1])
    na, nb = tensor(n, i2, i2), tensor(i2, n, i2)
    d = na-nb
    nt = na+nb
    zg = tensor(i2, i2, z)
    xg = tensor(i2, i2, x)
    ex4 = np.zeros((4, 4), complex)
    ex4[1, 2] = ex4[2, 1] = 1
    ex = tensor(ex4, i2)
    i8 = np.eye(8)
    rho = tensor(np.diag([.75, .25]), np.diag([.6, .4]),
                 np.array([[.5, .5], [.5, .5]]))
    rho2 = tensor(rho, rho)
    w = tensor(d, zg)+tensor(zg, d)
    checks = {}
    for c in (0, 1):
        h = nt+c*ex-xg
        h2 = tensor(h, i8)+tensor(i8, h)+w/2
        na1, nb1 = tensor(na, i8), tensor(nb, i8)
        ja = 1j*(h2@na1-na1@h2)
        jb = 1j*(h2@nb1-nb1@h2)
        checks[str(c)] = {
            "hermitian_defect": float(np.max(np.abs(h2-h2.conj().T))),
            "local_total_excitation_commutator": float(np.max(np.abs(ja+jb))),
            "initial_total_energy": float(np.trace(rho2@h2).real),
            "initial_exchange_energy": float(np.trace(rho@(c*ex)).real),
            "initial_reciprocal_energy": float(np.trace(rho2@(w/2)).real),
            "current_A_operator_nonzero": bool(np.max(np.abs(ja)) > 0),
        }
        assert checks[str(c)]["hermitian_defect"] == 0
        assert checks[str(c)]["local_total_excitation_commutator"] == 0
        assert abs(checks[str(c)]["initial_total_energy"]+.7) < 1e-14
        assert abs(checks[str(c)]["initial_exchange_energy"]) < 1e-14
        assert abs(checks[str(c)]["initial_reciprocal_energy"]) < 1e-14
        assert checks[str(c)]["current_A_operator_nonzero"] == bool(c)
    return checks

def calculate():
    all_a = [series(c) for c in (0, 1)]
    tail = F(1, 2**ORDER)
    int_tail = T*tail
    ig = [integral(a[5]) for a in all_a]
    tau = [bounds(T+EPS*i, EPS*int_tail) for i in ig]
    xi = [bounds(T+NU*i, NU*int_tail) for i in ig]
    rel = bounds((NU-EPS)*(ig[1]-ig[0]), 2*(NU-EPS)*int_tail)
    assert rel[1] < -F(1, 10**7)
    assert all_a[1][5][4]-all_a[0][5][4] == -F(1, 10)
    assert all(all_a[0][2][n] == 0 for n in range(1, ORDER+1))
    ea = bounds(R0*(evaluate(all_a[1][2])+1)/2, R0*tail/2)
    assert ea[0] > F(21, 10000)
    # Independent signs and the entire power-series conservation certificates.
    assert all_a[1][2][2] == 2
    # Enormous but finite K is handled by integer exponent inequalities only.
    w = 2*V
    lp = V*(2*EPS*F(4, 5)+18*NU/ELL)
    assert 4*lp*T < 2**22
    exponent = 6*w*T
    ceil_exponent = -((-exponent.numerator)//exponent.denominator)
    numerator_power = 22+2*ceil_exponent  # e<4
    sqrt_denominator_power = 900000       # sqrt(10^600000)>2^900000
    eta_power = sqrt_denominator_power-numerator_power
    assert eta_power > 40
    assert F(1, 2**40) < F(1, 10**12)
    # Actual one-cluster marginal: sqrt(alpha)+eta_P plus endpoint input errors.
    assert ea[0]-F(1, 10**12)-2*F(1, 10**10) > F(1, 500)
    # Exact finite-lattice momentum normalization; no continuum replacement.
    a5, a239 = atan_interval(F(1, 5)), atan_interval(F(1, 239))
    pi = (16*a5[0]-4*a239[1], 16*a5[1]-4*a239[0])
    q = (2*pi[0]/65, 2*pi[1]/65)
    sq = positive_sin_interval(q)
    kl = sq[0]/ELL, sq[1]/ELL
    assert F(3, 10) < kl[0] < kl[1] < F(2, 5)
    pd, px = [], []
    for c in (0, 1):
        sd = positive_sin_interval((tau[c][0]/5, tau[c][1]/5))
        sx = positive_sin_interval((kl[0]*xi[c][0], kl[1]*xi[c][1]))
        pd.append((F(1, 2)+sd[0]/2, F(1, 2)+sd[1]/2))
        px.append((F(1, 2)+sx[0]/4, F(1, 2)+sx[1]/4))
    delta = F(1, 10**10)+F(1, 10**12)
    for c in (0, 1):
        assert 0 < 2*(pd[c][0]-delta)-1
        assert 2*(pd[c][1]+delta)-1 < F(29, 1000)
        assert 0 < 4*(px[c][0]-delta)-2
        assert 4*(px[c][1]+delta)-2 < F(68, 1000)
    assert F(100)/(1-F(29, 1000)**2) < 11**2
    assert F(16)/(F(3, 10)**2*(1-F(68, 1000)**2)) < 14**2
    transported = 2*(11+14)*delta
    actual_gap = -rel[1]-transported
    assert actual_gap > F(9, 10**8)
    slope = F(27, 17)
    curve_factor = NU-slope*EPS
    assert curve_factor == F(19, 170)
    curve_ideal = curve_factor*(ig[0]-ig[1]-2*int_tail)
    curve_error = 2*(14+11*slope)*delta
    curve_actual = curve_ideal-curve_error
    assert curve_actual > F(4, 10**8)
    assert EPS*(ig[0]-ig[1]-2*int_tail) > 22*delta
    serbytes = json.dumps([[[str(x) for x in row] for row in a] for a in all_a],
                         separators=(",", ":")).encode()
    return {
        "verdict": "PASS_AUTHOR_FINITE_COMMON_RESPONSE",
        "author_groups": 1,
        "parameters": {"T": str(T), "couplings": "Delta=lambda=kappa=Omega=1",
                       "epsilon": str(EPS), "nu": str(NU), "V": V,
                       "ell": str(ELL), "copies": "10^600000", "order": ORDER},
        "exact_taylor": {
            "coefficient_sha256": hashlib.sha256(serbytes).hexdigest(),
            "majorant_coefficient_bound": "abs(a_n)<=4^n",
            "endpoint_tail_bound": "2^-64",
            "integral_tail_bound": "2^-67",
            "norm_and_energy_coefficients_through_order": ORDER,
            "exchange_induced_gz_t4_coefficient": "-1/10",
            "I_gz": [outward(bounds(i, int_tail)) for i in ig],
            "tau": [outward(i) for i in tau],
            "xi": [outward(i) for i in xi],
            "relative_transport_difference": outward(rel),
            "ideal_relative_gap_strict_lower": "1/10000000",
            "mean_field_energy_A_gain": outward(ea),
            "actual_energy_A_gain_strict_lower": "1/500",
        },
        "finite_resources": {
            "w": w, "L_probe": str(lp), "six_w_T": str(exponent),
            "sqrt_K_lower": "2^900000",
            "eta_P_upper_power_of_two": -eta_power,
            "used_eta_P": "1/1000000000000",
            "new_uniform_preparation_instrument_contract": "1/10000000000",
            "full_BG_kept": True, "unknown_probe_R_uniform": True,
        },
        "finite_ports": {
            "k_lattice": outward(kl), "p_clock": [outward(i) for i in pd],
            "p_position": [outward(i) for i in px],
            "inverse_Lipschitz_clock": 11, "inverse_Lipschitz_position": 14,
            "pair_transport_error": str(transported),
            "actual_relative_gap_lower_from_certificate": outward((actual_gap, actual_gap)),
            "published_actual_relative_gap_strict_lower": "9/100000000",
            "continuum_lattice_error_used": False,
            "baseline_dxi_dtau_upper": str(slope),
            "same_clock_curve_gap_ideal_lower": outward((curve_ideal, curve_ideal)),
            "same_clock_curve_pair_error": str(curve_error),
            "same_clock_curve_gap_actual_lower": outward((curve_actual, curve_actual)),
            "published_same_clock_curve_gap_strict_lower": "1/25000000",
        },
        "matrix_diagnostic": matrices(),
        "scope": {
            "new_quantitative_material_input": True,
            "different_registered_program_values": True,
            "same_other_initial_resources": True,
            "instrument_contract_constructed": False,
            "full_Hilbert_GR_source_certified": False,
            "Jacobian_or_horizon_thermal_identity": False,
            "stable_reference_recovery_certified": False,
            "old_1175_or_1272_numerical_certificate_reused": False,
        },
    }

def build():
    payload = calculate()
    names = ("fixed_spec.md", "proof.md", "admission.md", "source_baseline.json", "check.py")
    payload["source_bindings"] = {n: digest(HERE/n) for n in names}
    return payload

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    result = build()
    target = HERE/"results.json"
    if args.write:
        with target.open("x", encoding="utf-8", newline="\n") as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write("\n")
        mode = "WRITE_ONCE"
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        if saved != result:
            raise SystemExit("READ_ONLY_MISMATCH")
        mode = "READ_ONLY_MATCH"
    print(mode, result["verdict"])
    print("relative_geometry_gap_lower =", result["finite_ports"]["published_actual_relative_gap_strict_lower"])
    print("actual_local_energy_gain_lower = 1/500")
    print("author_groups = 1")

if __name__ == "__main__":
    main()

