"""1267: finite terminal instruments; one author group.

No ODE, seven-pair wavefunction propagation or device test.
Default: recompute and compare saved JSON without writes.
--write: exclusive first creation of instrument_results.json.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESULT = HERE / "instrument_results.json"

def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def fstr(q):
    q = Q(q)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"

def expminus_interval(x, even_order=20):
    x = Q(x)
    assert 0 <= x <= 1 and even_order % 2 == 0
    term, accum = Q(1), Q(1)
    lower = None
    for j in range(1, even_order + 1):
        term *= -x / j
        accum += term
        if j == even_order - 1:
            lower = accum
    assert lower <= accum
    return lower, accum

def exact_fourier_tests():
    # The continuous theorem is the Fourier integral proof, not this finite menu.
    menus = [(Q(3,5), Q(1)), (Q(-2,7), Q(4)), (Q(1,11), Q(1,2))]
    F2 = sum(abs(c) * k2 for c,k2 in menus)
    out = []
    for eps in (Q(1,8), Q(3,10), Q(1,32)):
        plain = aw = weighted = Q(0)
        character = []
        for coeff,k2 in menus:
            loq, hiq = expminus_interval(eps*k2/4)
            loh, hih = expminus_interval(eps*k2/2)
            defect2_upper = 2*(1-loq)
            assert defect2_upper <= eps*k2/2
            character.append({
                "norm_xi_squared": fstr(k2),
                "defect_squared_upper": fstr(defect2_upper),
                "epsilon_xi_squared_over_2": fstr(eps*k2/2)})
            plain += abs(coeff)*(1-loh)
            aw += abs(coeff)*hiq*(1-loh)
            weighted += abs(coeff)*(1-loq)
        assert plain <= eps*F2/2
        assert aw <= eps*F2/2
        assert weighted <= eps*F2/4
        out.append({
            "epsilon": fstr(eps), "F2": fstr(F2),
            "weyl_heat_difference_upper": fstr(plain),
            "aw_heat_difference_upper": fstr(aw),
            "record_weyl_heat_upper": fstr(weighted),
            "nonselective_bound": fstr(eps*F2/2),
            "weighted_bound": fstr(eps*F2/4),
            "character_defects": character})
    # 6/pi < 2 follows from pi > 3. Factor constants use fourth powers.
    assert Q(2)**4 * Q(2)**3 < Q(4)**4  # 2 * 2^(3/4) < 4
    assert Q(4) + Q(1) + Q(1,4) < 8
    return out

def gaussian_tests():
    # GH integrates a 14D product Gaussian without making a 48^14 mesh.
    nodes, weights = np.polynomial.hermite.hermgauss(48)
    weights = weights / np.sqrt(np.pi)
    mu = np.arange(14,dtype=float)/11 - 0.5
    xis = [np.arange(14,dtype=float)/23,
           np.cos(np.arange(14,dtype=float))/3,
           np.concatenate(([0.7,-0.4],np.zeros(12)))]
    max_error = 0.
    for eps in (0.125, 0.3):
        z = mu[:,None] + np.sqrt(2*eps)*nodes[None,:]
        for xi in xis:
            qmean = np.prod(np.sum(weights[None,:] *
                                   np.exp(1j*xi[:,None]*z), axis=1))
            expected = np.exp(1j*np.dot(xi,mu)-eps*np.dot(xi,xi)/2)
            post = np.exp(-eps*np.dot(xi,xi)/4)*qmean
            post_expected = np.exp(1j*np.dot(xi,mu)-3*eps*np.dot(xi,xi)/4)
            max_error = max(max_error, abs(qmean-expected),
                            abs(post-post_expected))
    one = np.sum(weights*np.exp(-nodes**2))
    returned = one**14
    assert max_error < 2e-13
    assert abs(returned - 1/128) < 2e-13
    return {
        "quadrature": "48-point Gauss-Hermite, factorized 14 dimensions",
        "maximum_character_error": f"{max_error:.3e}",
        "measure_prepare_return_probability": f"{returned:.14g}",
        "analytic_return_probability": "1/128",
        "analytic_half_trace_lower_bound": "127/128",
        "scope": "numerical diagnostic; MP return is not a Luders claim"}

def opnorm(a):
    return float(np.linalg.norm(a,2))

def root_psd(a):
    e,v = np.linalg.eigh((a+a.conj().T)/2)
    assert e.min() > -1e-13
    return (v*np.sqrt(np.maximum(e,0))) @ v.conj().T

def finite_frame_tests():
    # A six-point Parseval frame, not a truncation of the original field.
    row = np.arange(6)[:,None]
    col = np.arange(3)[None,:]
    W = np.exp(2j*np.pi*row*col/6)/np.sqrt(6)
    b = np.array([.2,-.3,.7,-.5,.4,.1])
    Mb = np.diag(b)
    B = np.array([[.3,.08j,.12],[-.08j,-.2,.05j],
                  [.12,-.05j,.1]], dtype=complex)
    T = Mb @ W - W @ B
    E = []
    P = []
    for c in range(3):
        mask = np.zeros(6)
        mask[2*c:2*c+2] = 1
        Pc = np.diag(mask)
        P.append(Pc)
        E.append(W.conj().T @ Pc @ W)
    roots = [root_psd(e) for e in E]
    V = np.vstack(roots)
    residuals = [opnorm(W.conj().T@W-np.eye(3)),
                 opnorm(sum(E)-np.eye(3)),
                 opnorm(V.conj().T@V-np.eye(3))]
    for e,p,r in zip(E,P,roots):
        comm = e@B-B@e
        rhs = -W.conj().T@p@T+T.conj().T@p@W
        residuals.append(opnorm(comm-rhs))
        delta = opnorm(comm)
        assert delta <= 2*opnorm(T)+1e-13
        assert opnorm(r@B-B@r) <= 2*np.sqrt(opnorm(B)*delta)+1e-13
    Lstar = sum(r@B@r for r in roots)
    qbound = 2*len(E)*np.sqrt(opnorm(B)*2*opnorm(T))
    assert opnorm(Lstar-B) <= qbound+1e-13
    rw = np.array([1.,-.4,.25])
    mr = np.repeat(rw,2)
    target = W.conj().T@np.diag(mr*b)@W
    lhs = sum(w*r@B@r for w,r in zip(rw,roots))-target
    rhs = sum(w*(r@B@r-e@B) for w,r,e in zip(rw,roots,E))
    rhs -= W.conj().T@np.diag(mr)@T
    residuals.append(opnorm(lhs-rhs))
    assert opnorm(lhs) <= qbound+opnorm(T)+1e-13

    # Full reference check on a maximally entangled 3x3 input.
    psi = np.eye(3,dtype=complex).reshape(-1)/np.sqrt(3)
    rho = np.outer(psi,psi.conj())
    outputs = []
    for r in roots:
        kr = np.kron(r,np.eye(3))
        outputs.append(kr@rho@kr.conj().T)
    out = sum(outputs)
    ref_before = np.trace(rho.reshape(3,3,3,3),axis1=0,axis2=2)
    ref_after = np.trace(out.reshape(3,3,3,3),axis1=0,axis2=2)
    residuals.append(opnorm(ref_before-ref_after))
    residuals.append(abs(np.trace(out)-1))
    for o in outputs:
        assert np.linalg.eigvalsh(o).min() > -1e-13
    assert max(residuals) < 2e-13
    return {
        "frame_shape": [6,3], "outcomes": 3,
        "max_algebra_or_reference_residual": f"{max(residuals):.3e}",
        "source_commutator_bound": f"{qbound:.12g}",
        "actual_nonselective_source_difference": f"{opnorm(Lstar-B):.12g}",
        "weighted_source_difference": f"{opnorm(lhs):.12g}",
        "scope": "finite identity diagnostic, not continuous propagation or physical realization"}

def exact_resources():
    old = read_json(STAGE/"1249/quantum_bounds_results.json")
    ledger = old["group_constructive_finite_resource_bound"]["ledger"]
    current = read_json(STAGE/"1266/results.json")
    inherited = current["inherited_constants"]
    M = 2**3059
    EH, E0 = 3093, 14008
    EF = 95928+100*M
    ceildiv = lambda a,b: (a+b-1)//b
    N = 1+max(ceildiv(EF+42,3),ceildiv(EH+EF+42,6))
    assert N == int(inherited["N_decimal"])
    assert N == int(ledger["integer_n_power_exponent_decimal"])
    assert EH == inherited["EH"] == ledger["Fourier_H_exponent"]
    assert E0 == inherited["E0"] == ledger["Fourier_menu_base_exponent"]
    assert EF == int(ledger["Fourier_transport_exponent_decimal"])
    cbeta = factorial(81)//factorial(40)**2
    assert cbeta*factorial(40)**2 == factorial(81)
    Bplat = 324*cbeta*2**40
    Lc = 56*Bplat
    assert Bplat == int(ledger["B_decimal"])
    assert Lc == int(ledger["Lcut_decimal"])
    As = 2**98*(1+Lc)**2
    D = (2*As*10**10+1)**3
    assert D < 2**4096
    es = E0+EH+15*M+22-2*N
    assert es == int(current["authoritative_negative_binary_exponents"]
                    ["bounded_source_error_upper_power2"])
    eL = E0+4098-N
    eR = E0+4099-N
    eMP = E0-1-4*N
    assert max(es,eL,eR,eMP) < -82
    assert Q(1,2**82)+Q(1,2**82) == Q(1,2**81)
    return {
        "Bplat_decimal": str(Bplat), "Lcut_decimal": str(Lc),
        "A_s_decimal": str(As),
        "D_upper_decimal": str(D),
        "D_bit_length": D.bit_length(),
        "D_less_than_2_power_4096": True,
        "N_decimal": str(N), "E0": E0, "EH": EH,
        "negative_binary_exponents": {
            "original_kappa3_pre_source": str(es),
            "Luders_nonselective": str(eL),
            "Luders_record_weighted_Weyl_or_AW": str(eR),
            "measure_prepare_nonselective": str(eMP)},
        "all_exponents_strictly_below_minus_82": True,
        "original_kappa3_total_post_source_strict_bound": "2^-81 < 2^-80",
        "other_parameter_families": "their own pre-error plus universal instrument error"}

def result():
    bindings = {
        "1267/instrument_proof.md": HERE/"instrument_proof.md",
        "1267/instrument_check.py": Path(__file__).resolve(),
        "1267/admission.md": HERE/"admission.md",
        "1249/quantum_bridge.md": STAGE/"1249/quantum_bridge.md",
        "1249/quantum_bounds_results.json": STAGE/"1249/quantum_bounds_results.json",
        "1251/proof.md": STAGE/"1251/proof.md",
        "1251/results.json": STAGE/"1251/results.json",
        "1266/proof.md": STAGE/"1266/proof.md",
        "1266/results.json": STAGE/"1266/results.json",
        "_admission/finite_common_task_after1266/README.md":
            STAGE/"_admission/finite_common_task_after1266/README.md",
        "_admission/finite_common_task_after1266/integration_checks.json":
            STAGE/"_admission/finite_common_task_after1266/integration_checks.json",
    }
    return {
        "round":1267, "part":"finite_terminal_instrument",
        "status":"PASS_AUTHOR_WITH_SCOPE", "instrument_author_groups":1,
        "primary_instrument":"finite-outcome Luders; two finite D label registers",
        "comparison_instrument":"coherent measure-and-prepare; separate post-state",
        "exact_fourier_checks": exact_fourier_tests(),
        "gaussian_diagnostic": gaussian_tests(),
        "finite_frame_diagnostic": finite_frame_tests(),
        "exact_resources": exact_resources(),
        "scope": [
            "normal CP and reference extension exact; physical implementation not certified",
            "same finite effects; marginal Luders and fine Luders post-states differ",
            "universal bounded-source operator bounds, no arbitrary input recovery",
            "2^-80 total source only for original kappa=3 two-program node states",
            "kappa family uses fixed a marginal and its own pre-measurement error",
            "new finite register physical energy, complete sources and CC remain inputs",
            "MP return 1/128 is not asserted for the Luders instrument",
            "no unbounded stress, psi^2, quantum field simulation or device performance claim"
        ],
        "source_bindings": {name:digest(path) for name,path in bindings.items()}
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    data = result()
    if args.write:
        with RESULT.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2,sort_keys=True)
            stream.write("\n")
        print("PASS_AUTHOR_WITH_SCOPE FIRST_WRITE_EXCLUSIVE")
    else:
        saved=read_json(RESULT)
        if data != saved:
            keys=[k for k in set(data)|set(saved) if data.get(k)!=saved.get(k)]
            raise SystemExit("MISMATCH: "+", ".join(sorted(keys)))
        print("PASS_AUTHOR_WITH_SCOPE READ_ONLY_MATCH")
    print("instrument_author_groups=1")
    print("D_bit_length="+str(data["exact_resources"]["D_bit_length"]))
    print("same_resource_post_source<2^-80; finite_CP_environment_not_physical_energy_certificate")

if __name__=="__main__":
    main()
