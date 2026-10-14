"""Independent audit of the finite-head algebra and quadrupole normalization.

No author code is imported.  Uses orthonormal head matrices and a Laurent
polynomial time integral, independent of the author's rational-Gram route.
The added exact arithmetic checks the continuous-parent norm and Dyson budget.
It is not a field simulation or an apparatus matching certificate.
"""
from fractions import Fraction as F
from math import comb, pi
import json
import numpy as np


def exact(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def comm(a, b):
    return a @ b - b @ a


def parent_certificate(ic):
    # Half-beta recurrence, independent of the author's polynomial expansion.
    def radial(r, n):
        result = F(1, 2*r+1)
        for j in range(1, n+1):
            result *= F(2*j, 2*r+2*j+1)
        return result

    radius = F(1, 4)
    norm = radial(1, 4)
    i0_pi = radial(1, 8)/(4*radius**3*norm**2)
    i2_pi = 16*radial(2, 6)/(radius**5*norm**2)
    assert i0_pi == F(776160, 4199)
    assert i2_pi == F(10644480, 221)
    i0, i2 = i0_pi/3, i2_pi/3
    im2 = F(2, 9)+i0/16
    l2 = 72**2*(F(4, 9)+i0/4)
    s2 = 72**2*im2
    w = 98*im2+388*i0+2*i2
    assert l2 < 300**2 and s2 < 150**2 and w < 60000

    # Direct exact 3x3 TT projector, including density/current in the
    # covariant numerator. These are checks of the 940 identity, not an
    # assertion that the TT piece alone has causal support.
    k = [F(2), F(3), F(5)]
    k2 = sum(x*x for x in k)
    P = [[F(i == j)-k[i]*k[j]/k2 for j in range(3)] for i in range(3)]
    def mul(a, b):
        return [[sum((a[i][l]*b[l][j] for l in range(3)), F())
                 for j in range(3)] for i in range(3)]
    tt_cases = []
    for omega in (F(1), F(7), F(11, 3)):
        e = [F(1), F(-1), F(0)]
        rho = -sum(e[i]*k[i]**2 for i in range(3))
        current = [-omega*e[i]*k[i] for i in range(3)]
        tensor = [[-omega**2*e[i] if i == j else F()
                   for j in range(3)] for i in range(3)]
        assert omega*rho == sum(k[i]*current[i] for i in range(3))
        assert all(omega*current[j] == sum(k[i]*tensor[i][j] for i in range(3))
                   for j in range(3))
        ptp = mul(mul(P, tensor), P)
        traceptp = sum(ptp[i][i] for i in range(3))
        tt = [[ptp[i][j]-P[i][j]*traceptp/2 for j in range(3)] for i in range(3)]
        normtt = sum(tt[i][j]**2 for i in range(3) for j in range(3))
        normv = sum(x*x for x in current)
        vl2 = sum(k[i]*current[i] for i in range(3))**2/k2
        sigma = sum(tensor[i][i] for i in range(3))
        constraint = rho*rho/2+rho*sigma-2*normv+vl2/2
        covariant = (rho*rho-2*normv+
                     sum(tensor[i][j]**2 for i in range(3) for j in range(3))
                     -(rho-sigma)**2/2)
        defect = covariant/(k2-omega**2)-normtt/(k2-omega**2)-constraint/k2
        assert defect == 0
        tt_cases.append({"omega": exact(omega), "complete_kernel_defect": exact(defect)})

    # The l linear / k bounded-constraint words have vacuum norm
    # <= u^l/sqrt(l!) * v^k/k!. Use a geometric majorant for the tail.
    g, time_upper = F(1, 10**15), F(22, 7)
    u, v = 600*g*time_upper, 60000*g*g*time_upper
    d = u*u+v
    tail = (u**3/(1-u)+u*v/((1-u)*(1-v))+v*v/(1-v))
    single = 2*tail*(1+u+d)+tail*tail+2*u*d+d*d
    pair = 2*single
    leading = g*g*abs(ic)
    lower = leading-pair
    assert lower > F(7, 10**34)
    positivity = 2-g*g*(2*150**2+60000)
    assert positivity > 1
    return {
        "radial_method": "exact half-beta recurrence",
        "pi_times_I0": exact(i0_pi), "pi_times_I2": exact(i2_pi),
        "I_minus2_upper": exact(im2), "L_squared_upper": exact(l2),
        "S_squared_upper": exact(s2), "W_norm_upper": exact(w),
        "TT_constraint_exact_cases": tt_cases,
        "single_probability_tail_exact": exact(single),
        "pair_probability_tail_exact": exact(pair),
        "normal_parent_difference_lower_exact": exact(lower),
        "normal_parent_difference_lower_display": float(lower),
        "positivity_lower_exact": exact(positivity),
        "actual_instrument_matching": False,
        "nonlinear_GR_or_physical_xi_selection": False,
    }


def run():
    # Fourier expansion of a(t)=5 sin(t)^6/16. All frequencies are even.
    coeff = {6 - 2*j: -F(5*comb(6, j)*(-1)**j, 1024) for j in range(7)}
    a2_pi = sum((c*coeff[-k] for k, c in coeff.items()), F(0))
    adot2_pi = sum((k*k*c*coeff[-k] for k, c in coeff.items()), F(0))
    radius, distance = F(1, 4), F(4)
    m2 = 2*F(3, 13)*radius**2
    c0_pi = F(3, 2)/distance**5
    c1_pi = -(3*distance**2+m2)/(4*distance**5)
    ic = c0_pi*a2_pi+c1_pi*adot2_pi
    assert a2_pi == F(5775, 262144)
    assert adot2_pi == F(4725, 65536)
    assert ic == -F(22700475, 27917287424)

    h = 3*np.eye(7, dtype=complex)
    for j in range(6):
        h[j, j+1] = h[j+1, j] = -np.sqrt((j+1)*(6-j))/2
    p = np.zeros((7, 7), dtype=complex)
    p[3, 3] = 1
    ad = [p]
    for _ in range(3):
        ad.append(comm(h, ad[-1]))
    endpoint_residual = max(np.linalg.norm(ad[j][:, k])
                            for j in range(3) for k in (0, 6))
    third_endpoint_norm = min(np.linalg.norm(ad[3][:, k]) for k in (0, 6))
    assert endpoint_residual < 1e-14 and third_endpoint_norm > 1
    dh = 1j*ad[1]
    dh_norm_squared = np.linalg.norm(dh, 2)**2
    assert abs(dh_norm_squared-6) < 1e-13
    w, v = np.linalg.eigh(h)
    upi = (v*np.exp(-1j*pi*w)) @ v.conj().T
    mirror_error = np.linalg.norm(upi-np.eye(7)[:, ::-1], 2)
    mean_error = 0.0
    for t in (0.0, pi/7, pi/4, pi/2, 3*pi/4, pi):
        psi = v @ (np.exp(-1j*t*w)*v[0, :].conj())
        mean_error = max(mean_error, abs(abs(psi[3])**2-5*np.sin(t)**6/16))
    assert mirror_error < 2e-14 and mean_error < 1e-14

    # Complete two-head/two-data finite operator. This is ONLY the instantaneous
    # source block, not a proposed replacement for TT + constraint propagation.
    z = np.diag([1, -1]).astype(complex)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    ey = (np.eye(2)+y)/2
    pauli_error = np.linalg.norm(comm(z, ey)+1j*x, 2)
    cmat = float(c0_pi)/pi*np.kron(p, p)+float(c1_pi)/pi*np.kron(dh, dh)
    complete_c = np.kron(cmat, np.kron(z, z))
    effect = np.kron(np.eye(49), np.kron(np.eye(2), ey))
    observable_derivative = 1j*comm(-complete_c, effect)
    assert np.linalg.norm(observable_derivative-observable_derivative.conj().T) < 1e-14
    norm_bound = (float(c0_pi)+6*abs(float(c1_pi)))/pi
    complete_norm = np.linalg.norm(complete_c, 2)
    assert complete_norm <= norm_bound + 1e-14

    # Integrate the exact free witness by independent Gauss-Legendre quadrature.
    # This is a numerical cross-check; the Laurent coefficients above are exact.
    nodes, weights = np.polynomial.legendre.leggauss(64)
    plusx = np.array([1, 1], complex)/np.sqrt(2)
    probability_derivative = 0.0
    for node, weight in zip(nodes, weights):
        t = (node+1)*pi/2
        hc = v @ (np.exp(-1j*t*w)*v[0, :].conj())
        heads = np.kron(hc, hc)
        vp = np.kron(heads, np.kron([1, 0], plusx))
        vm = np.kron(heads, np.kron([0, 1], plusx))
        difference = np.vdot(vp, observable_derivative@vp)-np.vdot(vm, observable_derivative@vm)
        assert abs(difference.imag) < 1e-13
        probability_derivative += weight*pi/2*difference.real
    coefficient_exact = -2*ic  # p(A+) - p(A-), H_inst,AB = -g^2 xi C_AB.
    assert abs(probability_derivative-float(coefficient_exact)) < 1e-14

    return {
        "status": "PASS_WITH_SCOPE",
        "source_algebra": {
            "head_spectrum_max_error": float(np.max(np.abs(w-np.arange(7)))),
            "endpoint_ad_0_1_2_residual": float(endpoint_residual),
            "third_ad_endpoint_norm": float(third_endpoint_norm),
            "head_mirror_operator_error": float(mirror_error),
            "mean_P_max_error": float(mean_error),
            "derivative_norm_squared": float(dh_norm_squared),
            "Pauli_commutator_error": float(pauli_error),
        },
        "normalization": {
            "integral_a2_over_pi": exact(a2_pi),
            "integral_adot2_over_pi": exact(adot2_pi),
            "I_C": exact(ic),
            "leading_probability_coefficient": exact(coefficient_exact),
            "independent_full_operator_quadrature": float(probability_derivative),
            "quadrature_difference": float(abs(probability_derivative-float(coefficient_exact))),
            "complete_C_dimension": int(complete_c.shape[0]),
            "complete_C_norm_diagnostic": float(complete_norm),
            "analytic_C_norm_upper": float(norm_bound),
        },
        "continuous_parent_budget": parent_certificate(ic),
        "scope": {
            "field_parent_certified_by_this_script": False,
            "COM_or_pointer_implementation_certified": False,
            "finite_coupling_arithmetic_checked": True,
            "analytic_field_domain_and_retarded_support_require_proof": True,
            "no_author_code_imported": True,
        },
    }


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2, sort_keys=True))
