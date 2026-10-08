"""1031: a common leading classical parent family, not a full QFT construction.

Full SU(3) matrix Euler/Gauss/Bianchi and Hilbert stress are evaluated before
the one-function reduction is compared. Fractions certify the finite gap;
RK4 is a numerical cross-check, not an EFT or integration error certificate.
Default execution only recomputes and compares. --write creates once.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "common_parent_coupling_freedom_results.json"
TOL = 2e-11


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def generators():
    z = np.zeros((3, 3), dtype=complex)
    result = []
    for i, j in ((0, 1),):
        a, b = z.copy(), z.copy()
        a[i, j] = a[j, i] = .5
        b[i, j], b[j, i] = -.5j, .5j
        result.extend((a, b))
    result.append(np.diag([.5, -.5, 0]).astype(complex))
    for i, j in ((0, 2), (1, 2)):
        a, b = z.copy(), z.copy()
        a[i, j] = a[j, i] = .5
        b[i, j], b[j, i] = -.5j, .5j
        result.extend((a, b))
    result.append(np.diag([1, 1, -2]).astype(complex)/(2*math.sqrt(3)))
    return np.array(result)


def comm(a, b):
    return a @ b-b @ a


def inner(a, b):
    return float(np.real(2*np.trace(a @ b)))


def maxabs(a):
    return float(np.max(np.abs(a)))


def matrix_structure_certificate():
    ts = generators()
    gram = np.array([[inner(a, b) for b in ts] for a in ts])
    structure = np.zeros((8, 8, 8))
    reconstruction = 0.0
    for a in range(8):
        for b in range(8):
            structure[a, b] = [inner(-1j*comm(ts[a], ts[b]), t) for t in ts]
            reconstruction = max(reconstruction, maxabs(
                comm(ts[a], ts[b])-1j*np.einsum("c,cij->ij", structure[a, b], ts)))
    outside = maxabs(structure[:3, :3, 3:])
    jacobi = np.einsum("abe,ecd->abcd", structure, structure)
    jacobi += np.einsum("bce,ead->abcd", structure, structure)
    jacobi += np.einsum("cae,ebd->abcd", structure, structure)
    assert maxabs(gram-np.eye(8)) < TOL
    assert reconstruction < TOL and outside < TOL and maxabs(jacobi) < TOL
    assert abs(structure[0, 1, 2]-1) < TOL
    cases = []
    for g in (1., 1.25, 1.5, 1.75, 2.):
        for label, rep in (("fundamental", ts), ("antifundamental", -ts.conj())):
            coupled = g*rep
            err = 0.0
            for a in range(8):
                for b in range(8):
                    err = max(err, maxabs(comm(coupled[a], coupled[b])
                        -1j*g*np.einsum("c,cij->ij", structure[a, b], coupled)))
            assert err < TOL
            cases.append(dict(g=g, representation=label, covariant_algebra_residual=err,
                              cubic_coefficient=g, quartic_coefficient=g*g))
    return dict(normalization="2 Tr(T_a T_b)=delta_ab; F=dA-i*g[A,A]",
                embedded_generators=[1, 2, 3], outside_subalgebra_residual=outside,
                matrix_commutator_residual=reconstruction, jacobi_residual=maxabs(jacobi),
                coupled_representation_cases=cases)


def connection_and_field(g, f, z, w, ts=None):
    ts = generators() if ts is None else ts
    aa = np.zeros((4, 3, 3), complex)
    da = aa.copy()
    aa[1:], da[1:] = f*ts[:3], z*ts[:3]
    ff = np.zeros((4, 4, 3, 3), complex)
    dff = ff.copy()
    for mu in range(4):
        for nu in range(4):
            ff[mu, nu] = ((da[nu] if mu == 0 else 0)
                         -(da[mu] if nu == 0 else 0)-1j*g*comm(aa[mu], aa[nu]))
            dff[mu, nu] = ((w*ts[nu-1] if mu == 0 and nu else 0)
                         -(w*ts[mu-1] if nu == 0 and mu else 0)
                         -1j*g*(comm(da[mu], aa[nu])+comm(aa[mu], da[nu])))
    return aa, ff, dff


def full_tensor_jet(g, f, z, a, adot, mp2, bare_lambda, vacuum_energy, bad_w=0.):
    """No reduced Euler or stress formula is used to construct these tensors."""
    ts = generators()
    signs = np.array([-1., 1., 1., 1.])
    eta = np.diag(signs)
    hc = adot/a
    w = -2*g*g*f**3+bad_w
    aa, ff, dff = connection_and_field(g, f, z, w, ts)
    raised = signs[:, None, None, None]*signs[None, :, None, None]*ff/a**4
    draised = signs[:, None, None, None]*signs[None, :, None, None]*(dff-4*hc*ff)/a**4
    gamma = np.zeros((4, 4, 4))
    for r in range(4):
        for mu in range(4):
            for nu in range(4):
                gamma[r, mu, nu] = hc*((r == mu)*(nu == 0)
                    +(r == nu)*(mu == 0)-eta[mu, nu]*eta[r, 0])
    euler = np.zeros((4, 3, 3), complex)
    for nu in range(4):
        euler[nu] += draised[0, nu]
        for mu in range(4):
            euler[nu] += -1j*g*comm(aa[mu], raised[mu, nu])
            for lam in range(4):
                euler[nu] += (gamma[mu, mu, lam]*raised[lam, nu]
                              +gamma[nu, mu, lam]*raised[mu, lam])
    bianchi = np.zeros((4, 4, 4, 3, 3), complex)
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for i, j, k in ((mu, nu, rho), (nu, rho, mu), (rho, mu, nu)):
                    bianchi[mu, nu, rho] += ((dff[j, k] if i == 0 else 0)
                                              -1j*g*comm(aa[i], ff[j, k]))
    f2_flat, df2_flat = 0., 0.
    for mu in range(4):
        for nu in range(4):
            f2_flat += signs[mu]*signs[nu]*inner(ff[mu, nu], ff[mu, nu])
            df2_flat += 2*signs[mu]*signs[nu]*inner(ff[mu, nu], dff[mu, nu])
    tflat, dtflat = np.zeros((4, 4)), np.zeros((4, 4))
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                tflat[mu, nu] += signs[rho]*inner(ff[mu, rho], ff[nu, rho])
                dtflat[mu, nu] += signs[rho]*(inner(dff[mu, rho], ff[nu, rho])
                                                     +inner(ff[mu, rho], dff[nu, rho]))
            tflat[mu, nu] -= eta[mu, nu]*f2_flat/4
            dtflat[mu, nu] -= eta[mu, nu]*df2_flat/4
    tcov = tflat/a**2
    tmix = signs[:, None]*tflat/a**4
    dtmix = signs[:, None]*(dtflat-4*hc*tflat)/a**4
    div = np.zeros(4)
    for nu in range(4):
        div[nu] = dtmix[0, nu]
        for mu in range(4):
            for lam in range(4):
                div[nu] += gamma[mu, mu, lam]*tmix[lam, nu]
                div[nu] -= gamma[lam, mu, nu]*tmix[mu, lam]
    h = z*z+g*g*f**4
    expected = np.diag([1.5*h, .5*h, .5*h, .5*h])/a**2
    rho, pressure = tcov[0, 0]/a**2, sum(tcov[i, i] for i in (1, 2, 3))/(3*a*a)
    electric = sum(inner(ff[0, i], ff[0, i]) for i in (1, 2, 3))/(2*a**4)
    magnetic = sum(inner(ff[i, j], ff[i, j]) for i in (1, 2, 3) for j in (1, 2, 3))/(4*a**4)
    lam_eff = bare_lambda+vacuum_energy/mp2
    addot = 2*lam_eff*a**3/3
    einstein = np.diag([3*hc*hc]+[hc*hc-2*addot/a]*3)
    metric = a*a*eta
    total_t = tcov-vacuum_energy*metric
    einstein_residual = einstein+bare_lambda*metric-total_t/mp2
    # Every color projection is checked, including five omitted by the ansatz.
    euler_components = np.array([[inner(ts[c], euler[mu]) for c in range(8)] for mu in range(4)])
    # A constant gauge rotation is independently reconstructed at connection level.
    mat = sum((c+1)*ts[c]/17 for c in range(8))
    ev, u = np.linalg.eigh(mat)
    rot = (u*np.exp(1j*ev)) @ u.conj().T
    aa2 = np.array([rot @ x @ rot.conj().T for x in aa])
    da2 = np.zeros_like(aa2)
    da2[1:] = np.array([rot @ (z*x) @ rot.conj().T for x in ts[:3]])
    ff2 = np.zeros_like(ff)
    for mu in range(4):
        for nu in range(4):
            ff2[mu, nu] = ((da2[nu] if mu == 0 else 0)
                          -(da2[mu] if nu == 0 else 0)-1j*g*comm(aa2[mu], aa2[nu]))
    cov_err = maxabs(ff2-np.array([[rot @ ff[mu, nu] @ rot.conj().T for nu in range(4)] for mu in range(4)]))
    f2_rot = sum(signs[mu]*signs[nu]*inner(ff2[mu, nu], ff2[mu, nu])
                 for mu in range(4) for nu in range(4))/a**4
    scalar_m = .5*(1+(f2_flat/a**4)/(4*rho))
    return dict(g=g, f=f, f_prime=z, scale_factor=a,
        color_matrix_euler_max=maxabs(euler),
        color_euler_max=maxabs(euler_components), gauss_max=maxabs(euler_components[0]),
        inactive_color_euler_max=maxabs(euler_components[:, 3:]),
        bianchi_max=maxabs(bianchi), full_stress_formula_max=maxabs(tcov-expected),
        stress_divergence_max=maxabs(div), pressure_radiation_residual=abs(pressure-rho/3),
        einstein_equation_max=maxabs(einstein_residual), momentum_source_max=maxabs(tcov[0, 1:]),
        source_density=rho, electric_density=electric, magnetic_density=magnetic,
        magnetic_fraction=magnetic/rho, scalar_fraction_residual=abs(scalar_m-magnetic/rho),
        gauge_rotation_F_max=cov_err, gauge_rotation_F_squared_residual=abs(f2_rot-f2_flat/a**4),
        energy_invariant=h, vacuum_energy=vacuum_energy, bare_lambda=bare_lambda,
        lambda_effective=lam_eff)


def tensor_certificate():
    cases = []
    for index, (g, f, z, a, mp2, lam, v0) in enumerate([
        (1., 0., 1., 1., 2., 0., 0.),
        (2., .2, .99, 1.125, 2., 0., 0.),
        (1.5, -.3, .8, 1.4, 3., .02, .03),
        (1.25, .35, -.7, .9, 2.5, .01, .04),
        (1.75, -.2, -.9, 1.2, 4., .03, .02),
    ]):
        energy = z*z+g*g*f**4
        adot = math.sqrt(energy/(2*mp2)+(lam+v0/mp2)*a**4/3)
        data = full_tensor_jet(g, f, z, a, adot, mp2, lam, v0)
        for key, value in data.items():
            if key.endswith(("_max", "_residual")):
                assert value < TOL, (index, key, value)
        cases.append(data)
    # Negative control violates the actual eight-color Euler equations.
    g, f, z, a, mp2 = 1.5, .2, .8, 1.1, 2.
    adot = math.sqrt((z*z+g*g*f**4)/(2*mp2))
    bad = full_tensor_jet(g, f, z, a, adot, mp2, 0., 0., bad_w=.125)
    assert bad["color_euler_max"] > .08
    assert bad["stress_divergence_max"] > .1
    assert bad["gauss_max"] < TOL and bad["bianchi_max"] < TOL
    return dict(samples=len(cases), full_matrix_color_equations_per_sample=32,
                cases=cases, wrong_acceleration_control=bad,
                calibration_scope="All SU(3) components, curvature Bianchi, covariant Hilbert stress, its divergence, Einstein tensor and vacuum bookkeeping at exact ansatz jets; not arbitrary QFT states.")


def rational_bounds_certificate():
    t, g1, g2, p = F(1, 4), F(1), F(2), F(1)
    epsilon = g2*g2*p*p*t**4
    assert epsilon <= 1
    f2_lower = p*t*(1-epsilon/10)
    m2_lower = g2*g2*f2_lower**4/(p*p)
    m1_upper = g1*g1*p*p*t**4
    bernoulli_gap = g2*g2*p*p*t**4*(1-4*epsilon/10)-m1_upper
    assert bernoulli_gap == F(119, 10240)
    assert m2_lower-m1_upper >= bernoulli_gap > 0
    lower_derivative = F(1, 128)*F(637, 640)*F(127, 128)
    upper_derivative = F(1, 64)
    assert lower_derivative == F(80899, 10485760)
    center, constraint_bound = F(3, 2), F(1, 64)
    center_upper = center*center*t**4
    strict_margin = constraint_bound-center_upper
    half_margin_radius = strict_margin/(2*upper_derivative)
    assert strict_margin == F(7, 1024) and half_margin_radius == F(7, 32)
    assert 1 < center-half_margin_radius < center+half_margin_radius < 2
    eps_each = F(1, 1024)
    transported_gap = bernoulli_gap-2*eps_each
    assert transported_gap == F(99, 10240) > 0
    return dict(common_p=p, conformal_endpoint=t, g_interval=[g1, g2],
        bootstrap_condition_g_squared_p_squared_T_fourth=epsilon,
        endpoint_M2_exact_lower=m2_lower, endpoint_M1_exact_upper=m1_upper,
        exact_gap_bound_before_Bernoulli=m2_lower-m1_upper,
        simple_strict_gap=bernoulli_gap,
        derivative_lower_at_T=lower_derivative, derivative_upper_on_full_interval=upper_derivative,
        derivative_formula="d_g M=2*x*y(x)^3*y'(x)/g; x=sqrt(p*g)*eta",
        finite_constraint_example=dict(name="declared classical M(T)<=1/64", center=center,
            bound=constraint_bound, center_prediction_upper=center_upper,
            certified_margin=strict_margin, uniform_Lipschitz_constant=upper_derivative,
            half_margin_radius=half_margin_radius,
            surviving_interval=[center-half_margin_radius, center+half_margin_radius]),
        conditional_error_transport=dict(hypothetical_each_prediction_error=eps_each,
            surviving_prediction_gap=transported_gap,
            actual_EFT_error_certified=False,
            meaning="Only IF independent actual errors obey these bounds; not a measured or derived QFT error."))


def rk4(g, p=1., endpoint=.25, steps=256):
    def rhs(v):
        return np.array([v[1], -2*g*g*v[0]**3])
    step = endpoint/steps
    state = np.array([0., p])
    max_energy = 0.
    for _ in range(steps):
        k1 = rhs(state)
        k2 = rhs(state+step*k1/2)
        k3 = rhs(state+step*k2/2)
        k4 = rhs(state+step*k3)
        state = state+step*(k1+2*k2+2*k3+k4)/6
        max_energy = max(max_energy, abs(state[1]**2+g*g*state[0]**4-p*p))
    return state, max_energy


def trajectory_certificate():
    rows = []
    for g in np.linspace(1, 2, 9):
        coarse, ce = rk4(g, steps=128)
        fine, fe = rk4(g, steps=256)
        yy, _ = rk4(1., endpoint=math.sqrt(g)/4, steps=256)
        scaled = np.array([yy[0]/math.sqrt(g), yy[1]])
        f, z = fine
        m = g*g*f**4
        derivative = 2*(math.sqrt(g)/4)*yy[0]**3*yy[1]/g
        lower = float(F(80899, 10485760))
        assert lower <= derivative <= 1/64
        assert maxabs(fine-coarse) < 5e-11 and maxabs(fine-scaled) < 5e-12
        assert fe < 5e-12
        # Shared exact geometry: p=1, M_P^2=2, Lambda_eff=0, a=1+eta/2.
        a, adot = 9/8, .5
        full = full_tensor_jet(g, f, z, a, adot, 2., 0., 0.)
        assert full["einstein_equation_max"] < 5e-12
        assert abs(full["source_density"]-3/(2*a**4)) < 5e-12
        assert abs(full["magnetic_fraction"]-m) < 5e-12
        rows.append(dict(g=float(g), f=float(f), f_prime=float(z), magnetic_fraction=float(m),
            d_g_M_from_scaling=float(derivative), common_a=a, common_a_prime=adot,
            common_rho_exact=3/(2*a**4),
            coarse_fine_difference=maxabs(fine-coarse), scaling_difference=maxabs(fine-scaled),
            numerical_energy_residual=fe, full_equations_at_numerical_endpoint=full))
    actual_gap = rows[-1]["magnetic_fraction"]-rows[0]["magnetic_fraction"]
    assert actual_gap > float(F(119, 10240))
    return dict(samples=len(rows), coarse_steps=128, fine_steps=256, rows=rows,
        endpoint_numerical_gap=actual_gap,
        numerical_integration_is_rigorous_error_bound=False,
        rigorous_bound_source="rational_bounds_certificate and analytic bootstrap/scaling proof, not RK convergence",
        common_geometry="a=1+eta/2; M_P^2=2; Lambda=Vvac=0; proper time tau=eta+eta^2/4",
        exact_common_proper_time_at_endpoint=F(17, 64),
        exact_common_density_at_endpoint=F(3, 2)/F(9, 8)**4)


def block_diag(*blocks):
    sizes = [x.shape[0] for x in blocks]
    out = np.zeros((sum(sizes), sum(sizes)), complex)
    i = 0
    for b in blocks:
        n = b.shape[0]
        out[i:i+n, i:i+n] = b
        i += n
    return out


def hermitian_basis(n):
    mats = []
    for i in range(n):
        x = np.zeros((n, n), complex)
        x[i, i] = 1
        mats.append(x)
    for i in range(n):
        for j in range(i+1, n):
            x, y = np.zeros((n, n), complex), np.zeros((n, n), complex)
            x[i, j] = x[j, i] = 1
            y[i, j], y[j, i] = 1j, -1j
            mats.extend((x, y))
    return mats


def anomaly_and_mass_certificate():
    ys = [F(1, 6), F(-2, 3), F(1, 3), F(-1, 2), F(1)]
    multiplicities = [6, 3, 3, 2, 1]
    traces = dict(mixed_gravity_U1=sum((n*y for n, y in zip(multiplicities, ys)), F(0)),
        U1_cubed=sum((n*y**3 for n, y in zip(multiplicities, ys)), F(0)),
        SU3_squared_U1_without_common_Dynkin_index=2*ys[0]+ys[1]+ys[2],
        SU2_squared_U1_without_common_Dynkin_index=3*ys[0]+ys[3],
        SU3_cubed_fundamental_units=F(2-1-1))
    assert all(v == 0 for v in traces.values())
    yu = np.diag([1/5, 2/5, 4/5]).astype(complex)
    axis = [1, 2, 3]
    mixing_exact = [[F(i == j)-F(axis[i]*axis[j], 7) for j in range(3)] for i in range(3)]
    orthogonal = [[sum((mixing_exact[i][k]*mixing_exact[j][k] for k in range(3)), F(0))
                   for j in range(3)] for i in range(3)]
    assert orthogonal == [[F(i == j) for j in range(3)] for i in range(3)]
    hd_exact = [[sum((mixing_exact[i][k]*F((k+1)**2, 49)*mixing_exact[j][k]
                      for k in range(3)), F(0)) for j in range(3)] for i in range(3)]
    off_diagonal = [hd_exact[0][1], hd_exact[0][2], hd_exact[1][2]]
    assert off_diagonal == [F(18, 343), F(12, 343), F(6, 343)]
    mixing = np.array([[float(v) for v in row] for row in mixing_exact])
    yd = mixing @ np.diag([1/7, 2/7, 3/7])
    ye = np.diag([1/11, 2/11, 3/11])
    c5 = np.diag([1/13, 2/13, 4/13])
    hu, hd = yu @ yu.conj().T, yd @ yd.conj().T
    columns = []
    for x in hermitian_basis(3):
        v = np.concatenate([comm(x, hu).ravel(), comm(x, hd).ravel()])
        columns.append(np.concatenate([v.real, v.imag]))
    rank = int(np.linalg.matrix_rank(np.array(columns).T, tol=1e-11))
    assert rank == 8
    def dirac(m):
        zz = np.zeros_like(m)
        return np.block([[zz, m], [m.T, zz]])
    v_h, scale5 = 1., 1.
    mu = np.kron(np.eye(3), yu)*v_h/math.sqrt(2)
    md = np.kron(np.eye(3), yd)*v_h/math.sqrt(2)
    me = ye*v_h/math.sqrt(2)
    mn = c5*v_h*v_h/(2*scale5)
    matrix = block_diag(dirac(mu), dirac(md), dirac(me), mn)
    q = np.diag([2/3]*9+[-2/3]*9+[-1/3]*9+[1/3]*9+[-1]*3+[1]*3+[0]*3)
    masses = np.linalg.svd(matrix, compute_uv=False)
    assert matrix.shape == (45, 45) and np.linalg.matrix_rank(matrix) == 45
    assert maxabs(q.T @ matrix+matrix @ q) < TOL
    neutral = np.linalg.svd(mn, compute_uv=False)
    assert min(abs(x-y) for i, x in enumerate(neutral) for y in neutral[i+1:]) > .03
    return dict(left_Weyl_species_per_generation=["q_L", "u_R^c", "d_R^c", "l_L", "e_R^c"],
        hypercharges=ys, multiplicities=multiplicities, exact_anomaly_coefficients=traces,
        total_SU2_doublets=12, SU2_global_parity=0,
        anomaly_dependence_on_positive_g3="Each zero representation coefficient remains zero; no assertion of completed all-sector quantum Ward construction.",
        fixed_Higgs_v=v_h, fixed_Weinberg_scale=scale5,
        exact_flavor_certificate=dict(mixing_matrix=mixing_exact,
            Hu_diagonal=[F(1, 25), F(4, 25), F(16, 25)], Hd=hd_exact,
            Hd_upper_off_diagonal=off_diagonal,
            argument="Distinct Hu eigenvalues force its commutant diagonal; every nonzero Hd off-diagonal forces all three diagonal entries equal. SVD below only cross-checks this proof."),
        quark_common_Hermitian_commutant_dimension=9-rank,
        full_tree_Weyl_mass_rank=45, mass_singular_values=masses,
        neutral_Takagi_masses=neutral, residual_charge_mass_Ward=maxabs(q.T @ matrix+matrix @ q),
        variation_with_g3="All displayed Yukawa/Weinberg matrices and v are fixed. Tree/matching masses and flavor identities are exactly independent of g3.",
        not_identical=["all-order pole masses", "QCD bound-state masses and confinement scale", "running matching coefficients at other scales", "quantum spectra in the time-dependent color background"])


def common_ledger():
    return [
        dict(id="representation_anomaly", status="exact_preserved", condition="Fixed P981 representations; standard chiral anomaly coefficients vanish for every positive g3.", source="1010-1014,1024", not_claimed="anomaly cancellation generates the representations or the complete quantum process"),
        dict(id="vector_vertices", status="exact_preserved", condition="Canonical F3=dA3-i*g3[A3,A3], cubic g3 and quartic g3^2 varied together.", source="1015", not_claimed="g3 value selected"),
        dict(id="matter_Noether", status="exact_preserved", condition="Every colored covariant derivative and related vertex uses the same g3; all remaining couplings fixed.", source="1016,1027", not_claimed="all EFT derivative operators classified"),
        dict(id="universal_gravity", status="exact_preserved", condition="One fixed Einstein coupling and total Hilbert source of the same full leading action; actual nonzero color source and evolving metric included.", source="1027-1029", not_claimed="Einstein action or microscopic common metric generated"),
        dict(id="tree_flavor_mass", status="exact_preserved", condition="v,Yu,Yd,Ye,C5 fixed; full-rank and nondegenerate neutral masses and quark common commutant retain old conditional contracts.", source="1010-1014", not_claimed="exact pole or bound-state masses invariant"),
        dict(id="other_classical_fields", status="exact_solution_branch", condition="Color-singlet Higgs at potential minimum; Aweak=B=0, classical Weyl=0; optional even-Z2 sigma=0; nonzero portal allowed.", source="P981,B993,1028", not_claimed="quantum vacuum has zero fluctuation source"),
        dict(id="finite_same_geometry", status="exact_solution_branch", condition="Same canonical initial electric data, p,v,Mp,a0,Lambda_eff and torus units; identical radiation source and geometry for all g3.", source="1031 analytic branch", not_claimed="same all-field history or realistic cosmology"),
        dict(id="linear_Higgs_coefficients", status="exact_preserved_only", condition="The adopted linear doublet has a=b=1 for all g3.", source="1016,1030", not_claimed="1030 stable physical channel map or actual delta established"),
        dict(id="high_order_matching", status="not_included", condition="Leading minimal action only; unspecified higher-curvature, nonminimal or quantum terms need separate matching and errors.", source="P981 effective scope", not_claimed="all such terms vanish in the real theory"),
        dict(id="quantum_state_source", status="not_proved", condition="No common interacting quantum state, renormalized quantum source or semiclassical error derived.", source="1024 common-object ledger", not_claimed="full QFT completion"),
        dict(id="QCD_running_confinement", status="not_proved", condition="Classical g3 is a fixed matched parameter; no claim that its running, confinement or physical pole spectra are identical.", source="P981 effective scope", not_claimed="real strong-interaction phenomenology recovered"),
        dict(id="FUCP_and_instruments", status="not_proved", condition="The invariant field scalar is a predicted quantity; its autonomous preparation, measurement and record have not been constructed.", source="1009,1024", not_claimed="a countermodel to the full cognitive programme or a complete FUCP physical implementation"),
    ]


def run():
    sources = [
        "archive_956_989/981/drafts/common_parent_contract_v1.md",
        "archive_990_1008/993/common_candidate_v1.md",
        "archive_1009_/1009/input_dependency_ledger_v0_1.md",
        "archive_1009_/research_note_1009.md",
        "archive_1009_/research_note_1015.md",
        "archive_1009_/research_note_1024.md",
        "archive_1009_/research_note_1027.md",
        "archive_1009_/research_note_1029.md",
        "archive_1009_/research_note_1030.md",
        "archive_554_584/research_note_572.md",
        "archive_554_584/research_note_573.md",
        "archive_742_763/research_note_753.md",
    ]
    return encode(dict(round=1031, status="scientific_calibration_verified",
        new_calibration_groups=1, cumulative_test_groups=3808, new_cognitive_axioms=0,
        goal_complete=False, all_scientific_calibrations_passed=True,
        scope="An explicit common leading classical P981/B993 parameter family preserves joined identities and has identical nonzero gravitational backreaction but distinct finite invariant field predictions. Not a full quantum/cognitive countermodel.",
        unit_and_effective_domain="p=1 means the common unit mu=sqrt(p), not retuning a physical cutoff. Ratios v/mu, Mp/mu, matched scales/mu, torus periods and the proper-time map stay fixed across g. Actual higher-order or quantum error control in any proposed common domain remains a separate premise.",
        common_constraint_ledger=common_ledger(),
        su3_structure=matrix_structure_certificate(), full_field_tensors=tensor_certificate(),
        rational_finite_certificate=rational_bounds_certificate(),
        numerical_trajectories=trajectory_certificate(),
        anomaly_and_fixed_tree_mass=anomaly_and_mass_certificate(),
        primary_sources=[dict(url="https://arxiv.org/pdf/1503.05222", sections="II, equations (5)-(8)",
            role="Mature isotropic Yang-Mills ansatz, stress and conformal oscillator; the embedding, common-parent ledger and finite rational certificate are explicitly checked here.")],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        historical_source_sha256={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in sources}))


def compare(actual, expected, path="root"):
    if isinstance(actual, dict):
        assert isinstance(expected, dict) and actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(actual, list):
        assert isinstance(expected, list) and len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(actual, float):
        assert isinstance(expected, (float, int)) and math.isclose(actual, expected, rel_tol=1e-10, abs_tol=TOL), (path, actual, expected)
    else:
        assert actual == expected, (path, actual, expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="passed", round=1031,
        mode="exclusive_first_write" if args.write else "read_only_recompute_compare",
        full_tensor_samples=result["full_field_tensors"]["samples"],
        trajectory_samples=result["numerical_trajectories"]["samples"],
        exact_finite_gap=result["rational_finite_certificate"]["simple_strict_gap"],
        numerical_gap=result["numerical_trajectories"]["endpoint_numerical_gap"],
        historical_sources=len(result["historical_source_sha256"]), output=str(OUT)), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
