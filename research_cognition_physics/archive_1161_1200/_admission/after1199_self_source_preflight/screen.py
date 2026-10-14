"""Unnumbered preflight: a retained ADM mixed term and the original 941 clock jet.

Exact rational bounds certify formulas in proof.md. Fock matrices and quadrature
are diagnostics on the stated polynomial support, NOT a finite-time truncation
certificate or an error bound relative to the full nonlinear ADM process.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def run():
    # Symmetry of G combines the two cross terms of (E+rho)G(-7A+B-rho)/8.
    mixed_a = F(-7 - 1, 8)
    mixed_b = F(1 - 1, 8)
    self_rho = F(-1, 8)
    assert (mixed_a, mixed_b, self_rho) == (F(-1), F(0), F(-1, 8))

    # Load the actual 925 model directly; 941 imports it through its historical
    # runtime layout. No historical source or saved result is modified here.
    old = ROOT / 'research_cognition_physics/archive_923_934/925/record_resource_exchange.py'
    spec = importlib.util.spec_from_file_location('preflight_original925', old)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    z = module.model()
    ids = [z['idx'](0, 0, 3), z['idx'](1, 0, 3)]
    embed = np.eye(z['d'])[:, ids]
    hsmall = embed.T @ z['h'] @ embed
    leakage = np.linalg.norm(z['h'] @ embed - embed @ hsmall)
    sigma_x = np.array([[0., 1.], [1., 0.]])
    m0, g, eps, bw = F(526, 5), F(2, 5), F(1, 1000), F(1, 12)
    mass = 100 * np.eye(2) + hsmall
    assert leakage < 1e-13
    assert np.linalg.norm(mass - float(m0) * np.eye(2) - float(g) * sigma_x) < 1e-13

    pi = np.pi
    omega = np.array([pi, 3 * pi])
    W = np.array([[1 / 12 + 1 / (4*pi*pi), -3 / (16*pi*pi)],
                  [-3 / (16*pi*pi), 1 / 12 + 1 / (36*pi*pi)]])
    nodes, weights = np.polynomial.legendre.leggauss(96)
    zz, weights = (nodes + 1) / 2, weights / 2
    modes = np.sqrt(2) * np.sin(np.outer(zz, omega))
    wint = modes.T @ ((weights * zz * (1-zz)/2)[:, None] * modes)
    quad_error = float(np.max(np.abs(wint-W)))
    assert quad_error < 2e-14 and np.linalg.eigvalsh(W)[0] > 0

    # L=5 retains all intermediate states needed for X|psi>, whose total
    # occupation is <=3. Norm(X psi)^2 computes the full polynomial moment.
    L = 5
    annihilator = np.diag(np.sqrt(np.arange(1, L)), 1)
    single_p = 1j * (annihilator.T - annihilator) / np.sqrt(2)
    ps = [np.sqrt(omega[0]) * np.kron(single_p, np.eye(L)),
          np.sqrt(omega[1]) * np.kron(np.eye(L), single_p)]
    xgeom = sum(W[i, j] * ps[i] @ ps[j] for i in range(2) for j in range(2))
    hgeo = np.diag([omega[0]*(n+.5) + omega[1]*(k+.5)
                    for n in range(L) for k in range(L)])
    e = np.eye(L*L)
    states = [(e[:, L] + sign * e[:, 1]) / np.sqrt(2) for sign in (1, -1)]
    moments = [(float(np.vdot(s, xgeom@s).real),
                float(np.vdot(xgeom@s, xgeom@s).real)) for s in states]
    dx = moments[0][0] - moments[1][0]
    dx2 = moments[0][1] - moments[1][1]
    dx_formula = -3*np.sqrt(3)/(8*pi)
    assert abs(dx - dx_formula) < 1e-13
    assert abs(dx2 - (pi+1/pi)*dx_formula) < 1e-12
    high_support = [n*L+k for n in range(L) for k in range(L) if n+k > 3]
    assert max(np.linalg.norm((xgeom@s)[high_support]) for s in states) < 1e-13

    alpha = 1-eps*bw*m0/4
    # Independently compute p_B''(0)=2||Pi_B H |A,psi>||^2 in the full matrix.
    H = (np.kron(mass, e) + np.kron(np.eye(2), hgeo)
         - float(eps/2)*np.kron(mass, xgeom)
         - float(eps*bw/8)*np.kron(mass@mass, e))
    jets = []
    for s, (x1, x2) in zip(states, moments):
        state = np.concatenate([s, np.zeros_like(s)])
        hb = (H@state)[L*L:]
        actual = float(2*np.vdot(hb, hb).real)
        formula = 2*float(g*g)*(float(alpha*alpha)-float(alpha*eps)*x1+float(eps*eps/4)*x2)
        assert abs(actual-formula) < 1e-13
        jets.append(actual)
    delta_jet_formula = 2*float(g*g*eps)*(-dx_formula)*(float(alpha)-float(eps)*(pi+1/pi)/4)
    assert abs(jets[0]-jets[1]-delta_jet_formula) < 1e-13

    # Rigorous inequalities: 3<pi<22/7 and sqrt(3)>17/10.
    minus_dx_lower = F(357, 1760)
    beta_lower = alpha - eps * (F(22, 7) + F(1, 3)) / 4
    jet_lower = 2*g*g*eps*minus_dx_lower*beta_lower
    target = F(64, 10**6)
    assert beta_lower > 0 and jet_lower > target
    # Positivity of retained quadratic full-Fock model, not full ADM stability.
    mass_max = m0 + g
    trace_w_upper = F(16, 81)
    kinetic_margin_lower = 1-eps*mass_max*trace_w_upper
    mass_factor_lower = 1-eps*bw*mass_max/8
    assert kinetic_margin_lower > 0 and mass_factor_lower > 0
    return {
        'schema': 'after1199_adm_original_clock_initial_jet_preflight_v1',
        'status': 'PASS_SCOPED_PREFLIGHT',
        'numbered_round': None,
        'exact_mixed_coefficients': {'rho_G_A': str(mixed_a), 'rho_G_B': str(mixed_b), 'rho_G_rho': str(self_rho)},
        'original_model_path': old.relative_to(ROOT).as_posix(),
        'original_model_sha256': hashlib.sha256(old.read_bytes()).hexdigest(),
        'clock': {'mass_center': str(m0), 'mass_x': str(g), 'initial': 'original A,q=0,n=3', 'effect': 'original B,q=0,n=3 projection', 'invariant_pair_leakage_diagnostic': float(leakage)},
        'inputs': {'epsilon': str(eps), 'coordinate_profile': 'w(z)=1 on [0,1]', 'boundary': 'Dirichlet Green; unit transverse torus', 'modes': [1, 3]},
        'exact_certificate': {'minus_delta_X_strict_lower': str(minus_dx_lower), 'beta_strict_lower': str(beta_lower), 'delta_p_second_derivative_strict_lower': str(jet_lower), 'target': str(target), 'margin_over_target': str(jet_lower-target), 'retained_full_fock_kinetic_margin_strict_lower': str(kinetic_margin_lower), 'retained_background_mass_factor_lower': str(mass_factor_lower)},
        'numerical_diagnostics': {'quadrature_max_error': quad_error, 'X_moments_plus_minus': moments, 'initial_p_second_derivatives_plus_minus': jets, 'delta_jet': float(jets[0]-jets[1]), 'fock_levels_per_mode': L, 'maximum_required_total_occupation_for_X_action': 3},
        'not_certified': ['finite_time_probability_difference_or_remainder', 'full_unknown_qutrit_and_reference_coupled_channel', 'stationary_support_or_complete_boundary_resource', 'full_ADM_quantum_dynamics_or_its_remainder', 'GR_generation', 'dark_matter_or_dark_energy_identity'],
    }


if __name__ == '__main__':
    print(json.dumps(run(), ensure_ascii=False, indent=2))
