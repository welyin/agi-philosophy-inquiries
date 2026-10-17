"""Check one normal outgoing pulse by independent source/flux integrals.

Analytic identities and scope are in proof.md. Floating-point quadrature is
a diagnostic, not an error certificate for a coupled Einstein-matter system.
The default command is read-only; --write refuses to overwrite saved results.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
H, A, B, L = 1.0, 1.0, 2.0, 3.0
M2 = 2.0 * H * H


def profile(s):
    s = np.asarray(s, dtype=float)
    y = 2.0 * (s - A) / (B - A) - 1.0
    f, fp = np.zeros_like(s), np.zeros_like(s)
    inside = np.abs(y) < 1.0
    q = y[inside]
    d = 1.0 - q*q
    f[inside] = np.exp(1.0 - 1.0/d)
    fp[inside] = -4.0*q/(B-A)/d**2 * f[inside]
    return f, fp


def nodes(lo, hi, n):
    z, w = np.polynomial.legendre.leggauss(n)
    return lo+(hi-lo)*(z+1.0)/2.0, (hi-lo)*w/2.0


def integrals(n):
    # Reduced wave energy, integrated in tortoise coordinate.
    x, wx = nodes(A, B, n)
    f, fp = profile(x)
    q_wave = 4.0*math.pi*float(wx @ (fp*fp))

    # Independently integrate physical Hilbert T_tt/F in areal radius.
    r, wr = nodes(math.tanh(H*A)/H, math.tanh(H*B)/H, n)
    F = 1.0-H*H*r*r
    x = np.arctanh(H*r)/H
    f, fp = profile(x)
    phi, phi_t = f/r, -fp/r
    phi_r = fp/(r*F)-f/(r*r)
    density = 2.0*math.pi*(r*r/F*phi_t**2
        + r*r*F*phi_r**2 + M2*r*r*phi**2)
    q_initial = float(wr @ density)

    # Physical outward current -4 pi R^2 T^r_t at a finite static sphere.
    R = math.tanh(H*L)/H
    FR = 1.0-H*H*R*R
    t, wt = nodes(L-B, L-A, n)
    f, fp = profile(L-t)
    phi_t = -fp/R
    phi_r = fp/(R*FR)-f/(R*R)
    power = -4.0*math.pi*R*R*FR*phi_r*phi_t
    q_sphere = float(wt @ power)
    cross_integral = float(wt @ (-4.0*math.pi*FR/R*f*fp))

    # Horizon integral in its own future affine coordinate, not x.
    lower, upper = math.exp(-H*B)/H, math.exp(-H*A)/H
    lam, wl = nodes(lower, upper, n)
    v = np.log(H*lam)/H
    f, fp = profile(-v)
    phi_lam = -fp/lam
    Tkk = phi_lam**2
    area = 4.0*math.pi/(H*H)
    q_horizon = area*float(wl @ (H*lam*Tkk))
    modular_flux = 2.0*math.pi*area*float(wl @ (lam*Tkk))
    beta0 = 2.0*math.pi/H
    N = math.sqrt(FR)
    # Geometric normalization only; no detector or thermostat is simulated.
    beta_local, local_energy = beta0*N, q_sphere/N
    return {
        'nodes_per_integral': n,
        'initial_physical_energy_r_integral': q_initial,
        'reduced_energy_x_integral': q_wave,
        'finite_sphere_integrated_outward_energy': q_sphere,
        'finite_sphere_cross_term_integral': cross_integral,
        'horizon_killing_energy_lambda_integral': q_horizon,
        'horizon_modular_flux': modular_flux,
        'relative_entropy_from_mature_theorem': beta0*q_initial,
        'local_beta_times_local_integrated_energy': beta_local*local_energy,
        'max_energy_identity_residual': max(abs(q_initial-q_wave),
            abs(q_initial-q_sphere), abs(q_initial-q_horizon)),
        'entropy_flux_identity_residual': abs(modular_flux-beta0*q_initial),
    }


def calculate():
    rows = [integrals(n) for n in (96, 192, 384)]
    last = rows[-1]
    energy_lower = math.pi*math.exp(-2.0/3.0)
    cs_energy_lower = 16.0*energy_lower
    assert last['initial_physical_energy_r_integral'] > cs_energy_lower > 0
    assert last['max_energy_identity_residual'] < 1e-9
    assert last['entropy_flux_identity_residual'] < 1e-8
    qkeys = ('initial_physical_energy_r_integral', 'reduced_energy_x_integral',
        'finite_sphere_integrated_outward_energy',
        'horizon_killing_energy_lambda_integral')
    refinement = {key: abs(rows[-1][key]-rows[-2][key]) for key in qkeys}
    assert max(refinement.values()) < 1e-8
    return {
        'status': 'PASS_FINITE_SOURCE_TRANSPORT_DIAGNOSTIC_WITH_SCOPE',
        'parameters': {'H': H, 'fixed_minimal_mass_squared': M2,
            'xi': 0, 'pulse_support_x': [A,B], 'sphere_tortoise_radius': L,
            'sphere_areal_radius': math.tanh(H*L)/H},
        'integrals': rows,
        'analytic_energy_lower_bound': energy_lower,
        'stronger_cauchy_schwarz_energy_lower_bound': cs_energy_lower,
        'analytic_relative_entropy_lower_bound': 2.0*math.pi/H*energy_lower,
        'quadrature_refinement_differences': refinement,
        'finite_sphere_static_time_window': [L-B,L-A],
        'finite_sphere_proper_duration': (B-A)/math.cosh(H*L),
        'horizon_future_affine_window': [math.exp(-H*B)/H,math.exp(-H*A)/H],
        'new_numbered_rounds': 0,
        'scope': [
            'Normal compact initial-state input with all field energy and correlations retained.',
            'A new minimal fixed-mass action, not the previous conformal Hilbert source.',
            'Same-state source propagation to a finite sphere and finite affine horizon segment.',
            'Coherent excitation energy is not a certified thermal-contact heat flow.',
            'Finite affine horizon support is not finite central-clock arrival or readout.',
            'Exact identities follow from proof.md; numerical quadrature is diagnostic only.',
            'No coupled metric solution, independent area reading or backreaction remainder.',
            'No new Einstein selection, adopted axiom, or repeated formal check count.']
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write:
        assert not path.exists(), 'Do not overwrite saved results.'
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        status = 'SAVED'
    else:
        expected = json.loads(path.read_text(encoding='utf-8'))
        assert result == expected, 'Saved result differs from read-only recomputation.'
        status = 'READ_ONLY_MATCH'
    row = result['integrals'][-1]
    print(json.dumps({'status':status,
        'energy':row['initial_physical_energy_r_integral'],
        'relative_entropy':row['relative_entropy_from_mature_theorem'],
        'max_energy_residual':row['max_energy_identity_residual'],
        'affine_window':result['horizon_future_affine_window']}))


if __name__ == '__main__':
    main()
