"""Adjoint return amplitude, its complementary channel, mass-ratio and mixing tests."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PRIOR = runpy.run_path(str(HERE.parent / "1668/check.py"))
CFG = json.loads((HERE / "inputs.json").read_bytes())


def spectrum(a, b, modulus, phase, R):
    legs = np.array([1., 1., modulus])
    u, v = legs*a, legs*b
    aa, bb, ab = u@u, v@v, u@v
    determinant = aa*bb-ab*ab
    coefficient = determinant*(R + 1/R) - 2*np.cos(phase)*ab**2
    disc = coefficient**2 - 4*(aa*bb)**2
    roots = [] if disc < 0 or coefficient <= 0 else [
        2*aa**2/(coefficient + np.sqrt(disc)),
        (coefficient + np.sqrt(disc))/(2*bb**2)]
    for r in roots:
        mass = np.outer(u,u) + r*np.exp(1j*phase)*np.outer(v,v)
        eig = np.linalg.svd(mass, compute_uv=False)
        np.testing.assert_allclose(eig[0]/eig[1], R, rtol=1e-12)
    return [aa,bb,ab,determinant], disc, roots


def run():
    adj = np.diag(CFG["adjoint_vev"])
    phi = np.pi * CFG["phi_rad_over_pi"]
    omega = np.exp(4j*phi)
    w = np.diag(np.exp(2j*np.pi*np.diag(adj)/3))
    np.testing.assert_allclose(w, np.diag([omega,omega,omega,1,1]), atol=1e-14)
    np.testing.assert_allclose(np.linalg.det(w), 1, atol=1e-14)
    direction, identity = adj/np.sqrt(30), np.eye(5)/np.sqrt(5)
    moved = w @ direction
    z, complementary = np.vdot(direction,moved), np.vdot(identity,moved)
    np.testing.assert_allclose(z, (3+2*omega)/5, atol=1e-14)
    np.testing.assert_allclose(complementary, np.sqrt(6)/5*(omega-1), atol=1e-14)
    np.testing.assert_allclose(moved, z*direction + complementary*identity, atol=1e-14)
    np.testing.assert_allclose([abs(z)**2,abs(complementary)**2], [7/25,18/25])
    np.testing.assert_allclose(2-omega, -omega*(3+2*omega), atol=1e-14)
    couplings = np.array([1., np.exp(2j*phi), np.exp(3j*phi)*z])
    common_pulse = np.pi/(1+abs(z))
    output, residual = [], []
    for value in couplings:
        hamiltonian = np.array([[0j, value.conjugate()], [value, 0j]])
        eig, modes = np.linalg.eigh(hamiltonian)
        evolution = (modes*np.exp(-1j*common_pulse*eig)) @ modes.conj().T
        output.append(evolution[1, 0])
        residual.append(evolution[0, 0])
    output, residual = np.array(output), np.array(residual)
    np.testing.assert_allclose(output, -1j*np.sin(common_pulse)*couplings/abs(couplings))
    np.testing.assert_allclose(residual, np.cos(abs(z)*common_pulse)*np.array([-1.,-1.,1.]))
    closed = PRIOR["gap_box"](2)
    Rmax = closed[1]
    a,b = np.array(CFG["a"]), np.array(CFG["b"])
    phase = np.pi * CFG["eta_rad_over_pi"][1]
    gram, maximum_disc, roots = spectrum(a, b, abs(z), phase, Rmax)
    assert roots == [] and maximum_disc < -1
    for endpoint in closed:
        for sign in (-1,1):
            assert not spectrum(a,b,abs(z),sign*phase,endpoint)[2]
    aa, bb, ab, determinant = gram
    critical_sum = (2*aa*bb+2*np.cos(phase)*ab**2)/determinant
    critical_R = (critical_sum + np.sqrt(critical_sum**2-4))/2
    # Saturating r=aa/bb verifies this bound independently by direct SVD.
    u, v = np.array([1.,1.,abs(z)])*a, np.array([1.,1.,abs(z)])*b
    sat = np.linalg.svd(np.outer(u,u) + (aa/bb)*np.exp(1j*phase)*np.outer(v,v),
                        compute_uv=False)
    np.testing.assert_allclose(sat[0]/sat[1], critical_R)
    assert critical_R > Rmax
    # Necessary raw-leg modulus at the upper edge of the two-frequency rectangle.
    h = Rmax+1/Rmax-2
    xroots = np.sort(np.roots([1, 6-5*h, 9-h]))
    assert 0 < xroots[0] < 1 < xroots[1]
    min_modulus = np.sqrt(xroots[0])
    cfg = PRIOR["read"](PRIOR["INPUT"]["charged"])
    charged,_ = PRIOR["OLD"]["charged_basis"](cfg)
    angles = np.array(json.loads((HERE/CFG["angle_windows"]).read_bytes())["joint_windows"])
    n = np.array([2.,-1.,1/abs(z)])
    terms = abs(n * charged[:,0])/np.linalg.norm(n)
    zero_range = [max(0., 2*terms.max()-terms.sum())**2, terms.sum()**2]
    # For each reactor solution p, F=|U_e1|²=(1-p)(1-sin²theta12).
    lo, hi = angles[0]
    pl, ph = (1-np.sqrt(1-np.array([lo,hi])))/2
    slo, shi = angles[1]
    required_F = [[pl*(1-shi),ph*(1-slo)],[(1-ph)*(1-shi),(1-pl)*(1-slo)]]
    assert required_F[0][1] < zero_range[0] <= zero_range[1] < required_F[1][0]
    # The unit-phase realization has a physical witness, with both roots recomputed.
    ratio = 5.74
    _, _, norm_roots = spectrum(a,b,1.,phase,ratio)
    normalised = np.diag([1.,np.exp(2j*phi),np.exp(3j*phi)*z/abs(z)])
    mass = normalised @ (np.outer(a,a)+norm_roots[0]*np.exp(1j*phase)*np.outer(b,b)) @ normalised.T
    _, states = PRIOR["OLD"]["takagi"](mass)
    readouts = PRIOR["OLD"]["angles"](charged.conj().T@states)
    assert np.all((angles[:,0] < readouts) & (readouts < angles[:,1]))
    equivalent_mass = np.diag(output) @ (np.outer(a,a)+norm_roots[0]*np.exp(1j*phase)*np.outer(b,b)) @ np.diag(output)
    np.testing.assert_allclose(equivalent_mass, -np.sin(common_pulse)**2*mass)
    return {"return_complex":[z.real,z.imag], "return_modulus":abs(z),
            "return_phase_degrees":np.degrees(np.angle(z)),
            "singlet_prob":abs(complementary)**2, "return_prob":abs(z)**2,
            "raw_gram_aa_bb_ab_det":gram, "JUNO_two_width_R":closed.tolist(),
            "max_discriminant":maximum_disc, "positive_roots_anywhere":0,
            "critical_R":critical_R, "saturating_r":aa/bb,
            "frequency_required_modulus_min":min_modulus,
            "minimum_gain_over_raw_for_frequencies":min_modulus/abs(z),
            "zero_electron_all_phases":zero_range, "required_zero_electron_two_branches":required_F,
            "common_pulse":{"g_times_t":common_pulse, "each_success_probability":abs(output[0])**2,
                            "each_residual_probability":abs(residual[0])**2},
            "unit_phase_at_R_5_74":{"two_roots":norm_roots,"angles_small_root":readouts.tolist()}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = PRIOR["tidy"](run())
    path = HERE/"results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf8")
    else:
        assert json.loads(path.read_bytes()) == result
    print(json.dumps(result, indent=2))
