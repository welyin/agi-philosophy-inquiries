"""Count0 phase-convention checks; no apparatus simulation or empirical fit."""
from fractions import Fraction as F
from pathlib import Path
import cmath
import json
import math
import sys


def joint(phi_g, phi_e, analysis_phase):
    return [[(1 + q * math.cos(phi + analysis_phase)) / 4
             for q in (1, -1)] for phi in (phi_g, phi_e)]


def reconstructed_unit_phase(phi_g, phi_e, state):
    p0 = joint(phi_g, phi_e, 0)[state]
    p90 = joint(phi_g, phi_e, math.pi / 2)[state]
    return complex(2*(p0[0]-p0[1]), -2*(p90[0]-p90[1]))


def calculate():
    # These rational phase data are an algebra diagnostic, not an atom design.
    phi_g, delta_base, rate = F(1, 5), F(3, 7), F(1, 16)
    ti, tip = F(1, 3), F(5, 6)
    exact_double = rate*(tip-ti)
    assert exact_double == F(1, 32)
    theta = F(2, 9)
    prob_residuals, normalization, coarse_residuals = [], [], []
    reconstructed, probability_tables = [], []
    for time in (ti, tip):
        pg, pe = float(phi_g), float(phi_g+delta_base+rate*time)
        tables = []
        for analysis_phase in (0, math.pi/2):
            probs = joint(pg, pe, analysis_phase)
            assert min(value for row in probs for value in row) >= 0
            normalization.append(abs(sum(map(sum, probs))-1))
            # The four amplitudes span the same complete (internal, exit) labels.
            amps = [[cmath.exp(1j*(0 if s == 0 else float(theta))) *
                     (1+q*cmath.exp(1j*(phi+analysis_phase))) /
                     (2*math.sqrt(2)) for q in (1, -1)]
                    for s, phi in enumerate((pg, pe))]
            for s in range(2):
                for q in range(2):
                    prob_residuals.append(abs(abs(amps[s][q])**2-probs[s][q]))
            # Removing all g/e cross terms cannot change these diagonal effects.
            flat = amps[0]+amps[1]
            rho = [[a*b.conjugate() for b in flat] for a in flat]
            dephased = [[entry if i//2 == j//2 else 0j
                         for j, entry in enumerate(row)]
                        for i, row in enumerate(rho)]
            assert all(rho[i][i] == dephased[i][i] for i in range(4))
            coarse = probs[0][0]+probs[1][0]
            half_angle = .5 + .5*math.cos((pe-pg)/2)*math.cos((pe+pg)/2+analysis_phase)
            coarse_residuals.append(abs(coarse-half_angle))
            tables.append(probs)
        probability_tables.append(tables)
        ze = reconstructed_unit_phase(pg, pe, 1)
        zg = reconstructed_unit_phase(pg, pe, 0)
        reconstructed.append(ze*zg.conjugate())
    double_phase = cmath.phase(reconstructed[1]*reconstructed[0].conjugate())
    double_residual = abs(double_phase-float(exact_double))
    wrong_half_residual = abs(double_phase-float(exact_double)/2)
    assert double_residual < 2e-15
    assert wrong_half_residual > .01
    assert max(prob_residuals+normalization+coarse_residuals) < 2e-15

    # Eq.(48) retained-order identity: common position offset cancels first.
    # k*dv*T from laser phases equals m*vrec*dv*T/hbar from kinetic phases.
    mass, hbar, wave_number, duration = F(7, 3), F(2, 5), F(11, 7), F(5, 9)
    dv, offset = F(-1, 13), F(2, 17)
    recoil = hbar*wave_number/mass
    laser_change = wave_number*((offset+dv*duration)-offset)
    kinetic_change = -mass*((recoil+dv)**2-dv**2-recoil**2)*duration/(2*hbar)
    assert laser_change+kinetic_change == 0

    # Eq.(4) units, dimension vector order: length, mass, time.
    energy, hbar_dim, light, grav, length, time_dim = (
        (2, 1, -2), (2, 1, -1), (1, 0, -1),
        (1, 0, -2), (1, 0, 0), (0, 0, 1))
    dimension = tuple(e-h-2*c+g+z+t for e,h,c,g,z,t in
                      zip(energy,hbar_dim,light,grav,length,time_dim))
    assert dimension == (0, 0, 0)
    # Public proposal-scale illustration; no measured g or uncertainty budget.
    speed = 299792458
    wavelength, gravity, separation, interval = F(698, 10**9), F(981, 100), F(1, 100), F(1)
    phase_over_two_pi = gravity*separation*interval/(wavelength*speed)
    phase = 2*math.pi*float(phase_over_two_pi)
    return {
        "scope": "count0 fixed algebra and proposal scale; no measured data, full Q/R instrument, or finite-pulse error certificate",
        "fixed_phase_inputs": {"phi_g": str(phi_g), "delta_base": str(delta_base),
                               "rate": str(rate), "initialization_times": [str(ti), str(tip)]},
        "probability_order": "time, analysis phase (0, pi/2), internal (g,e), port (+,-)",
        "joint_probability_tables": probability_tables,
        "max_amplitude_probability_residual": max(prob_residuals),
        "max_normalization_residual": max(normalization),
        "max_coarse_half_angle_residual": max(coarse_residuals),
        "double_phase_exact": str(exact_double),
        "double_phase_from_joint_probabilities": double_phase,
        "double_phase_residual": double_residual,
        "wrong_half_coefficient_residual": wrong_half_residual,
        "energy_basis_dephasing_leaves_all_four_probabilities_unchanged": True,
        "retained_recoil_identity": {"laser_increment": str(laser_change),
                                     "kinetic_increment": str(kinetic_change), "sum": "0"},
        "phase_dimensions_length_mass_time": list(dimension),
        "proposal_scale_only": {"wavelength_m": str(wavelength), "g_m_s2": str(gravity),
                                "height_m": str(separation), "interval_s": str(interval),
                                "c_m_s": speed, "D_over_2pi_exact": str(phase_over_two_pi),
                                "D_radian": phase, "D_milliradian": 1000*phase},
        "all_assertions_passed": True,
    }


if __name__ == "__main__":
    result = calculate()
    target = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("No arguments: read-only check; --save-exclusive: create once.")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert saved == result
    print("PASS: four outcomes, full/half phase distinction, retained recoil identity and proposal units.")
