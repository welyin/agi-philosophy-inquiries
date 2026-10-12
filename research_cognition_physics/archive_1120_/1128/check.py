"""1128 diagnostics. No continuum or Gaussian-clock cutoff is certified here.

The fixed 91,125-site result follows from exact rational bounds plus proof.md.
Small matrix models only cross-check independent algebra and mechanisms.
Default: read and compare results.json. --write-results explicitly registers it.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import itertools
import json
from math import factorial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 2e-10


def norm(a):
    return float(np.linalg.norm(a))


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def evolution(h, t=1.0):
    e, q = np.linalg.eigh(h)
    return (q * np.exp(-1j * t * e)) @ q.conj().T


def annihilators(n):
    ops = []
    for j in range(n):
        a = np.zeros((2**n, 2**n), complex)
        for bits in range(2**n):
            if (bits >> j) & 1:
                a[bits ^ (1 << j), bits] = (-1) ** ((bits & ((1 << j) - 1)).bit_count())
        ops.append(a)
    return ops


def creation(v, ann):
    return sum(v[j] * ann[j].conj().T for j in range(len(ann)))


def lift(u):
    n = len(u)
    g = np.zeros((2**n, 2**n), complex)
    g[0, 0] = 1
    for k in range(1, n + 1):
        subsets = list(itertools.combinations(range(n), k))
        for rows in subsets:
            i = sum(1 << j for j in rows)
            for cols in subsets:
                j = sum(1 << q for q in cols)
                g[i, j] = np.linalg.det(u[np.ix_(rows, cols)])
    return g


def quadratic(h, ann):
    return sum(h[i, j] * ann[i].conj().T @ ann[j]
               for i in range(len(h)) for j in range(len(h)))


def budget_group():
    radius, side = 20, 45
    sites = (4*radius**3 + 6*radius**2 + 8*radius + 3)//3
    counted = sum(abs(x)+abs(y)+abs(z) <= radius
                  for x in range(-radius, radius+1)
                  for y in range(-radius, radius+1)
                  for z in range(-radius, radius+1))
    assert sites == counted == 11521
    tail = F(6**21, factorial(21))/(1-F(6, 22))
    pi_hi, sqrt2_hi = F(22, 7), F(99, 70)
    assert sqrt2_hi**2 > 2
    platform = (pi_hi + 2*sqrt2_hi)*tail
    sigma, mass = F(1, 10**11), 10**24
    sigma_hi = sigma*F(80001, 80000)
    assert sigma_hi**2 >= sigma**2 + 1/(F(4*mass**2)*sigma**2)
    lam_hi, lip = pi_hi/2, 512
    factor = lam_hi*8*sites*lip
    clock = factor*sigma_hi
    force = lam_hi*8*lip
    endpoint = factor*(sigma_hi+force/F(2*mass))
    total = platform+clock
    controller_energy = sites*(F(mass, 2)+1/(F(8*mass)*sigma**2))
    field_bound = 16*side**3*10
    interaction_bound = lam_hi*8*sites*16
    propagation_amplitude = F(1732, 1000)/100-F(1, 200)
    assert F(1732, 1000)**2 < 3
    propagation_signal = propagation_amplitude**2-2*clock/100
    assert total < F(1, 100)
    assert endpoint < F(741567, 10**9)
    assert propagation_signal > F(136, 10**6)
    # The omitted geometric tail is bounded analytically, not by truncation.
    finite_sum = sum((F(6**n, factorial(n)) for n in range(21, 100)), F(0))
    assert finite_sum < tail
    return dict(box_side=side, box_sites=side**3, platform_sites=sites,
                full_car_modes=16*side**3, target_half_diamond=0.01,
                tail_fraction=str(tail), tail_upper=float(tail),
                platform_error_upper=float(platform), clock_error_upper=float(clock),
                total_error_fraction=str(total), total_error_upper=float(total),
                endpoint_interaction_energy_upper=float(endpoint),
                per_clock_force_upper=float(force),
                initial_controller_energy_exact=str(controller_energy),
                initial_controller_energy=float(controller_energy),
                field_hamiltonian_norm_upper=field_bound,
                interaction_operator_norm_upper=float(interaction_bound),
                actual_early_spatial_signal_lower=float(propagation_signal),
                scope="Exact rational certificate for the fixed 3D instance; no giant Fock matrix.")


def wilson_group():
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1, -1]).astype(complex)
    alphas = [np.kron(sx, s) for s in (sx, sy, sz)]
    beta, gamma = np.kron(sz, np.eye(2)), np.kron(sy, np.eye(2))
    coords = list(itertools.product(range(-1, 2), repeat=3))
    index = {x: i for i, x in enumerate(coords)}
    n = len(coords)
    shifts = []
    for axis in range(3):
        shift = np.zeros((n, n), complex)
        for x, i in index.items():
            dest = list(x)
            dest[axis] += 1
            if tuple(dest) in index:
                shift[index[tuple(dest)], i] = 1
        shifts.append(shift)
    h = np.kron(np.eye(n), beta)
    for a, t in zip(alphas, shifts):
        h += np.kron((t-t.conj().T)/(2j), a)
        h += np.kron(np.eye(n)-(t+t.conj().T)/2, beta)
    gam = np.kron(np.eye(n), gamma)
    bet = np.kron(np.eye(n), beta)
    e, q = np.linalg.eigh(h)
    s = q[:, e < 0] @ q[:, e < 0].conj().T
    f = np.zeros(4*n, complex)
    f[4*index[(0, 0, 0)]:4*index[(0, 0, 0)]+4] = np.array([1, 0, 1j, 0])/np.sqrt(2)
    occupation = float(np.vdot(f, s @ f).real)
    gap = float(np.min(np.abs(e)))
    mass_form = (bet @ h+h @ bet)/2
    neighbors = np.repeat([sum(abs(v) for v in x) == 1 for x in coords], 4)
    hop_norm = norm((h @ f)[neighbors])
    early = evolution(h, 0.01) @ f
    early_norm = norm(early[neighbors])
    axis_masses = []
    for axis in range(3):
        mask = np.repeat([x[axis] != 0 for x in coords], 4)
        axis_masses.append(norm(early[mask])**2)
    # A genuinely 3D small-box spectral diagnostic with four drive slices.
    t_end = 0.05
    eta = np.diag([float(sum(abs(v) for v in x) <= 1) for x in coords])
    hfull = np.kron(h, np.eye(2))
    coupling = np.kron(np.kron(eta, np.eye(4)), sx)
    local = np.eye(8*n, dtype=complex)
    for fraction in (0.15, 0.35, 0.35, 0.15):
        local = evolution(t_end*hfull/4+(np.pi/2)*fraction*coupling) @ local
    source = np.kron(f, [1, 0])
    free = evolution(h, t_end) @ f
    global_target = -1j*np.kron(free, [0, 1])
    plate_error = norm(local @ source-global_target)
    drive_global = np.kron(evolution(h, t_end), evolution((np.pi/2)*sx))
    global_error = norm(drive_global @ source-global_target)
    small_tail = (6*t_end)**2/2/(1-6*t_end/3)
    assert norm(h-h.conj().T) < TOL and norm(gam @ h+h @ gam) < TOL
    assert gap >= 1-TOL and np.linalg.eigvalsh(mass_form)[0] >= 1-TOL
    assert abs(occupation-0.5) < TOL and abs(hop_norm-np.sqrt(3)) < TOL
    assert early_norm > 0.01232 and min(axis_masses) > 0
    assert global_error < TOL and plate_error < (np.pi/2)*small_tail
    return dict(scope="27-site open 3D Wilson diagnostic, not the 91,125-site simulation.",
                sites=n, spinor_dimension=4*n, positive_levels=int(sum(e > 0)),
                negative_levels=int(sum(e < 0)), gap=gap,
                hermiticity_error=norm(h-h.conj().T),
                chiral_anticommutator_error=norm(gam @ h+h @ gam),
                source_negative_weight=occupation, source_empty_probability=(1-occupation)**2,
                six_neighbor_hf_norm=hop_norm, early_neighbor_amplitude=early_norm,
                early_axis_probabilities=axis_masses,
                global_rotation_error=global_error,
                local_platform_mode_error=plate_error,
                local_platform_mode_error_bound=(np.pi/2)*small_tail)


def car_rotation_diagnostic():
    ann = annihilators(6)
    rng = np.random.default_rng(1128)
    q, _ = np.linalg.qr(rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6)))
    s = q[:, :3] @ q[:, :3].conj().T
    vacuum = np.zeros(64, complex)
    vacuum[0] = 1
    for v in q[:, :3].T[::-1]:
        vacuum = creation(v, ann) @ vacuum
    source_indices = (0, 2)
    mask = np.array([float(not(bits & 1) and not(bits & 4)) for bits in range(64)])
    raw = mask*vacuum
    probability = norm(raw)**2
    base = raw/np.sqrt(probability)
    code = np.column_stack([ann[i].conj().T @ base for i in source_indices])
    u = np.eye(6, dtype=complex)
    single_errors = []
    distances = []
    for pair, theta, phase in (((0, 1), 0.12, 0.08), ((2, 3), 0.17, -0.06)):
        a = np.cos(theta)*np.exp(1j*phase)
        b = np.sin(theta)
        block = np.array([[a, -b], [b, a.conjugate()]])
        one = np.eye(6, dtype=complex)
        one[np.ix_(pair, pair)] = block
        distance = norm(block[:, 0]-np.array([1, 0]))
        single_errors.append(abs(opnorm(lift(one)-np.eye(64))-distance))
        distances.append(distance)
        u = one @ u
    # Pure phase (b=0) must use its conjugate companion mode.
    pure = np.eye(6, dtype=complex)
    pure[0, 0], pure[1, 1] = np.exp(0.19j), np.exp(-0.19j)
    phase_error = abs(opnorm(lift(pure)-np.eye(64))-abs(np.exp(0.19j)-1))
    lifted = lift(u)
    out = lifted @ code
    units_error = 0.0
    for i, source_i in enumerate(source_indices):
        for j, source_j in enumerate(source_indices):
            ci = creation(u[:, source_i], ann)
            cj = creation(u[:, source_j], ann)
            expected = np.zeros((2, 2), complex)
            expected[i, j] = 1
            units_error = max(units_error, norm(out.conj().T @ ci @ cj.conj().T @ out-expected))
    rotated_distance = opnorm(lifted-np.eye(64))
    f = np.eye(6)[:, source_indices]
    probability_error = abs(probability-np.linalg.det(np.eye(2)-f.conj().T @ s @ f).real)
    correction_error = norm(lifted.conj().T @ out-code)
    assert max(single_errors+[phase_error, units_error, probability_error, correction_error]) < TOL
    assert rotated_distance <= sum(distances)+TOL
    assert norm(s[np.ix_(source_indices, (1, 3, 4, 5))]) > 0.1
    return dict(modes=6, full_fock_dimension=64,
                complex_overlap_single_rotation_error=max(single_errors),
                pure_phase_pair_error=phase_error,
                two_rail_full_fock_distance=rotated_distance,
                two_rail_distance_upper=sum(distances),
                all_code_matrix_units_error=units_error,
                conditional_empty_probability=probability,
                determinant_probability_error=probability_error,
                full_output_correction_error=correction_error)


def clock_reaction_diagnostic():
    """Two small clock matrices: exact finite algebra, NOT Gaussian/CCR proof."""
    ann = annihilators(4)
    number = [a.conj().T @ a for a in ann]
    hsite = np.array([[0.4, 0.6], [0.6, -0.4]])
    hone = np.kron(hsite, np.eye(2))
    eigen = np.linalg.eigvalsh(hone)
    offset = -float(sum(eigen[eigen < 0]))
    bond = 0.6*sum(ann[i].conj().T @ ann[i+2]+ann[i+2].conj().T @ ann[i] for i in (0, 1))
    d0 = 0.4*(number[0]+number[1])+offset*np.eye(16)/2
    d1 = -0.4*(number[2]+number[3])+offset*np.eye(16)/2
    hf = d0+d1+bond
    currents = [ann[i].conj().T @ ann[i+1]+ann[i+1].conj().T @ ann[i] for i in (0, 2)]
    nc, dx, mc, g, end = 5, 0.75, 1.3, 0.8, 0.35
    positions = dx*(np.arange(nc)-(nc-1)/2)
    momenta = 2*np.pi*np.fft.fftfreq(nc, d=dx)
    fourier = np.exp(1j*np.outer(positions, momenta))/np.sqrt(nc)
    p = (fourier*momenta) @ fourier.conj().T
    x = np.diag(positions)
    hc = p @ p/(2*mc)
    gate = np.diag(np.exp(-(positions-0.1)**2/(2*0.45**2)))
    i5, i16, i25 = np.eye(nc), np.eye(16), np.eye(nc**2)
    control0 = np.kron(i16, np.kron(hc, i5))
    control1 = np.kron(i16, np.kron(i5, hc))
    coupling0 = g*np.kron(currents[0], np.kron(gate, i5))
    coupling1 = g*np.kron(currents[1], np.kron(i5, gate))
    field = np.kron(hf, i25)
    control = control0+control1
    interaction = coupling0+coupling1
    total = field+control+interaction
    e0 = np.kron(d0+bond/2, i25)+control0+coupling0
    e1 = np.kron(d1+bond/2, i25)+control1+coupling1
    initial_field = np.zeros(16, complex)
    initial_field[1] = initial_field[2] = 1/np.sqrt(2)
    phi0 = np.exp(-(positions+0.45)**2/(4*0.35**2)+0.7j*positions)
    phi1 = np.exp(-(positions-0.35)**2/(4*0.40**2)-0.4j*positions)
    phi0 /= norm(phi0)
    phi1 /= norm(phi1)
    initial = np.kron(initial_field, np.kron(phi0, phi1))
    final = evolution(total, end) @ initial
    mean = lambda a, state: float(np.vdot(state, a @ state).real)
    delta_field = mean(field, final)-mean(field, initial)
    delta_control = mean(control, final)-mean(control, initial)
    delta_interaction = mean(interaction, final)-mean(interaction, initial)
    p0 = np.kron(i16, np.kron(p, i5))
    delta_pointer = mean(p0, final)-mean(p0, initial)
    force = 1j*(total @ p0-p0 @ total)
    local_force = 1j*(coupling0 @ p0-p0 @ coupling0)
    flow = 1j*(e1 @ e0-e0 @ e1)
    derivative0 = 1j*(total @ e0-e0 @ total)
    derivative1 = 1j*(total @ e1-e1 @ total)
    balance_error = abs(delta_field+delta_control+delta_interaction)
    energy_error = abs(mean(total, final)-mean(total, initial))
    local_flow_error = max(norm(derivative0-flow), norm(derivative1+flow))
    force_error = norm(force-local_force)
    ccr_defect = norm(x @ p-p @ x-1j*np.eye(nc))
    assert max(balance_error, energy_error, local_flow_error, force_error) < TOL
    assert abs(delta_pointer) > 1e-5 and abs(delta_interaction) > 1e-5
    assert ccr_defect > 1
    return dict(scope="Two 5-dimensional clock matrices and full 4-mode CAR; not canonical Gaussian clocks.",
                total_dimension=len(total), duration=end,
                field_energy_change=delta_field, controller_energy_change=delta_control,
                interaction_energy_change=delta_interaction,
                total_energy_conservation_error=energy_error,
                sum_of_sector_energy_changes_error=balance_error,
                deleting_interaction_balance_error=abs(delta_field+delta_control),
                actual_pointer_momentum_change=delta_pointer,
                exact_matrix_force_error=force_error,
                local_energy_flow_error=local_flow_error,
                finite_clock_ccr_defect=ccr_defect)


def compare(actual, expected, path="results"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for k in actual:
            compare(actual[k], expected[k], path+"."+k)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, path+f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isclose(actual, expected, atol=1e-10, rtol=1e-10), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = dict(round=1128, schema=1, status="PASS",
                  budget=budget_group(), wilson=wilson_group(),
                  car_and_clock=dict(car=car_rotation_diagnostic(), clock=clock_reaction_diagnostic()))
    path = HERE/"results.json"
    if args.write_results:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print("PASS: three groups; results.json registered explicitly.")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        print("PASS: three groups; recorded results match; no files written.")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
