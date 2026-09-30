"""Thermal/vacuum shared-reference difference records on a supplied lattice.

Only a conditional interface: geometry, thermal preparation, reference gradient,
zero-mode preparation and finite-noise final instruments are inputs.  No dimension
is selected uniquely.  Fourier sums diagnose the separately proved IR bounds.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "thermal_reference_dimension_results.json"
RATIO, MU2, NU = .5, 1.25, .25


def thermal_q(ell, temperature):
    """Position covariance for hbar=kB=mass=1; strictly positive ell only."""
    ell = np.asarray(ell, dtype=float)
    assert np.all(ell > 0) and temperature >= 0
    omega = np.sqrt(ell)
    if temperature == 0:
        return .5 / omega
    return .5 / (omega * np.tanh(omega / (2 * temperature)))


def fourier(n, d, temperature, separations):
    """Axis-separated pairs; O(n**(d-1)) memory, exact finite-mode summation."""
    one = 4 * np.sin(np.pi * np.arange(n) / n)**2
    tail = np.zeros(1)
    for _ in range(d-1):
        tail = (tail[:, None] + one[None, :]).ravel()
    qs_r, qs_m, energies_r, energies_m, resistances = [], [], [], [], []
    for value in one:
        ell = tail + value
        active = ell > 0
        cr = np.zeros_like(ell)
        cr[active] = thermal_q(ell[active], temperature)
        cm = thermal_q(ell + MU2, temperature)
        qs_r.append(float(np.sum(cr)))
        qs_m.append(float(np.sum(cm)))
        energies_r.append(float(np.sum(ell * cr)))
        energies_m.append(float(np.sum((ell + MU2) * cm)))
        resistances.append(float(np.sum(1 / ell[active])))
    phase = 2 * np.pi * np.arange(n) / n
    volume = n**d
    out = []
    for sep in separations:
        factor = 2 * (1 - np.cos(phase * sep)) / volume
        dr, dm = float(factor @ qs_r), float(factor @ qs_m)
        total = dr + dm / RATIO**2 + 2*(1+RATIO**2)*NU**2/RATIO**2
        assert dr >= 0 and dm >= 0
        assert dm <= 2*float(thermal_q(MU2, temperature)) + 1e-12
        out.append(dict(separation=sep, common_variance=dr,
                        relative_variance=dm/RATIO**2,
                        readout_variance=2*(1+RATIO**2)*NU**2/RATIO**2,
                        actual_difference_variance=total,
                        common_classical_variance=temperature*float(factor @ resistances)))
    return dict(n=n, d=d, temperature=temperature,
                nonzero_common_site_covariance=float(np.sum(qs_r)/volume),
                common_fluctuation_energy_density=float(np.sum(energies_r)/volume),
                relative_fluctuation_energy_density=float(np.sum(energies_m)/volume),
                pairs=out)


def dense_check(n, d, temperature, zero_variance):
    sites = list(itertools.product(range(n), repeat=d))
    index = {p: i for i, p in enumerate(sites)}
    size = len(sites)
    lap = 2*d*np.eye(size)
    for p, row in index.items():
        for axis in range(d):
            for sign in (-1, 1):
                q = list(p)
                q[axis] = (q[axis]+sign) % n
                lap[row, index[tuple(q)]] -= 1
    ell, u = np.linalg.eigh(lap)
    active = ell > 1e-10
    qspec = np.full(size, zero_variance)
    pspec = np.full(size, 1/(4*zero_variance))
    qspec[active] = thermal_q(ell[active], temperature)
    pspec[active] = ell[active]*qspec[active]
    qr, pr = (u*qspec)@u.T, (u*pspec)@u.T
    qm = (u*thermal_q(ell+MU2, temperature))@u.T
    pm = (u*((ell+MU2)*thermal_q(ell+MU2, temperature)))@u.T
    omega = np.block([[np.zeros_like(lap), np.eye(size)],
                      [-np.eye(size), np.zeros_like(lap)]])
    uncertainty = np.linalg.eigvalsh(np.block([[qr,np.zeros_like(lap)],
                                            [np.zeros_like(lap),pr]]) + .5j*omega)[0]
    assert uncertainty > -1e-12
    # Full phi,psi covariance is obtained by the actual canonical rotation.
    rot = np.array([[1., -RATIO], [RATIO, 1.]])/np.sqrt(1+RATIO**2)
    bigrot = np.kron(rot, np.eye(size))
    normal_q = np.block([[qr,np.zeros_like(lap)],[np.zeros_like(lap),qm]])
    original_q = bigrot.T @ normal_q @ bigrot
    record_q = (1+RATIO**2)/RATIO**2*(original_q[size:,size:]+NU**2*np.eye(size))
    contrast = np.zeros(size)
    contrast[0], contrast[index[(2,)+(0,)*(d-1)]] = 1., -1.
    direct = float(contrast @ record_q @ contrast)
    predicted = fourier(n,d,temperature,[2])["pairs"][0]["actual_difference_variance"]
    assert abs(direct-predicted) < 5e-12
    # Two local positive-noise instruments preserve Q and add P covariance.
    local_noise = np.diag([0., 1/(4*NU**2)])
    normal_noise = rot @ local_noise @ rot.T
    expected = np.array([[RATIO**2,-RATIO],[-RATIO,1.]])/(4*NU**2*(1+RATIO**2))
    assert np.max(abs(normal_noise-expected)) < 1e-13
    assert abs(np.trace(normal_noise)-1/(4*NU**2)) < 1e-13
    return dict(n=n,d=d,temperature=temperature,zero_mode_variance=zero_variance,
                actual_difference_variance=direct,formula_residual=abs(direct-predicted),
                common_uncertainty_minimum=float(uncertainty),
                two_read_field_energy_increment=float(np.trace(normal_noise)))


def run():
    # 1. Thermal oscillator covariance, independently from the occupation law.
    oscillator_error = 0.
    for freq in (.3, 1., 3.):
        for temp in (.1, 1.):
            prob = (1-math.exp(-freq/temp))*np.exp(-freq/temp*np.arange(512))
            direct = float(prob @ ((np.arange(512)+.5)/freq))
            oscillator_error = max(oscillator_error,abs(direct-float(thermal_q(freq*freq,temp))))
    assert oscillator_error < 1e-13
    # 2. Independent dense joint field/record covariance including physical zero mode.
    dense = [dense_check(5,d,temp,zero) for d in (1,2)
             for temp in (0.,1.) for zero in (.2,17.)]
    for row_a,row_b in zip(dense[::2],dense[1::2]):
        assert abs(row_a["actual_difference_variance"]-row_b["actual_difference_variance"]) < 5e-12
    # 3. Exact finite-cycle resistance identity and complete pair-average identity.
    cycle_error, average_error = 0., 0.
    for n in (8,17,32):
        stats=fourier(n,1,1.,range(n))
        for row in stats["pairs"]:
            sep=row["separation"]
            cycle_error=max(cycle_error,abs(row["common_classical_variance"]-sep*(n-sep)/n))
        average=float(np.mean([p["common_variance"] for p in stats["pairs"]]))
        average_error=max(average_error,abs(average-2*stats["nonzero_common_site_covariance"]))
    assert cycle_error < 2e-13 and average_error < 2e-13
    # 4. Finite-size diagnostics only; threshold is proved by spectral shells.
    sizes={1:(64,256,1024),2:(16,64,256),3:(8,32,64),4:(8,16,32)}
    rows=[fourier(n,d,temp,[1,n//2]) for d,ns in sizes.items()
          for temp in (0.,1.) for n in ns]
    for row in rows:
        d,temp=row["d"],row["temperature"]
        assert row["pairs"][0]["common_variance"] <= (temp+math.sqrt(d))/d + 1e-12
        assert row["common_fluctuation_energy_density"] <= temp+math.sqrt(d)+1e-12
        assert row["relative_fluctuation_energy_density"] <= temp+.5*math.sqrt(4*d+MU2)+1e-12
        g=row["nonzero_common_site_covariance"]
        assert row["pairs"][1]["common_variance"] <= 4*g+1e-12
    # 5. Positive readout noise rules out fixed tolerance for infinitely many pairs.
    noise_variance=2*(1+RATIO**2)*NU**2/RATIO**2
    tolerance=1.
    pair_cap=math.erf(tolerance/math.sqrt(2*noise_variance))
    simultaneous=[dict(disjoint_pairs=k,success_upper_bound=pair_cap**k) for k in (1,10,100)]
    assert 0 < pair_cap < 1 and simultaneous[-1]["success_upper_bound"] < 1e-8
    # 6. A low-dimensional relative-accuracy task survives the absolute-variance failure.
    low=[r for r in rows if r["d"]==1 and r["temperature"]==1.]
    relative=[dict(n=r["n"],separation=r["n"]//2,
                   mean_square_fractional_error=r["pairs"][1]["actual_difference_variance"]/(r["n"]//2)**2,
                   absolute_variance=r["pairs"][1]["actual_difference_variance"])
              for r in low]
    assert all(a["absolute_variance"] < b["absolute_variance"] for a,b in zip(relative,relative[1:]))
    assert all(a["mean_square_fractional_error"] > b["mean_square_fractional_error"] for a,b in zip(relative,relative[1:]))
    return dict(date="2026-09-30",diagnostic_tests=6,failures=0,errors=0,
        numbered_round_created=False,numbered_test_increment=0,
        inputs=dict(hbar=1.,boltzmann=1.,ratio=RATIO,massive_squared=MU2,readout_sigma=NU,
                    zero_mode="finite Gaussian kept; no massless full Gibbs state asserted",
                    geometry="unit nearest-neighbor periodic fluctuations; affine twist supplied for coordinate mean"),
        diagnostics=dict(oscillator_covariance_residual=oscillator_error,dense_joint_checks=dense,
             cycle_identity_residual=cycle_error,pair_average_residual=average_error,
             finite_size_rows=rows,infinite_joint_obstruction=dict(tolerance=tolerance,
                 meter_only_pair_success_cap=pair_cap,finite_disjoint_pair_bounds=simultaneous),
             relative_accuracy_counterexample=relative),
        provenance_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
             for name in ("reference_alignment_completion.py","reference_alignment_completion_results.json")},
        scope=dict(fixed_positive_temperature_threshold="d > 2",vacuum_threshold="d > 1",
             proof_not_numerical_fit=True,uniform_bound_is_per_pair=True,
             infinite_simultaneous_fixed_tolerance=False,spatial_dimension_three_selected=False,
             thermalization_derived=False,quantum_einstein_backreaction_solved=False,
             arbitrary_unknown_quantum_input_bound=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    if args.check:
        assert json.loads(TARGET.read_text("utf8")) == json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
