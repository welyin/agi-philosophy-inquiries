"""Positive local reference alignment: mode, stress and noisy-record interfaces.

This is an explicitly chosen quadratic potential. It does not select dimension,
derive gravity, prepare the aligned state, or establish a continuum quantum limit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import material_reference_geometry as bg

HERE = Path(__file__).resolve().parent
TARGET = HERE / "reference_alignment_completion_results.json"


def rotation(r):
    return np.array([[1., -r], [r, 1.]]) / np.sqrt(1+r*r)


def laplacian_min(n, d):
    return float(4*d*np.sin(np.pi/(2*(n+1)))**2)


def stress(gradients, inverse_metric):
    cov = np.linalg.inv(inverse_metric)
    outer = gradients.T @ gradients
    return outer - .5*cov*np.trace(inverse_metric @ outer)


def run():
    m, r, nu = 1., .25, .4
    lam, mass2 = r*m*m, m*m*(1+r*r)
    rot = rotation(r)
    mass = m*m*np.array([[r*r, r], [r, 1.]])
    assert np.max(abs(rot @ rot.T-np.eye(2))) < 1e-14
    assert np.max(abs(rot @ mass @ rot.T-np.diag([0., mass2]))) < 1e-14
    # Large-region obstruction, evaluated only in the exact lowest spatial mode.
    stability = []
    for n in (3, 127):
        ell = laplacian_min(n, 3)
        bare = np.array([[ell, lam], [lam, ell+m*m]])
        completed = ell*np.eye(2)+mass
        low_old, low_new = map(float, (np.linalg.eigvalsh(bare)[0],
                                      np.linalg.eigvalsh(completed)[0]))
        assert abs(low_new-ell) < 1e-13
        stability.append(dict(d=3, n=n, laplacian_min=ell,
                              bare_min_eigenvalue=low_old,
                              completed_min_eigenvalue=low_new))
    assert stability[0]["bare_min_eigenvalue"] > 0
    assert stability[1]["bare_min_eigenvalue"] < 0
    # Canonical kinetic, gradients, potential and total stress transform together.
    rng = np.random.default_rng(731)
    transformation_residual = 0.
    for _ in range(12):
        z, p = rng.normal(size=(2, 7)), rng.normal(size=(2, 7))
        normal, pn = rot@z, rot@p
        original = .5*np.sum(p*p)+.5*np.sum(z*(mass@z))
        transformed = .5*np.sum(pn*pn)+.5*mass2*np.sum(normal[1]**2)
        transformation_residual = max(transformation_residual, abs(original-transformed))
    assert transformation_residual < 1e-12
    # Same classical Einstein-reference solution, embedded by a constant rotation.
    geometric_rows = []
    for d in (2, 3, 4):
        a, q, b = 1.7, .6, .4
        inv_g = np.diag([-1.] + [a**-2]*d)
        refs = np.diag([q/a**d] + [b]*d)
        shared = [refs[0]]
        for i in range(1, d+1):
            shared.extend([refs[i]/np.sqrt(1+r*r),
                           -r*refs[i]/np.sqrt(1+r*r)])
        difference = float(np.max(abs(stress(np.array(shared), inv_g)-stress(refs, inv_g))))
        rho = q*q/(2*a**(2*d)) + d*b*b/(2*a*a)
        h = float(bg.adot(a,d,q,b)/a)
        friedmann = abs(d*(d-1)*h*h/2-rho)
        assert difference < 1e-14 and friedmann < 1e-14
        geometric_rows.append(dict(d=d, stress_residual=difference,
                                   friedmann_residual=friedmann))
    # Finite Dirichlet quantum patch: retain both normal-mode vacuum fluctuations.
    n = 5
    lap = 2*np.eye(n)-np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1)
    values, vec = np.linalg.eigh(lap)
    vr = (vec*(.5/np.sqrt(values)))@vec.T
    vm = (vec*(.5/np.sqrt(values+mass2)))@vec.T
    spatial_mass = np.kron(np.eye(2),lap)+np.kron(mass,np.eye(n))
    full_values, full_vec = np.linalg.eigh(spatial_mass)
    full_q = (full_vec*(.5/np.sqrt(full_values)))@full_vec.T
    inferred = (1+r*r)/r**2 * (full_q[n:,n:]+nu**2*np.eye(n))
    predicted = vr + vm/r**2 + (1+r*r)/r**2*nu**2*np.eye(n)
    covariance_error = float(np.max(abs(inferred-predicted)))
    assert covariance_error < 1e-12
    assert np.min(np.linalg.eigvalsh(predicted-vr)) > 0
    # A finite-resolution Gaussian read of psi adds momentum noise, shared by modes.
    noise = np.diag([0., 1/(4*nu*nu)])
    normal_noise = rot@noise@rot.T
    predicted_noise = np.array([[r*r, -r],[-r,1.]])/(4*nu*nu*(1+r*r))
    assert np.max(abs(normal_noise-predicted_noise)) < 1e-14
    assert abs(.5*np.trace(normal_noise)-1/(8*nu*nu)) < 1e-14
    # The aligned mean is a correlated source, not a blank independently prepared probe.
    x = np.arange(1,n+1,dtype=float)
    aligned = rot.T @ np.vstack([x,np.zeros(n)])
    assert np.max(abs(rot@aligned-np.vstack([x,np.zeros(n)]))) < 1e-14
    assert np.linalg.norm(aligned[1]) > 0
    # With fixed lambda, making the bare probe heavy amplifies calibrated massive noise.
    heavy_rows = []
    for bare_m in (1., 4.):
        rr=lam/(bare_m*bare_m)
        mm2=bare_m*bare_m*(1+rr*rr)
        variance=.5/(rr*rr*np.sqrt(1.+mm2))
        heavy_rows.append(dict(bare_mass=bare_m,fixed_lambda=lam,
                               calibrated_massive_mode_variance_at_laplacian_one=float(variance)))
    assert heavy_rows[1]["calibrated_massive_mode_variance_at_laplacian_one"] > heavy_rows[0]["calibrated_massive_mode_variance_at_laplacian_one"]
    return dict(date="2026-09-30",diagnostic_tests=6,failures=0,errors=0,
        numbered_round_created=False,numbered_test_increment=0,
        dependency_sha256={"material_reference_geometry.py":hashlib.sha256((HERE/"material_reference_geometry.py").read_bytes()).hexdigest()},
        inputs=dict(m=m,r=r,lam=lam,massive_mode_mass_squared=mass2,readout_sigma=nu),
        diagnostics=dict(stability=stability,canonical_energy_residual=float(transformation_residual),
            classical_geometry=geometric_rows,
            quantum_record=dict(covariance_identity_residual=covariance_error,
                reference_variance=float(vr[2,2]),massive_mode_calibration_noise=float(vm[2,2]/r**2),
                readout_calibration_noise=(1+r*r)/r**2*nu**2,
                total_calibrated_variance=float(predicted[2,2])),
            readout_recoil=dict(normal_momentum_covariance_increment=normal_noise.tolist(),
                kinetic_energy_increment_per_read=float(.5*np.trace(normal_noise))),
            source_and_heavy_limit=dict(aligned_probe_mean_norm=float(np.linalg.norm(aligned[1])),
                fixed_lambda_examples=heavy_rows)),
        scope=dict(positive_potential_is_added_model_choice=True,
            exact_classical_reference_geometry_embedding=True,
            same_blank_probe_protocol_inherited=False,
            finite_patch_quantum_record_covariance=True,
            quantum_vacuum_stress_solved_with_einstein=False,
            quantum_clock_control_included=False,
            spatial_dimension_selected=False,complete_unification=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    if args.check:
        assert json.loads(TARGET.read_text("utf8"))==json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
