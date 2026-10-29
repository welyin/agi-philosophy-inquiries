"""1333: finite diagnostics for a uniform one-step Yukawa grammar.

The all-n statements are proved in proof.md. This script checks identities,
boundary cases and finite perturbations; it does not simulate a formation device.
No dependencies on earlier scientific implementations; NumPy only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 2e-9


def shift(n):
    return np.diag(np.ones(n-1), -1).astype(complex)


def norm(x):
    return float(np.linalg.norm(x, 2))


def cpair(z):
    return [float(np.real(z)), float(np.imag(z))]


def matrices(n, a, b, c, d):
    s = shift(n)
    yu = a*np.eye(n)+b*s
    yd = c*np.eye(n)+d*s.conj().T
    return yu, yd, yu@yu.conj().T, yd@yd.conj().T


def spectrum_checks():
    rows = []
    coeffs = [(0, 0), (1, 0), (0, 1), (1, 1), (1, 1j),
              (.001, 1), (1, .001), (.6+.8j, -1.2+.5j), (2, -3j)]
    for n in (3, 4, 5, 8, 16):
        s = shift(n)
        for a, b in coeffs:
            for reverse in (False, True):
                t = s.conj().T if reverse else s
                y = a*np.eye(n)+b*t
                h = y@y.conj().T
                endpoint = n-1 if reverse else 0
                p = np.zeros((n, n)); p[endpoint, endpoint] = 1
                A, B, r = abs(a)**2, abs(b)**2, abs(a*b)
                expected = (A+B)*np.eye(n)-B*p+np.conj(a)*b*t+a*np.conj(b)*t.conj().T
                idx = [0, 1] if reverse else [n-2, n-1]
                compression = np.linalg.eigvalsh(h[np.ix_(idx, idx)])
                sv = np.linalg.svd(y, compute_uv=False)
                phasefree = np.linalg.svd(abs(a)*np.eye(n)+abs(b)*t, compute_uv=False)
                errs = [norm(h-expected),
                        float(np.max(np.abs(compression-np.array([A+B-r, A+B+r])))),
                        float(np.max(np.abs(sv-phasefree)))]
                passed = (max(errs) < TOL and sv[1]**2+TOL >= A+B-r
                          and sv[0]**2 <= A+B+2*r+TOL and sv[1]+TOL >= .5*sv[0])
                rows.append(dict(n=n, a=cpair(a), b=cpair(b), reverse=reverse,
                    ratio=None if sv[0] == 0 else float(sv[1]/sv[0]),
                    identity_error=errs[0], compression_error=errs[1], phase_spectrum_error=errs[2],
                    homogeneous_bound_margin=float(sv[1]-.5*sv[0]), passed=bool(passed)))
    return rows


def boundary_checks():
    rows = []
    for a, b, expected in [(0, 1, 0.), (1, 1, (3-math.sqrt(5))/2), (1, 0, 1.)]:
        sv = np.linalg.svd(a*np.eye(2)+b*shift(2), compute_uv=False)
        ratio = float(sv[1]/sv[0])
        rows.append(dict(n=2, a=a, b=b, ratio=ratio, expected_ratio=expected,
                         violates_n_ge_3_bound=ratio < .5,
                         passed=bool(abs(ratio-expected) < TOL)))
    rows.append(dict(n=1, second_singular_value_defined=False,
                     zero_matrix_ratio_defined=False, passed=True))
    return rows


def perturbation_checks():
    rng = np.random.default_rng(1333)
    rows = []
    for n in (3, 5, 8):
        y = (.8+.3j)*np.eye(n)+(1.1-.2j)*shift(n)
        u, sv, vh = np.linalg.svd(y)
        base = sv[0]
        for epsilon in (0., .05, .2, .45):
            z = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
            z = epsilon*base*z/norm(z)
            stretched = np.maximum(sv-epsilon*base, 0.)
            stretched[0] = sv[0]+epsilon*base
            perturbations = [('complex_dense', z),
                             ('singular_spread', u@np.diag(stretched-sv)@vh)]
            for label, delta_y in perturbations:
                actual = norm(delta_y)/base
                target_sv = np.linalg.svd(y+delta_y, compute_uv=False)
                ratio = float(target_sv[1]/target_sv[0])
                lower = (.5-epsilon)/(1+epsilon)
                weyl = float(np.max(np.abs(target_sv-sv)))
                rows.append(dict(n=n, kind=label, epsilon=epsilon, actual_relative_change=actual,
                    ratio=ratio, guaranteed_lower_ratio=lower, weyl_max_shift=weyl,
                    passed=bool(actual <= epsilon+TOL and weyl <= norm(delta_y)+TOL
                                and ratio+TOL >= lower)))
        # Explicit target with top two singular values 1 and rho, without fitting the grammar.
        for rho in (0., .01, .1, .3):
            wanted = np.linspace(rho, 0., n-1)
            target = u@np.diag(np.r_[1., wanted])*base@vh
            delta = norm(target-y)/base
            necessary = (.5-rho)/(1+rho)
            target_sv = np.linalg.svd(target, compute_uv=False)
            ratio = float(target_sv[1]/target_sv[0])
            rows.append(dict(n=n, kind='explicit_small_ratio_target', rho=rho,
                ratio=ratio, actual_relative_change=delta, necessary_change=necessary,
                passed=bool(abs(ratio-rho) < TOL and delta+TOL >= necessary)))
    return rows


def invariant_checks():
    rows = []
    examples = [(1, 1, 2, 1j), (1.1, .7, 1.6, .9*np.exp(.83j)),
                (.4+.8j, .9-.3j, 1.1+.2j, -.5+.7j), (1, 1, 1, 1),
                (0, 1, 2, 1j), (1, 0, 2, 1j)]
    for n in (3, 4, 5, 8):
        s = shift(n)
        for a, b, c, d in examples:
            yu, yd, hu, hd = matrices(n, a, b, c, d)
            A, B, C, D = abs(a)**2, abs(b)**2, abs(c)**2, abs(d)**2
            z = np.conj(a)*b*np.conj(c)*d
            mixed = n*(A+B)*(C+D)-B*(C+D)-D*(A+B)+2*(n-1)*z.real
            comm = hu@hd-hd@hu
            u, v = np.conj(a)*b, c*np.conj(d)
            w = u*np.conj(v)-np.conj(u)*v
            predicted_comm = np.zeros((n, n), complex)
            predicted_comm[0, 0] = -w; predicted_comm[-1, -1] = w
            predicted_comm[1, 0] += B*v; predicted_comm[0, 1] -= B*np.conj(v)
            predicted_comm[-1, -2] += D*u; predicted_comm[-2, -1] -= D*np.conj(u)
            predicted_cube = 6j*z.imag*abs(b*d)**2*(abs(b*c)**2-abs(a*d)**2)
            actual_cube = np.trace(comm@comm@comm)
            trace_error = max(abs(np.trace(hu)-n*A-(n-1)*B),
                              abs(np.trace(hd)-n*C-(n-1)*D))
            determinant_error = max(abs(np.linalg.det(hu)-A**n)/(1+A**n),
                                    abs(np.linalg.det(hd)-C**n)/(1+C**n))
            mixed_error = abs(np.trace(hu@hd)-mixed)
            comm_error = norm(comm-predicted_comm)
            cube_error = abs(actual_cube-predicted_cube)
            # Common left rephasing and independent right rephasings: same physical invariants.
            q = np.diag(np.exp(.37j*np.arange(n)))
            yup = np.exp(.23j)*q@yu@q.conj().T
            ydp = np.exp(-.51j)*q@yd@q.conj().T
            hp, kp = yup@yup.conj().T, ydp@ydp.conj().T
            rephase_error = max(norm(hp-q@hu@q.conj().T), norm(kp-q@hd@q.conj().T),
                                abs(np.trace(hp@kp)-np.trace(hu@hd)))
            errors = [trace_error, determinant_error, mixed_error, comm_error, cube_error, rephase_error]
            rows.append(dict(n=n, coefficients=[cpair(x) for x in (a,b,c,d)],
                trace_error=float(trace_error), relative_determinant_error=float(determinant_error),
                mixed_trace_error=float(mixed_error), commutator_error=float(comm_error),
                CP_trace_cube=cpair(actual_cube), predicted_CP_trace_cube=cpair(predicted_cube),
                CP_formula_error=float(cube_error), rephasing_error=float(rephase_error),
                passed=bool(max(errors) < TOL)))
    # Reconstruct five coordinates from weak-basis invariants and check local rank numerically.
    def recover(params, n=3):
        A, B, C, D, phi = params
        _, _, hu, hd = matrices(n, math.sqrt(A), math.sqrt(B), math.sqrt(C),
                               math.sqrt(D)*np.exp(1j*phi))
        ra = float(np.linalg.det(hu).real)**(1/n)
        rc = float(np.linalg.det(hd).real)**(1/n)
        rb = (np.trace(hu).real-n*ra)/(n-1)
        rd = (np.trace(hd).real-n*rc)/(n-1)
        baseline = n*(ra+rb)*(rc+rd)-rb*(rc+rd)-rd*(ra+rb)
        cosphi = (np.trace(hu@hd).real-baseline)/(2*(n-1)*math.sqrt(ra*rb*rc*rd))
        return np.array([ra, rb, rc, rd, cosphi])
    point = np.array([1.21, .49, 2.56, .81, .83])
    recovered = recover(point)
    target = np.r_[point[:4], math.cos(point[4])]
    step = 1e-5
    jac = np.column_stack([(recover(point+step*np.eye(5)[i])-recover(point-step*np.eye(5)[i]))/(2*step)
                           for i in range(5)])
    expected_jac = np.diag([1, 1, 1, 1, -math.sin(point[4])])
    rank = int(np.linalg.matrix_rank(jac, tol=1e-6))
    return rows, dict(parameters=point.tolist(), recovered_invariants=recovered.tolist(),
        reconstruction_error=float(np.max(np.abs(recovered-target))), jacobian_rank=rank,
        jacobian_error=float(np.max(np.abs(jac-expected_jac))),
        strong_CP_and_fixed_source_phases_excluded=True,
        passed=bool(rank == 5 and np.max(np.abs(jac-expected_jac)) < 1e-7))


def physical_example():
    upper_up = (.619+.084)/(171.7-3.)
    upper_down = (.055+.016)/(2.89-.09)
    rho = .1
    return dict(source='https://arxiv.org/html/0712.1419', table='II', historical_year=2008,
        scheme='MSbar', scale='M_Z', direct_generation_identification=3,
        upper_ratio_rectangle=dict(charm_top=upper_up, strange_bottom=upper_down),
        rectangle_is_not_joint_confidence_region=True, broad_target_rho=rho,
        necessary_relative_change=(.5-rho)/(1+rho), is_latest_fit=False,
        extra_unknown_generations_excluded_by_data=False,
        passed=bool(max(upper_up, upper_down) < rho and math.isclose((.5-rho)/(1+rho), 4/11)))


def compute():
    spectra=spectrum_checks(); boundaries=boundary_checks(); perturbations=perturbation_checks()
    invariants, freedom=invariant_checks(); physical=physical_example()
    passed=all(row['passed'] for row in spectra+boundaries+perturbations+invariants+[freedom,physical])
    paths=[HERE/'check.py', HERE.parent/'_admission'/'round1333_admission.md']
    provenance={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in paths}
    return dict(round=1333, status='PASS' if passed else 'FAIL', provenance=provenance,
        spectrum_checks=spectra, boundary_checks=boundaries, perturbation_checks=perturbations,
        invariant_checks=invariants, five_freedom_check=freedom, historical_physical_example=physical,
        counts=dict(spectra=len(spectra), boundaries=len(boundaries), perturbations=len(perturbations),
                    invariant_cases=len(invariants), other=2),
        numerical_tolerance=TOL, all_n_proof_is_analytic=True, floating_interval_certification=False,
        quantum_formation_device_constructed=False, new_axioms=0, empirical_experiments=0,
        stage_complete=False, long_goal_complete=False, next_round_admitted=False)


def compare(a, b, path='root'):
    if isinstance(a, dict):
        assert set(a)==set(b), (path, 'keys')
        for k in a: compare(a[k], b[k], path+'.'+k)
    elif isinstance(a, list):
        assert len(a)==len(b), (path, 'length')
        for i, (x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a, float):
        assert math.isclose(a,b,rel_tol=1e-7,abs_tol=2e-9), (path,a,b)
    else:
        assert a==b, (path,a,b)


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write-new',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=compute()
    assert result['status']=='PASS', 'A finite diagnostic failed; no result is written.'
    output=HERE/'results.json'
    if args.write_new:
        with output.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
    else:
        compare(result,json.loads(output.read_text(encoding='utf-8')))
    print(json.dumps(dict(round=1333,status='PASS',counts=result['counts'],results_written=args.write_new)))


if __name__=='__main__':
    main()
