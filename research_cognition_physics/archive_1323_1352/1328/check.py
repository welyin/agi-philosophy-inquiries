"""1328: original W/J parent, all-dictionary spectral gap and prepared-reference process.

NumPy only. Default/--check is read-only; --write-new creates results.json exclusively.
The main finite experiment is built in the original D/U/L parent, not in the target.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def norm(a):
    return float(np.linalg.norm(a, 2))


def evolve(h, t):
    w, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * w)) @ v.conj().T


def kernel(m, p):
    n = len(m)
    return np.block([[p * np.eye(n), m.conj().T],
                     [m, -p * np.eye(n)]])


def parent(N, md, mu, p=0.0):
    n = len(md)
    w = np.zeros((2 * N, N + 1), complex)
    j = np.zeros((2 * N, N - 1), complex)
    # Boson occupation basis uses h = number in mode 1.
    for h in range(N + 1):
        if h < N:
            w[h, h] = math.sqrt((N - h) / N)
        if h:
            w[N + h - 1, h] = math.sqrt(h / N)
    for h in range(N - 1):
        j[h + 1, h] = math.sqrt((h + 1) / N)
        j[N + h, h] = -math.sqrt((N - 1 - h) / N)
    c = np.column_stack([w, j])
    t = np.column_stack([np.kron(w, md), np.kron(j, mu)])
    half = 2 * N * n
    v = np.block([[np.zeros((half, half)), t.conj().T],
                  [t, np.zeros((half, half))]])
    h0 = np.diag(np.r_[np.full(half, p), np.full(half, -p)])
    # S maps the fixed CG dictionary back into the original physical parent.
    s = np.zeros_like(v)
    s[:half, :half] = np.kron(c.conj().T, np.eye(n))
    s[half:, half:] = np.eye(half)
    eside = np.kron(np.kron(np.eye(2), np.eye(N)[:, :1]), np.eye(n))
    edict = np.zeros((2 * half, 4 * n), complex)
    edict[:half, :2 * n] = eside
    edict[half:, 2 * n:] = eside
    e = s @ edict
    md_mu = np.block([[md, np.zeros_like(md)], [np.zeros_like(mu), mu]])
    target = kernel(md_mu, p)
    nd = np.diag(np.r_[np.ones((N + 1) * n), np.zeros((3 * N - 1) * n)])
    nu = np.diag(np.r_[np.zeros((N + 1) * n), np.ones((N - 1) * n),
                       np.zeros(half)])
    nl = np.diag(np.r_[np.zeros(half), np.ones(half)])
    nb = N * nd + (N - 2) * nu + (N - 1) * nl
    return dict(W=w, J=j, C=c, S=s, E=e, H=v+h0, V=v, H0=h0,
                target=target, nd=nd, nu=nu, nl=nl, nb=nb, n=n, N=N)


def spectral_cases():
    out = []
    d, u = 1.3, .4
    for N in [2, 3, 8, 17]:
        a = parent(N, np.array([[d]]), np.array([[u]]))
        eig, vec = np.linalg.eigh(a['H'])
        groups = eig.reshape(4, N)
        centers = np.array([(r[0] + r[-1]) / 2 for r in groups])
        optimal = (vec * np.repeat(centers, N)) @ vec.conj().T
        fixed = (vec * np.repeat([-d, -u, u, d], N)) @ vec.conj().T
        row = dict(N=N, parent_dim=len(eig),
                   best_arbitrary_h=norm(a['H']-optimal),
                   expected_best=(d-u)/2,
                   best_fixed_spectrum=norm(a['H']-fixed),
                   expected_fixed=d-u,
                   projector_residual=norm(a['C'].conj().T@a['C']-np.eye(2*N)))
        row['passed'] = (abs(row['best_arbitrary_h']-(d-u)/2)<1e-11 and
                         abs(row['best_fixed_spectrum']-(d-u))<1e-11 and
                         row['projector_residual']<1e-11)
        out.append(row)
    return out


def scalar_time_cases():
    out = []
    d, u = 1.3, .4
    for N in [2, 3, 8, 17]:
        a = parent(N, np.array([[d]]), np.array([[u]]))
        e, h, target = a['E'], a['H'], a['target']
        sd = parent(N, np.ones((1, 1)), np.zeros((1, 1)))['V']
        su = parent(N, np.zeros((1, 1)), np.ones((1, 1)))['V']
        td = kernel(np.diag([1., 0.]), 0.)
        tu = kernel(np.diag([0., 1.]), 0.)
        initial = e[:, 3]  # Natural physical left, lower direction.
        rng = np.random.default_rng(1328)
        logical = rng.normal(size=4)+1j*rng.normal(size=4)
        logical /= np.linalg.norm(logical)
        source_bias = float(np.real(logical.conj() @ tu @ logical))/N
        for time in [.7, 2.4, 6.7]:
            uh, ut = evolve(h, time), evolve(target, time)
            actual = uh @ initial
            pd = float(np.real(actual.conj()@a['nd']@actual))
            pu = float(np.real(actual.conj()@a['nu']@actual))
            flow = float(np.real(actual.conj()@a['nb']@actual)-(N-1))
            raw = uh @ e @ logical
            effective = ut @ logical
            delta_d = float(np.real(raw.conj()@sd@raw-effective.conj()@td@effective))
            delta_u = float(np.real(raw.conj()@su@raw-effective.conj()@tu@effective))
            expected_error = 2*abs(math.sin((d-u)*time/2))/math.sqrt(N)
            row = dict(N=N, time=time, complete_error=norm(uh@e-e@ut),
                       predicted_error=expected_error, probability_D=pd, probability_U=pu,
                       reference_transfer=flow, source_d_bias=delta_d, source_u_bias=delta_u)
            residuals = [abs(row['complete_error']-expected_error),
                         abs(pd-math.sin(d*time)**2/N),
                         abs(pu-(N-1)*math.sin(u*time)**2/N),
                         abs(flow-pd+pu), abs(delta_d-source_bias), abs(delta_u+source_bias)]
            row['max_formula_residual'] = max(residuals)
            row['passed'] = max(residuals)<1e-10
            out.append(row)
    return out


def matrix_cases():
    md0 = np.array([[.9, .15+.07j], [-.03j, .55]])
    mu0 = np.array([[.4, .08j], [.12, .65]])
    n, N0, p0 = 2, 3, .27
    calibration = parent(N0, md0, mu0, p0)
    compression = calibration['E'].conj().T @ calibration['H'] @ calibration['E']
    forward = compression[2*n:, :2*n]
    md = forward[:n, :n]
    mu = (N0*forward[n:, n:]-md)/(N0-1)
    calibration_res = max(norm(md-md0), norm(mu-mu0))
    out = []
    for N in [3, 8, 17]:
        for p in [0., p0]:
            a = parent(N, md, mu, p)
            e, h, target = a['E'], a['H'], a['target']
            half = len(h)//2
            nl, nu, nd, nb = (a[k] for k in ['nl', 'nu', 'nd', 'nb'])
            # Original independent source matrices are differentiated before work conditions.
            perturb = np.zeros((n, n), complex); perturb[0, 1] = 1j
            source = parent(N, perturb, np.zeros((n, n)))['V']
            target_source = kernel(np.block([[perturb, np.zeros((n,n))],
                                            [np.zeros((n,n)), np.zeros((n,n))]]), 0.)
            source_res = norm(source@e-e@target_source)
            psource = np.diag(np.r_[np.ones(half), -np.ones(half)])
            ptarget = np.diag(np.r_[np.ones(2*n), -np.ones(2*n)])
            k_res = norm(nb+nl+2*nu-N*np.eye(len(h)))
            # Real preparation of every encoded right input from the natural left family.
            equal = parent(N, .7*np.eye(n), .7*np.eye(n), 0.)
            el, er = e[:, 2*n:], e[:, :2*n]
            prep = evolve(equal['H'], math.pi/(2*.7)) @ el
            prep_res = norm(prep+1j*er)
            for time in [.7, 2.4, 6.7]:
                uh, ut = evolve(h, time), evolve(target, time)
                hd, hu = kernel(md,p), kernel(mu,p)
                expected = norm(evolve(hd,time)-evolve(hu,time))/math.sqrt(N)
                error = norm(uh@e-e@ut)
                # Both physically qualified input families; the full 4n identity is conditional.
                families = []
                for name, initial in [('natural_left', el), ('prepared_right', 1j*prep)]:
                    final = uh@initial
                    m = initial.shape[1]
                    # Turn off V first. Then store all three results in distinct true records.
                    recorded = np.vstack([nd@final, nu@final, nl@final])
                    gram_res = norm(recorded.conj().T@recorded-np.eye(m))
                    partner_res = norm(recorded.conj().T@recorded/m-np.eye(m)/m)
                    record_energy = sum((proj@final).conj().T@a['H0']@(proj@final)
                                        for proj in [nd,nu,nl])
                    energy_res = norm(final.conj().T@h@final-initial.conj().T@h@initial)
                    record_energy_res = norm(record_energy-final.conj().T@a['H0']@final)
                    # Woff = -<V> is included as an operator identity on every unknown input.
                    off_work = -final.conj().T@a['V']@final
                    work_res = norm(record_energy-initial.conj().T@h@initial-off_work)
                    flow_res = norm(final.conj().T@(nb+nl+2*nu)@final-
                                    initial.conj().T@(nb+nl+2*nu)@initial)
                    rng = np.random.default_rng(804)
                    probe = rng.normal(size=m)+1j*rng.normal(size=m)
                    probe /= np.linalg.norm(probe)
                    psi = final@probe
                    probs = [float(np.real(psi.conj()@proj@psi)) for proj in [nd,nu,nl]]
                    transfer = float(np.real(psi.conj()@nb@psi-
                                             (initial@probe).conj()@nb@(initial@probe)))
                    families.append(dict(family=name, probabilities_D_U_L=probs,
                                         reference_transfer=transfer,
                                         max_record_resource_residual=max(gram_res,partner_res,energy_res,
                                             record_energy_res,work_res,flow_res,abs(sum(probs)-1))))
                residuals = [abs(error-expected), abs(source_res-1/math.sqrt(N)),
                             norm(psource@e-e@ptarget), k_res, prep_res,
                             abs(norm(h@e-e@target)-norm(md-mu)/math.sqrt(N))]
                residuals += [r['max_record_resource_residual'] for r in families]
                out.append(dict(N=N,p=p,time=time,complete_error=error,predicted_error=expected,
                                maximum_all_time_bound=2/math.sqrt(N),
                                source_intertwining_residual=source_res,
                                equal_coupling_preparation_residual=prep_res,
                                qualified_families=families,max_identity_residual=max(residuals),
                                passed=max(residuals)<1e-10))
    return dict(calibration_N=N0,calibration_error=calibration_res,
                p_is_registered_level_dictionary=True,parameter_refits_after_calibration=0,
                cases=out,passed=calibration_res<1e-12 and all(r['passed'] for r in out))


def compute():
    paths = [HERE/'proof.md', HERE/'check.py',
             HERE.parent/'_admission'/'round1328_admission.md',
             HERE.parent/'1325'/'proof.md', HERE.parent/'1327'/'proof.md']
    provenance = {str(p.relative_to(HERE.parent)).replace('\\','/'):
                  hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    spectrum = spectral_cases()
    scalar = scalar_time_cases()
    matrix = matrix_cases()
    passed = all(r['passed'] for r in spectrum+scalar) and matrix['passed']
    return dict(round=1328,status='PASS' if passed else 'FAIL',new_axioms=0,
                physical_scope='natural left plus actual equal-coupling prepared right; full-role code conditional',
                provenance=provenance,spectral_cases=spectrum,scalar_cases=scalar,
                matrix_calibration_and_holdouts=matrix,
                formal_stage_complete=False,long_goal_complete=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),(path,'keys')
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-10),(path,a,b)
    else:
        assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    result=compute()
    assert result['status']=='PASS',result
    destination=HERE/'results.json'
    if args.write_new:
        with destination.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
    else:
        compare(result,json.loads(destination.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',round=1328,spectral_cases=4,scalar_holdouts=12,
                          matrix_holdouts=18,qualified_preparation_families=2,
                          results_written=args.write_new),ensure_ascii=False))


if __name__=='__main__':
    main()
