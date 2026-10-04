"""718: finite-time recorded histories and transport of the actual instrument.

The full-H statements are proved in the note.  Matrix tests use an exact
conditional neutral sector of the original 64-mode mass/hopping coefficients.
They do not simulate the dynamic bosonic or interacting physical Gibbs state.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import joint_smooth_mode_contract as shared

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'round718_drafts'))
from reflected_history_entry import quadratic_matrix
TARGET = HERE/'joint_record_history_transport_results.json'
IDS = np.array([26, 27, 30, 31, 58, 59, 62, 63])
DIM = 256
I = np.eye(DIM)


def comm(a, b):
    return a@b-b@a


def fixture(lam=0.17):
    points = np.array([[.55, 0, 0, 0, .32], [.31, 0, 0, 0, -.27]])
    h = np.zeros((64, 64), complex)
    d = np.zeros_like(h)
    for j, p in enumerate(points):
        s = slice(32*j, 32*j+32)
        h[s, s], d[s, s] = shared.matter.mass_matrices(p)
    hop = np.zeros_like(h)
    hop[:32, 32:] = .29*np.exp(.37*lam)*np.eye(32)
    hop[32:, :32] = hop[:32, 32:].conj().T
    other = np.array([j for j in range(64) if j not in IDS])
    residual = max(np.linalg.norm((h+hop)[np.ix_(other, IDS)]),
                   np.linalg.norm(d[np.ix_(other, IDS)]))
    assert residual < 1e-14
    block = np.ix_(IDS, IDS)
    H = quadratic_matrix((h+hop)[block], d[block])
    G = quadratic_matrix(.37*hop[block], np.zeros((8, 8)))
    G2 = .37*G
    w = np.array([.8, 1.2])*np.exp(np.array([1.1, -.7])*lam)
    spinor = np.array([[1., 0.], [.7*np.exp(.3j), 0.]])
    full_f, _ = shared.profile(w, spinor)
    f = full_f[IDS]
    slopes = np.repeat([1.1, -.7], 4)
    mean = float(np.sum(abs(f)**2*slopes))
    var = float(np.sum(abs(f)**2*(slopes-mean)**2))
    fd = .5*(slopes-mean)*f
    fdd = (.25*(slopes-mean)**2-.5*var)*f
    P = np.outer(f, f.conj())
    Pd = np.outer(fd, f.conj())+np.outer(f, fd.conj())
    Pdd = (np.outer(fdd, f.conj())+np.outer(f, fdd.conj())
           +2*np.outer(fd, fd.conj()))
    zero = np.zeros((8, 8), complex)
    R = I-2*quadratic_matrix(P, zero)
    Rd = -2*quadratic_matrix(Pd, zero)
    Rdd = -2*quadratic_matrix(Pdd, zero)
    m = np.diag([float(bool(i & (1 << 7))) for i in range(DIM)])
    psi = np.zeros(DIM, complex)
    psi[0] = np.sqrt(.4)
    psi[12] = np.sqrt(.35)*np.exp(.6j)
    psi[192] = np.sqrt(.25)*np.exp(-.2j)
    assert np.linalg.norm(comm(m, R)) < 1e-13
    return dict(H=H, G=G, G2=G2, R=R, Rd=Rd, Rdd=Rdd, m=m,
                psi=psi, sector_residual=float(residual), f=f)


def spectrum(H):
    return np.linalg.eigh(H)


def evolution(sp, t):
    e, v = sp
    return (v*np.exp(-1j*t*e))@v.conj().T


def evolve_columns(sp, t, columns):
    e, v = sp
    return v@(np.exp(-1j*t*e)[:, None]*(v.conj().T@columns))


def evolution_derivative(sp, G, t):
    e, v = sp
    diff = e[:, None]-e[None, :]
    mid = .5*(e[:, None]+e[None, :])
    divided = -1j*t*np.exp(-1j*t*mid)*np.sinc(t*diff/(2*np.pi))
    return v@(divided*(v.conj().T@G@v))@v.conj().T


def thermal(sp, beta, G=None):
    e, v = sp
    x = e-e.min()
    weights = np.exp(-beta*x)
    z = weights.sum()
    rho = (v*(weights/z))@v.conj().T
    if G is None:
        return rho, (v*np.sqrt(weights/z))@v.conj().T
    diff = x[:, None]-x[None, :]
    mid = .5*(x[:, None]+x[None, :])
    u = beta*diff/2
    ratio = np.ones_like(u)
    np.divide(np.sinh(u), u, out=ratio, where=abs(u)>1e-12)
    divided = -beta*np.exp(-beta*mid)*ratio
    rawd = v@(divided*(v.conj().T@G@v))@v.conj().T
    rd = rawd/z-rho*np.trace(rawd).real/z
    return rho, rd


def trace_norm_hermitian(a):
    return float(np.sum(abs(np.linalg.eigvalsh((a+a.conj().T)/2))))


def one_wait_check(data):
    H, R = data['H'], data['R']
    D = (R@H@R-H)/2
    sp = spectrum(H)
    rho, root = thermal(sp, .73)
    resource = float(np.linalg.norm(D@root))
    rows = []
    for t in (.15, .65, 1.7):
        U = evolution(sp, t)
        Uprime = R@U@R
        branch = float(np.linalg.norm((Uprime-U)@root))
        cq = 0.
        iso2 = 0.
        for sign in (1, -1):
            K = (I+sign*R)/2
            early, late = U@K@root, K@U@root
            cq += trace_norm_hermitian(early@early.conj().T-late@late.conj().T)
            iso2 += np.linalg.norm(early-late)**2
        bound = min(2., 2*np.sqrt(2)*t*resource)
        identity = abs(iso2-branch**2/2)
        assert branch <= 2*t*resource+1e-12
        assert cq <= bound+1e-12 and identity < 1e-12
        rows.append(dict(time=t,retained_record_and_poststate_trace_norm=cq,
                         bound=float(bound),comparison_HS_norm=branch,
                         exact_isometry_identity_error=float(identity)))
    return dict(conditional_thermal_D_second_moment=resource**2,rows=rows)


def pure_difference(u, v):
    # Exact trace norm of |u><u|-|v><v|, including unnormalised records.
    a, b = np.vdot(u, u).real, np.vdot(v, v).real
    return float(np.sqrt(max(0., (a+b)**2-4*abs(np.vdot(u, v))**2)))


def multi_history_check(data):
    H, R, m, psi = (data[k] for k in ('H', 'R', 'm', 'psi'))
    D = (R@H@R-H)/2
    sp = spectrum(H)
    # Original sterile occupation instruments commuting with R, used only for
    # the finite-CAR calibration. Analytic physical result includes sin(s).
    ls = [np.diag(np.sqrt(np.diag(.5*I+sign*.25*m)))
          for sign in (1, -1)]
    waits = [.43, .77]
    prefix = psi[:, None]
    integrals = []
    for t in waits:
        absc, weights = np.polynomial.legendre.leggauss(32)
        integral = 0.
        for x, w in zip(absc, weights):
            rays = evolve_columns(sp, .5*t*(x+1), prefix)
            integral += .5*t*w*np.linalg.norm(D@rays)
        integrals.append(float(integral))
        end = evolve_columns(sp, t, prefix)
        prefix = np.concatenate([l@end for l in ls], axis=1)
    records = [(np.eye(DIM), ())]
    for t in waits:
        U = evolution(sp, t)
        records = [(l@U@v, labels+(r,)) for v, labels in records
                   for r, l in enumerate(ls)]
    distance, iso2, norm_early, norm_late = 0., 0., 0., 0.
    for V, _ in records:
        for sign in (1, -1):
            K = (I+sign*R)/2
            early, late = V@K@psi, K@V@psi
            distance += pure_difference(early, late)
            iso2 += np.linalg.norm(early-late)**2
            norm_early += np.vdot(early, early).real
            norm_late += np.vdot(late, late).real
    bound = min(2., 2*np.sqrt(2)*sum(integrals))
    assert abs(norm_early-1)<1e-12 and abs(norm_late-1)<1e-12
    assert distance <= bound+1e-11 and np.sqrt(iso2) <= np.sqrt(2)*sum(integrals)+1e-11
    return dict(wait_times=waits,unconditional_prefix_resource_integrals=integrals,
                all_record_blocks=8,trace_norm=distance,bound=float(bound),
                isometry_norm=float(np.sqrt(iso2)),normalization_error=float(
                    max(abs(norm_early-1),abs(norm_late-1))),
                rare_outcome_probabilities_not_divided=True)


def probability(data, sign, t=1.31):
    H, R, m = (data[k] for k in ('H', 'R', 'm'))
    sp = spectrum(H)
    rho, _ = thermal(sp, .73)
    K = (I+sign*R)/2
    V = m@evolution(sp, t)
    return float(np.trace(V@K@rho@K@V.conj().T).real)


def geometric_history_check(data):
    H, G, G2, R, Rd, Rdd, m = (data[k] for k in
                              ('H', 'G', 'G2', 'R', 'Rd', 'Rdd', 'm'))
    Hp = R@H@R
    gamma = Rd@R
    gammad = Rdd@R+Rd@Rd
    Gp = R@G@R+comm(gamma, Hp)
    Gp2 = (R@G2@R+2*comm(gamma, R@G@R)+comm(gammad, Hp)
           +comm(gamma, comm(gamma, Hp)))
    direct2 = (Rdd@H@R+R@H@Rdd+2*Rd@H@Rd
               +2*Rd@G@R+2*R@G@Rd+R@G2@R)
    second_error = float(np.linalg.norm(Gp2-direct2))
    sp, spp = spectrum(H), spectrum(Hp)
    t = 1.31
    U, Ud = evolution(sp, t), evolution_derivative(sp, G, t)
    Up = R@U@R
    Upd = R@Ud@R+comm(gamma, Up)
    independently = evolution_derivative(spp, Gp, t)
    rho, rhod = thermal(sp, .73, G)
    V, Vd, Vp, Vpd = m@U, m@Ud, m@Up, m@Upd

    def bilinear(v, vd, w, wd):
        a = v@rho@w.conj().T
        ad = vd@rho@w.conj().T+v@rhod@w.conj().T+v@rho@wd.conj().T
        return a, ad

    a, ad = bilinear(V, Vd, V, Vd)
    b, bd = bilinear(Vp, Vpd, Vp, Vpd)
    c, cd = bilinear(V, Vd, Vp, Vpd)
    endpoint = float(np.trace(c@Rd).real)
    crossd = np.trace(cd@R+c@Rd)
    rows = []
    for sign in (1, -1):
        K, Kd = (I+sign*R)/2, sign*Rd/2
        actual, actuald = bilinear(V@K, Vd@K+V@Kd, V@K, Vd@K+V@Kd)
        derivative = float(np.trace(actuald).real)
        combined = float((np.trace(ad).real+np.trace(bd).real
                          +2*sign*crossd.real)/4)
        noendpoint = combined-sign*endpoint/2
        nocross = float((np.trace(ad).real+np.trace(bd).real)/4)
        no_thermal = derivative-float(np.trace(V@K@rhod@K@V.conj().T).real)
        finite = []
        for eps in (2e-3, 5e-4, 1.25e-4):
            fd = (probability(fixture(.17+eps),sign)
                  -probability(fixture(.17-eps),sign))/(2*eps)
            finite.append(dict(step=eps,derivative=fd,error=abs(fd-derivative)))
        err = abs(combined-derivative)
        assert err<1e-12 and finite[-1]['error']<2e-8
        rows.append(dict(first_outcome_sign=sign,actual_joint_probability=float(np.trace(actual).real),
                         derivative=derivative,transport_identity_error=err,finite_differences=finite,
                         omit_endpoint_error=abs(noendpoint-derivative),
                         omit_cross_history_error=abs(nocross-derivative),
                         omit_prepared_thermal_response_error=abs(no_thermal-derivative)))
    assert max(r['omit_endpoint_error'] for r in rows)>1e-5
    assert max(r['omit_prepared_thermal_response_error'] for r in rows)>1e-5
    evolution_error = float(np.linalg.norm(Upd-independently))
    assert second_error<1e-12 and evolution_error<1e-11
    return dict(second_source_connection_identity_error=second_error,
                independently_differentiated_evolution_error=evolution_error,
                cross_endpoint_derivative=endpoint,rows=rows)


def run():
    data = fixture()
    one = one_wait_check(data)
    multi = multi_history_check(data)
    geo = geometric_history_check(data)
    dependencies = ('research_note_598.md','research_note_623.md','research_note_624.md',
                    'research_note_625.md','research_note_633.md','research_note_634.md',
                    'research_note_704.md','research_note_716.md','research_note_717.md',
                    'joint_smooth_mode_contract.py','round718_drafts/reflected_history_entry.py')
    return dict(round=718,tests_run=3,failures=0,errors=0,
                original_coefficient_modes=64,conditional_neutral_modes=8,
                conditional_Fock_dimension=DIM,neutral_sector_residual=data['sector_residual'],
                finite_time_record=one,multiple_record_history=multi,geometry_source=geo,
                scope=dict(full_H_finite_time_bound_proved_analytically=True,
                           first_classical_record_and_poststate_retained=True,
                           shared_domain_and_second_source_framework_reused=True,
                           numerical_fixture_is_conditional_not_full_boson_history=True,
                           numerical_thermal_state_is_not_original_full_Gibbs=True,
                           no_reset_of_actual_postrecord_state=True,
                           no_uniform_spatial_limit_or_causal_device_claim=True,
                           no_higher_moment_counterexample_scan=True,
                           goal_not_completed=True),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                                   for n in dependencies})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
