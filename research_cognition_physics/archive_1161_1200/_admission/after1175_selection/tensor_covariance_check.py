"""Small independent Clifford/covariance audit; default read-only, --record once."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.])
GAMMA = np.array([np.kron(Z, I2), np.kron(X, X),
                  np.kron(X, Y), np.kron(X, Z)])
PAIRS = [(a, b) for a in range(4) for b in range(a + 1, 4)]
LGEN = np.array([1j * GAMMA[a] @ GAMMA[b] for a, b in PAIRS])
BASIS = np.concatenate([GAMMA, LGEN])
EPS, NU, COUP, OM, MASS = .15, .35, .12, .5, 2.


def reconstruct(c):
    return .5 * (np.eye(4) + np.einsum("...a,aij->...ij", c, BASIS))


def coefficients(p):
    return .5 * np.einsum("aij,...ji->...a", BASIS, p).real


def rhs10(c, d):
    s = np.zeros(c.shape[:-1] + (4, 4))
    for j, (a, b) in enumerate(PAIRS):
        s[..., a, b] = c[..., 4 + j]
        s[..., b, a] = -c[..., 4 + j]
    r = c[..., :4]
    dr = -2 * np.einsum("...ab,...b->...a", s, d)
    ds = np.stack([2 * (d[..., b] * r[..., a] - d[..., a] * r[..., b])
                   for a, b in PAIRS], axis=-1)
    return np.concatenate([dr, ds], axis=-1)


def sources(c, k):
    return np.r_[2 * MASS * np.mean(c[:, 0]),
                 2 * np.mean(k * c[:, 1:4], axis=0)]


def fields(gs, gt, k):
    vel = np.array([1 + NU * gs[2] + COUP * gt[2],
                    1 + NU * gs[2] - COUP * gt[2], 1 + NU * gs[2]])
    return np.column_stack([np.full(len(k), MASS * (1 + EPS * gs[2])),
                            k * vel])


def rhs(state, k, matrix=False):
    q, gs, gt = state
    c = coefficients(q) if matrix else q
    em, tx, ty, tz = sources(c, k)
    d = fields(gs, gt, k)
    if matrix:
        h = np.einsum("na,aij->nij", d, GAMMA)
        dq = -1j * (h @ q - q @ h)
    else:
        dq = rhs10(c, d)
    return (dq, 2 * np.cross([-OM, 0, EPS * em + NU * (tx + ty + tz)], gs),
            2 * np.cross([-OM, 0, COUP * (tx - ty)], gt))


def rk4(y, k, dt, matrix):
    def add(a, b, f):
        return tuple(x + f * z for x, z in zip(a, b))
    k1 = rhs(y, k, matrix)
    k2 = rhs(add(y, k1, dt / 2), k, matrix)
    k3 = rhs(add(y, k2, dt / 2), k, matrix)
    k4 = rhs(add(y, k3, dt), k, matrix)
    return tuple(a + dt * (b + 2*c + 2*d + e) / 6
                 for a, b, c, d, e in zip(y, k1, k2, k3, k4))


def energy(y, k, matrix=False):
    q, gs, gt = y
    c = coefficients(q) if matrix else q
    return float(2 * np.mean(np.sum(fields(gs, gt, k) * c[:, :4], axis=1))
                 - OM * (gs[0] + gt[0]))


def calculate():
    rng = np.random.default_rng(1175001)
    cliff = max(np.linalg.norm(GAMMA[a] @ GAMMA[b] + GAMMA[b] @ GAMMA[a]
                              - 2*(a == b)*np.eye(4)) for a in range(4) for b in range(4))
    c = rng.normal(size=(37, 10)) * .02
    d = rng.normal(size=(37, 4))
    h = np.einsum("na,aij->nij", d, GAMMA)
    p = reconstruct(c)
    direct = -1j * (h @ p - p @ h)
    compact = .5 * np.einsum("na,aij->nij", rhs10(c, d), BASIS)
    algebra_error = float(np.max(np.abs(direct - compact)))
    closure_error = float(np.max(np.abs(direct - .5*np.einsum(
        "na,aij->nij", coefficients(direct), BASIS))))
    n, ell = 5, .8
    labels = np.array([(a,b,c) for a in range(-2,3) for b in range(-2,3)
                       for c in range(-2,3)], dtype=int)
    k = np.sin(2*np.pi*labels/n) / ell
    c0 = np.zeros((len(k), 10))
    c0[:, 1] = np.sign(labels[:, 0])
    c0[labels[:, 0] == 0, 0] = 1
    cb0 = np.zeros_like(c0)
    cb0[:, 2] = np.sign(labels[:, 1])
    cb0[labels[:, 1] == 0, 0] = 1
    x0 = np.array([1.,0.,0.])
    y = (c0.copy(), x0.copy(), x0.copy())
    ym = (reconstruct(c0), x0.copy(), x0.copy())
    yb = (cb0, x0.copy(), x0.copy())
    e0 = energy(y, k)
    steps, dt = 160, .4/160
    for _ in range(steps):
        y = rk4(y, k, dt, False)
        ym = rk4(ym, k, dt, True)
        yb = rk4(yb, k, dt, False)
    unit = (np.eye(4) - GAMMA[1] @ GAMMA[2]) / np.sqrt(2.)
    lookup = {tuple(v): j for j,v in enumerate(labels)}
    invrot = [lookup[(b,-a,c)] for a,b,c in labels]
    pred_b = unit @ reconstruct(y[0])[invrot] @ unit.conj().T
    joint_error = max(float(np.max(np.abs(reconstruct(y[0])-ym[0]))),
                      float(np.max(np.abs(y[1]-ym[1]))),
                      float(np.max(np.abs(y[2]-ym[2]))))
    proj = reconstruct(y[0])
    symmetry_error = max(float(np.max(np.abs(reconstruct(yb[0])-pred_b))),
                         float(np.max(np.abs(yb[1]-y[1]))),
                         float(np.max(np.abs(yb[2]-y[2]*[1,-1,-1]))))
    vals = np.linalg.eigvalsh(proj)
    e_sym = abs(energy(yb,k)-energy(y,k))
    # Full Fock norm of one Hermitian quadratic equals max(sum positive,
    # -sum negative). Use a tiny central site with its +/-x and +/-y neighbours.
    def ti(axis):
        out = np.zeros((20,20), complex)
        for site, sign in ((1+2*axis,1), (2+2*axis,-1)):
            out[:4,4*site:4*(site+1)] = -1j*sign*GAMMA[axis+1]/(4*ell)
            out[4*site:4*(site+1),:4] = out[:4,4*site:4*(site+1)].conj().T
        return out
    ev = np.linalg.eigvalsh(ti(0)-ti(1))
    source_norm = max(float(ev[ev>0].sum()),float(-ev[ev<0].sum()))
    source_bound = 4/ell
    probe_norm = 3*source_norm
    result = {
        "status": "PASS", "scope": "un-numbered small-grid algebra/calibration, not candidate signal certification",
        "clifford": {"maximum_anticommutator_error":float(cliff),
                     "rhs10_vs_full_commutator_max_error":algebra_error,
                     "ten_component_closure_max_error":closure_error},
        "self_consistent_small_grid": {
            "grid":n,"modes":len(k),"ell":ell,"T":steps*dt,"RK4_steps":steps,
            "epsilon":EPS,"nu":NU,"tensor_g":COUP,"Omega_s_t":OM,"source_mass":MASS,
            "compact_vs_matrix_max_error":joint_error,
            "projector_max_error":float(np.max(np.abs(proj@proj-proj))),
            "occupancy_trace_max_error":float(np.max(np.abs(np.trace(proj,axis1=1,axis2=2)-2))),
            "coefficient_casimir_max_error":float(np.max(np.abs(np.sum(y[0]**2,axis=1)-1))),
            "covariance_eigenvalue_min":float(vals.min()),"covariance_eigenvalue_max":float(vals.max()),
            "energy_initial":e0,"energy_final":energy(y,k),"energy_drift":abs(energy(y,k)-e0),
            "A_B_rotation_max_error":symmetry_error,"A_B_energy_difference":e_sym,
            "g_s_final":y[1].tolist(),"g_t_A_final":y[2].tolist(),
            "g_t_B_final":yb[2].tolist(),"source_components_A_final":sources(y[0],k).tolist()},
        "full_Fock_local_bound": {
            "source_tix_minus_tiy_exact_norm":source_norm,"source_upper_bound":source_bound,
            "three_flavour_probe_exact_norm":probe_norm,"three_flavour_probe_upper_bound":12/ell},
    }
    assert algebra_error<1e-13 and closure_error<1e-13 and joint_error<1e-12
    assert symmetry_error<1e-12 and e_sym<1e-12
    assert result["self_consistent_small_grid"]["projector_max_error"]<1e-8
    assert result["self_consistent_small_grid"]["energy_drift"]<1e-8
    assert source_norm<=source_bound and probe_norm<=12/ell
    return result


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        assert abs(a-b)<=2e-12*(1+abs(a)),(a,b)
    else: assert a==b,(a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--record",action="store_true")
    args=parser.parse_args()
    result=calculate()
    result["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    path=HERE/"tensor_covariance_checks.json"
    if args.record:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
        print("RECORDED",path.name)
    else:
        compare(result,json.loads(path.read_text(encoding="utf-8")))
        print("READ_ONLY_MATCH")
    print(json.dumps(result,ensure_ascii=False))
