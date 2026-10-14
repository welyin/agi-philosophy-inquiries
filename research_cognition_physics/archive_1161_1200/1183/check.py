"""Read-only verification unless --write is explicitly supplied.

Classical initial-constraint witness, not a physical quantum implementation.
NumPy checks the full nonlinear constraint along its actual Hamiltonian tangent;
the proof, not finite sampling, supplies the global parameter/lapse statements.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def initial(lam, x, amplitude=1.0):
    z = amplitude * np.sin(x)
    zp = amplitude * np.cos(x)
    zpp = -z
    pi = np.stack((2 + 0*z, 2 + z, 2 + 0*z), axis=1)
    trace = 6 + z
    hg = (pi*pi).sum(axis=1) - lam*trace**2
    hgz = 4 - 12*lam + 2*(1-lam)*z
    p = np.sqrt(-2*hg)
    px = -hgz*zp/p
    qdot = 2*(pi-lam*trace[:, None])
    qdx = np.stack((-2*lam*zp, 2*(1-lam)*zp, -2*lam*zp), axis=1)
    qdxx = np.stack((-2*lam*zpp, 2*(1-lam)*zpp, -2*lam*zpp), axis=1)
    pidot = -2*pi*pi + 2*lam*trace[:, None]*pi
    pix = np.stack((0*z, zp, 0*z), axis=1)
    pidotx = -4*pi*pix + 2*lam*(zp[:, None]*pi+trace[:, None]*pix)
    return hg, p, px, pi, pix, qdot, qdx, qdxx, pidot, pidotx


def scalar_curvature(q, qx, qxx):
    # Exact scalar curvature of diag(A(x),B(x),C(x)); no target bracket used.
    A, B, C = q.T
    Ap, Bp, Cp = qx.T
    _, Bpp, Cpp = qxx.T
    return (-Bpp/(A*B)-Cpp/(A*C)
            + Ap*Bp/(2*A*A*B)+Ap*Cp/(2*A*A*C)
            + Bp*Bp/(2*A*B*B)+Cp*Cp/(2*A*C*C)
            - Bp*Cp/(2*A*B*C))


def tangent_constraints(lam, x, t, amplitude=1.0):
    hg, p, px, pi, pix, qd, qdx, qdxx, pid, pidx = initial(lam, x, amplitude)
    q = 1+t*qd
    assert q.min() > 0
    qx, qxx = t*qdx, t*qdxx
    mom, momx = pi+t*pid, pix+t*pidx
    rootq = np.sqrt(q.prod(axis=1))
    tr = (q*mom).sum(axis=1)
    grav = (((q*mom)**2).sum(axis=1)-lam*tr**2)/rootq
    grav -= rootq*scalar_curvature(q, qx, qxx)
    # phi=t*p, p_phi=p. Its gradient and q response must both be retained.
    matter = p*p/(2*rootq) + rootq*(t*px)**2/(2*q[:, 0])
    d1 = (-2*(qx[:, 0]*mom[:, 0]+q[:, 0]*momx[:, 0])
          +(mom*qx).sum(axis=1)+p*t*px)
    return grav+matter, d1, scalar_curvature(q, qx, qxx)


def compute():
    lo, hi = F(2, 5), F(3, 5)
    def hg(lam, z):
        return 12-36*lam+(4-12*lam)*z+(1-lam)*z*z
    assert hg(lo, F(-1)) == -1
    assert hg(hi, F(1)) == F(-62, 5)
    # dHg/dlambda=-(6+z)^2<0. At lo the z-polynomial is convex;
    # at hi its derivative (4z-16)/5 is negative throughout [-1,1].
    exact = {str(lam): str(2*lam-1) for lam in (lo, F(1, 2), hi)}
    x = np.arange(256)*2*np.pi/256
    rows = []
    largest_initial = 0.0
    for lam in (0.4, 0.5, 0.6):
        h0, d0, _ = tangent_constraints(lam, x, 0)
        largest_initial = max(largest_initial, float(np.max(np.abs(h0))), float(np.max(np.abs(d0))))
        M = 2+np.sin(x)
        for dt in (2e-6, 1e-6):
            hp, dp, rp = tangent_constraints(lam, x, dt)
            hm, dm, rm = tangent_constraints(lam, x, -dt)
            hd = (hp-hm)/(2*dt)
            dd = (dp-dm)/(2*dt)
            rd = (rp-rm)/(2*dt)
            target = -2*(1-2*lam)*np.sin(x)
            error = float(np.max(np.abs(hd-target)))
            derror = float(np.max(np.abs(dd)))
            rerror = float(np.max(np.abs(rd+target)))
            assert error < 2e-8 and derror < 2e-8 and rerror < 2e-8
            smeared = float(np.mean(M*hd))
            assert abs(smeared-(2*lam-1)) < 2e-8
            rows.append(dict(lam=lam, dt=dt, constraint_density_derivative_error=error,
                             momentum_derivative_error=derror, curvature_derivative_error=rerror,
                             smeared_derivative=smeared))
    assert largest_initial < 1e-13
    # Smooth positive lapse of genuine spatial dependence. No lapse sampling
    # is substituted for the integration-by-parts proof in proof.md.
    lapses = []
    for a in (0.0, 0.25, 0.5):
        N = 1+a*np.sin(x)
        Nx = a*np.cos(x)
        f = N*(-np.sin(x))+2*Nx*np.cos(x)
        lhs = float(np.mean((6+np.sin(x))*N*f))
        rhs = -float(np.mean(N*N*np.cos(x)**2))
        assert abs(lhs-rhs) < 1e-13
        assert -rhs >= (1-a)**2/2-1e-14
        lapses.append(dict(a=a, minimum_lapse=float(N.min()),
                           integrated_positive_obstruction=-rhs, integration_error=abs(lhs-rhs)))
    # CMC deletion: amplitude zero. Same positive scalar source remains.
    cmc = []
    for lam in (0.4, 0.6):
        hp, _, _ = tangent_constraints(lam, x, 1e-5, amplitude=0)
        hm, _, _ = tangent_constraints(lam, x, -1e-5, amplitude=0)
        val = float(np.max(np.abs((hp-hm)/(2e-5))))
        assert val < 2e-7
        cmc.append(dict(lam=lam, first_tangent_residual=val))
    sources = ['../archive_342_369/research_note_351.md',
               '../archive_342_369/research_note_352.md',
               '../archive_935_955/research_note_937.md',
               '../archive_935_955/research_note_938.md',
               '_shared/goal.md', '_admission/after1178_reference/local_boundary_mechanism.md']
    return dict(round=1183, status='PASS', groups=4,
                scope='classical Hamiltonian and momentum initial-constraint surface only; no quantum implementation or finite instrument certificate',
                source_sha256={name:hashlib.sha256((STAGE/name).read_bytes()).hexdigest() for name in sources},
                exact=dict(parameter_interval=['2/5','3/5'], gravity_density_bounds=['-62/5','-1'],
                           positive_matter_density_bounds=['1','62/5'], scalar_momentum_squared_bounds=['2','124/5'],
                           total_initial_constraints='H(x)=0, D_i(x)=0',
                           smeared_unit_lapse_derivative=exact,
                           positive_lapse_obstruction='integral N^2 |grad pi|^2 > 0 for every smooth strictly positive N',
                           no_all_close_lambda_finite_error_exclusion=True),
                initial_max_residual=largest_initial, nonlinear_tangent_checks=rows,
                positive_lapse_checks=lapses, cmc_deletion_checks=cmc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    data = compute()
    out = HERE/'results.json'
    if args.write:
        out.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        status = 'WRITTEN'
    else:
        assert data == json.loads(out.read_text(encoding='utf-8-sig')), 'saved results differ'
        status = 'READ_ONLY_MATCH'
    print(json.dumps(dict(status=status, groups=4, scope=data['scope']), ensure_ascii=False))


if __name__ == '__main__':
    main()
