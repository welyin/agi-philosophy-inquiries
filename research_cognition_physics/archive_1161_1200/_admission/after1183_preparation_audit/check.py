"""Read-only by default: exact-flow diagnostics for a rejected preparation shortcut.

No finite-dimensional CCR truncation, field-equation evolution or physical
quantum preparation is certified by this classical calculation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]


def calc():
    x = np.arange(256) * (2 * np.pi / 256)
    z, dz, ddz = np.sin(x), np.cos(x), -np.sin(x)
    rows = []
    for lam in (0.4, 0.5, 0.6):
        pi = np.stack((2 + 0*z, 2 + z, 2 + 0*z))
        dpi = np.stack((0*z, dz, 0*z))
        trace = pi.sum(axis=0)
        a = pi - lam * trace
        az = np.array([-lam, 1-lam, -lam])[:, None]
        w = -2 * ((pi*pi).sum(axis=0) - lam*trace**2)
        wz = (24*lam-8) + (4*lam-4)*z
        wzz = 4*lam-4
        f = np.sqrt(w)
        ratio = a/f
        rz = az/f - a*wz/(2*f**3)
        rzz = -az*wz/f**3 - a*wzz/(2*f**3) + 3*a*wz**2/(4*f**5)
        for eps in (0.01, 0.1, 0.5, 1.0):
            # Full spatially covariant canonical flow: q_i*pi_i conserved.
            u = 2*eps*ratio
            du = 2*eps*rz*dz
            ddu = 2*eps*(rzz*dz**2+rz*ddz)
            q = np.exp(u)
            pif = pi/q
            dpif = (dpi-pi*du)/q
            dq = q*du
            determinant_root = np.sqrt(q.prod(axis=0))
            ricci = (-ddu[1]-ddu[2] - (du[1]**2+du[2]**2+du[1]*du[2])/2
                     +du[0]*(du[1]+du[2])/2)/q[0]
            # D_1 = sum pi_i q_i' - 2 (q_1*pi_1)' with constant phi.
            d1 = (pif*dq).sum(axis=0) - 2*(dq[0]*pif[0]+q[0]*dpif[0])
            afinal = ((q*pif)**2).sum(axis=0) - lam*(q*pif).sum(axis=0)**2
            kin_source = (afinal+f*f/2)/determinant_root
            h = kin_source-determinant_root*ricci
            assert np.max(np.abs(d1)) < 3e-14
            assert np.max(np.abs(kin_source)) < 8e-14
            row = dict(lam=lam, phi=eps, min_metric=float(q.min()),
                       momentum_max=float(np.max(np.abs(d1))),
                       algebraic_energy_max=float(np.max(np.abs(kin_source))),
                       full_H_at_x0=float(h[0]))
            if lam == 0.5:
                exact = -math.exp(-eps/(2*math.sqrt(3)))*(eps/(3*math.sqrt(3))-2*eps**2/27)
                assert exact < 0
                assert abs(h[0]-exact) < 2e-14
                row['closed_form_H_x0'] = exact
            # Coordinate-frozen Pi-only flow fails D, even at lambda=1/2.
            fpi = -2*ratio
            dfpi = -2*rz*dz
            qbad = 1-eps*fpi
            dqbad = -eps*dfpi
            dbad = (pi*dqbad).sum(axis=0)-2*(dqbad[0]*pi[0]+qbad[0]*dpi[0])
            expected = 4*eps*dfpi[0]
            assert np.max(np.abs(dbad-expected)) < 2e-14
            row['Pi_only_D_max'] = float(np.max(np.abs(dbad)))
            rows.append(row)
    sources = ['1183/proof.md', '1183/acceptance.json',
               '_admission/after1183_preparation_screen.md']
    return dict(status='PASS', kind='unnumbered_route_screen',
                new_scientific_rounds=0, new_counted_groups=0,
                source_sha256={n:hashlib.sha256((STAGE/n).read_bytes()).hexdigest() for n in sources},
                exact_H_x0='-exp(-phi/(2*sqrt(3)))*(phi/(3*sqrt(3))-2*phi^2/27)',
                scope='classical full canonical flow; not actual quantum preparation; approximate preparation not excluded',
                diagnostic_rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    actual = calc()
    output = HERE/'results.json'
    if args.write:
        if output.exists():
            raise SystemExit('Refusing to overwrite saved evidence')
        output.write_text(json.dumps(actual, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print('WRITE_ONCE_PASS')
    else:
        assert actual == json.loads(output.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH; unnumbered preparation screen; science counts unchanged')
