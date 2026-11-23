"""An independent 100 TeV matching and a path from the electroweak vacuum."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def gluino_response(g, mz, m_squark, m_gluino):
    # Masses in TeV. Run directly in inverse alpha.
    inv = 4 * np.pi / g**2 + (7*np.log(m_gluino/mz) + 5*np.log(m_squark/m_gluino))/(2*np.pi)
    t, w = np.polynomial.legendre.leggauss(200)
    x = (1 + t)/2
    # Equal squarks: a one-dimensional Feynman parameter in the exact kernel.
    group = x / (x*m_squark**2 + (1-x)*m_gluino**2)
    integral = w @ group / 2
    z = (m_gluino/m_squark)**2
    exact = (1-z+z*np.log(z))/(1-z)**2/m_squark**2
    np.testing.assert_allclose(integral, exact, rtol=5e-7)
    return 2*m_gluino*exact/(3*np.pi*inv), 1/inv


def gaussian_bound(a, yuk, vd, vu, mu, m, force_numeric):
    # Neutral Higgs quadratic, including mu^2, and both squark soft masses
    # each equal m^2. On the explicit path the D term does not vary.
    extra = 2*yuk*yuk*vd*vd + 2*abs(yuk*mu*vu)
    mass2 = 3*m*m + extra
    trilinear = abs(a*yuk)
    den = 4*trilinear**2/(81*mass2) - 3*yuk**2/16
    if den <= 0:
        return None
    amplitude = 9*mass2/(4*trilinear)
    action = 9*np.pi**2/(4*den)
    radius = np.sqrt(1.5/(amplitude**2 * den))
    out = dict(Y=yuk, soft_a_TeV=a, physical_T_TeV=trilinear,
               field_at_center_TeV=amplitude, radius_TeV_inverse=radius, action_upper=action)
    if force_numeric:
        # Quadrature along the actual curved Higgs/squark path. Use the worst
        # possible sign of the mu cross term and the vacuum's Bmu > 0.
        mz = 0.0911876
        cos2 = (vd**2-vu**2)/(vd**2+vu**2)
        Bmu = (m*m+0.5*mz**2*cos2)*vd/vu
        xs, ws = np.polynomial.legendre.leggauss(140)
        r, weight = 5*(1+xs), 5*ws
        ph = amplitude*np.exp(-r*r)
        hd = np.hypot(vd, ph)
        grad = 2+ph**2/(vd*vd+ph**2)
        kinetic = 2*np.pi**2*(weight @ (r**3*grad*(2*r*ph)**2))
        pot = (mass2*ph**2 - 2*trilinear*hd*ph**2 + 3*yuk**2*ph**4
               - 2*Bmu*vu*ph**2/(hd+vd))
        vol = 2*np.pi**2*(weight @ (r**3*pot))
        assert vol < 0 and Bmu > 0
        upper_exact_path = -kinetic**2/(4*vol)
        assert upper_exact_path <= action*(1+1e-9)
        out["curved_path_action_upper_with_worst_mu_sign"] = float(upper_exact_path)
    return out


def calculate():
    m, gluino, mu, tb = 100., 1., 1., 10.
    g = load("../1369/inputs.json")
    eta, alpha = gluino_response(g["g3"][0], g["mu_ref_GeV"]/1000, m, gluino)
    M0 = 3*m*m
    full_vev = 0.174
    vd = full_vev/np.sqrt(1+tb*tb)
    vu = tb*vd
    matching = load("../1402/inputs.json")["scales"]["1e5"]
    d_sm,s_sm = (matching[i][0] for i in (0,1))
    X = np.array([d_sm,s_sm])*np.sqrt(1+tb*tb)
    gs_base = 1+eta*mu*tb
    records={}
    for tag in ("release_down_A","release_down_B"):
        vals=load("../1556/results.json")["candidates"][tag]["values"]
        high_ratio=vals[3]/vals[4]
        source_ratio=load("../1557/independent_results.json")["cases"][tag[-1]]["full"]["low_ratios"][0]
        assert abs(source_ratio/high_ratio-1)<0.001
        R=X[0]/X[1]/source_ratio
        R_low=(d_sm-2*matching[0][1])/(s_sm+2*matching[1][1])/source_ratio
        t=(R-1)/(eta*(R/X[0]+1/X[1]))
        t_box=(R_low-1)/(eta*(R_low/((d_sm-2*matching[0][1])*np.sqrt(101))
                                  +1/((s_sm+2*matching[1][1])*np.sqrt(101))))
        a_star=gs_base*(R-1)/((R+1)*eta)
        # For G=1+eta(mu tan beta-a), a sharp solution has T_d<0<T_s.
        gs_equalT=gs_base/(1+eta*t/X[1])
        candidates=[("strange_fixed",1.),
                    ("minimax_normalized_soft",gs_base-eta*a_star),
                    ("minimax_physical_trilinear",gs_equalT)]
        samples={}
        for name,gs in candidates:
            gains=np.array([R*gs,gs])
            Y=X/gains
            a=(gs_base-gains)/eta
            np.testing.assert_allclose(Y[0]/Y[1],source_ratio,rtol=3e-12)
            np.testing.assert_allclose(Y*(gs_base-eta*a),X)
            if name=="minimax_physical_trilinear":
                np.testing.assert_allclose(Y*a,[-t,t],rtol=3e-12)
            if name=="minimax_normalized_soft":
                np.testing.assert_allclose(a,[-a_star,a_star],rtol=1e-10)
            bounds=[gaussian_bound(aa,yy,vd,vu,mu,m,True) for aa,yy in zip(a,Y)]
            action=min(b["action_upper"] for b in bounds if b is not None)
            samples[name] = dict(gains=gains.tolist(),Y=Y.tolist(),soft_a_TeV=a.tolist(),
                                 physical_T_TeV=(Y*a).tolist(),ray_bounds=bounds,
                                 fastest_action_upper=action)
        # Exact lower bound over all complex G with |Gd/Gs| fixed. Rewriting
        # tau=(Xd-Xs*r)/(eta*(1+r)) also proves monotonicity over the box.
        np.testing.assert_allclose(t, (X[0]-source_ratio*X[1])/(eta*(1+source_ratio)))
        records[tag]=dict(upstream_ds_ratio=source_ratio, gain_central=R, gain_two_widths_lower=R_low,
                          minimum_max_physical_T_central_TeV=t,minimum_max_physical_T_two_width_TeV=t_box,
                          minimum_common_soft_ratio=a_star/np.sqrt(3*M0),scenarios=samples)
    return dict(alpha_s_100TeV=alpha, eta_inverse_TeV=eta, cases=records)


def compare(expected, computed):
    if isinstance(computed, dict):
        assert set(computed)==set(expected)
        for k,v in computed.items():
            compare(expected[k],v)
    elif isinstance(computed, list):
        assert len(expected)==len(computed)
        for a,b in zip(expected,computed):
            compare(a,b)
    elif computed is None:
        assert expected is None
    else:
        np.testing.assert_allclose(computed, expected, rtol=2e-8, atol=1e-8)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    result=calculate()
    if args.write_new:
        with (HERE/"independent_results.json").open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    else:
        compare(load("independent_results.json"),result)
    print("1560 independent: matched trilinears and electroweak-origin tunneling paths verified")
