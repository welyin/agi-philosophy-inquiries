"""Harrison--Scott pair: five CP-even inputs and a Run 3 charged-B holdout."""
import argparse
import itertools
import json
from math import cos, sin, pi
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))


def model(parameters, orientation=-1):
    scale, b, coupling, up_light, down_light = parameters
    up = orientation * 1j * scale * sin(pi/8)
    down = scale * cos(pi/8)
    systems = []
    for lam, light in [(up, up_light), (down, down_light)]:
        matrix = np.array([[light*abs(lam)**4, b*lam**3, 0],
                           [b*lam.conjugate()**3, b*abs(lam)**2, coupling*lam**2],
                           [0, coupling*lam.conjugate()**2, 1]])
        values, axes = np.linalg.eigh(matrix)
        order = np.argsort(abs(values))
        systems.append((values[order], axes[:,order]))
    (u, left_up), (d, left_down) = systems
    v = left_up.conj().T @ left_down
    angles = np.angle([-v[2,0]*v[2,2].conjugate()/(v[0,0]*v[0,2].conjugate()),
                       -v[1,0]*v[1,2].conjugate()/(v[2,0]*v[2,2].conjugate()),
                       -v[0,0]*v[0,2].conjugate()/(v[1,0]*v[1,2].conjugate())],deg=True)
    fitted = np.array([u[0]/u[1], d[0]/d[1], d[1]/d[2], abs(v[0,1]), abs(v[1,2])])
    return fitted, np.array([abs(u[1]/u[2]), abs(v[0,2]/v[1,2]),
                             *angles, (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag])


def calibration():
    ratios = DATA["mass_ratios"]
    raw = [ratios[name] for name in ("mu_mc","md_ms","ms_mb")]
    p, ep = DATA["K_product"]
    f, ef = DATA["K_form"]
    cb, *unc = DATA["BelleII_Vcb"]
    high, low = DATA["run_cb_ratio"]
    extra_us, extra_cb = DATA["relative_transfer_extra"]
    values = [row[0] for row in raw]+[p/f, cb*high/low]
    lower = [val-err for val,err in raw]+[(p-ep)/(f+ef)*(1-extra_us),
                                           (cb-sum(unc))*high/low*(1-extra_cb)]
    upper = [val+err for val,err in raw]+[(p+ep)/(f-ef)*(1+extra_us),
                                           (cb+sum(unc))*high/low*(1+extra_cb)]
    return np.array(values), np.array(lower), np.array(upper)


def jacobian(fn, point):
    result = []
    for k in range(len(point)):
        step = np.zeros_like(point)
        step[k] = 2.e-6*(1+abs(point[k]))
        result.append((fn(point+step)-fn(point-step))/(2*step[k]))
    return np.column_stack(result)


def solve(target, seed):
    def residual(params):
        return model(params)[0]/target-1
    params = np.array(seed, dtype=float)
    for _ in range(24):
        err = residual(params)
        if max(abs(err)) < 4e-10:
            return params
        delta = np.linalg.solve(jacobian(residual,params), err)
        for scale in (1., .5, .25, .125, .0625, .03125, .015625):
            proposal = params-scale*delta
            if np.linalg.norm(residual(proposal)) < np.linalg.norm(err):
                params = proposal
                break
        else:
            raise ValueError("Local branch did not converge")
    raise ValueError("Maximum Newton steps reached")


def verify_orientation(params, response):
    scale, b, coupling, cu, cd = params
    u, d = scale*sin(pi/8), scale*cos(pi/8)
    a, bb, dd, e = cu*u**4, b*u*u, cd*d**4, b*d*d
    p, q, r, s = b*u**3, -coupling*u*u, b*d**3, coupling*d*d
    up = np.array([[a,1j*p,0],[-1j*p,bb,q],[0,q,1]])
    down = np.array([[dd,r,0],[r,e,s],[0,s,1]])
    direct = np.linalg.det(up@down-down@up)/(2j)
    bracket = (s*(bb-1)+q*(1-e))*(s*(a-1)+q*(1-dd))-p*p*s*s-r*r*q*q
    def gaps(matrix):
        vals = np.linalg.eigvalsh(matrix)
        first, second, third = vals[np.argsort(abs(vals))]
        return (second-first)*(third-first)*(third-second)
    assert np.isclose(direct, p*r*bracket, rtol=1e-9, atol=1e-22)
    assert np.isclose(direct, gaps(up)*gaps(down)*response[-1], rtol=1e-9, atol=1e-22)
    assert .18 <= scale <= .28 and 0 < abs(b) <= .70
    assert max(abs(cu),abs(cd)) <= 2 and coupling > 0
    bound = .28*cos(pi/8)
    h = max(2*bound**4, .70*bound**2)
    ell = .70*bound**3
    assert bracket >= ((1-h)**2-ell**2)*(s+abs(q))**2


def run():
    central, lo, hi = calibration()
    branches = []
    for middle, up_sign, down_sign in itertools.product((-1,1),repeat=3):
        signs = np.array([up_sign,down_sign,middle,1,1])
        seed = [.224, .47*middle, .86,
                .47*middle*(1+up_sign*.002/(.224*sin(pi/8))**2),
                .47*middle*(1+down_sign*.0497/(.224*cos(pi/8))**2)]
        params = solve(central*signs,seed)
        fitted, response = model(params)
        verify_orientation(params, response)
        otherfit, other = model(params,1)
        assert np.allclose(fitted,otherfit)
        assert np.allclose(response[:2],other[:2])
        assert np.allclose(response[2:],-other[2:])
        # Jacobian in actual coordinates, using relative input variations.
        sing = np.linalg.svd(jacobian(lambda x:model(x)[0]/(central*signs),params),
                             compute_uv=False)
        samples = []
        for mask in itertools.product((0,1),repeat=5):
            target = np.where(mask,hi,lo)*signs
            trial = solve(target,params)
            samples.append(model(trial)[1])
        samples = np.array(samples)
        branches.append({
            "signs_middle_down_and_light_up_over_charm_down_over_strange": [middle,up_sign,down_sign],
            "parameters_lambda0_b_A0_cu_cd": params.tolist(), "fit_max_relative_error": max(abs(fitted/(central*signs)-1)),
            "smallest_jacobian_singular_value":min(sing),
            "response_mc_mt_abs_ub_over_cb_alpha_beta_gamma_J":response.tolist(),
            "CP_conjugate_gamma_mod180":float(other[4] % 180),
            "corner_response_min":samples.min(axis=0).tolist(),"corner_response_max":samples.max(axis=0).tolist(),
        })
    return {"phase_blind_inputs":central.tolist(),"low":lo.tolist(),"high":hi.tolist(),
            "branches":branches, "Run3_gamma_two_width_mod180_with_transfer":[
                DATA["Run3_gamma_deg"][0]-2*DATA["Run3_gamma_deg"][1]-DATA["gamma_transfer_deg"],
                DATA["Run3_gamma_deg"][0]+2*DATA["Run3_gamma_deg"][1]+DATA["gamma_transfer_deg"]]}


def rounded(values):
    if isinstance(values,dict):
        return {k:rounded(v) for k,v in values.items()}
    if isinstance(values,list):
        return [rounded(i) for i in values]
    if isinstance(values,(float,np.floating)):
        return round(float(values),9)
    return values


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = rounded(run())
    path = HERE / "results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    elif args.check:
        assert result == json.loads(path.read_text(encoding="utf-8-sig"))
        print("1445 eight local signed solutions, their conjugates and corner sensitivities verified")
    else:
        print(json.dumps(result,indent=2))
