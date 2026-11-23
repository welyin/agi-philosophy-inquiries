"""Rebuild the type-I mass matrices and verify the registered data windows."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def observations(config):
    q = read(HERE/config["up_quark_CKM"])["high_scale"]
    lepton = read(HERE/config["charged_leptons"])["scales"]["1e16"][3:]
    osc = read(HERE/config["oscillation"])
    no = osc["NO"]
    rows = q["up"] + q["down"] + lepton + q["angles_theta12_theta13_theta23_delta"]
    rows += [no["solar_gap"], no["atmospheric_gap"], no["solar_sine_squared"],
             osc["DayaBay_amplitude"], no["atmospheric_sine_squared"]]
    central, lower, upper = zip(*[(r[0],r[1],r[-1]) for r in rows])
    return np.array(central), np.array(lower), np.array(upper)


def ckm(angles):
    a,b,c,phase = angles
    sa,sb,sc = np.sin([a,b,c])
    ca,cb,cc = np.cos([a,b,c])
    z = np.exp(1j*phase)
    return np.array([[ca*cb, sa*cb, sb/z],
        [-sa*cc-ca*sc*sb*z, ca*cc-sa*sc*sb*z, sc*cb],
        [sa*sc-ca*cc*sb*z, -ca*sc-sa*cc*sb*z, cc*cb]])


def mixing(v):
    p = abs(v)**2
    c13 = np.sqrt(1-p[0,2])
    s12,s13,s23 = abs(v[0,1])/c13, abs(v[0,2]), abs(v[1,2])/c13
    c12,c23 = abs(v[0,0])/c13, abs(v[2,2])/c13
    j = np.imag(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj())
    sin_delta = j/(c12*s12*c23*s23*c13**2*s13)
    cos_delta = (s12**2*s23**2+c12**2*c23**2*s13**2-p[2,0])/(2*c12*s12*c23*s23*s13)
    return s12,s13,s23,np.arctan2(sin_delta,cos_delta)


def takagi(m):
    left, spec, right_h = np.linalg.svd(m)
    matrix = (left*np.sqrt(np.diag(right_h@left.conj())))[:,::-1]
    assert np.linalg.norm(m-(matrix*spec[::-1])@matrix.T) < 1e-10*spec[0]
    return matrix, spec[::-1]


def from_charge_chart(point, central, widths):
    z = np.array(point)
    charged = central[:13]+z[7:20]*widths[:13]
    up = np.diag(charged[:3]*np.exp(1j*np.r_[z[:2],0]))
    axes = ckm(charged[9:13]).conj()
    down = (axes*(charged[3:6]*np.exp(1j*np.r_[z[2:4],0])))@axes.T
    q = (z[4]+1j*z[5])*charged[5]/charged[2]
    e0 = down+q*up
    c = np.exp(1j*z[6])*charged[8]/np.linalg.norm(e0,2)
    j = c*q
    e = c*e0
    x = down-e
    h,f = (3*down+e)/4, x/4
    w = up+(3+c)/j*x
    raw_n = w@np.linalg.solve(x,w.T)
    mass = np.linalg.svd(raw_n,compute_uv=False)[::-1]
    scale = np.sqrt((central[14]+z[20]*widths[14])/(mass[2]**2-mass[1]**2))
    rot,hdiag = takagi(h)
    upper = (rot.conj().T@f@rot.conj())[np.triu_indices(3)]
    r = abs((1-c)/j)
    s = (-3-c)/(1-c)
    t = scale*r**2/4
    return dict(H_diag=hdiag.tolist(),F_real_upper=upper.real.tolist(),F_imag_upper=upper.imag.tolist(),
                r=float(r),s_real=float(s.real),s_imag=float(s.imag),t_eV=float(t))


def observables_from_tensors(p):
    h = np.diag(p["H_diag"])
    inds = np.triu_indices(3)
    f = np.zeros((3,3),dtype=complex)
    f[inds] = np.array(p["F_real_upper"])+1j*np.array(p["F_imag_upper"])
    f[inds[1],inds[0]] = f[inds]
    r,s,t = p["r"],complex(p["s_real"],p["s_imag"]),p["t_eV"]
    w = h-3*s*f
    n = t*w@np.linalg.solve(f,w.T)
    matrices = [r*(h+s*f),h+f,h-3*f,n]
    specs,rv = [],[]
    for m in matrices:
        _,sv,vh = np.linalg.svd(m)
        specs.append(sv[::-1])
        rv.append(vh.conj().T[:,::-1])
    u,d,e,m = specs
    a,b,c,delta_q = mixing(rv[0].conj().T@rv[1])
    pmns = rv[2].conj().T@rv[3]
    al,bl,cl,delta_l = mixing(pmns)
    obs = np.r_[u,d,e,np.arcsin([a,b,c]),delta_q,m[1]**2-m[0]**2,m[2]**2-m[1]**2,
                al**2,4*bl**2*(1-bl**2),cl**2]
    evec = rv[2][:,0]
    mdhat = (d[0]+u[0])/2
    qmass = np.sqrt((d[1]**2-mdhat**2)/(d[0]**2-u[0]**2))
    return obs, dict(delta_lepton_rad=float(delta_l),m0_meV=float(1000*m[0]),
        m_beta_meV=float(1000*np.linalg.norm(pmns[0]*m)),
        m_bb_meV=float(1000*abs(evec.T@n@evec)),eta_Q_G=float(qmass))


def test_point(p,cfg,central,lower,upper):
    vals, predictions = observables_from_tensors(p)
    pulls = (vals-central)/np.where(vals<central,lower,upper)
    arc_cfg = read(HERE/cfg["phase_holdout"])
    lo,hi = arc_cfg["phase_credible_arcs_pi"]["NO_99.73"]
    guard = arc_cfg["joint_published_rounding_guard_pi"]
    angle = predictions["delta_lepton_rad"]/np.pi
    inside = any(lo-guard<=angle+2*i<=hi+guard for i in [-1,0,1])
    q0,qw = cfg["eta_Q"]
    qg = predictions["eta_Q_G"]
    allowed_t = []
    for widths in (1,2,3.05):
        a = max(math.sqrt((central[k]-widths*lower[k])/vals[k]) for k in (13,14))
        b = min(math.sqrt((central[k]+widths*upper[k])/vals[k]) for k in (13,14))
        allowed_t.append([float(a),float(b)] if a<=b else None)
    return dict(tensors=p,values=vals.tolist(),pulls=pulls.tolist(),max_abs_pull=float(max(abs(pulls))),
                max_other_than_d=float(max(abs(np.delete(pulls,3)))),
                max_other_than_atmospheric=float(max(abs(pulls[:-1]))),
                predictions=predictions,delta_in_rounded_NO_arc=inside,
                Q_G_over_two_width_Q_low_energy_interval=[qg/(q0+2*qw),qg/(q0-2*qw)],
                t_relative_intervals_from_two_mass_gaps_only=allowed_t,
                down_transfer_needed_for_central=central[3]/vals[3])


def local_refit(cfg,central,lower):
    """Optional NumPy-only refit: withhold y_d and scan an independently chosen source phase."""
    first = np.array(cfg["models"]["release_down_A"])
    x = np.r_[first[:7],math.log((central[3]+first[10]*lower[3])/1e-6)]
    # Start slightly off the saved branch. Hold two source coordinates fixed.
    x[[2,3,4,5,6,7]] += np.array([.001,-.002,.0003,-.0004,.0005,.001])
    movable = [2,3,4,5,6,7]
    target_ids = [6,7,13,15,16,17]
    def residual(v):
        coords=np.zeros(21)
        coords[:7]=v[:7]
        coords[10]=(1e-6*np.exp(v[-1])-central[3])/lower[3]
        obs,_=observables_from_tensors(from_charge_chart(coords,central,lower))
        return np.log(obs[target_ids]/central[target_ids])
    for iteration in range(60):
        f=residual(x)
        if max(abs(f))<2e-10:
            return dict(steps=iteration,down=float(1e-6*np.exp(x[-1])),
                        max_log_residual=float(max(abs(f))))
        jac=[]
        for k in movable:
            step=np.zeros(8)
            step[k]=1e-6
            jac.append((residual(x+step)-residual(x-step))/2e-6)
        dx=np.linalg.solve(np.column_stack(jac),-f)
        for i in range(20):
            trial=x.copy()
            trial[movable]+=dx*2**(-i)
            if np.linalg.norm(residual(trial))<np.linalg.norm(f):
                x=trial
                break
        else:
            raise ArithmeticError("local correction failed")
    raise ArithmeticError("iteration limit")


def evaluate():
    cfg=read(HERE/"inputs.json")
    center,lower,upper=observations(cfg)
    refs=runpy.run_path(str((HERE/cfg["historical_source"]).resolve()))
    x=np.array(list(map(float,refs["fractions"](refs["POINTS"]["MN_no_RGE_eq15"]))))
    parent=dict(H_diag=x[:3].tolist(),F_real_upper=x[3:9].tolist(),F_imag_upper=x[9:15].tolist(),
                r=x[15],s_real=x[16],s_imag=x[17],t_eV=x[18])
    points={"published_seed":test_point(parent,cfg,center,lower,upper)}
    for key, point in cfg["models"].items():
        points[key]=test_point(from_charge_chart(point,center,lower),cfg,center,lower,upper)
    assert points["historical_search_joint"]["max_abs_pull"] < 3.047
    assert points["profile_down"]["max_other_than_d"] < 1.349
    assert points["profile_atmospheric"]["max_other_than_atmospheric"] < 1.949
    for key in ("release_down_A","release_down_B"):
        assert points[key]["max_other_than_d"]<1e-6
        assert abs(points[key]["pulls"][3])>30
    return {"labels":["yu","yc","yt","yd","ys","yb","ye","ymu","ytau","theta12_q","theta13_q",
        "theta23_q","delta_q","solar_gap","atmospheric_gap_32","s12_l_squared","sin2_2theta13_l",
        "s23_l_squared"], "center":center.tolist(),"lower_width":lower.tolist(),"upper_width":upper.tolist(),
        "candidates":points,"local_refit":local_refit(cfg,center,lower)}


def compare(actual,expected):
    if isinstance(actual,dict):
        assert actual.keys()==expected.keys()
        for k,v in actual.items():
            compare(v,expected[k])
    elif isinstance(actual,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):
            compare(a,b)
    elif isinstance(actual,(float,np.floating)):
        assert np.isclose(actual,expected,rtol=3e-8,atol=5e-8)
    else:
        assert actual==expected


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    options=parser.add_mutually_exclusive_group(required=True)
    options.add_argument("--check",action="store_true")
    options.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    output=evaluate()
    path=HERE/"results.json"
    if args.check:
        compare(output,read(path))
        print("1556: charged, oscillation, held-out entries and local refit verified")
    else:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(output,stream,indent=2)
            stream.write("\n")
