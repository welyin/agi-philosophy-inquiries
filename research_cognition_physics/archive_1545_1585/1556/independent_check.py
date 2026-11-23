"""Reconstruct the matched poles and evaluate independent normalization budgets."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_json(p):
    return json.loads(p.read_text(encoding="utf-8"))


def takagi(m):
    _, s, vh = np.linalg.svd(m)
    axes = vh.conj().T[:,::-1]
    axes = axes*np.exp(-.5j*np.angle(np.diag(axes.T@m@axes)))
    np.testing.assert_allclose(axes.T@m@axes,np.diag(s[::-1]),atol=1e-10*max(s))
    return s[::-1],axes


def mixing(v):
    p = abs(v)**2
    s13 = p[0,2]
    s12, s23 = p[0,1]/(1-s13), p[1,2]/(1-s13)
    k = np.sqrt(s12*(1-s12)*s23*(1-s23)*s13)
    j = np.imag(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj())
    cosine = (p[1,0]-s12*(1-s23)-(1-s12)*s23*s13)/(2*k)
    phase = np.arctan2(j/(k*(1-s13)),cosine)
    return np.array([s12,s13,s23,phase])


def scale_free_required_width(reads, targets):
    """Exact minimax residual after permitting one common positive gain."""
    a = np.array(targets)
    x = a[:,0]/reads
    e = a[:,1]/reads
    costs = abs(x[:,None]-x)/(e[:,None]+e)
    i,j = np.unravel_index(np.argmax(costs),costs.shape)
    needed = costs[i,j]
    lo,hi = (x-needed*e).max(),(x+needed*e).min()
    np.testing.assert_allclose(lo,hi,atol=2e-12)
    return {"minimax_width":float(needed), "gain":float((lo+hi)/2),
            "limiting_pair":[int(i),int(j)]}

def finite_spectrum(up_dirac, majorana, light_ev, light_masses):
    # Canonical blocks in GeV, leading three masses in eV. All bounds are spectral norms.
    d = up_dirac
    right = majorana
    q = float(np.linalg.norm(np.linalg.solve(right,d.T),ord=2))
    c = np.linalg.norm(d,ord=2)/np.linalg.svd(right,compute_uv=False)[-1]
    assert 4*c*q<1
    rad = 2*q
    delta_ev = 1e9*np.linalg.norm(d,ord=2)*(4*c*q*q+rad**3)
    leading = -1e9*d@np.linalg.solve(right,d.T)
    np.testing.assert_allclose(leading,light_ev,rtol=3e-8,atol=1e-12)
    delta_h = (2*max(light_masses)+delta_ev)*delta_ev
    gap = np.min(np.diff(np.sort(np.square(light_masses))))
    assert gap>2*delta_h
    return {"right_singular_GeV":np.sort(np.linalg.svd(right,compute_uv=False)).tolist(),
            "active_weak_norm_loss_upper":rad**2,"light_mass_absolute_error_upper_eV":delta_ev,
            "light_direction_sine_upper":float(delta_h/(gap-delta_h)),
            "physical_heavy_absolute_shift_upper_GeV":float(np.linalg.norm(d,ord=2))}

def double_ratio(down, charged):
    return down[0]*charged[1]/(down[1]*charged[0])


def first_family_check(down, charged, down_inputs, charged_inputs):
    dv = np.array(down_inputs)
    ev = np.array(charged_inputs)
    lo_d,hi_d = dv[:,0]-2*dv[:,1],dv[:,0]+2*dv[:,1]
    lo_e,hi_e = ev[:,0]-2*ev[:,1],ev[:,0]+2*ev[:,1]
    low = lo_d[0]*lo_e[1]/(hi_d[1]*hi_e[0])
    high = hi_d[0]*hi_e[1]/(lo_d[1]*lo_e[0])
    val = double_ratio(down,charged)
    factor = max(1.,low/val,val/high)
    return {"first_family_double_ratio":float(val), "two_width_interval":[float(low),float(high)],
            "minimum_equal_relative_shape_budget":float((np.sqrt(factor)-1)/(np.sqrt(factor)+1))}

def reconstruct(point):
    h, re, im, r, s, t = point
    f = np.zeros((3,3),complex)
    upper = np.triu_indices(3)
    f[upper] = np.array(re)+1j*np.array(im)
    f[upper[::-1]] = f[upper]
    h = np.diag(h)
    up,down,lepton = r*(h+s*f),h+f,h-3*f
    n = h-3*s*f
    light = t*n@np.linalg.solve(f,n.T)
    u,bu = takagi(up)
    d,bd = takagi(down)
    e,be = takagi(lepton)
    masses,bn = takagi(light)
    angles = mixing(be.conj().T@bn)
    # A second read of the coherent electron pair uses the light matrix directly.
    mbb = abs((be.T@light@be)[0,0])
    np.testing.assert_allclose(mbb,abs(np.sum((be.conj().T@bn)[0]**2*masses)))
    high = load_json(HERE.parent/"1397"/"inputs.json")["high_scale"]
    charged = load_json(HERE.parent/"1402"/"inputs.json")["scales"]["1e16"][3:]
    low = load_json(HERE.parent/"1386"/"inputs.json")["yukawa"]
    ratio = first_family_check(d,e,high["down"],charged)
    low_window = first_family_check(d,e,low["down"],low["charged_leptons"])["two_width_interval"]
    # Worst signs/cross-family running. No claim of a full RGE integration.
    long_log = np.log(1e16/91.1876)
    quark_budget = long_log*(.02**2+.02**2+1.2**2*.01**2)
    variation = max(0.,np.log(low_window[0]/ratio["first_family_double_ratio"]))
    required = max(0.,32*np.pi**2/3*variation-quark_budget)
    scale = 174.  # common *example* single-Higgs normalization, in GeV
    finite = finite_spectrum(scale*r*n,-(scale*r)**2*1e9*f/t,light,masses)
    return {"first_family_comparison":ratio,"weak_scale_double_ratio_interval":low_window,
            "down_sector_best_common_gain":scale_free_required_width(d,high["down"]),
            "neutral_spread_integral_required":float(required),
            "one_loop_peak_active_yukawa_required_whole_interval":float(np.sqrt(required/long_log)),
            "one_loop_peak_if_lowest_threshold_1e9_GeV":float(np.sqrt(required/np.log(1e16/1e9))),
            "ckm_angles_and_delta":(np.r_[np.arcsin(np.sqrt(mixing(bu.conj().T@bd)[:3])),
                                           mixing(bu.conj().T@bd)[3]]).tolist(),
            "pmns_sine_squared_and_delta":angles.tolist(),"light_masses_eV":masses.tolist(),
            "effective_majorana_eV":float(mbb),"finite_check":finite}


def verify_all_exported_models():
    data = load_json(HERE/"results.json")
    high = load_json(HERE.parent/"1397"/"inputs.json")["high_scale"]
    ec = load_json(HERE.parent/"1402"/"inputs.json")["scales"]["1e16"][3:]
    osc = load_json(HERE.parent/"1386"/"inputs.json")
    no = osc["NO"]
    measured = (high["up"]+high["down"]+ec+high["angles_theta12_theta13_theta23_delta"]+
                [no["solar_gap"],no["atmospheric_gap"],no["solar_sine_squared"],
                 osc["DayaBay_amplitude"],no["atmospheric_sine_squared"]])
    central = np.array([entry[0] for entry in measured])
    for name,row in data["candidates"].items():
        par = row["tensors"]
        f = np.zeros((3,3),complex)
        tri = np.triu_indices(3)
        f[tri] = np.array(par["F_real_upper"])+1j*np.array(par["F_imag_upper"])
        f[tri[::-1]] = f[tri]
        h = np.diag(par["H_diag"])
        s = complex(par["s_real"],par["s_imag"])
        n = h-3*s*f
        m, bn = takagi(par["t_eV"]*n@np.linalg.solve(f,n.T))
        u,bu = takagi(par["r"]*(h+s*f))
        d,bd = takagi(h+f)
        e,be = takagi(h-3*f)
        q = mixing(bu.conj().T@bd)
        p = mixing(be.conj().T@bn)
        vals = np.r_[u,d,e,np.arcsin(np.sqrt(q[:3])),q[3],np.diff(m**2),
                     p[0],4*p[1]*(1-p[1]),p[2]]
        np.testing.assert_allclose(vals,row["values"],atol=1e-11,rtol=3e-8)
        # The input lists the lower and then the upper marginal half-width.
        width = np.array([entry[1] if val<entry[0] else entry[-1]
                          for val,entry in zip(vals,measured)])
        pull = (vals-central)/width
        np.testing.assert_allclose(pull,row["pulls"],atol=2e-7,rtol=3e-8)
        if name=="historical_search_joint":
            assert 3.0461 < max(abs(pull)) < 3.0462
        if name.startswith("release_down"):
            assert max(abs(np.delete(pull,3))) < 2e-7


def evaluate():
    verify_all_exported_models()
    a = ([5.605448712486393e-7,5.145801393492981e-5,.00662964860481175],
         [-1.0083742649221273e-6,1.3217056207479311e-5,8.291913809274802e-5,
          -.00013042360786900102,-.00019530504214036796,-.0007167805874705575],
         [-2.423993647875405e-6,-3.493304967527641e-6,5.982783540856674e-5,
          8.062516051781105e-5,.0003338404276463674,-.0012020945848813963],
         68.45732325605479,.2846300737073716+.05598839915933038j,106.32429347477591)
    import runpy
    old = runpy.run_path(str(HERE.parent/"1554"/"independent_check.py"))
    x = np.array(old["fractions"](old["POINTS"]["MN_no_RGE_eq15"]),dtype=float)
    historic = (x[:3],x[3:9],x[9:15],x[15],x[16]+1j*x[17],x[18])
    result = {"A_charged_first_family_held_out":reconstruct(a),
              "original_eq15":reconstruct(historic)}
    branch = result["A_charged_first_family_held_out"]
    assert branch["first_family_comparison"]["first_family_double_ratio"] < 6
    assert branch["one_loop_peak_active_yukawa_required_whole_interval"] > 1
    assert branch["finite_check"]["active_weak_norm_loss_upper"] < 1e-12
    return result


def compare(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for key in x:
            compare(x[key],y[key])
    else:
        np.testing.assert_allclose(x,y,atol=0,rtol=3e-8)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    actual = evaluate()
    path = HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(actual,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        compare(actual,load_json(path))
        print("1556: finite six pole errors and one-loop first-family transfer budgets verified")
