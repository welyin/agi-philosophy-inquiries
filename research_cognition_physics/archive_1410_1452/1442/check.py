"""Analytic all-mass bounds with independent positive-mass SVD spot checks."""

from argparse import ArgumentParser
from pathlib import Path
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def interval(values, widths):
    return values[0] - widths*values[1], values[0] + widths*values[-1]


def pmns(s12sq, s13sq, s23sq, delta):
    a,b,c = np.sqrt([s12sq,s13sq,s23sq])
    x,y,z = np.sqrt([1-s12sq,1-s13sq,1-s23sq])
    phase = np.exp(1j*delta)
    return np.array([[x*y,a*y,b/phase],
                     [-a*z-x*c*b*phase,x*z-a*c*b*phase,c*y],
                     [a*c-x*z*b*phase,-x*c-a*z*b*phase,z*y]])


def calculate():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf8"))
    old = json.loads((HERE/cfg["historical_input"]).read_text(encoding="utf8"))
    widths = cfg["half_widths"]
    y = {sector:[interval(row,widths) for row in rows] for sector, rows in old["yukawa"].items()}
    electron,mu,tau = y["charged_leptons"]
    f = mu[0]**2 / (electron[1]*tau[1])
    a13 = interval(old["DayaBay_amplitude"],widths)
    s13lo,s13hi = [(1-math.sqrt(1-v))/2 for v in a13]
    quark = old["quark_sines_and_phase"]
    s13q = interval(quark[1],widths)[1]
    s23q = interval(quark[2],widths)[1]
    q_upper = y["up"][1][1]*y["down"][1][1] / (
        y["up"][2][0]*y["down"][2][0]*(1-s13q**2)*(1-s23q**2))
    out = {"charge_hierarchy_lower":f,"quark_middle_to_largest_upper":q_upper}
    for ordering in ["NO","IO"]:
        osc = old[ordering]
        sol = interval(osc["solar_gap"],widths)
        atm = interval(osc["atmospheric_gap"],widths)
        c12sq = 1-interval(osc["solar_sine_squared"],widths)[1]
        s23hi = interval(osc["atmospheric_sine_squared"],widths)[1]
        if ordering == "NO":
            a = sol[0]/(sol[0]+atm[1])
            electron_light = math.sqrt(c12sq*(1-s13hi))
            r_floor = mu[0]/tau[1]*electron_light*math.sqrt(a)
            geom = f * electron_light**3 * 2*math.sqrt(a*(1-a))
            details = {"solar_to_largest_gap_lower":a}
        else:
            r2 = atm[0]/(atm[0]+sol[1])
            r_floor = mu[0]/tau[1]*math.sqrt(s13lo*r2)
            u = cfg["u_split_IO"]
            tau_light_lo = (1-s23hi)*(1-s13hi)
            norm_scale = 1-tau_light_lo*(1-u*u)+(mu[1]**2+electron[1]**2)/tau[0]**2
            assert norm_scale < 1
            amin = s13lo
            turning = math.sqrt(amin/(2*(1-amin)))
            assert u >= turning
            small = f*r2*(1.5*math.sqrt(3)*amin*math.sqrt(1-amin))/norm_scale**1.5
            large = f*r2*(amin+(1-amin)*u*u)**1.5/u
            geom = min(small,large)
            details = {"r_squared_floor":r2,"small_u_bound":small,"large_u_bound":large,
                       "normalized_tau_Frobenius_upper_small_u":norm_scale}
        # These diagnostics are separate from the analytic continuum certificate.
        smallest = math.inf
        for m in (1e-7, 0.0005, 0.005, 0.016, 0.15, 1.5):
            A, S = osc["atmospheric_gap"][0], osc["solar_gap"][0]
            offsets = [0,S,S+A] if ordering=="NO" else [A,A+S,0]
            for delta in np.linspace(0,2*np.pi,25):
                V=pmns(osc["solar_sine_squared"][0],(1-math.sqrt(1-old["DayaBay_amplitude"][0]))/2,
                       osc["atmospheric_sine_squared"][0],delta)
                K=np.diag([a[0] for a in old["yukawa"]["charged_leptons"]]) @ V @ np.diag(np.sqrt(m*m+np.array(offsets)))
                high,mid,low = np.linalg.svd(K,compute_uv=False)
                measured=mid**2/(high*low)
                assert measured >= geom and mid/high >= r_floor
                smallest=min(smallest,measured)
        out[ordering] = {**details,"geometric_ratio_lower":geom,"middle_to_largest_lower":r_floor,
                         "shape_gap_factor":r_floor/q_upper,"spot_minimum_only":smallest,
                         "necessary_log_change_per_singular_value":math.log(geom)/4}
        assert geom>1 and r_floor>q_upper
    return {k:({kk:round(vv,11) for kk,vv in v.items()} if isinstance(v,dict) else round(v,11))
            for k,v in out.items()}


if __name__=="__main__":
    parser=ArgumentParser()
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument("--write",action="store_true")
    modes.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    dest=HERE/"results.json"
    if args.write:
        with dest.open("x",encoding="utf-8") as file:
            json.dump(result,file,indent=2)
            file.write("\n")
    else:
        actual=json.loads(dest.read_text(encoding="utf8"))
        for key,value in result.items():
            if isinstance(value,dict):
                assert set(value)==set(actual[key])
                assert all(math.isclose(v,actual[key][k],abs_tol=1e-8) for k,v in value.items())
            else:
                assert math.isclose(value,actual[key],abs_tol=1e-10)
    print(json.dumps(result))
