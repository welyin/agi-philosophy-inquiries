"""LEP backward-muon bins and the scalar source blind to muonium conversion."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ALPHA = 1 / 128
SW2 = 0.23122
MZ = 91.1876
GZ = 2.4952
PB = 389379366.0
MH = 125.
HEAVY = 600.
LIGHT_G2 = 8 * math.pi * .00407 * .000044 / MH
ANGLE = .95
# Energy, z endpoints, measured bin average, experimental error, printed SM.
LEP = [(189, -1., -.8, .614, .080, .532), (189, -.8, -.6, .420, .065, .514),
       (207, -1., -.8, .518, .087, .481), (207, -.8, -.6, .496, .075, .456)]


def moments(row, nodes=24):
    energy, zl, zu, _, _, _ = row
    zc, zh = (zl+zu)/2, (zu-zl)/2
    t, w = np.polynomial.legendre.leggauss(nodes)
    z = zc + zh*t
    s = energy**2
    q = s*(1-z)/2
    clr = (4*math.pi*ALPHA/s + 4*math.pi*ALPHA*(SW2-.5) /
           ((1-SW2)*(s-MZ**2+1j*MZ*GZ)))
    f = s*PB*(1-z)**2/(128*math.pi)
    p = 2/(HEAVY**2+q)
    a = np.dot(w/2, f*p**2/4)
    b = np.dot(w/2, f*p*clr.real)
    e = 2*LIGHT_G2*(MH**-2-HEAVY**-2)
    fmax = s*PB*(1-zl)**2/(128*math.pi)
    # Rigorous truncation for integrating (1+u)^-1 and (1+u)^-2.
    radius = s*zh/2/(HEAVY**2+s*(1-zc)/2)
    assert radius < .006 and nodes >= 10
    # degree-10 Taylor times the quadratic angular numerator fits Gauss degree 19.
    # relative error incl. quadrature and exact positive integral, round outward.
    assert 2*13*radius**11*(1+radius)**2/(1-radius)**2 < 1e-18
    return dict(a=float(a), b=float(b), e0=fmax*e*(abs(clr.real)+e/4),
                e1=fmax*e/HEAVY**2, g=float(clr.real), q=q, f=f, weights=w/2)


def kernel(q, x):
    s2 = LIGHT_G2/x
    ai = HEAVY**-2+s2*(MH**-2-HEAVY**-2)
    return LIGHT_G2/(MH**2+q)+(x-LIGHT_G2)/(HEAVY**2+q)+x/(1/ai+q)


def shift(m, x):
    k = kernel(m["q"], x)
    return float(np.dot(m["weights"], m["f"]*(-m["g"]*k+k*k/4)))


def separate_response(ms, size, gain=(.8,1.2)):
    """Allow independent positive gains on the square and on the interference."""
    l, h = gain
    intervals = [(LIGHT_G2,math.inf)]
    for row,m in zip(LEP,ms):
        obs,err,sm = row[3:]
        low,high = obs-sm-size*err, obs-sm+size*err
        # e0 bounds the sum of the constant pieces from both terms.
        less = polynomial_intervals(l*m["a"]*(1-1e-12),h*m["b"]*(1+1e-12),
                                    high+h*m["e0"],True)
        more = polynomial_intervals(h*m["a"]*(1+1e-12),l*m["b"]*(1-1e-12)-h*m["e1"],
                                    low-h*m["e0"],False)
        intervals = intersect(intervals,intersect(less,more))
    witnesses = []
    for lo,hi in intervals:
        x = (lo+hi)/2
        details = []
        for row,m in zip(LEP,ms):
            k = kernel(m["q"],x)
            mag_i = float(np.dot(m["weights"],m["f"]*m["g"]*k))
            sq = float(np.dot(m["weights"],m["f"]*k**2/4))
            lower,upper = l*sq-h*mag_i,h*sq-l*mag_i
            t = max(0.,min(1.,(row[3]-row[5]-lower)/(upper-lower)))
            k_sq,k_int = l+(h-l)*t,h-(h-l)*t
            residual = (row[5]+k_sq*sq-k_int*mag_i-row[3])/row[4]
            assert abs(residual) < size
            details.append(dict(k_interference=k_int,k_square=k_sq,residual_sigma=residual))
        witnesses.append(dict(X=x,bins=details))
    return dict(X_necessary=intervals,
                misalignment_sq_necessary=[(LIGHT_G2/hh,LIGHT_G2/ll) for ll,hh in reversed(intervals)],
                witnesses=witnesses)


def intersect(left, right):
    return [(max(a,c),min(b,d)) for a,b in left for c,d in right if max(a,c) <= min(b,d)]


def polynomial_intervals(a, b, value, less):
    # Set on the real line where a*x^2-b*x <= value (or >= value).
    discr = b*b+4*a*value
    if discr < 0:
        return [] if less else [(-math.inf,math.inf)]
    high = (b+math.sqrt(discr))/(2*a)
    low = -value/(a*high) if high != 0 else 0.
    return [(low,high)] if less else [(-math.inf,low),(high,math.inf)]


def scan(inputs, response, size):
    intervals = [(LIGHT_G2,math.inf)]
    for row,m in zip(LEP,inputs):
        obs,err,sm = row[3:]
        adjusted = [(obs-sm+d*size*err)/r for d in (-1,1) for r in response]
        low, high = min(adjusted),max(adjusted)
        # small relative margins also enclose floating quadrature roundoff
        left = polynomial_intervals(m["a"]*(1-1e-12),m["b"]*(1+1e-12)+m["e1"],
                                    high+m["e0"],True)
        right = polynomial_intervals(m["a"]*(1+1e-12),m["b"]*(1-1e-12)-m["e1"],
                                     low-m["e0"],False)
        intervals = intersect(intervals,intersect(left,right))
    points = []
    for l,h in intervals:
        x = (l+h)/2
        deltas = [shift(m,x) for m in inputs]
        # Choose permitted positive gain in each bin to approach its residual.
        gains = [min(response[1],max(response[0],(row[3]-row[5])/delta)) if delta != 0 else 1.
                 for row,delta in zip(LEP,deltas)]
        std = [(row[5]+gain*delta-row[3])/row[4] for row,gain,delta in zip(LEP,gains,deltas)]
        assert max(map(abs,std)) < size
        points.append(dict(X=x, residuals_in_one_exp_error=std, gain=gains))
    return dict(X_necessary=intervals, y_abs_necessary=[(math.sqrt(l),math.sqrt(h)) for l,h in intervals],
                misalignment_sq_necessary=[(LIGHT_G2/h,LIGHT_G2/l) for l,h in reversed(intervals)],
                witnesses=points)


def run():
    ms = [moments(row) for row in LEP]
    x0 = LIGHT_G2/(1-ANGLE**2)
    ai = HEAVY**-2+(1-ANGLE**2)*(MH**-2-HEAVY**-2)
    rows = []
    for row,m in zip(LEP,ms):
        v = shift(m,x0)
        assert abs(v-shift(moments(row,nodes=48),x0)) < 1e-17
        rows.append(dict(energy_GeV=row[0],cos_theta_minus=list(row[1:3]),
                         measured_pb=row[3],sigma_pb=row[4],published_SM_pb=row[5],
                         born_shift_pb=v, shift_in_experimental_errors=v/row[4],
                         polynomial_A_pb=m["a"],polynomial_B_pb=m["b"],
                         extra_error_constant_pb=m["e0"],extra_error_linear_pb=m["e1"]))
    assert all(r["born_shift_pb"] < 0 for r in rows)
    for x in np.geomspace(LIGHT_G2,10,55):
        for m in ms:
            diff = kernel(m["q"],x)-2*x/(HEAVY**2+m["q"])
            assert min(diff)>-1e-18 and max(diff)<2*LIGHT_G2*(MH**-2-HEAVY**-2)+1e-18
            assert abs(shift(m,x)-(m["a"]*x**2-m["b"]*x)) < m["e0"]+m["e1"]*x+1e-12
    return dict(inputs=dict(alpha=ALPHA,sin2theta_eff=SW2,MZ=MZ,GammaZ=GZ,GeV_minus2_to_pb=PB,
                            gh2_reference=LIGHT_G2,vev_projection=.95,reference_MA_GeV=ai**-.5,
                            R_light_to_measured_upper_reference=True),
                quantities="X=|gh|^2/(1-c^2); necessary rectangles, not confidence limits",
                CP_even_odd_contact_residual_GeV_minus2=LIGHT_G2/MH**2+(x0-LIGHT_G2)/HEAVY**2-x0*ai,
                bins=rows, two_marginal_errors={
                    "born_increment":scan(ms,(1,1),2),
                    "independent_bin_gains_on_net_0.8_to_1.2":scan(ms,(.8,1.2),2),
                    "separate_interference_and_square_gains_0.8_to_1.2":separate_response(ms,2)})


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    out=json.dumps(run(),indent=2,ensure_ascii=False,allow_nan=False)+"\n"
    args=ap.parse_args()
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8"))==json.loads(out)
        print("1794 independent OK",hashlib.sha256(out.encode()).hexdigest()[:12])
    else:
        path.write_text(out,encoding="utf-8")
