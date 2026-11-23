"""Independent one-loop shooting check in the top-dominance SM, ordinary gY."""
import json
import math
import sys
from pathlib import Path

OUT = Path(__file__).with_name("independent_results.json")
M_TOP = 173.34
LAMBDA_REF = 0.12604
TOP_REF = 0.93690
GAUGE = (0.35830, 0.64779, 1.16660)
BETA = (41/6, -19/6, -7)
LOOP = 16 * math.pi**2

def rhs(t, state):
    u, lam, sens = state  # u = y_t squared; sens = derivative wrt low-scale quartic
    gy, gw, gs = (a*a/(1-2*b*a*a*t/LOOP) for a,b in zip(GAUGE, BETA))
    du = u * (9*u - (17/6)*gy - (9/2)*gw - 16*gs)
    coeff = 12*u - 9*gw - 3*gy
    dl = 24*lam**2 + coeff*lam - 6*u**2 + 3/8*(2*gw**2+(gw+gy)**2)
    return (du/LOOP, dl/LOOP, (48*lam+coeff)*sens/LOOP)

def plus(y, k, h):
    return tuple(a+h*b for a,b in zip(y,k))

def run(scale, yt=TOP_REF, lam=LAMBDA_REF, steps=128):
    width = math.log(scale/M_TOP)/steps
    y = (yt**2, lam, 1.)
    for i in range(steps):
        t = i*width
        a = rhs(t, y)
        b = rhs(t+width/2, plus(y,a,width/2))
        c = rhs(t+width/2, plus(y,b,width/2))
        d = rhs(t+width, plus(y,c,width))
        y = tuple(v + width*(p+2*q+2*r+s)/6 for v,p,q,r,s in zip(y,a,b,c,d))
    return y

def bisect(func, lo, hi):
    left = func(lo)
    assert left * func(hi) < 0
    for _ in range(45):
        middle = (lo+hi)/2
        if func(middle)*left > 0:
            lo = middle
        else:
            hi = middle
    return (lo+hi)/2

def evidence():
    rows = {}
    for scale in (1e6, 3e6, 1e7):
        _, lam, variation = run(scale)
        _, fine, fine_variation = run(scale, steps=256)
        assert abs(fine-lam)<3e-10 and abs(fine_variation-variation)<3e-9
        step = 1e-5
        diff = (run(scale,lam=LAMBDA_REF+step)[1]-run(scale,lam=LAMBDA_REF-step)[1])/(2*step)
        assert abs(diff-variation)<3e-8
        base = bisect(lambda q: run(scale,lam=q)[1], 0.04, 0.15)
        top = bisect(lambda y: run(scale,yt=y)[1], TOP_REF, 1.2)
        assert abs(run(scale,lam=base,steps=256)[1])<3e-10
        assert abs(run(scale,yt=top,steps=256)[1])<8e-10
        rows[str(int(scale))] = dict(high_quartic=round(fine,9),
            d_high_d_low_quartic=round(fine_variation,8),
            zero_boundary_low_quartic=round(base,9),
            zero_boundary_low_y_top=round(top,9),
            required_relative_top_increase=round(top/TOP_REF-1,8))
    return dict(boundary_scale_GeV=M_TOP, input_lambda=LAMBDA_REF,
        input_y_top=TOP_REF, gauge_ordinary_Y_2_3=GAUGE, matched_scale_GeV=rows,
        source="https://arxiv.org/abs/1307.3536",
        scope="Historical MS-bar matched benchmark; one-loop SM, only top Yukawa retained. These are finite-scale response tests, not precision pole predictions.")

if __name__ == "__main__":
    data = json.dumps(evidence(), indent=2, ensure_ascii=False) + "\n"
    if sys.argv[1:] == ["--write-new"]:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(data)
    elif sys.argv[1:] in ([], ["--check"]):
        assert json.loads(OUT.read_text(encoding="utf-8")) == json.loads(data)
        print("1493 independent quartic running, inverse constraints and Jacobian verified.")
    else:
        raise SystemExit("Expected --check or --write-new.")
