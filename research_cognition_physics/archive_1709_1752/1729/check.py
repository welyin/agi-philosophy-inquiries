"""Transport the actual three-family matrix pencil and compare its simple root with Vud."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2


def read(path):
    return json.loads((HERE/path).read_bytes())


def pencil(y):
    return np.linalg.solve(y[1,:,:3].T, y[2,:,:3]).T


def projectors(y,dy=None):
    data = []
    for i in (0,1):
        Y = y[i]
        val,vec = np.linalg.eigh(Y @ Y.conj().T)
        light = vec[:,0]
        P = np.outer(light,light.conj())
        dP = np.zeros_like(P)
        if dy is not None:
            dH = dy[i]@Y.conj().T + Y@dy[i].conj().T
            dv = vec[:,1:]@((vec[:,1:].conj().T@dH@light)/(val[0]-val[1:]))
            tmp = np.outer(dv,light.conj())
            dP = tmp+tmp.conj().T
        data.append((P,dP))
    (u,du),(d,dd) = data
    p = np.trace(u @ d).real
    dp = np.trace(du @ d + u @ dd).real
    return math.sqrt(p),dp/(2*p)


def radial(y,z):
    roots = np.linalg.eigvals(pencil(y))
    j = np.argmin(abs(roots-z))
    return roots[j],np.delete(roots,j)


def advance(flow,y,g,dt):
    # Single MSSM segment with the gauge spectators active at 10 PeV.
    bs = np.array([33/5,1.,-3.])+sum(w for _,w in flow["physics"]["thresholds"])
    gg = [1/(1/g-k*dt*bs/LOOP) for k in (0,1,2)]
    h = dt/LOOP
    beta = flow["beta"]
    k1 = beta(y,gg[0],True)
    k2 = beta(y+h*k1/2,gg[1],True)
    k3 = beta(y+h*k2/2,gg[1],True)
    k4 = beta(y+h*k3,gg[2],True)
    return y+h*(k1+2*k2+2*k3+k4)/6,gg[2]


def transport(flow,y,g,dt,step):
    count = max(1,math.ceil(abs(dt)/step))
    y,g = y.copy(),g.copy()
    for _ in range(count):
        y,g=advance(flow,y,g,dt/count)
    return y


def run():
    cfg = read("inputs.json")
    mod = runpy.run_path(str(HERE/cfg["frame"]))
    flow,rep = mod["charged_background"](.1)
    unpack = mod["unpack"]
    y = np.zeros((4,3,4),complex)
    y[0,:,:3] = np.diag(rep["up_texture_diag_10PeV"])
    y[1,:,:3] = unpack(rep["down_texture_10PeV"])
    y[2,:,:3] = unpack(rep["charged_texture_10PeV"])
    y[3] = unpack(rep["neutral_columns_10PeV"])
    g = flow["segments"](.1)[-1][3][-1]
    high = cfg["high_scale_GeV"]
    assert high==flow["END"] and max(m for m,_ in flow["physics"]["thresholds"])<1e6
    E,D = y[2,:,:3].T,y[1,:,:3]
    z,z0 = E[2,2]/D[2,2],E[0,1]/D[0,1]
    ref_phase=z/abs(z)
    B=pencil(y)
    dy=flow["beta"](y,g,True)/LOOP
    dB=np.linalg.solve(D.T,(dy[2,:,:3].T-B@dy[1,:,:3]).T).T
    def polynomial_rate(w):
        mat=w*np.eye(3)-B
        total=0j
        for j in range(3):
            changed=mat.copy()
            changed[:,j]=-dB[:,j]
            total+=np.linalg.det(changed)
        return total
    zrate = -polynomial_rate(z)/(z-z0)**2
    opening = -polynomial_rate(z0)/(z0-z)
    vud,dlog_vud=projectors(y,dy)
    h=.001
    split_test = {}
    for sign in (-1,1):
        yy=advance(flow,y,g,sign*h)[0]
        root,pair=radial(yy,z)
        np.testing.assert_allclose((root-z)/(sign*h),zrate,rtol=.001)
        # (lambda_1-lambda_2)^2/4 cancels the common first-order branch motion.
        split_test[str(sign)] = ((pair[0]-pair[1])**2/(4*sign*h))/opening
    assert abs(opening)>1e-13
    def sample(dt,step=.008):
        yy=transport(flow,y,g,dt,step)
        root,pair=radial(yy,z)
        v=projectors(yy)[0]
        error=root/(v*ref_phase)
        ee,dd=yy[2,:,:3].T,yy[1,:,:3]
        matrix=ee-v*ref_phase*dd
        left,sv,right=np.linalg.svd(matrix)
        closest=matrix-sv[-1]*np.outer(left[:,-1],right[-1])
        assert np.linalg.svd(closest,compute_uv=False)[-1]<1e-15
        return dict(scale_GeV=high*math.exp(dt),Vud=v,root_abs=abs(root),
                    relative_modulus_residual=abs(error)-1,
                    phase_residual_degrees=math.degrees(np.angle(error)),
                    minimum_additive_change_on_E_for_exact_two_arms=sv[-1],
                    minimum_relative_change_over_lepton_operator_norm=sv[-1]/np.linalg.norm(ee,2),
                    root_pair_separation=abs(pair[0]-pair[1]))
    lo,hi=0.,math.log(cfg["nearby_upper_GeV"]/high)
    assert sample(lo)["relative_modulus_residual"]<0<sample(hi)["relative_modulus_residual"]
    for _ in range(44):
        mid=(lo+hi)/2
        if sample(mid)["relative_modulus_residual"]<0:
            lo=mid
        else:
            hi=mid
    dt=(lo+hi)/2
    match=sample(dt,.002)
    assert abs(match["relative_modulus_residual"])<3e-11
    def comp(x):
        return [float(x.real),float(x.imag)]
    return dict(high_scale_GeV=high,Vud_high=vud,relative_root_high=comp(z/(vud*ref_phase)),
                simple_root_dlog_abs_per_log_scale=(zrate/z).real,
                Vud_dlog_abs_per_log_scale=dlog_vud,
                simple_root_phase_degrees_per_log_scale=math.degrees((zrate/z).imag),
                simple_root_universal_gauge_dlog=(16*g[2]/3-4*g[0]/3)/LOOP,
                double_root_opening_coefficient=comp(opening),
                double_root_first_order_ratio_at_plus_minus_one_per_mille={k:comp(v) for k,v in split_test.items()},
                original_scale=sample(0),down_one_decade=sample(-math.log(10)),nearby_modulus_match=match,
                old_flow_at_first_alternative_weak_pole=sample(math.log(
                    read('../1720/results.json')['central_conditional_combined_point']['four_weak_poles_GeV'][0]/high)),
                close_step_stability_relative_modulus=abs(sample(dt)["relative_modulus_residual"]-match["relative_modulus_residual"]))


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda x:float(f"{float(x):.10g}"))
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
        print("1729 verified covariant pencil derivatives, CKM projector running and a finite modulus crossing.")
    else:
        with path.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
