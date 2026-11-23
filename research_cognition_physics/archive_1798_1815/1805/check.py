"""Common-weight intertwining, canonical Higgs moments and one-scale SM comparisons."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
L = 16 * math.pi**2
BETAS = np.array([4.1, -19/6, -7.])
GAUGE_DAMPING = np.array([[.85, 2.25, 8.], [.25, 2.25, 8.], [2.25, 2.25, 0.]])
MULT = np.repeat([3., 3., 1.], 3)
PAD = 1e-12


def read(path):
    return json.loads((HERE/path).read_bytes())


def mixing(angles, phase):
    a,b,c = np.sin(angles)
    A,B,C = np.cos(angles)
    p = np.exp(1j*phase)
    return np.array([[A*B, a*B, b/p],
                     [-a*C-A*c*b*p, A*C-a*c*b*p, c*B],
                     [a*c-A*C*b*p, -A*c-a*C*b*p, C*B]])


def stability(matrix):
    P = abs(matrix)**2
    np.testing.assert_allclose(P.sum(axis=0),1.,atol=1e-12)
    gap = 1-np.linalg.svd(P, compute_uv=False)[1]**2
    return dict(minimum_squared_overlap=float(P.min()), dephasing_squared_gap=float(gap),
                deviation_over_readout_error=float(1/math.sqrt(gap)))


def gauges(g, t):
    return 1 / (1/g**2 - 2*BETAS*t/L)


def central_run(tmax, g, charged, v, lam, steps):
    # Gram matrices directly evolve the physical left-port observables.
    qu = np.diag(charged[:3]**2)
    qd = (v * charged[3:6]**2) @ v.conj().T
    le = np.diag(charged[6:]**2)
    z = np.r_[qu.ravel(), qd.ravel(), le.ravel(), lam].astype(complex)
    def beta(t,z):
        matrices = z[:27].reshape(3,3,3)
        sq = matrices @ matrices
        tr = np.trace(matrices,axis1=1,axis2=2).real
        tr4 = np.trace(sq,axis1=1,axis2=2).real
        T, B = np.dot([3,3,1],tr), np.dot([3,3,1],tr4)
        gs = gauges(g,t)
        damp = GAUGE_DAMPING @ gs
        cross = matrices[0] @ matrices[1] + matrices[1] @ matrices[0]
        dgram = 3*sq + 2*(T-damp)[:,None,None]*matrices
        dgram[:2] -= 1.5*cross
        val = z[-1].real
        weak, hyper = gs[1], .6*gs[0]
        dl = 24*val**2 + (4*T-9*weak-3*hyper)*val - 2*B
        dl += 3/8*(2*weak**2+(weak+hyper)**2)
        return np.r_[dgram.ravel(),dl]/L
    h = tmax/steps
    for j in range(steps):
        t = h*j
        a = beta(t,z)
        b = beta(t+h/2,z+h*a/2)
        c = beta(t+h/2,z+h*b/2)
        d = beta(t+h,z+h*c)
        z += h*(a+2*b+2*c+d)/6
    squared = np.linalg.eigvalsh(z[:27].reshape(3,3,3)).ravel()
    return squared, float(z[-1].real)


def gauge_interval(glo, ghi, lo, hi):
    ends = np.array([gauges(g,t) for g in (glo,ghi) for t in (lo,hi)])
    return ends.min(axis=0)-PAD, ends.max(axis=0)+PAD


def times_product(l1,h1,l2,h2):
    vals = [l1*l2,l1*h2,h1*l2,h1*h2]
    return min(vals),max(vals)


def box(y,dy,g,dg,lam,dlam,cfg):
    """Scalar differential inequalities for all flavour orientations; no corner interpolation."""
    gl, gh = g-dg, g+dg
    xl, xh = (y-dy)**2, (y+dy)**2
    ll, lh = lam-dlam, lam+dlam
    lo_t, hi_t = cfg["gauge_time_outer"]
    ratio_lo, ratio_hi = cfg["gauge_ratio_outer"]
    h = cfg["log_step"]
    count = round(hi_t/h)
    assert abs(count*h-hi_t)<1e-13
    intervals = []
    for j in range(count):
        ql, qh = gauge_interval(gl,gh,j*h,(j+1)*h)
        # Bootstrap |d ln y_i^2 /dt| <= 1 and |d lambda/dt| <= .16 on every cell.
        low, high = xl*math.exp(-h)*(1-PAD), xh*math.exp(h)*(1+PAD)
        Tlo, Thi = MULT@low, MULT@high
        max_peer = np.repeat([high[3:6].max(),high[:3].max(),0.],3)
        min_peer = np.repeat([low[3:6].min(),low[:3].min(),0.],3)
        rateslo = 2*(Tlo+1.5*(low-max_peer)-np.repeat(GAUGE_DAMPING@qh,3))/L-PAD
        rateshi = 2*(Thi+1.5*(high-min_peer)-np.repeat(GAUGE_DAMPING@ql,3))/L+PAD
        assert min(rateslo)>-.99 and max(rateshi)<.99
        tlo, thi = ll-.16*h-PAD, lh+.16*h+PAD
        assert thi < .25  # The same cells enter the Weinberg-operator comparison.
        wlo, whi, ulo, uhi = ql[1],qh[1],.6*ql[0],.6*qh[0]
        linear_lo,linear_hi = times_product(4*Tlo-9*whi-3*uhi,4*Thi-9*wlo-3*ulo,tlo,thi)
        min_sq = 0 if tlo<=0<=thi else min(tlo*tlo,thi*thi)
        beta_lo = (24*min_sq+linear_lo + .375*(2*wlo**2+(wlo+ulo)**2)-2*MULT@(high**2))/L-PAD
        beta_hi = (24*max(tlo*tlo,thi*thi)+linear_hi + .375*(2*whi**2+(whi+uhi)**2)-2*MULT@(low**2))/L+PAD
        assert max(abs(beta_lo), abs(beta_hi)) < .159
        if (j+1)*h>=lo_t-PAD:
            A_lo, A_hi = 3*sum(low[:6])+ratio_lo*sum(low[6:]),3*sum(high[:6])+ratio_hi*sum(high[6:])
            B_lo = 3*sum(low[:6]**2)+ratio_lo*sum(low[6:]**2)
            S_lo, S_hi = (3+ratio_lo)*wlo,(3+ratio_hi)*whi
            deficit = S_lo-A_hi
            assert deficit>0
            min_lambda = (B_lo+deficit**2/(3*ratio_hi))/S_hi
            intervals.append([tlo,thi,A_lo,A_hi,B_lo,S_lo,S_hi,deficit,min_lambda,min_lambda-thi,
                              B_lo/S_hi-thi])
        xl, xh = xl*np.exp(h*rateslo)*(1-PAD), xh*np.exp(h*rateshi)*(1+PAD)
        ll, lh = ll+h*beta_lo-PAD, lh+h*beta_hi+PAD
    values=np.array(intervals)
    return dict(lambda_data=[float(values[:,0].min()),float(values[:,1].max())],
                A_ch=[float(values[:,2].min()),float(values[:,3].max())],
                B_ch_lower=float(values[:,4].min()),
                S=[float(values[:,5].min()),float(values[:,6].max())],
                DeltaA_lower=float(values[:,7].min()),
                rank3_lambda_lower=float(values[:,8].min()),
                needed_lambda_upward_shift_if_other_matching_fixed=float(values[:,9].min()),
                charged_only_fourth_moment_gap_even_for_unbounded_new_rank=float(values[:,10].min()),
                time_cells_with_possible_matching=len(intervals), all_transport_steps=count)


def calculate():
    cfg = read("inputs.json")
    weak, gauge, quartic = (read(cfg[k]) for k in ("charged","gauge","quartic"))
    match=read(cfg["gauge_match"])
    for narrow,outer in [(match["one_loop_corner_bounds"]["t"],cfg["gauge_time_outer"]),
                         (match["one_loop_corner_bounds"]["wl_wq"],cfg["gauge_ratio_outer"])]:
        assert outer[0] <= narrow[0] and narrow[1] <= outer[1]
    g,dg = np.array([gauge[k] for k in ("g1_GUT","g2","g3")]).T
    y,dy = np.concatenate([weak["yukawa"][k] for k in ("up","down","charged_leptons")]).T
    lam,dlam = np.array(quartic["lambda_paper"])/quartic["lambda_paper_factor"]
    ang = np.array(weak["quark_sines_and_phase"])[:,0]
    v = mixing(np.arcsin(ang[:3]),ang[3])
    t = match["one_loop_central"]["t"]
    omega = match["one_loop_central"]["wl_wq"]
    sq, lambda_data = central_run(t,g,y,v,lam,cfg["central_steps"])
    fine, lf = central_run(t,g,y,v,lam,2*cfg["central_steps"])
    np.testing.assert_allclose(fine,sq,atol=1e-10,rtol=1e-9)
    assert abs(lf-lambda_data)<1e-10
    g2 = gauges(g,t)[1]
    weights = np.repeat([3.,3.,omega],3)
    A, B, S = weights@sq, weights@(sq**2),(3+omega)*g2
    missing = S-A
    central = dict(omega=omega,g2_squared=float(g2),S=float(S),A_ch=float(A),B_ch=float(B),
                   DeltaA=float(missing),lambda_data=lambda_data,
                   critical_lambda_rank={str(k):float((B+missing**2/(k*omega))/S) for k in [1,2,3]},
                   even_unlimited_rank_lower_lambda=float(B/S),
                   pure_charged_canonical_ratio=float(A/S),low_to_high_yukawa_squared=sq.tolist(),
                   rk4_step_halving_lambda_difference=abs(lf-lambda_data))
    irrep={"quark_MZ":stability(v)}
    match10=read(cfg["high_flavour"])
    gamma=math.radians(match10["gamma_10PeV_deg"])
    angles=np.array(match10["ckm_angles_12_13_23_rad"])
    delta=gamma+math.asin(math.sin(gamma)*math.cos(angles[0])*math.tan(angles[2])*math.sin(angles[1])/math.sin(angles[0]))
    irrep["quark_10PeV_registered_two_doublet_endpoint"]=stability(mixing(angles,delta))
    for key, data in read(cfg["mixed_points"])["same_two_tensor_mixed_type_I_and_II"].items():
        _,a13,s12,s23=data["low_ratio_A13_solar_atmospheric"]
        s13=a13/(2*(1+math.sqrt(1-a13)))
        irrep[key+"_low_lepton"]=stability(mixing(np.arcsin(np.sqrt([s12,s13,s23])),math.radians(data["low"]["delta_deg"])))
    return dict(commutation=irrep,central=central,outer=box(y,dy,g,dg,lam,dlam,cfg))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    raw=json.dumps(result,indent=2,allow_nan=False)+"\n"
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
    else:
        path.write_bytes(raw.encode())
    print("1805",hashlib.sha256(raw.encode()).hexdigest()[:12])

