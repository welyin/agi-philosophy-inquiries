"""On-shell W and baryonic flavour inputs identify dispersive B mixing magnitudes."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


CFG = read("inputs.json")
PORT = read(CFG["angles_and_relative_tree"])
MIX = read(CFG["mixing"])["data"]
ABS = read(CFG["absolute"])["mixing_coefficient_ps_inv"]
RAD = math.pi/180


def envelope(x,n):
    e = math.hypot(*x[1:])
    return np.array([x[0]-n*e,x[0]+n*e])


def angle_inputs(n,mirror):
    s = envelope(MIX["golden_S"],n)
    beta = .5*np.arcsin(s)
    if mirror:
        beta = math.pi/2-beta[::-1]
    tree = PORT["time_integrated_gamma_deg"]
    gamma = np.array([tree["central"]]*2 if n==0 else tree["68" if n==1 else "95"])*RAD
    return beta,gamma


def flavour(w,r,beta,gamma):
    c,x = w*w,(r*w)**2
    t = 1-c-x
    a = (math.sin(beta)/math.sin(gamma))**2
    b = (math.sin(beta+gamma)/math.sin(gamma))**2
    f = (1+r*r)*(a/r**2+b)-1
    y = t*f                       # |Vts/Vtd|^2
    total = t*c*(1+r*r)
    pd, ps = total/(1+y), total*y/(1+y)
    np.testing.assert_allclose(a*pd/x+b*pd/c+pd/t,1.,rtol=4e-13)
    # Rebuild both complex orthogonal columns independently of F.
    products = math.sqrt(pd)*np.array([math.sqrt(a)*np.exp(1j*gamma),
                                       -math.sqrt(b),np.exp(-1j*beta)])
    colb = np.sqrt([x,c,t])
    cold = products/colb
    cols = np.cross(cold.conj(),colb)
    matrix = np.column_stack([cold,cols,colb])
    np.testing.assert_allclose(matrix.conj().T@matrix,np.eye(3),atol=2e-13)
    np.testing.assert_allclose(abs(matrix[2,1]*matrix[2,2])**2,ps,rtol=2e-13)
    return dict(top_bottom=t,ts_over_td_squared=y,pd_squared=pd,ps_squared=ps,
                Vus=abs(matrix[0,1]),Vub=r*w,Vcb=w)


def observed(n,xi_name):
    md,ms = (envelope(MIX[k],n) for k in ("delta_md","delta_ms"))
    bd,bs = (envelope(MIX[k],n) for k in ("mass_d","mass_s"))
    xi = envelope(MIX[xi_name],n)
    q2 = xi**2*md/ms[::-1]*bs/bd[::-1]
    return q2,ms


def ratio_port(name,n):
    center,neg,pos = PORT["ratio_ports"][name]
    return np.array([center-n*math.hypot(*neg),center+n*math.hypot(*pos)])


def root_r(w,b,g,q2):
    lo,hi=.025,.45
    for _ in range(64):
        mid=(lo+hi)/2
        y=flavour(w,mid,b,g)["ts_over_td_squared"]*q2
        if y>1:
            lo=mid
        else:
            hi=mid
    val=(lo+hi)/2
    assert abs(flavour(w,val,b,g)["ts_over_td_squared"]*q2-1)<1.e-11
    return val


def evaluation(name,n,mirror,ratio="Lambda_b_original",xi_name="HPQCD",stress=0.):
    wc,wl,wh = CFG["real_W_cb"][name]
    w = np.array([(wc-n*wl)*(1-stress),(wc+n*wh)*(1+stress)])
    r = ratio_port(ratio,n)
    b,g = angle_inputs(n,mirror)
    q2,ms = observed(n,xi_name)
    cs = envelope(ABS,n)
    # Rigorous tF extrema for the observed range (see derivative certificate in proof.md).
    low = flavour(w[1],r[1],b[0],g[1])
    high = flavour(w[0],r[0],b[1],g[0])
    y = np.array([low["ts_over_td_squared"],high["ts_over_td_squared"]])
    k = w**2*(1+r**2)
    assert 0<k[0]<=k[1]<.5
    total = k*(1-k)
    ps = total*y/(1+y)
    hs = ms/cs[::-1]/ps[::-1]    # an outer box, since ps numerator and y reuse w and r
    xq = q2*y
    row = dict(hd_over_hs_sharp=xq,hs_outer=hs,VtsVtb_squared_outer=ps,
               hd_eq_hs_in_ratio_band=bool(xq[0]<=1<=xq[1]))
    # Z = [sin²(beta+gamma) + (w²/t) sin²(gamma)] / sin²(beta).
    # Work with a rigorous trigonometric envelope even if beta+gamma crosses pi/2.
    end_sines = np.sin(b+g)**2
    bsum_max = 1 if b[0]+g[0]<=math.pi/2<=b[1]+g[1] else max(end_sines)
    sine_sum = np.array([min(end_sines),bsum_max])
    zz = (sine_sum+(w**2/(1-k))*(np.sin(g)**2))/(np.sin(b[::-1])**2)
    ru = r**2*zz
    us = np.sqrt(np.maximum(0,ru/(1+ru)-(r*w)[::-1]**2))
    row["Vus_outer_without_either_clock_or_absolute_lattice"]=us
    kaon = envelope(read(CFG["kaon"])["holdout_Kl3"],n)
    row["kaon_Vus_measurement"]=kaon
    row["kaon_intersects_Vus_outer"]=bool(max(us[0],kaon[0])<=min(us[1],kaon[1]))
    if ratio=="Lambda_b_original":
        # r grows with q2 and beta and decreases with gamma and w.
        row["required_r_if_hd_eq_hs"]=[root_r(w[1],b[0],g[1],q2[0]),
                                        root_r(w[0],b[1],g[0],q2[1])]
        h_common = ms/cs[::-1]*(1+q2)/total[::-1]
        row["common_h_if_absent_phase_used_outer"]=h_common
    return row


def calculate():
    primary=next(iter(CFG["real_W_cb"]))
    data = {}
    for experiment, vals in CFG["real_W_cb"].items():
        cases = {}
        for mirror in (False,True):
            cases["reflected_beta" if mirror else "small_beta"] = {
                str(n):evaluation(experiment,n,mirror) for n in (0,1,2)}
        data[experiment]=cases
    ports = {name:evaluation(primary,2,False,name)["hd_over_hs_sharp"]
             for name in PORT["ratio_ports"]}
    variation = {label:evaluation(primary,2,mirror,xi_name=label,stress=CFG["finite_range_transfer"])
                 for label,mirror in [("HPQCD",True),("FNAL_MILC",True)]}
    material={}
    for tag,row in read(CFG["six"])["refits"].items():
        m=row["moduli"]
        r=m["Vub"]/m["Vcb"]
        w=m["Vcb"]
        pred=flavour(w,r,.5*math.asin(row["S"]),row["gamma"]*RAD)
        np.testing.assert_allclose([pred["ps_squared"],1/math.sqrt(pred["ts_over_td_squared"])],
                                   [row["tbs_squared"],row["td_ts"]],rtol=2e-8)
        material[tag]=dict(w=w,r=r,
            Vus=m["Vus"],
            kaon_signed_pull=(m["Vus"]-read(CFG["kaon"])["holdout_Kl3"][0])/
                math.hypot(*read(CFG["kaon"])["holdout_Kl3"][1:]),
            W_signed_pulls={k:(w-v[0])/(v[2] if w>v[0] else v[1])
                            for k,v in CFG["real_W_cb"].items()},
            baryon_signed_pull=(r-.083)/math.hypot(.004,.004))
    ratio=.0363
    # On a fixed CKM shape and q, h scales as 1/[w^2*(1-k)].
    # Require symmetric fractional intervals on w for two hypotheses to be disjoint.
    w0=.04252
    r=.083
    k0=w0*w0*(1+r*r)
    k1=2*k0*(1-k0)/(1+ratio)/(1+math.sqrt(1-4*k0*(1-k0)/(1+ratio)))
    w1=math.sqrt(k1/(1+r*r))
    separation=(w0-w1)/(w0+w1)
    # Is there an actual small-beta weak matrix satisfying the independent kaon
    # central value and both B clocks with common mixing? Fix w and the angles
    # in advance, use r to calibrate the kaon, and record which lattice xi is needed.
    beta,gamma=(a[0] for a in angle_inputs(0,False))
    kaon=read(CFG["kaon"])["holdout_Kl3"][0]
    low,high=.07,.10
    for _ in range(64):
        mid=.5*(low+high)
        if flavour(w0,mid,beta,gamma)["Vus"]<kaon:
            low=mid
        else:
            high=mid
    joint=flavour(w0,.5*(low+high),beta,gamma)
    assert abs(joint["Vus"]-kaon)<1.e-13
    q2,ms=observed(0,"HPQCD")
    xi=MIX["HPQCD"][0]/math.sqrt(q2[0]*joint["ts_over_td_squared"])
    hs=ms[0]/(ABS[0]*joint["ps_squared"])
    joint_point=dict(true_weak_matrix=joint,r=.5*(low+high),required_xi=xi,
                     required_common_h=hs,
                     xi_pulls_each_lattice={label:(xi-MIX[label][0])/math.hypot(*MIX[label][1:])
                                            for label in ("HPQCD","FNAL_MILC")})
    return dict(on_shell_and_baryon=data,source_sensitivity_small_beta_two_widths=ports,
                alternate_xi_and_one_percent_W_transfer_reflected_two_widths=variation,
                six_registered_charged_materials=material,kaon_calibrated_weak_witness=joint_point,
                central_half_width_fraction_for_disjoint_hypotheses_1_and_1p0363=separation,
                width_floor_fraction_from_Cs_alone=ABS[1]/ABS[0])


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):
        return [clean(v) for v in x]
    if isinstance(x,bool):
        return x
    return float(f"{x:.10g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    result=clean(calculate())
    output=HERE/"results.json"
    if check:
        assert result==json.loads(output.read_bytes())
        print("1772 OK: independent W and baryonic inputs, exact two-B flavour ratio and six CKM matrices")
    else:
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")

