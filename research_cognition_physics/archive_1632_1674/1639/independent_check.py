"""Independent one-dimensional outer projection of the fixed-texture two-width task."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(n, fn="inputs.json"):
    return json.loads((HERE.parent/str(n)/fn).read_text(encoding="utf-8-sig"))


def interval(a, n=2):
    return a[0]-n*sum(a[1:]),a[0]+n*sum(a[1:])


def from_matrix(raw):
    m=np.asarray(raw)
    u,s,vh=np.linalg.svd(m[:,:,0]+1j*m[:,:,1])
    order=np.argsort(s)
    u=u[:,order]
    prod=lambda row: row[1].conjugate()*row[2]
    return abs(prod(u[2])),np.angle((prod(u[2])/prod(u[1]))**2)


def solve(md, label, bins=40000):
    weak=read(1404); nu=read(1406); mix=read(1468); hp=read(1407); rare=read(1464)
    prod,phi0=from_matrix(read(1635,"results.json")[label]["weak_basis_mass_GeV"])
    # Uniform over |F| <= 0.0027; derived in the review from K=V sqrt(G).
    poles=read(1633,"results.json")
    rb=math.sqrt(poles["right_doublet_probability"])
    sb2=2*poles["left_difference"]/poles["k"]-2*(4.8*rb/md)**2
    # First use measured charged-current caps and column closure, independently of the SVD.
    caps=read(1632)["charged_current_caps"]
    limits_nu=nu["BelleII"]; rnu=rare["HPQCD_tableI_eq15"]
    cbefore=rnu["Xt_ref"]/(math.pi*rnu["sine_squared_ref"]*rnu["alpha_inverse_ref"])
    abs_n_before=cbefore*rnu["lambda_ref"]*math.sqrt(
        (limits_nu["mu"]+2*limits_nu["width"])/(1-2*rnu["hadronic_rate_width"]))
    f_initial=abs_n_before+cbefore*(math.hypot(*caps.values())+.1*4.8/md**2)
    assert f_initial < .0027
    f0=.00081; fall=.0027
    delta_G=(.1/(md*rb))**2*(fall**2+f0**2)+.1*4.8/md**2*(fall+f0)
    gram_min=1-sb2-(.1*fall/(rb*md))**2-(4.8*rb/md)**2
    delta_K=delta_G/(2*math.sqrt(gram_min))
    assert 2*delta_K < 1e-8
    dl=1e-8
    dphase=2e-6
    lr=[prod-dl,prod+dl]
    c0=rare["HPQCD_tableI_eq15"]
    c=c0["Xt_ref"]/(math.pi*c0["sine_squared_ref"]*c0["alpha_inverse_ref"])
    ell=np.multiply(lr,c); nref=c*c0["lambda_ref"]
    belle=nu["BelleII"]; bdata=interval([belle["mu"],belle["width"]])
    nmin=nref*math.sqrt(bdata[0]/(1+2*c0["hadronic_rate_width"]))/ell[1]
    nmax=nref*math.sqrt(bdata[1]/(1-2*c0["hadronic_rate_width"]))/ell[0]
    cons=weak["const_GeV"]; gf,mw,mes=cons["GF"],cons["MW"],cons["M_Bs"]
    pre=2*mes/(3*hp["hbar_GeV_ps"]); nloop=gf**2*mw**2/(4*math.pi**2)
    frange=interval(rare["fBs_GeV"]); th=interval(weak["tau_ps"]["s_heavy"])
    fbhat=interval(hp["HPQCD_inputs"]["fBs_sqrtB_GeV"]); bag=interval(weak["bag"]["s"])
    hb=hp["hbar_GeV_ps"]; mu=cons["m_mu"]
    h=gf**2*mu**2*mes*math.sqrt(1-4*(mu/mes)**2)/(16*math.pi*hb)
    hmu=[h*frange[0]**2*th[0]*(1-rare["mu_rate_transfer_width"]),
         h*frange[1]**2*th[1]*(1+rare["mu_rate_transfer_width"])]
    cms=nu["CMS"]; cmsmin=cms["branching_in_1e-9"]-2*sum(cms["lower_widths"])
    r3=weak["LHCb_run3_2024_25"]["s_in_1e9"]; r3max=r3[0]+2*(r3[2]+r3[4])
    residual_h=(mes/125.1)**2*.0027
    pmin=math.sqrt((cmsmin*1e-9/hmu[1])-residual_h**2)/ell[1]
    pmax=math.sqrt(r3max*1e-9/(hmu[0]*rare["run3_working_light_response"]))/ell[0]
    sr=weak["short_distance"]; eb=interval(sr["etaB"]); s0=interval(sr["S0_published"])
    xt=sr["xt"][0]; mt=sr["mtMS"][0]; mz=91.1876; az=mix["alpha_s_ref"]
    at=az/(1+az*(23/3)*math.log(mt/mz)/(2*math.pi)); am=at/(1+at*7*math.log(md/mt)/(2*math.pi))
    qcd=(at/az)**(6/23)*(am/at)**(2/7)
    ms,mb=read(1635)["witness"]["light_masses_GeV"][1:]
    t=2*ms*mb/md**2; Lf=(1-ms**2/md**2)*(1-mb**2/md**2)
    weakv=read(1635)["witness"]["Higgs_vev_GeV"]
    h1=-(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)*math.log(xt)/(2*(xt-1)**2)
    kernel=xt*(math.log(md/mw)+h1)*Lf
    tree=gf*(1+t*t)/math.sqrt(2)
    box=qcd*(2*md**2/weakv**2*Lf)**2/(64*math.pi**2*md**2)
    zeta=mix["eta_Z_RGI_range"]
    b4=interval([.0534,.0032]);b5=interval([.0493,.0037])
    a=[pre*fbhat[i]**2*nloop*lr[i]**2*eb[i]*s0[i] for i in (0,1)]
    b=[pre*fbhat[i]**2*(tree+box)*zeta[i]*ell[i]**2 for i in (0,1)]
    d=[pre*nloop*lr[i]*ell[i]*kernel*(1.18*b4[i]+1.42*b5[i]) for i in (0,1)]
    clock=interval(mix["omega_ps_inverse"]); phases=mix["phase_inputs_rad"]
    kk=[interval(phases[key]) for key in ("LHCb_KK_2015_18","CMS_KK_13TeV_final")]
    ptargets=[max(row[0] for row in kk)-mix["extra_decay_phase_cap_rad"]-phi0-dphase,
              min(row[1] for row in kk)+mix["extra_decay_phase_cap_rad"]-phi0+dphase]
    imcap=clock[1]*max(abs(math.sin(q)) for q in ptargets)
    kmin,kmax=nu["sm_coefficient_ratio"]["range"]
    # Physical odd vertices yield an initially loose but independent compact domain.
    edge=nmax+2
    ext=np.linspace(-edge,edge,bins+1)
    xl,xh=ext[:-1],ext[1:]
    def square_lo(x,y):
        return np.where(x*y<=0,0,np.minimum(x*x,y*y))
    def square_hi(x,y):
        return np.maximum(x*x,y*y)
    tN=1+t; tP=1-t
    ylo=np.maximum.reduce([np.zeros(bins),
        (nmin**2-square_hi(1+tN*xl,1+tN*xh))/tN**2,
        (pmin**2-square_hi(kmin-tP*xh,kmax-tP*xl))/tP**2])
    yhi=(pmax**2-square_lo(kmin-tP*xh,kmax-tP*xl))/tP**2
    # Im C=(d+2bx)y, with C the complex frequency before the tiny charm-phase rotation.
    dlo=d[0]+2*np.minimum(b[0]*xl,b[1]*xl); dhi=d[1]+2*np.maximum(b[0]*xh,b[1]*xh)
    yhi=np.minimum(yhi, imcap**2/np.maximum(square_lo(dlo,dhi),1e-40))
    # Re C <= omega_high, allowing maximal independent amplitude variations.
    min_real=a[0]+np.minimum(d[0]*xl,d[1]*xl)
    q=min_real-clock[1]
    ylo=np.maximum(ylo,square_lo(xl,xh)+np.minimum(q/b[0],q/b[1]))
    # Six coarse, analytic pieces. They also certify the result without the fine partition.
    x_min=(nmin**2/tN**2-pmax**2/tP**2-1/tN**2+kmin**2/tP**2)/(2*(1/tN+kmax/tP))
    assert x_min > .17
    def segments(extra):
        gaps={}
        for lo,hi in ((.17,.30),(.30,.33),(.33,.50),(.50,.95)):
            cap=((imcap+extra)/(d[0]+2*b[0]*lo))**2
            from_n=(nmin**2-(1+tN*hi)**2)/tN**2
            from_p=(pmin**2-max(abs(kmin-tP*hi),abs(kmax-tP*lo))**2)/tP**2
            gaps[f"{lo:.2f} to {hi:.2f}"]=max(from_n,from_p)-cap
        lo=.95
        cap=((imcap+extra)/(d[0]+2*b[0]*lo))**2
        gaps["0.95 to infinity (ps_inverse)"]=a[0]+d[0]*lo+b[0]*lo**2-b[1]*cap-clock[1]-extra
        return gaps
    # Least additive complex change certified by these coarse inequalities.
    lo,hi=0.,2.
    for _ in range(55):
        mid=(lo+hi)/2
        if min(segments(mid).values())>0:
            lo=mid
        else:
            hi=mid
    # Preserve an outward-rounded simple certificate for the text.
    epsilon_coarse=math.floor(lo*1000)/1000
    assert min(segments(epsilon_coarse).values()) > 0
    # Four necessary inequalities for all positions within each narrow x slab.
    # The interval enclosures below are enlarged by 1e-10 in dimensionless y^2.
    base_lo=np.maximum.reduce([np.zeros(bins),
        (nmin**2-square_hi(1+tN*xl,1+tN*xh))/tN**2,
        (pmin**2-square_hi(kmin-tP*xh,kmax-tP*xl))/tP**2])
    base_hi=(pmax**2-square_lo(kmin-tP*xh,kmax-tP*xl))/tP**2
    def minimum_gap(extra):
        q=min_real-clock[1]-extra
        lower=np.maximum(base_lo,square_lo(xl,xh)+np.minimum(q/b[0],q/b[1]))
        upper=np.minimum(base_hi,(imcap+extra)**2/np.maximum(square_lo(dlo,dhi),1e-40))
        return float(np.min(lower-upper))-1e-10
    lo,hi=0.,2.
    for _ in range(48):
        mid=(lo+hi)/2
        if minimum_gap(mid)>0:
            lo=mid
        else:
            hi=mid
    epsilon_fine=math.floor(lo*1000)/1000
    assert minimum_gap(epsilon_fine)>1e-8
    assert min(segments(0).values())>0
    margin=ylo-yhi
    possible=np.nonzero(margin<=0)[0]
    return {"M_GeV":md,"nmin":nmin,"nmax":nmax,"pmin":pmin,"pmax":pmax,
            "ell_range":ell.tolist(),"reference_top_modulus":prod,"phi0":phi0,
            "max_imag_clock":imcap,"a":a,"b":b,"d":d,
            "CMS_min_1e9":cmsmin,"R3_max_1e9":r3max,"num_bins":bins,
            "F_abs_outer_for_belle":(nmax+1)*ell[1]/(1+t),
            "smallest_y_squared_separation":float(np.min(margin)),
            "minimizer_x":[xl[margin.argmin()],xh[margin.argmin()]],
            "possible_cells":len(possible),"six_pieces_gaps":segments(0),
            "coarse_complex_residual_lower_ps_inv":epsilon_coarse,
            "fine_complex_residual_lower_ps_inv":epsilon_fine,
            "fine_post_residual_y2_gap":minimum_gap(epsilon_fine),
            "common_sqrt_Gram_floor":gram_min,"initial_F_bound":f_initial,
            "K_operator_norm_perturbation":delta_K,"analytic_x_floor":x_min,
            "possible_range":[] if len(possible)==0 else [xl[possible[0]],xh[possible[-1]]],
            "alone_f_sqrt_B":fbhat,"f":frange,"Bhat":bag,
            "f_sqrt_B_from_marginals":[frange[i]*math.sqrt(bag[i]) for i in (0,1)]}


def whole_interval(points):
    """One common enclosure for all M in [1600,2600], without interpolating data."""
    low,high=points["common_witness"],points["high_mass_witness"]
    poles=read(1633,"results.json")
    rb2=poles["right_doublet_probability"]
    # The only new Gram entry sensitive to M at fixed bottom charges is G_bb.
    # Over the whole interval add its worst change to the previous operator bound.
    gram_cross=(4.8**2*rb2)*(1600**-2-2600**-2)
    dk=low["K_operator_norm_perturbation"]+gram_cross/(2*math.sqrt(.988))
    assert 2*dk<5e-7
    dl=5e-7
    sratio=dl/min(row["reference_top_modulus"] for row in (low,high))
    assert (1+sratio)**2 <1.0001
    tiny=5*(2*.1*4.8/1600**2)
    # |u z-z|, |v z-z| <= tiny. Use |z| < 5 throughout the Belle rectangle.
    assert max(row["nmax"] for row in (low,high))*(1+sratio)+1<5
    n=min(row["nmin"] for row in (low,high))/(1+sratio)-tiny
    p=min(row["pmin"] for row in (low,high))/(1+sratio)-tiny
    pmax=max(row["pmax"] for row in (low,high))/(1-sratio)+tiny
    A=min(row["a"][0] for row in (low,high))/1.0001
    D=low["d"][0]/1.0001
    B0=low["b"][0]/1.0001
    B1=high["b"][1]*1.0001
    clock_hi=interval(read(1468)["omega_ps_inverse"])[1]
    im_hi=max(row["max_imag_clock"] for row in (low,high))+clock_hi*1e-4
    # Endpoint factors in B: use the independent adverse endpoint for 1+t^2.
    # Relative drift of this near-unit factor is < 1.5e-13, far inside 1e-4.
    assert (1+(2*.1*4.8/1600**2)**2)*(1+sratio)**2 < 1.0001
    xfloor=(n*n-pmax*pmax-1+.6**2)/(2*(1+.7))
    assert xfloor>.17
    def gaps(extra):
        values=[]
        for a,b in ((.17,.30),(.30,.33),(.33,.50),(.50,.95)):
            yhi=((im_hi+extra)/(D+2*B0*a))**2
            ylo=max(n*n-(1+b)**2,p*p-max(abs(.6-b),abs(.7-a))**2)
            values.append(ylo-yhi)
        a=.95
        values.append(A+D*a+B0*a*a-B1*((im_hi+extra)/(D+2*B0*a))**2-clock_hi-extra)
        return values
    assert min(gaps(.075))>1e-5
    return {"mass_interval_GeV":[1600,2600],"global_2deltaK":2*dk,
            "x_floor":xfloor,"nmin":n,"pmin":p,"pmax":pmax,"A_lower":A,"D_lower":D,
            "B_lower_upper":[B0,B1],"I_upper":im_hi,"zero_error_six_piece_gaps":gaps(0),
            "additional_complex_residual_test_ps_inverse":.075,
            "all_six_piece_gaps_at_added_residual":gaps(.075),"possible_cells":0}



def round_numbers(data):
    if isinstance(data,dict):return {k:round_numbers(v) for k,v in data.items()}
    if isinstance(data,(list,tuple)):return [round_numbers(v) for v in data]
    if isinstance(data,(float,np.floating)):return float(f"{data:.10g}")
    return data


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--check",action="store_true")
    args=p.parse_args()
    data=round_numbers({lab:solve(md,lab) for md,lab in ((1600,"common_witness"),(2600,"high_mass_witness"))})
    data["whole_mass_interval"]=round_numbers(whole_interval(data))
    target=HERE/"independent_results.json"
    if args.check:
        assert data==json.loads(target.read_text(encoding="utf-8"))
        print("1639 independent projected intervals: PASS",[(k,v["possible_cells"]) for k,v in data.items()])
    else:
        target.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
        print(json.dumps(data))
