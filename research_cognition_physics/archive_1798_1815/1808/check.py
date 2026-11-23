"""Equal-trace Higgs and type-I tests on the full D,T threshold mass band."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())

CFG=read("inputs.json")
FLOW=runpy.run_path(str(HERE/CFG["higgs_flow"]))
FIELDS=runpy.run_path(str(HERE/CFG["threshold_solver"]))
SM=FIELDS["BASE_B"]
D,T=FIELDS["DOUBLET"],FIELDS["TRIPLET"]
TAU=2*math.pi
LOOP=16*math.pi**2
WEIGHTS=np.repeat([3.,3.,1.],3)
PAD=1e-12


def inputs(label):
    gauge=read(CFG["gauge_data"])
    weak=read(CFG["masses_and_ckm"])
    quart=read(CFG["higgs_data"])
    g,dg=np.array([gauge[k] for k in ["g1_GUT","g2","g3"]]).T
    a=4*math.pi/g**2
    al=4*math.pi/(g+2*dg)**2
    ah=4*math.pi/(g-2*dg)**2
    if label=="independent_lattice":
        external=read(CFG["lattice_gauge"])
        alpha,da=external["strong_lattice_5f"]["ALPHA_2026"]
        mt,dm=external["top_MS_matching_GeV"]
        def strong(x,y):
            return 1/x+math.log(y/gauge["mu_ref_GeV"])/(3*math.pi)
        a[2],al[2],ah[2]=strong(alpha,mt),strong(alpha+2*da,mt-2*dm),strong(alpha-2*da,mt+2*dm)
    ys,dy=np.concatenate([weak["yukawa"][k] for k in ["up","down","charged_leptons"]]).T
    lam,dl=np.array(quart["lambda_paper"])/quart["lambda_paper_factor"]
    weak_input=np.array(weak["quark_sines_and_phase"])[:,0]
    return dict(a=a,al=al,ah=ah,y=ys,dy=dy,lam=lam,dl=dl,
                v=FLOW["mixing"](np.arcsin(weak_input[:3]),weak_input[3]),MZ=gauge["mu_ref_GeV"])


def thresholds(point, mz):
    return [(math.log(m/mz),rep) for key,rep in [("M_D_GeV",D),("M_T_GeV",T)] for m in point[key]]


def couplings(data, fields, t):
    a=data["a"]-SM*t/TAU
    for birth,delta in fields:
        a-=delta*max(t-birth,0)/(3*math.pi)
    return 4*math.pi/a


def advance(x, start,stop, events, rhs, hmax):
    times=[start,*sorted(set(t for t,_ in events if start<t<stop)),stop]
    for a,b in zip(times,times[1:]):
        n=max(1,math.ceil((b-a)/hmax))
        h=(b-a)/n
        for i in range(n):
            t=a+i*h
            k1=rhs(t,x); k2=rhs(t+h/2,x+h*k1/2)
            k3=rhs(t+h/2,x+h*k2/2); k4=rhs(t+h,x+h*k3)
            x=x+h*(k1+2*k2+2*k3+k4)/6
    return x


def low_flow(data,fields,hi,hmax):
    ys=data["y"]
    u=np.diag(ys[:3]**2)
    d=(data["v"]*ys[3:6]**2)@data["v"].conj().T
    e=np.diag(ys[6:]**2)
    seed=np.r_[u.ravel(),d.ravel(),e.ravel(),data["lam"],np.zeros(3)]
    def rhs(t,z):
        gr=z[:27].reshape(3,3,3)
        g=couplings(data,fields,t)
        sq=gr@gr
        tr=float(np.dot([3,3,1],np.trace(gr,axis1=1,axis2=2).real))
        b=float(np.dot([3,3,1],np.trace(sq,axis1=1,axis2=2).real))
        quark=gr[0]@gr[1]+gr[1]@gr[0]
        deriv=3*sq+2*(tr-FLOW["GAUGE_DAMPING"]@g)[:,None,None]*gr
        deriv[:2]-=1.5*quark
        l=z[27].real
        w,gy=g[1],.6*g[0]
        dl=24*l*l+(4*tr-9*w-3*gy)*l-2*b+.375*(2*w*w+(w+gy)**2)
        wk=(4*l-3*w+2*tr)/2-1.5*np.diag(gr[2]).real
        return np.r_[deriv.ravel(),dl,wk]/LOOP
    z=advance(seed,0,hi,fields,rhs,hmax)
    sq=np.linalg.eigvalsh(z[:27].reshape(3,3,3)).ravel()
    return sq,z[27].real,np.exp(z[28:].real)


def log_ranges(box,n):
    t0,t1=box["log_scale"]
    if n==1:
        k0,k1=box["two_D_one_T_log_min_geometric_mean_D_over_T"]
        return [(0.,8*(t1-k0)/7),(k0,t1)]
    k0,k1=box["two_complete_pairs_log_mean_T_over_D"]
    return [(2*k0,2*t1),(0.,2*t1-2*k0)]


def cumulative(n,lo,hi,tau,t0,t1):
    delta=max(0.,t0-tau)
    lower=max(0.,lo-n*max(0.,t1-tau))
    k=min(n,int(hi/t1))
    up=k*max(0.,t1-delta)+max(0.,hi-k*t1-delta)
    return max(0.,lower-1e-11),min(n*tau,up)+1e-11


def interval_run(data,box,n):
    # Total mass logarithms need not fix the ordering of the individual poles.
    # Convexity provides cumulative active-threshold envelopes before Q.
    t0,t1=box["log_scale"]
    ranges=log_ranges(box,n)
    h=CFG["interval_step"]
    yl=(data["y"]-data["dy"])**2
    yh=(data["y"]+data["dy"])**2
    ll,lh=data["lam"]-data["dl"],data["lam"]+data["dl"]
    grow=0.
    rows=[]
    w_lows=[]
    for i in range(math.ceil(t1/h)):
        t=i*h
        ds=np.array([cumulative(j,*span,t,t0,t1) for j,span in zip((2,n),ranges)])
        de=np.array([cumulative(j,*span,t+h,t0,t1) for j,span in zip((2,n),ranges)])
        qlo=4*math.pi/(data["ah"]-np.minimum(SM*t,SM*(t+h))/TAU-(D*ds[0,0]+T*ds[1,0])/(3*math.pi))-PAD
        qhi=4*math.pi/(data["al"]-np.maximum(SM*t,SM*(t+h))/TAU-(D*de[0,1]+T*de[1,1])/(3*math.pi))+PAD
        low, high=yl*math.exp(-h)*(1-PAD),yh*math.exp(h)*(1+PAD)
        a,b=WEIGHTS@low,WEIGHTS@high
        peers=np.repeat([high[3:6].max(),high[:3].max(),0.],3)
        lower=2*(a+1.5*(low-peers)-np.repeat(FLOW["GAUGE_DAMPING"]@qhi,3))/LOOP-PAD
        upper=2*(b+1.5*high-np.repeat(FLOW["GAUGE_DAMPING"]@qlo,3))/LOOP+PAD
        assert min(lower)>-.99 and max(upper)<.99
        left,right=ll-.16*h-PAD,lh+.16*h+PAD
        c,d=FLOW["times_product"](4*a-9*qhi[1]-1.8*qhi[0],4*b-9*qlo[1]-1.8*qlo[0],left,right)
        def fourth(g):
            return .375*(2*g[1]**2+(g[1]+.6*g[0])**2)
        small=0 if left<=0<=right else min(left*left,right*right)
        rate_l=(24*small+c+fourth(qlo)-2*WEIGHTS@high**2)/LOOP-PAD
        rate_h=(24*max(left*left,right*right)+d+fourth(qhi)-2*WEIGHTS@low**2)/LOOP+PAD
        assert max(abs(rate_l),abs(rate_h))<.159
        growth=(4*right-3*qlo[1]+2*b)/LOOP+PAD
        grow+=h*max(0.,growth)
        if t+h>=t0:
            # Use the exact gauge trace C=(X-sum l_T/(3pi))/(6+n).
            stlo,sthi=ranges[1]
            cl=(box["X"][0]-sthi/(3*math.pi))/(6+n)
            ch=(box["X"][1]-stlo/(3*math.pi))/(6+n)
            slo,shi=(8*math.pi/(3*ch),8*math.pi/(3*cl))
            deficit=slo-b
            minb=WEIGHTS@low**2
            assert deficit>0
            rows.append([a,b,minb,slo,shi,deficit,left,right,
                         minb/shi-right,(minb+deficit**2/3)/shi-right])
            w_lows.append(deficit*CFG["ordinary_vev_GeV"]**2/
                          (2*data["MZ"]*math.exp(t1)*math.exp(grow))*1e9)
        yl*=np.exp(h*lower)*(1-PAD); yh*=np.exp(h*upper)*(1+PAD)
        ll+=h*rate_l-PAD; lh+=h*rate_h+PAD
    a=np.array(rows)
    lower_sum=min(w_lows)
    return dict(S_interval=[min(a[:,3]),max(a[:,4])],
                charged_squared_interval=[min(a[:,0]),max(a[:,1])],
                top_and_other_fourth_lower=min(a[:,2]),missing_kinetic_lower=min(a[:,5]),
                Higgs_from_data=[min(a[:,6]),max(a[:,7])],
                charged_alone_min_extra_quartic=min(a[:,8]),rank3_min_extra_quartic=min(a[:,9]),
                Weinberg_nuclear_norm_upward_gain_upper=math.exp(grow),
                positive_typeI_light_sum_lower_eV=lower_sum,
                positive_typeI_exceeds_cosmology=bool(lower_sum>CFG["comparison_cosmological_upper_eV"]))


def neutral_check(Q,S,missing,b,lam,left):
    """Common physical threshold: (bare, condensate)=(Q/2,Q/2) in each N."""
    nu=read(CFG["light_neutrinos"])["NO"]
    v=CFG["ordinary_vev_GeV"]
    u=FLOW["mixing"](np.arcsin(np.sqrt([nu["s12"][0],nu["s13"][0],.55])),math.radians(210))
    # Opposite electron double-beta contributions, selected independently of Higgs.
    u=u@np.diag([1,1j,1j*np.exp(1j*math.radians(210))])
    jac=np.array(left)
    target_B=S*lam-b+2*missing**2/3
    normalized=target_B/missing**2
    if normalized<.5:
        t=math.sqrt(max(0.,(3*normalized-1)/2))
        determinant_lower=missing**3*(1+t)**2*(1-2*t)/27
        mass_lower=3*(determinant_lower/math.prod(jac)**2)**(1/3)*v*v/(2*Q)*1e9
    else:
        mass_lower=0.
    u2=missing*(1-math.sqrt((3*normalized-1)/2))/3
    max_light_lower=v*v*1e9*u2/(2*Q*max(jac)**2)
    all_nu=read(CFG["light_neutrinos"])
    spread=max(all_nu[o]["solar"][2]+all_nu[o]["atm_abs"][2] for o in ["NO","IO"])
    osc_sum_bound=max_light_lower+2*math.sqrt(max(0.,max_light_lower**2-spread))
    def sources(m0,t):
        m=np.sqrt(m0*m0+np.array([0.,nu["solar"][0],nu["atm_abs"][0]]))
        low=(u*m)@u.T
        root=((jac[:,None]*u)*np.sqrt(m))*(math.sqrt(2*Q*1e-9)/v)
        co,si=math.cosh(t),math.sinh(t)
        R=np.eye(3,dtype=complex)
        R[1:,1:]=[[co,1j*si],[-1j*si,co]]
        y=root@R
        A=y.conj().T@y
        return A,low,y,m
    def root(fn,lo,hi):
        a,b=fn(lo),fn(hi)
        assert a*b<=0
        for _ in range(55):
            mid=(lo+hi)/2
            if fn(mid)*a>0:
                lo=mid
            else:
                hi=mid
        return (lo+hi)/2
    def kinetic(m0,t):
        return sources(m0,t)[0].trace().real-missing
    ceiling=root(lambda mass:kinetic(mass,0),0.,1.)
    def boosted(mass):
        z=0. if abs(kinetic(mass,0))<1e-12 else root(lambda x:kinetic(mass,x),0,5)
        return z,sources(mass,z)
    def potential(mass):
        a=boosted(mass)[1][0]
        return np.trace(a@a).real-target_B
    m0=root(potential,0.,ceiling)
    z,(A,massmat,y,light)=boosted(m0)
    np.testing.assert_allclose(v*v/(2*Q)*1e9*(y@y.T)/(jac[:,None]*jac[None,:]),massmat,atol=1e-12)
    np.testing.assert_allclose(A.trace().real,missing,atol=1e-12)
    lh=(b+np.trace(A@A).real)/S
    ls=2*S/3
    lp=2*missing/3
    B0=math.sqrt(S/3)
    w=Q/(2*B0)
    np.testing.assert_allclose(2*ls*w*w,Q*Q,rtol=1e-14)
    np.testing.assert_allclose(lh-lp**2/ls,lam,atol=3e-14)
    sum_floor=max(mass_lower,osc_sum_bound)
    beta=math.sqrt(float(abs(u[0])**2@(light*light)))
    return dict(neutral_fourth_over_squared_sum=normalized,
                determinant_light_sum_lower_eV=mass_lower,
                second_singular_light_max_lower_eV=max_light_lower,
                any_CP_low_sum_nu_necessary_eV=sum_floor,
                all_CP_exceed_comparison_cosmology=bool(sum_floor>CFG["comparison_cosmological_upper_eV"]),
                masses_low_eV=light.tolist(),sum_light_eV=sum(light),
                beta_endpoint_eV=beta,within_direct_beta_upper=bool(beta<all_nu["katrin"]["m_beta_upper_eV"]),
                double_beta_eV=abs(massmat[0,0]),imaginary_orthogonal_angle=z,
                Yukawa_squared_eigenvalues=np.linalg.eigvalsh(A).tolist(),
                triplet_of_equal_heavy_masses_GeV=Q,sigma_curvature_mass_GeV=Q,
                vev_GeV=w,bare_Majorana_over_physical=.5,
                quartics_h_hsigma_sigma=[lh,lp,ls],quartic_threshold_shift=lp**2/ls,
                scalar_to_doublet_partial_width_over_mass=(2*lp*w/Q)**2/(8*math.pi),
                soft_squared_ratio_required_before_weak_residual=lp/ls)


def examples(data,point):
    events=thresholds(point,data["MZ"])
    scale=point["matching_scale_GeV"]
    hi=math.log(scale/data["MZ"])
    sq,lam,le=low_flow(data,events,hi,.12)
    fine,l2,le2=low_flow(data,events,hi,.06)
    np.testing.assert_allclose(sq,fine,rtol=1e-8)
    assert abs(lam-l2)<1e-9
    S=8*math.pi/(3*point["common_trace_coefficient"])
    np.testing.assert_allclose(S,4*(6+2)/6*couplings(data,events,hi)[1],rtol=1e-12)
    a,b=WEIGHTS@fine,WEIGHTS@fine**2
    deficit=S-a
    return dict(Q_GeV=scale,physical_D_T_GeV=[point["M_D_GeV"],point["M_T_GeV"]],
                g_squared_at_Q=couplings(data,events,hi).tolist(),squared_yukawas_at_Q=fine.tolist(),
                lambda_observed=l2,S=S,charged_A=a,charged_B=b,missing_A=deficit,
                neutral_rank3_threshold=(b+deficit**2/3)/S-l2,
                Weinberg_left_gains=le2.tolist(),
                common_scalar_and_three_neutrino_threshold=neutral_check(scale,S,deficit,b,l2,le2))


def run():
    stored=read(CFG["threshold_examples"])
    out={}
    for label in ("common_MZ","independent_lattice"):
        data=inputs(label)
        cases={}
        for n,key in [(1,"two_doublets_one_triplet"),(2,"two_split_complete_pairs")]:
            cases[key]=dict(band=interval_run(data,stored[label]["two_marginal_box"],n),
                            point=examples(data,stored[label][key]))
        out[label]=cases
    return out


def clean(v):
    if isinstance(v,dict): return {k:clean(x) for k,x in v.items()}
    if isinstance(v,list): return [clean(x) for x in v]
    if isinstance(v,(float,np.floating)): return float(f"{v:.10g}")
    return v


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result=clean(run())
    raw=json.dumps(result,indent=2,allow_nan=False)+"\n"
    dst=HERE/"results.json"
    if p.parse_args().check: assert result==json.loads(dst.read_bytes())
    else: dst.write_bytes(raw.encode())
    print("1808",hashlib.sha256(raw.encode()).hexdigest()[:12])
