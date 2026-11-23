"""Profile the published same-charge mass covariance and bound a top-directed threshold."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
def read(p): return json.loads((HERE/p).read_bytes())
def api(p): return runpy.run_path(str(HERE/p))
CFG=read("inputs.json")
TABLE=CFG["table_IV"]
MU=np.array(TABLE["means"][:3]); SIG=np.array(TABLE["sigmas"][:3])
COV=np.outer(SIG,SIG)*np.array(TABLE["correlation"])[:3,:3]
PREC=np.linalg.inv(COV)
KDS,KUD,KBS=read(CFG["matching"])["MZ_over_QCD_ratio_ds_ud_bs_reference"]
SEM=read(CFG["source_weak"])["direct_semileptonic_moduli"]
Y=np.array(read(CFG["quark_inputs"])["yukawa"]["down"])
engine=api(CFG["joint"])
PHASE=api("../1736/independent_check.py")

def support(a,R):
    return a@MU+R*math.sqrt(a@COV@a)

def ratio_high(a,R):
    # Fieller support for numerator (a.x) divided by x_s in the ellipsoid.
    c=np.array([0.,1.,0.])
    lo,hi=-100.,100.
    assert support(a-hi*c,R)<0<support(a-lo*c,R)
    for _ in range(70):
        v=(lo+hi)/2
        if support(a-v*c,R)>0: lo=v
        else: hi=v
    return hi

def reference_shape(mz):
    d,s,b=mz
    return np.array([d/(s*KDS),1.,b/(s*KBS)])

def singular_after(mz,top,t,down_gain=1.):
    arm=reference_shape(mz)
    arm[0]*=down_gain
    w=top*arm
    G=np.diag(arm*arm)+math.expm1(2*t)*np.outer(w.conj(),w)
    return np.sqrt(np.linalg.eigvalsh(G))

def profile(arm):
    h=arm[::-1]/arm[1]
    norm=(h@PREC@MU)/(h@PREC@h)
    diff=norm*h-MU
    return math.sqrt(diff@PREC@diff)

def CK(row):
    f=PHASE["ckm"](row["weak_offsets"],row["gamma"])
    us,ub,cb=f["us_ub_cb"]
    c13=math.sqrt(1-ub**2); s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12*s12),math.sqrt(1-s23*s23)
    z=np.exp(1j*np.radians(f["pdg_delta_deg"]))
    return np.array([s12*s23-c12*c23*ub*z,-c12*s23-s12*c23*ub*z,c23*c13])

def table_row(row):
    mz=np.array(row["down_MZ"])
    top=CK(row)
    def calc(t): return profile(singular_after(mz,top,t))
    # The selected minimum is a numerical witness. Domain-wide exclusions below
    # use either the support function or analytic bounds, not this local search.
    grid=np.linspace(*CFG["reference_t_range"],61)
    j=np.argmin([calc(t) for t in grid])
    l,r=grid[max(0,j-1)],grid[min(len(grid)-1,j+1)]
    for _ in range(50):
        a,b=(2*l+r)/3,(l+2*r)/3
        if calc(a)<calc(b): r=b
        else: l=a
    t=(l+r)/2
    arm=singular_after(mz,top,t)
    target=MU[[2,0]]/MU[1]
    def residual(lt):
        s=singular_after(mz,top,lt[1],math.exp(lt[0]))
        return np.log((s[[0,2]]/s[1])/target)
    lt=np.zeros(2)
    for _ in range(6):
        v=residual(lt)
        if max(abs(v))<2e-10: break
        jac=np.column_stack([(residual(lt+e)-residual(lt-e))/(2e-4) for e in np.eye(2)*1e-4])
        lt-=np.linalg.solve(jac,v)
    assert np.linalg.norm(residual(lt))<1e-9
    amended=singular_after(mz,top,lt[1],math.exp(lt[0]))
    us0=row.get("up_MZ",engine["HAD"]["UP"][0])/mz[1]/(KUD*KDS)
    def fourth(s):
        full=np.array([s[2],s[1],s[0],us0])/s[1]
        mu=np.array(TABLE["means"]); sig=np.array(TABLE["sigmas"])
        fullprec=np.linalg.inv(np.outer(sig,sig)*TABLE["correlation"])
        n=(full@fullprec@mu)/(full@fullprec@full)
        d=n*full-mu
        S=2*s[1]/(us0+s[0])
        Q=math.sqrt((s[1]**2-((us0+s[0])/2)**2)/(s[0]**2-us0**2))
        return {"four_flavour_radius":math.sqrt(d@fullprec@d),"u_over_d":us0/s[0],"S":S,"eta_Q":Q}
    f=PHASE["ckm"](row["weak_offsets"],row["gamma"])
    p=np.array([1-f["us_ub_cb"][0]**2-f["us_ub_cb"][1]**2,f["us_ub_cb"][0]**2,f["us_ub_cb"][1]**2])
    F=(p[0]*mz[0]-p[1]*mz[1]+p[2]*mz[2])/mz[1]
    return {"weak":row["weak_offsets"],"gamma_deg":row["gamma"],"D_mz":mz[0]/mz[1],"B_mz":mz[2]/mz[1],
            "F_mz_at_actual_point":F,"t0_radius":calc(0),"profile_t":t,"profile_witness_radius":calc(t),
            "profiled_low_d_over_s":arm[0]/arm[1],"profiled_low_b_over_s":arm[2]/arm[1],
            "fourth_flavour_at_profile":fourth(arm),
            "additional_d_direction_at_central_same_charge_ratios":{"d_factor":math.exp(lt[0]),"t":lt[1],**fourth(amended)}}

def drift():
    # See proof, dimensionless time tau=log(mu/MZ)/(16*pi^2). Scalar traces
    # cancel in F/s. This is a sufficient corridor for the registered one-loop flow.
    duration=math.log(7171.08/91.1876)/(16*math.pi**2)
    u=np.array([1.,.23,.0048]); top=np.array([.046,.046,1.])
    nu=np.outer(top,top)+.005**2
    D=.056; B=59.; Bmin=48.
    M=np.array([[D/2,D/(1-D),D],[D/(1-D),.5,Bmin/(Bmin-1)],[D,Bmin/(Bmin-1),B/2]])
    term=3*np.sum(M*np.outer(u,u)*nu)
    yb,ys=.020,.00035
    self=1.5*(ys**2*(D**3+u[1]**2)+yb**2*B*u[2]**2)
    du=1.52*(yb**2*u[2]+ys**2)
    leg=2*du*np.linalg.norm(u*np.array([D,1.,B]))
    lam=3*(1e-5)**2+1.5*(ys**2+nu[1,1])
    # Close the assumed spectral and angle corridors uniformly, by continuation.
    eff=1.51
    ratios=[(Y[0,0]+2*Y[0,1])/(Y[1,0]-2*Y[1,1])*math.exp(eff*duration),
            (Y[2,0]-2*Y[2,1])/(Y[1,0]+2*Y[1,1])*math.exp(-eff*duration),
            (Y[2,0]+2*Y[2,1])/(Y[1,0]-2*Y[1,1])*math.exp(eff*duration)]
    assert ratios[0]<D and ratios[1]>Bmin and ratios[2]<B
    vb=[(SEM[k][0]+2*math.hypot(*SEM[k][1:])+.0014)*math.exp(eff*duration)-.0014 for k in ("Vub","Vcb")]
    assert vb[0]<.0048 and math.hypot(*vb)<.046
    usmax=SEM["Vus"][0]+2*math.hypot(*SEM["Vus"][1:])+duration*1.52*(yb**2+nu[0,1]+.0048*nu[1,2])
    assert usmax<.23
    g3sq_min=1/(1/1.20**2+14*duration)
    assert 1.5+3*(1+2*.005**2)+9*yb**2+.2<8*g3sq_min
    return {"duration":duration,"dot_sources":[term,self,leg],"linear_coefficient":lam,
            "closed_Dmax_Bmin_Bmax":ratios,"uniform_ub_cb_us_upper":[*vb,usmax],"g3_squared_lower":g3sq_min,
            "upper_F_mz_if_zero_at_upper":(term+self+leg)*math.expm1(lam*duration)/lam}

def threshold_bound(rad):
    # Bound the top orientation during the finite mass matching. The b column alone
    # has the displayed angle; the two small columns perturb its squared-mass projector.
    eta=sum((SEM[k][0]+2*math.hypot(*SEM[k][1:]))**2 for k in ("Vub","Vcb"))
    Bmin=(Y[2,0]-2*Y[2,1])/(Y[1,0]+2*Y[1,1])/KBS
    Bmax=(Y[2,0]+2*Y[2,1])/(Y[1,0]-2*Y[1,1])/KBS
    BLmin=-ratio_high(np.array([-1.,0.,0.]),rad); BLmax=ratio_high(np.array([1.,0.,0.]),rad)
    # Rule out disconnected solutions at arbitrary positive gain f=exp(t).
    def limit(k,sign):
        mean,*errors=SEM[k]
        return mean+sign*2*math.hypot(*errors)
    tdmin=limit("Vus",-1)*limit("Vcb",-1)-limit("Vub",1)
    D0max=(Y[0,0]+2*Y[0,1])/(Y[1,0]-2*Y[1,1])/KDS
    def small_gain_upper(x):
        numerator=Bmax*(eta+x)+(1+D0max**2)/Bmin
        return numerator/math.sqrt(x+(1-x)*tdmin**2)
    small_upper=max(small_gain_upper(0),small_gain_upper(.25))
    large_lower=2*Bmin*(1-eta)
    assert tdmin>0 and small_upper<BLmin and large_lower>BLmax
    def leakage(fmin,fmax):
        gap=fmin*fmin*(1-eta)
        perturb=max(1,fmax*fmax)/(Bmin*Bmin)
        return (math.sqrt(eta/gap)+perturb/(gap-perturb))**2
    leak=leakage(.5,2.)
    abs_log=max(abs(math.log(BLmin/Bmax)),abs(math.log(BLmax/Bmin)))
    tmax=abs_log/(1-2*leak)
    assert tmax<math.log(2)
    leak=leakage(math.exp(-tmax),math.exp(tmax))
    tmax=abs_log/(1-2*leak)
    Dmax=ratio_high(np.array([0.,0.,1.]),rad)
    return {"B0_bounds":[Bmin,Bmax],"low_B_bounds":[BLmin,BLmax],"max_abs_t":tmax,
            "all_real_t_tail_bounds":{"abs_Vtd_lower":tdmin,"f_below_half_B_upper":small_upper,
                "f_above_two_B_lower":large_lower},
            "max_low_top_probability_in_light_lines":leak,
            "abs_change_D_due_to_top":Dmax*math.expm1(leak*tmax),
            "abs_change_B_due_to_top":max(Bmax-BLmin,BLmax-Bmin)}

def weak_certificate(axis,vector,rad):
    box=engine["kaon_projection"](2.,axis,vector)
    if box is None: return {"empty_weak_projection":True}
    # Lower strange and upper beauty maximize the signed return at any fixed masses.
    ps=box[0]**2
    pb=(SEM["Vub"][0]+2*math.hypot(*SEM["Vub"][1:]))**2
    pd=1-ps-pb
    a=np.array([pb*KBS,-ps,pd*KDS])
    raw=ratio_high(a,rad)
    b=threshold_bound(rad)
    err=pd*KDS*b["abs_change_D_due_to_top"]+pb*KBS*b["abs_change_B_due_to_top"]
    prop=drift()["upper_F_mz_if_zero_at_upper"]
    return {"weak_Vus_projection":box,"highest_signed_return_no_top":raw,
            "top_response_compensation_upper":err,"max_F_after_top":raw+err,
            "remaining_gap_after_running":-raw-err-prop,
            "additional_fractional_d_response_necessary_if_only_d_changes":max(0.,-raw-err-prop)/(pd*KDS*(ratio_high(np.array([0.,0.,1.]),rad)+b["abs_change_D_due_to_top"]))}

def main():
    rr=read(CFG["points"])["witnesses"]
    out={k:table_row(v) for k,v in rr.items()}
    cp=engine["shoot"](CFG["B_K_CP_weak_offsets_from_independent_inverse"],0.,.25)
    out["B_and_K_CP_centres"]=table_row(cp)
    # Re-shoot an endpoint favouring d/s at the physical edge of the full weak box.
    weak={"us":(.223994-SEM["Vus"][0])/math.hypot(*SEM["Vus"][1:]),"ub":2.,"cb":2.}
    border=engine["shoot"](weak,-1.8,.25)
    out["RBC_weak_edge"]=table_row(border)
    out["RBC_weak_edge"]["external_residuals_in_marginal_widths"]=border["signed_widths"]
    d=read("../1732/results.json")["exploratory_k3"]["registered_with_lower_us"]
    out["1732_primary"]=table_row({"weak_offsets":{"us":-.95,"ub":-.93,"cb":.93},"gamma":d["gamma_MZ_deg"],"down_MZ":d["down_strange_bottom_MZ"]})
    return {"independent_profiled_points":out,"flow_bound":drift(),
            "top_bound":{str(r):threshold_bound(r) for r in CFG["radii"]},
            "domain":{axis+"_"+vector:{str(r):weak_certificate(axis,vector,r) for r in CFG["radii"]}
                      for axis in ("RM123_axial","RBC_axial") for vector in ("original","alternative")}}

if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    result=engine["HAD"]["tidy"](main())
    path=HERE/"results.json"
    if check:
        assert result==json.loads(path.read_bytes())
        print("1738 joint lattice profile, signed return and common top threshold verified")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
