"""Test the full charged GJ basis against a rank-two Majorana source and three mixing angles."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def charged_basis(cfg):
    old=runpy.run_path(str(HERE.parent/"1657"/"independent_check.py"))
    point=read(cfg["charged"])["extended_box_GJ_test"]
    m=old["mtx"](old["texture"](point["marginal_coordinates"]),[1,1,3],point["t"]).T
    m*=point["k_times_yb"]
    squared,Ue=np.linalg.eigh(m@m.conj().T)
    Yukawas=np.sqrt(squared)
    np.testing.assert_allclose(Yukawas,[point["predicted_e"],*point["marginal_coordinates"][-2:]],rtol=1e-8)
    return Ue,Yukawas


def takagi(M):
    _,U=np.linalg.eigh(M.conj().T@M)
    diag=np.diag(U.T@M@U)
    U=U@np.diag(np.exp(-.5j*np.angle(diag)))
    d=abs(diag)
    np.testing.assert_allclose(U.T@M@U,np.diag(d),atol=2e-14)
    return d,U


def rank_two(ds,da,eta,branch):
    R=math.sqrt(1+da/ds)
    B=6*(R+1/R)+16
    r=8/(B+math.sqrt(B*B-1936)) if branch=="small" else (B+math.sqrt(B*B-1936))/242
    a,b=np.array([0,1,1]),np.array([1,3,1])
    M=np.outer(a,a)+r*np.exp(1j*eta)*np.outer(b,b)
    d,U=takagi(M)
    scale=math.sqrt(ds)/d[1]
    np.testing.assert_allclose(np.diff((d*scale)**2),[ds,da],rtol=2e-12)
    return M*scale,U,r


def angles(V):
    p=abs(V)**2
    s13=p[0,2]
    return np.array([4*s13*(1-s13),p[0,1]/(1-s13),p[1,2]/(1-s13)])


def predict(Ue,M,alpha,beta):
    # The declared P phases enter the *mass*, so the Takagi modes transform with P*.
    P=np.diag(np.exp(1j*np.array([0.,alpha,beta])))
    physical=P@M@P
    d,U=takagi(physical)
    v=Ue.conj().T@U
    a=angles(v)
    J=(v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    mbeta=math.sqrt(sum(abs(v[0,:]*d)**2))
    mbb=abs(sum(v[0,:]**2*d))
    np.testing.assert_allclose(mbb,abs((Ue.T@physical@Ue)[0,0]),atol=2e-16)
    return {"A13_solar_atmospheric":a.tolist(),"J_dirac":J,"m_beta_eV":mbeta,"m_beta_beta_eV":mbb}


def search_point(Ue,M,U,centers,half,steps):
    grid=np.arange(steps)*2*math.pi/steps-math.pi
    phases=np.exp(-1j*grid)
    bestscore=math.inf
    best=None
    for i,p1 in enumerate(phases):
        def prob(l,j):
            z=Ue[:,l].conj()*U[:,j]
            return abs(z[0]+p1*z[1]+phases*z[2])**2
        e3=prob(0,2)
        readouts=np.array([4*e3*(1-e3),prob(0,1)/(1-e3),prob(1,2)/(1-e3)])
        score=np.max(abs(readouts-centers[:,None])/half[:,None],axis=0)
        j=int(np.argmin(score))
        if score[j]<bestscore:
            bestscore=score[j];best=(grid[i],grid[j])
    result=predict(Ue,M,*best)
    assert abs(max(abs(result["A13_solar_atmospheric"]-centers)/half)-bestscore)<1e-10
    return {**result,"mass_P_alpha_beta_rad":list(best),"maximum_fraction_of_enclosing_halfwidth":bestscore,
            "all_three_inside":bool(bestscore<1)}


def large_bound(Ue,osc):
    def ends(row):return row[0]-2*row[1],row[0]+2*row[2]
    ds,da=ends(osc["KamLAND_dm21_eV2"]),ends(osc["T2K_dm32_NO_eV2"])
    lo,hi=math.sqrt(1+da[0]/ds[1]),math.sqrt(1+da[1]/ds[0])
    def root(R):
        B=6*(R+1/R)+16
        return (B+math.sqrt(B*B-1936))/242
    # Electron entry of the heavy spectral projector = (11 r - 2/R)/(6(R-1/R)).
    xmax=math.sqrt((11*root(hi)-2/hi)/(6*(lo-1/lo)))
    xmin=math.sqrt((11*root(lo)-2/lo)/(6*(hi-1/hi)))
    a=abs(Ue[0,0]);b=math.sqrt(1-a*a)
    pmin=(max(0,a*xmin-b*math.sqrt(1-xmin*xmin)))**2
    pmax=(a*xmax+b)**2
    assert pmax<.5
    low=4*pmin*(1-pmin)
    assert low>ends(osc["DayaBay_sin2_2theta13"])[1]
    return {"R_range":[lo,hi],"bare_heavy_electron_modulus_range":[xmin,xmax],
            "physical_heavy_electron_weight_outer":[pmin,pmax],"A13_lower_bound":low}


def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items()}
    if isinstance(v,list):return [stable(x) for x in v]
    if isinstance(v,bool):return v
    return float(f"{v:.10g}")


def solve():
    cfg=read("inputs.json")
    osc=read(cfg["oscillation"])
    Ue,y=charged_basis(cfg)
    # The published upper-octant NOvA two-width row is asymmetric; .56 ± .07 encloses precisely [.49,.63].
    atm=cfg["atmospheric_NOvA_upper_octant"]
    lo=atm["central"]-atm["sigmas"]*atm["low"]
    hi=atm["central"]+atm["sigmas"]*atm["high"]
    target=np.array([osc["DayaBay_sin2_2theta13"][0],osc["JUNO_sin2_theta12_NO"][0],(hi+lo)/2])
    hh=np.array([2*osc["DayaBay_sin2_2theta13"][1],2*osc["JUNO_sin2_theta12_NO"][1],(hi-lo)/2])
    ds=osc["KamLAND_dm21_eV2"][0];start=osc["T2K_dm32_NO_eV2"]
    scenarios={}
    for label,da in [("central_gaps",start[0]),("common_alternative",cfg["alternative_dm32_eV2"])]:
        assert start[0]-2*start[1]<=da<=start[0]+2*start[2]
        rows={}
        for eta in (-120,120):
            M,U,r=rank_two(ds,da,math.radians(eta),"small")
            rows[str(eta)]={"r":r,**search_point(Ue,M,U,target,hh,cfg["grid_side"])}
        scenarios[label]={"dm21":ds,"dm32":da,"small_root_by_eta_deg":rows}
    assert all(p["all_three_inside"] for scenario in scenarios.values() for p in scenario["small_root_by_eta_deg"].values())
    return stable({"charged_yukawas":y.tolist(),"abs_Ue":abs(Ue).tolist(),
                   "joint_windows":np.array([target-hh,target+hh]).T.tolist(),
                   "scenarios":scenarios,"large_root_full_frequency_box":large_bound(Ue,osc)})


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true");parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args();result=solve();out=HERE/"results.json"
    if args.write_new:
        with out.open("x",encoding="utf8") as f:f.write(json.dumps(result,indent=2)+"\n")
    else:
        assert json.loads(out.read_bytes())==result
        print("1659 PASS: real charged spectrum, joint phases, two signed sources, Takagi and all-phase large-root bound")

