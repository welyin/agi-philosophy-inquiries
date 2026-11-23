"""Two-family chirality and threshold audits for a single active charged-lepton row."""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def gamma_prefactor(mh):
    return mh/(8*math.pi)


def matrices(q,s,cl):
    b=math.hypot(q*math.sqrt(1-s*s),s)
    a=q/b
    d=(1-q*q)*s*math.sqrt(1-s*s)/b
    ca,cb,cc=cl
    return np.array([[0.,ca*a],[cb*b,cc*d]])


def decompose(y):
    U,values,Vh=np.linalg.svd(y)
    return U[:,::-1],values[::-1],Vh.T[:,::-1]


def observables(g,mh,reference_width):
    norm=abs(g[0,1])**2+abs(g[1,0])**2
    # e mass and mu mass are kinematically negligible at 125 GeV.
    excess=abs(g[0,1])**2-abs(g[1,0])**2
    return {"h_em_and_me":g[[0,1],[1,0]].tolist(),"LFV_norm":norm,
            "partial_width_GeV":gamma_prefactor(mh)*norm,"BR_if_fixed_total":gamma_prefactor(mh)*norm/reference_width,
            "ratio_to_same_h_mumu":norm/abs(g[1,1])**2,"chirality_asymmetry":excess/norm}


def charged_point(q,s,cl,rat,xi,c,v,m_mu,mh,ref_width):
    Y=matrices(q,s,cl)
    U,y,W=decompose(Y)
    P=np.diag([1.,0.])
    g_before=U.T@(c*np.eye(2)+2*xi*P)@Y@W
    factor=U.T@P@U
    np.testing.assert_allclose(g_before,(c*np.eye(2)+2*xi*factor)@np.diag(y),atol=2e-15)
    assert abs(factor[0,1]**2-factor[0,0]*factor[1,1])<1e-16
    # t_mu = 1, t_e = high ratio / measured pole ratio. Both preserve the mass directions.
    transfer=np.diag([y[1]/y[0]/rat,1])
    norm=m_mu/(v*y[1])
    g_R=norm*g_before@transfer
    g_L=norm*transfer@g_before
    common=norm*g_before.copy()
    # Local term with F(0)=0 and chosen F'(0) can cancel both FC entries while leaving the matched mass intact.
    derivative=np.diag(np.diag(g_R))-g_R
    safe_mass=norm*np.diag(y)@transfer
    step=1e-5
    def low(t):
        return safe_mass+t*(g_R+derivative)
    np.testing.assert_allclose((low(step)-low(-step))/(2*step),np.diag(np.diag(g_R)),atol=1e-13)
    np.testing.assert_allclose(g_R[0,1]/g_R[1,0],rat,rtol=1e-9)
    a,b=g_L[0,1],g_R[0,1]
    assert abs(a/b-transfer[0,0])<1e-13
    saved={tag:observables(g,mh,ref_width) for tag,g in [("right",g_R),("left",g_L),("common_vertex",common)]}
    return {"rotation_sine":abs(U[0,1]),"high_mass_ratio":y[1]/y[0],"t_e_over_t_mu":transfer[0,0],
            "right":saved["right"],"left":saved["left"],"common_vertex":saved["common_vertex"]}


def solve():
    mass=read("../1480/inputs.json")
    charm=read("../1476/inputs.json")
    cfg=read("../1486/inputs.json")
    old=read("../1486/results.json")["textures"]
    test=read("../1653/inputs.json")
    q=mass["lattice"]["md"][0]/mass["lattice"]["ms"][0]
    s=next(row[1] for row in charm["edges"] if row[0]=="cd")
    rat=mass["muonium"]["pole_ratio"][0]
    v,st=test["v_GeV"],test["sin_theta"]
    xi=st*v/test["f_GeV"];c=math.sqrt(1-st*st)
    mh,m_mu,total,cms=125.,.1056583755,.00407,4.4e-5
    rows={}
    for name,cl in cfg["clebsch"].items():
        point=charged_point(q,s,cl,rat,xi,c,v,m_mu,mh,total)
        hist=old[name]
        assert abs(point["rotation_sine"]-math.sin(math.radians(hist["real_minus_witness"]["rotation_degrees"])))<1e-12
        lo,hi=hist["sine_rotation"]
        amin,amax=mass["muonium"]["pole_ratio"][0]-2*mass["muonium"]["pole_ratio"][1],mass["muonium"]["pole_ratio"][0]+2*mass["muonium"]["pole_ratio"][1]
        def dfrac(a):
            return a*a*(1-a*a)
        assert 0<lo<hi<1/math.sqrt(2)
        Ynorm=(2*xi*m_mu/v)**2
        # Reuse the certified high-scale angle and mass-ratio intervals as outer bounds. No statistical factorization.
        ylo,yhi=hist["lepton_ratio"]; tlo,thi=ylo/amax,yhi/amin
        pref=gamma_prefactor(mh)/total*Ynorm
        rr=[pref*dfrac(lo)*(1+amax**-2),pref*dfrac(hi)*(1+amin**-2)]
        ll=[pref*dfrac(lo)*(tlo*tlo+yhi**-2),pref*dfrac(hi)*(thi*thi+ylo**-2)]
        vv=[pref*dfrac(lo)*(1+yhi**-2),pref*dfrac(hi)*(1+ylo**-2)]
        ratio=[4*xi*xi*dfrac(lo)*(1+amax**-2)/(c+2*xi*hi*hi)**2,
               4*xi*xi*dfrac(hi)*(1+amin**-2)/(c+2*xi*lo*lo)**2]
        for tag,pred in [("right",rr),("left",ll),("common_vertex",vv)]:
            assert pred[0]<point[tag]["BR_if_fixed_total"]<pred[1]
        rows[name]={"central":point,"right_BR_window":rr,"left_BR_outer":ll,"common_BR_outer":vv,
                    "right_ratio_to_h_mumu_outer":ratio,"CMS_to_right_upper":cms/rr[1],
                    "old_mass_accepts_central":cfg["mass_matching_window"][0]<1/point["t_e_over_t_mu"]<cfg["mass_matching_window"][1],
                    "old_neutrino_accepts_a_branch":any(mode["overlap_possible"] for mode in hist["branches"].values())}
    for key in ("right_BR_window","left_BR_outer","common_BR_outer"):
        assert rows["small_left"][key][1]<rows["large_left"][key][0]
    return {"benchmark":{"xi":xi,"h_doublet_component":c,"Gamma_total_GeV_assumed":total,"CMS_BR_reference_upper":cms,
                        "physical_chirality_asymmetry_right":(rat**2-1)/(rat**2+1)},
            "placements":rows}


def stable(obj):
    if isinstance(obj,dict): return {k:stable(v) for k,v in obj.items()}
    if isinstance(obj,list): return [stable(v) for v in obj]
    if isinstance(obj,bool):return obj
    return float(f"{obj:.10g}")


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--check",action="store_true")
    parsed=args.parse_args()
    result=stable(solve())
    path=HERE/"independent_results.json"
    if parsed.check:
        assert result==json.loads(path.read_bytes())
        print("1655 independent PASS: row projector, chirality, three thresholds, two-family holdout widths")
    else:
        with path.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")

