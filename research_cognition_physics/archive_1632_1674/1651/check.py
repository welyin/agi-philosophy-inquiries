"""Renormalizable heavy-doublet matching of the two measured radial flavor contrasts."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def pair(z):
    return [float(z.real),float(z.imag)]


def texture(cfg):
    data=read(cfg["texture"])["reconstruction"]
    a,b,c,e,f=data["normalized_real_edges_A_B_C_E_F"]
    r,d=data["two_real_coefficients"]
    ph=d*cmath.exp(1j*math.radians(data["eta_example_deg"]))
    Y=np.array([[0,a,0],[b,c,r+ph],[0,e,f]],dtype=complex)
    yb=read(cfg["ordinary_scale"])["yukawa"]["down"][2][0]
    L,sv,RH=np.linalg.svd(Y*yb);o=np.argsort(sv)
    return Y*yb,yb,r,ph,L[:,o],sv[o],RH.conj().T[:,o]


def potential(cfg,M):
    v=cfg["v_GeV"];mh=cfg["mh_GeV"];v2=v*v
    cref=cfg["reference_cosine"];sref=math.sqrt(1-cref**2)
    mix6=(mh**2-cfg["reference_heavy_mass_GeV"]**2)*cref*sref/v2
    theta=.5*math.asin(2*(-mix6)*v2/(M*M-mh*mh))
    s,c=math.cos(theta),math.sin(theta)
    z1=(s*s*mh**2+c*c*M**2)/v2
    d22=c*c*mh**2+s*s*M**2
    z5=(d22-M**2)/v2;z4=z5
    z2,z3,z7=cfg["Z2"],cfg["Z3"],cfg["Z7"]
    y2=M**2-z3*v2/2
    mass=np.array([[z1*v2,mix6*v2],[mix6*v2,d22]])
    eigen=np.array([[s,c],[c,-s]])
    np.testing.assert_allclose(eigen@mass@eigen.T,np.diag([mh**2,M**2]),atol=1e-9)
    y3=-mix6*v2/2;y1=-z1*v2/2
    assert z7==0 and y2>0 and min(z1,z2)>0
    # Complete the square in x-v²/2 and Re z; use |z|² <= xy.
    margin=z3-max(0.,2*mix6**2/z1-z4-z5,-z4+z5)
    assert margin>0
    static=-mix6*v2/d22
    pole=-mix6*v2/(d22-mh**2)
    np.testing.assert_allclose(pole,c/s)
    np.testing.assert_allclose(pole/static,d22/(d22-mh**2))
    # Independently solve the real heavy-field stationary equation near the vacuum.
    z345=z3+z4+z5
    def stat(x):
        # x=sqrt(2)*Re H1^0, y=sqrt(2)*Re H2^0, with canonical kinetic terms.
        roots=np.roots([z2/2,0.,y2+z345*x*x/2,x*(y3+mix6*x*x/2)])
        rr=[float(i.real) for i in roots if abs(i.imag)<1e-12]
        assert len(rr)==1
        return rr[0]
    step=.0001*v
    np.testing.assert_allclose((stat(v+step)-stat(v-step))/(2*step),static,rtol=1e-7)
    return {"M_H_A_charged_GeV":M,"mh_GeV":mh,"s":s,"c":c,"Y1_Y2_Y3_GeV2":[y1,y2,y3],
       "Z1_to_Z7":[z1,z2,z3,z4,z5,mix6,z7],"M22_GeV2":d22,
       "global_vacuum_sufficient_margin":margin,"static_response":static,
       "pole_response_per_hVV":pole,"pole_over_static":pole/static,
       "relative_error_static_substituted_for_normalized_pole":1-static/pole,
       "H_to_VV_over_SM_same_mass":c*c}


def project(raw,L,R):
    return L.conj().T@raw@R


def observe(K,m):
    terms=np.diag(K)/m
    def width(i,j):
        return (abs(K[i,j])**2+abs(K[j,i])**2)/m[-1]**2
    return {"p":terms.imag.tolist(),"scalar":terms.real.tolist(),
            "ds_db_sb_over_SM_bb":[width(0,1),width(0,2),width(1,2)]}


def match(cfg,pot,model,reference_cosine=None):
    Y,yb,r,d,L,m,R=texture(cfg)
    s,c=pot["s"],pot["c"]
    cref=c if reference_cosine is None else reference_cosine
    sref=math.sqrt(1-cref**2)
    qr,qd=cfg["example_targets"][model]
    # Chosen aligned term also restores the ordinary CP-even diagonal at the reference pole.
    aligned=(1-sref)/cref
    extra=yb*(qr*r+qd*d)/cref
    rho=aligned*Y.copy()
    rho[1,2]+=extra
    h=s*Y+c*rho;H=c*Y-s*rho
    g,gh,gH,gA=[project(t,L,R) for t in (Y,h,H,1j*rho)]
    np.testing.assert_allclose(s*gh+c*gH,g,atol=1e-15)
    np.testing.assert_allclose(c*gh-s*gH,project(rho,L,R),atol=1e-15)
    hdata,Hdata=observe(gh,m),observe(gH,m)
    np.testing.assert_allclose(np.array(hdata["p"]),-c/s*np.array(Hdata["p"]),atol=2e-12)
    old=read(cfg["target"])["sample"][model]
    np.testing.assert_allclose(np.array(hdata["p"]),np.array(old["diagonal_p"])*c/cref,rtol=2e-9)
    if abs(c-cref)<1e-13:
        np.testing.assert_allclose(hdata["scalar"],old["diagonal_scalar"],rtol=2e-9)
    # Recover the nonaligned weak edge from the heavy coupling and a fixed signed common basis.
    reconstructed=L@(c*g-gH)@R.conj().T/s-aligned*Y
    edge=reconstructed[1,2]/yb
    ud=edge.imag/d.imag
    ur=(edge.real-ud*d.real)/r
    np.testing.assert_allclose([c*ur,c*ud],[qr*c/cref,qd*c/cref],atol=5e-13)
    # Distinct matching momenta. A has opposite imaginary signs for conjugate chiralities.
    mh=cfg["mh_GeV"];M=pot["M_H_A_charged_GeV"]
    same=1+(s*s-1)*mh**2/(c*c*M**2)
    opposite=1+(s*s+1)*mh**2/(c*c*M**2)
    for i,j in [(0,1),(0,2),(1,2)]:
        qsame=gh[i,j]**2/mh**2+(gH[i,j]**2+gA[i,j]**2)/M**2
        qlr=gh[i,j]*gh[j,i].conjugate()/mh**2
        qtotal=qlr+(gH[i,j]*gH[j,i].conjugate()+gA[i,j]*gA[j,i].conjugate())/M**2
        np.testing.assert_allclose(qsame/(gh[i,j]**2/mh**2),same,atol=2e-11)
        np.testing.assert_allclose(qtotal/qlr,opposite,atol=2e-11)
    return {"additional_weak_23_rho":pair(extra),"rho_flavor_norm":abs(extra),
            "aligned_ordinary_Y_coefficient":aligned,
            "reconstructed_physical_qR_qD":[c*ur,c*ud],"h":hdata,"H":Hdata,
            "same_chirality_neutral_tree_over_h_only":same,
            "LR_chirality_neutral_tree_over_h_only":opposite}


def degenerate_observables(cfg,pot):
    Y,yb,r,d,L,m,R=texture(cfg)
    s,c=pot["s"],pot["c"]
    aligned=(1-s)/c
    rows=[]
    for e in (2*yb*d/c,2*yb*abs(d)/c):
        rho=aligned*Y.copy();rho[1,2]+=e
        ph=observe(project(s*Y+c*rho,L,R),m)
        pH=observe(project(c*Y-s*rho,L,R),m)
        rows.append({"additional_23_abs":abs(e),"h":ph,"H":pH})
    for boson in ("h","H"):
        np.testing.assert_allclose(rows[0][boson]["ds_db_sb_over_SM_bb"],rows[1][boson]["ds_db_sb_over_SM_bb"],rtol=1e-11)
    assert abs(rows[0]["h"]["p"][0]-rows[1]["h"]["p"][0])>.1
    return rows


def normalize(x):
    if isinstance(x,dict):return {k:normalize(v) for k,v in x.items()}
    if isinstance(x,list):return [normalize(y) for y in x]
    if isinstance(x,float):return 0. if abs(x)<1e-15 else float(f"{x:.11g}")
    return x


def main():
    p=argparse.ArgumentParser();p.add_argument("--check",action="store_true");args=p.parse_args()
    cfg=read("inputs.json");obs,err=cfg["ATLAS_mass_GeV"]
    assert abs(cfg["mh_GeV"]-obs)<=err
    ref=potential(cfg,cfg["reference_heavy_mass_GeV"])
    high=potential(cfg,cfg["second_heavy_mass_GeV"])
    out={"600":ref,"1000_fixed_Z6":high,"models":{}}
    Y,yb,r,d,L,m,R=texture(cfg)
    budget=cfg["scale_budget"]
    for name,(qr,qd) in cfg["example_targets"].items():
        hit=yb*abs(qr*r+qd*d)
        # The budget applies to the additional nonaligned edge, independently of aligned Y.
        bound=math.sqrt(cfg["mh_GeV"]**2+cfg["v_GeV"]**2*budget["Z6_abs_max"]*
                        budget["rho_non_aligned_23_abs_max"]/(math.sqrt(budget["hVV_squared_min"])*hit))
        out["models"][name]={"reference":match(cfg,ref,name),
            "1000_same_Yukawas":match(cfg,high,name,reference_cosine=ref["c"]),
            "1000_recalibrated":match(cfg,high,name),"conditional_upper_M_GeV_for_same_h_response":bound}
    out["equal_untagged_FC_widths"]=degenerate_observables(cfg,ref)
    out=normalize(out)
    if args.check: assert out==read("results.json")
    else: (HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    for key in ("600","1000_fixed_Z6"):
        v=out[key];print(key,v["s"],v["c"],v["pole_over_static"],v["global_vacuum_sufficient_margin"])
    for name,x in out["models"].items():print(name,x["conditional_upper_M_GeV_for_same_h_response"])

if __name__=="__main__":
    main()
