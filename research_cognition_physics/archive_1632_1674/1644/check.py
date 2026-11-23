"""Neutral Bd -> D_CP h0 signed holdout on the physical round-1643 branches."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def width(a):
    return [a[0]-2*sum(a[1:]),a[0]+2*sum(a[1:])]


def phase_mod(z):
    return cmath.phase(z)


def short_distance(x,y):
    if x==y:
        return x*(4-11*x+x*x)/(4*(1-x)**2)-3*x**3*math.log(x)/(2*(1-x)**3)
    def term(a):
        return (4-8*a+a*a)*math.log(a)/(4*(1-a)**2)
    return x*y*((term(x)-term(y))/(x-y)-3/(4*(1-x)*(1-y)))


def budgets(theta,coeff,kweak):
    if math.sin(theta)>=coeff: return {"phase_rad":0.}
    amin=math.asin(coeff)
    distance=min(abs(cmath.phase(cmath.exp(1j*(theta-a)))) for a in (amin,math.pi-amin))
    # If the nearest phase ray is more than pi/2 away, the nearest point of that ray is the origin.
    rel=math.sin(min(math.pi/2,distance))
    r=math.sin(distance/2)
    return {"phase_rad":distance,"phase_deg":math.degrees(distance),"mixing_only_min_relative":rel,
            "decay_only_min_ratio_each_side_cap":r,"hadronic_ratio_if_only_wrong_B_tree":r/kweak}


def witness(label,M,F,angles,tree,loop,obs,build):
    p=dict(tree);p["ordinary_s12_s23_s13_delta_deg"]=angles;p["F"]=[F.real,F.imag]
    mass_matrix,m,K,_,_,_,_,_=build(p,M)
    v_left,s,vh=np.linalg.svd(mass_matrix)
    order=np.argsort(s)
    v_left=v_left[:,order];v_right=vh.conj().T[:,order]
    # Compare current magnitudes so SVD rephasing of individual columns is unrestricted.
    assert np.max(abs(abs(v_left[:3,:3])-abs(K[:,:3])))<1e-8
    K=v_left[:3,:3]
    assert max(abs(v_right[4,0]),abs(v_left[4,0]),abs(v_left[3,0]))<3e-9
    lam=K[:,0].conjugate()*K[:,2];lam_ds=K[:,0].conjugate()*K[:,1]
    assert max(abs(sum(lam)),abs(sum(lam_ds)))<2e-11
    B=K[1,2]*K[0,0].conjugate()
    D=K[1,1].conjugate()*K[0,1]
    dk=B*D
    top=lam[2]**2*loop[2,2]; full=lam@loop@lam
    beta2=phase_mod(top/lam[1]**2)
    lead=phase_mod(top/dk**2); total=phase_mod(full/dk**2)
    rw=abs(K[0,2]*K[1,0]/B)
    return {"source":label,"M_GeV":M,"masses_GeV":m.tolist(),"beta_top_deg":math.degrees(beta2/2),
            "arg_DD_over_charm_axis_shift_rad":phase_mod((lam[1]/dk)**2),
            "GIM_shift_phase_rad":phase_mod(full/top),"relative_GIM_box_modulus":abs(full/top-1),
            "theta_top_KK_rad":lead,"theta_LO_KK_rad":total,"S_reduced":math.sin(total),
            "direct_C_single_weak_amplitude":0.,"wrong_to_right_CKM":rw,
            "one_width_needed":budgets(total,width_1(obs),rw),
            "two_width_needed":budgets(total,width(obs)[0],rw)}


def width_1(row):
    return row[0]-sum(row[1:])


def domain(cfg,loop):
    inp=read(cfg["first_row"]); arch=read(cfg["tree_domains"]); rad=math.radians
    ud=width(inp["ud_beta"]); us=[x/y for x,y in zip(width(inp["K_product"]),reversed(width(inp["K_form"])))]
    rho=width(inp["ub_over_cb"])
    tree3=read(cfg["tree"]); tab=read(tree3["semileptonic"])
    cbmax=max(width(tab["BelleII_Vcb"])[1],width(arch["inclusive_holdout"]["Vcb"])[1])
    pol=read("../1633/results.json"); rb=math.sqrt(pol["right_doublet_probability"])
    # sqrt(G) differs from identity only in its s,b block.
    gmin=1-2*pol["left_difference"]/pol["k"]-(.1*.0027/(1600*rb))**2-(4.8*rb/1600)**2
    pert=1-math.sqrt(gmin)
    s13=cbmax*rho[1]+pert
    c13=math.sqrt(1-s13*s13)
    s12=(us[1]+pert)/c13; s23=(cbmax+pert)/c13
    td=s12*s23+s13
    ts=s23+s12*s13+pert
    cd=math.sqrt(1-ud[1]**2-td**2)
    cdmax=math.sqrt(1-ud[0]**2)
    rlow=ud[0]*rho[0]/cdmax;rhigh=ud[1]*rho[1]/cd
    wrongweak=rho[1]*cdmax/ud[0]
    tag=math.asin(td*ts/(ud[0]*us[0]))
    g=arch["gamma_Bpm"]
    gam=[rad(g["deg"]-2*g["down"]-2*g["transfer_cap_deg"])-tag,
         rad(g["deg"]+2*g["up"]+2*g["transfer_cap_deg"])+tag]
    # Unitless cap from the full LO GIM boxes after top normalization.
    ratios=[rhigh/(1-rhigh),1/(1-rhigh),1.]
    gcap=sum(ratios[i]*ratios[j]*abs(loop[i,j]/loop[2,2])
             for i,j in itertools.product(range(3),repeat=2) if (i,j)!=(2,2))
    boxes={}
    for name,shift in (("positive",0.),("mirror",math.pi)):
        gs=[a+shift for a in gam]; corners=[]
        for r in (rlow,rhigh):
            angles=gs.copy()
            angles.extend(a for a in (math.acos(r),2*math.pi-math.acos(r)) if gs[0]<=a<=gs[1])
            corners.extend(2*phase_mod(1-r*cmath.exp(-1j*a)) for a in angles)
        th=[min(corners)-2*tag-math.asin(gcap),max(corners)+2*tag+math.asin(gcap)]
        assert max(map(abs,th))<math.pi/2
        boxes[name]={"theta_LO_kk_bounds_rad":th,"reduced_S_LO_range":[math.sin(a) for a in th]}
        if shift:
            boxes[name]["two_width_needed"]=budgets(th[1],width(cfg["BaBar_Belle_2015"]["minus_eta_S"])[0],wrongweak)
    return {"r_bounds":[rlow,rhigh],"abs_Ktd_ts_caps":[td,ts],
            "gamma_s_minus_gamma_d_or_charm_axis_cap_rad":tag,
            "LO_GIM_relative_modulus_cap":gcap,"wrong_tree_CKM_cap":wrongweak,"branches":boxes}


def round_values(x):
    if isinstance(x,dict):return {k:round_values(v) for k,v in x.items()}
    if isinstance(x,list):return [round_values(v) for v in x]
    if isinstance(x,(float,np.floating)):return float(f"{x:.10g}")
    return x


def solve():
    cfg=read("inputs.json");tc=read(cfg["tree"]);tree=read(cfg["tree_domains"]);weak=read(cfg["weak"])
    masses=cfg["up_charm_mass_GeV"]+[weak["short_distance"]["mtMS"][0]]
    xs=[(m/weak["const_GeV"]["MW"])**2 for m in masses]
    kernel=np.array([[short_distance(x,y) for y in xs] for x in xs])
    build=runpy.run_path(str(HERE.parent/"1640"/"check.py"))["build"]
    obs=cfg["BaBar_Belle_2015"]["minus_eta_S"]; points={}
    for name,vals in tc["physical_witnesses"].items():
        p=dict(zip(tc["witness_keys"],vals))
        ang=[p["s12"],p["s23"],p["s13"],math.degrees(p["delta_rad"])]
        points[name]=witness(name,p["M_GeV"],complex(p["ReF"],p["ImF"]),ang,tree,kernel,obs,build)
    for row in read(cfg["archived_mirror"])["witnesses"]:
        key="mirror_"+str(int(row["M"]))
        points[key]=witness(key,row["M"],complex(*row["F_in_input_basis"]),row["U_input"],tree,kernel,obs,build)
    bound=domain(cfg,kernel)
    assert bound["branches"]["mirror"]["reduced_S_LO_range"][1]<0
    assert width(obs)[0] < points["positive_sine"]["S_reduced"] < width(obs)[1]
    assert points["negative_sine"]["S_reduced"] < 0
    assert width(cfg["BaBar_Belle_2015"]["C"])[0]<0<width(cfg["BaBar_Belle_2015"]["C"])[1]
    return round_values({"data_two_width_S":width(obs),"data_two_width_C":width(cfg["BaBar_Belle_2015"]["C"]),
        "GIM_loop":kernel.tolist(),"full_tree_domain":bound,"physical_points":points})


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true");parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args();data=solve();target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as out:out.write(json.dumps(data,indent=2)+"\n")
    else:
        assert json.loads(target.read_text(encoding="utf-8"))==data
    print("1644 PASS: five physical mass matrices, Bd closure, charm-tag and GIM phases, full mirror-domain budget")

