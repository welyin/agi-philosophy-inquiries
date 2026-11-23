"""Recover physical CKM and bound two general dispersive B-mixing corrections."""
import argparse
import cmath
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
def read(p): return json.loads((HERE/p).read_bytes())
OLD=runpy.run_path(str(HERE/"../1732/check.py"))
SEM=read("../1702/inputs.json")["direct_semileptonic_moduli"]
ANCHOR=runpy.run_path(str(HERE/"../1707/independent_check.py"))["DATA"]
DATA={"delta_md":[.5068,.0019],"delta_ms":[17.766,.006],"golden_S":[.710,.011],
      "phi_s_ccs":[-.053,.013],"mass_d":ANCHOR["mass_d_MeV"],"mass_s":ANCHOR["mass_s_MeV"],
      "FNAL_MILC":[1.206,.018,.006],"HPQCD":[1.216,.016]}

def ckm(offsets,gamma):
    us,ub,cb=[SEM[k][0]+offsets[n]*math.hypot(*SEM[k][1:]) for k,n in
              (("Vus","us"),("Vub","ub"),("Vcb","cb"))]
    c13=math.sqrt(1-ub**2)
    s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12*s12),math.sqrt(1-s23*s23)
    def matrix(delta):
        z=cmath.exp(1j*delta)
        return np.array([[c12*c13,s12*c13,ub*z.conjugate()],
                         [-s12*c23-c12*s23*ub*z,c12*c23-s12*s23*ub*z,s23*c13],
                         [s12*s23-c12*c23*ub*z,-c12*s23-s12*c23*ub*z,c23*c13]])
    def angle(V):
        return np.angle(-V[0,0]*V[0,2].conjugate()/(V[1,0]*V[1,2].conjugate()))
    # The physical gamma is a quartet; solve it directly (the fitted gamma is not delta).
    lo,hi=0.,np.pi
    for _ in range(60):
        mid=(lo+hi)/2
        if angle(matrix(mid))<math.radians(gamma): lo=mid
        else: hi=mid
    V=matrix((lo+hi)/2)
    assert abs(math.degrees(angle(V))-gamma)<1e-11 and np.linalg.norm(V@V.conj().T-np.eye(3))<1e-14
    beta=np.angle(-V[1,0]*V[1,2].conjugate()/(V[2,0]*V[2,2].conjugate()))
    beta_s=np.angle(-V[2,1]*V[2,2].conjugate()/(V[1,1]*V[1,2].conjugate()))
    return {"us_ub_cb":[us,ub,cb],"gamma_deg":gamma,"pdg_delta_deg":math.degrees((lo+hi)/2),
            "ratio":abs(V[2,0]/V[2,1]),"angle_d_rad":2*beta,"angle_s_rad":-2*beta_s,
            "golden_S":math.sin(2*beta)}

def interval(name,w):
    mean,*errors=DATA[name]
    return np.array([mean-w*math.hypot(*errors),mean+w*math.hypot(*errors)])

def moduli(f,w,group):
    fd,fs=interval("delta_md",w),interval("delta_ms",w)
    md,ms=interval("mass_d",w),interval("mass_s",w)
    kernel=np.sqrt(fd/fs[::-1]*ms/md[::-1])
    xi=f["ratio"]/kernel[::-1]
    # R=|Delta_d/Delta_s|=(xi_lattice*kernel/CKM)^2.
    R=(interval(group,w)*kernel/f["ratio"])**2
    return xi,R

def nearest(interval):
    a,b=interval
    return min(max(0.,a),b)

def correction(f,w,group,penguin_deg=0):
    need_x,R=moduli(f,w,group)
    pad=math.radians(penguin_deg)*np.array([-1.,1.])
    sin=interval("golden_S",w)
    angle=np.arcsin(sin)
    directions=[nearest(val-f["angle_d_rad"]+pad) for val in (angle,np.pi-angle[::-1])]
    d=min(directions,key=abs)
    s=nearest(interval("phi_s_ccs",w)-f["angle_s_rad"]+pad)
    def radii(e,phi):
        if e<abs(math.sin(phi)) or math.cos(phi)<0: return None
        a=math.sqrt(max(0,e*e-math.sin(phi)**2))
        return [math.cos(phi)-a,math.cos(phi)+a]
    def feasible(e):
        d1,d2=radii(e,d),radii(e,s)
        if d1 is None or d2 is None: return None
        lo=max(d1[0]/d2[1],R[0])
        hi=min(d1[1]/d2[0],R[1])
        if lo>hi: return None
        ratio=(lo+hi)/2
        joint=[max(d1[0]/ratio,d2[0]),min(d1[1]/ratio,d2[1])]
        rs=sum(joint)/2
        return ratio*rs*cmath.exp(1j*d),rs*cmath.exp(1j*s)
    lo,hi=0.,.9
    assert feasible(hi)
    for _ in range(65):
        mid=(lo+hi)/2
        if feasible(mid): hi=mid
        else: lo=mid
    D,S=feasible(hi)
    assert max(abs(D-1),abs(S-1))<=hi+2e-14
    return {"xi_required_if_zero_correction":need_x.tolist(),"ratio_new_moduli":R.tolist(),
            "min_max_abs_Delta_minus_one":hi,"witness_amplitudes":[[z.real,z.imag] for z in (D,S)],
            "weak_phase_changes_degrees":[math.degrees(d),math.degrees(s)]}

def point(name,offsets,gamma):
    f=ckm(offsets,gamma)
    return {"weak_offsets":offsets,"prediction":f,
            "relative_to_golden_LHCb_individual":(f["golden_S"]-ANCHOR["golden_S"][0])/math.hypot(*ANCHOR["golden_S"][1:]),
            "alternatives":{group:{f"width{w}_decay_deg{p}":correction(f,w,group,p) for w,p in
                                 ((0,0),(1,0),(1,1),(2,0),(2,1))} for group in ("FNAL_MILC","HPQCD")}}

def calculate():
    key="registered_with_lower_us"
    gamma=read("../1732/results.json")["exploratory_k3"][key]["gamma_MZ_deg"]
    out={"1732":point(key,OLD["WEAK"][key],gamma)}
    # Re-solve the independent nested charge-sector equations on every true
    # R1736 witness, including any alternative K_l3 input and shifted y_u.
    for name,row in read("results.json")["witnesses"].items():
        weak=row["weak_offsets"]
        OLD["WEAK"][name]=weak
        OLD["PARENT"]["FLOW"]["low"]["yukawa"]["up"][0][0]=row["up_MZ"]
        charged=OLD["solve"](3,name,OLD["CONFIG"]["log_step"])
        assert abs(charged["gamma_MZ_deg"]-row["gamma"])<2e-7
        out[name]=point(name,weak,charged["gamma_MZ_deg"])
        np.testing.assert_allclose(out[name]["prediction"]["us_ub_cb"],[row["moduli_MZ"][k]
                                   for k in ("Vus","Vub","Vcb")],rtol=5e-10)
    return OLD["PARENT"]["tidy"]({"data":DATA,"points":out})

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    res=calculate()
    path=HERE/"independent_results.json"
    if check:
        assert res==json.loads(path.read_bytes())
        print("1736 independent B-mixing complex disks and CKM reconstruction verified")
    else:
        path.write_text(json.dumps(res,indent=2,allow_nan=False)+"\n",encoding="utf-8")
