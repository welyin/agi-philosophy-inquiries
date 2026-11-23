"""Two tagged-jet flavour slots and beta-endpoint complement to 1457."""
import argparse
from itertools import product
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def check():
    cfg=read(HERE/"inputs.json"); prev=read(HERE/cfg["oscillation"])
    mass=read(HERE/cfg["solar_and_mass_gaps"])["juno_207_days"]
    solno=read(HERE/cfg["solar_no"])["NO_solar_sin2"]
    mb=read(HERE.parent/"1370"/"inputs.json")["katrin"]["m_beta_upper_eV"]
    zbox=lambda a,w:((1-math.sqrt(1-(a-2*w)))/2,(1-math.sqrt(1-(a+2*w)))/2)
    z=zbox(*prev["DayaBay_amplitude"])
    ma,maerr=mass["solar_21"]
    solar=(ma-2*maerr,ma+2*maerr)
    result={}
    for ordering,ref in [("NO",solno),("IO",prev["IO"]["solar_sine_squared"])]:
        center,da=ref
        atm=mass["normal_atmospheric_31" if ordering=="NO" else "inverted_abs_31"]
        # Published registration stores (central,minus,plus); use larger for a safe symmetric enclosure.
        abox=(atm[0]-2*max(atm[1:]),atm[0]+2*max(atm[1:]))
        deltas=[]
        for a,t,small,big in product((center-2*da,center+2*da),z,solar,abox):
            correction=small*(1-3*a*(1-t))
            deltas.append(correction+(1 if ordering=="NO" else -1)*big*(1-3*t))
        f=[1/(3+d/mb**2) for d in deltas]
        if ordering=="NO": assert max(f)<1/3
        else: assert min(f)>1/3
        result[ordering]={"signed_mass_delta_eV2":[min(deltas),max(deltas)],
            "one_sided_electron_fraction_edge_at_KATRIN_cap":max(f) if ordering=="NO" else min(f)}
    result["finite_forbidden_fraction_gap"]=result["IO"]["one_sided_electron_fraction_edge_at_KATRIN_cap"]-result["NO"]["one_sided_electron_fraction_edge_at_KATRIN_cap"]
    # For a column frame A and a single omitted b direction v, target A^dagger v=z.
    K=np.array([[.65, .3j],[.30,.5],[.04j,.10]])
    G=K.conj().T@K
    r=np.array([.14,.04])
    z1=np.sqrt(r)*[1,np.exp(-1j*np.angle(G[0,1]))]
    root=np.linalg.solve(G,z1)
    v0=K@root
    u=np.linalg.svd(K,full_matrices=True)[0][:,2]
    assert np.vdot(v0,v0).real<1
    v=v0+u*math.sqrt(1-np.vdot(v0,v0).real)
    assert np.isclose(np.linalg.norm(v),1) and np.allclose(abs(K.conj().T@v)**2,r)
    area=np.linalg.det(G).real
    def required(v):
        return (G[1,1].real*v[0]+G[0,0].real*v[1]-2*abs(G[0,1])*math.sqrt(v[0]*v[1]))/area
    assert np.isclose(required(r),np.vdot(v0,v0).real)
    corners=np.diag(G).real
    assert required(corners)>1
    result["two_bottom_flavour_rates"]={"feasible_target":r.tolist(),"required_norm_squared":float(required(r)),
        "individual_maxima":corners.tolist(),"joint_norm_squared_if_both_maxima":float(required(corners))}
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(); g=p.add_mutually_exclusive_group()
    g.add_argument("--write",action="store_true");g.add_argument("--check",action="store_true")
    opt=p.parse_args(); data=json.dumps(check(),indent=2)+"\n"; dest=HERE/"independent_results.json"
    if opt.write:
        with dest.open("x",encoding="utf-8") as f: f.write(data)
    elif data!=dest.read_text(encoding="utf-8"): raise SystemExit("1457 mismatch")
    else: print("1457 independent OK")
