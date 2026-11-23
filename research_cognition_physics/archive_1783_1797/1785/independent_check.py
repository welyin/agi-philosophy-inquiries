"""Project the full inclusive Z width into integer left and right quark-doublet ranks."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
SLO,SHI=.22,.24
MZ,SIGMZ=91.1876,.0021
GF,SIGGF=1.16637e-5,1.e-10   # LEP/SLD convention and muon-lifetime calibration
OBS=(1.7458-.0054,1.7458+.0054)
MATCH=(.97,1.10)  # entire pole-width transfer, with color counted explicitly


def scale(m,g):
    return 3*g*m**3/(3*math.sqrt(2)*math.pi)


def extrema(ranks):
    n=sum(ranks)
    k=sum(max(r-1,0) for r in ranks)
    free=sum(r in [1,2] for r in ranks)
    amin=(22-4*free)/9
    bmin=-(n+2*k)/3
    cmin=(n+k)/4
    lower_s=min(SHI,max(SLO,-bmin/(2*amin)))
    lower=amin*lower_s**2+bmin*lower_s+cmin
    def upper(s):
        return 22*s*s/9-(n+2*(k+free))*s/3+(n+k+free)/4
    upper_s=max([SLO,SHI],key=upper)
    return lower,upper(upper_s),lower_s,upper_s


def unitary(theta,phi):
    c,s=np.cos(theta),np.sin(theta)
    return np.array([[c,s*np.exp(1j*phi),0],[-s*np.exp(-1j*phi),c,0],[0,0,1]])


def project(rank,t,twist):
    v=np.array([math.sqrt(t),0,math.sqrt(1-t)])
    proj=np.zeros((3,3))
    if rank==3:
        proj=np.eye(3)
    elif rank in [1,2]:
        proj=np.outer(v,v)
        if rank==2:
            proj[1,1]=1
    return twist@proj@twist.conj().T


def trace_direct(ranks,s,t,t2):
    total=0.
    for r,x,ang in zip(ranks,[t,t2],[.38,.71]):
        pu=project(r,x,unitary(ang,.63))[:2,:2]
        pd=project(r,.37,unitary(ang,.27))
        up=.5*pu-(2*s/3)*np.eye(2)
        down=-.5*pd+s/3*np.eye(3)
        total+=np.vdot(up,up).real+np.vdot(down,down).real
    return total


def compute():
    c0=scale(MZ-2*SIGMZ,GF-2*SIGGF)
    c1=scale(MZ+2*SIGMZ,GF+2*SIGGF)
    groups={}
    details=[]
    for l in range(4):
        for r in range(4):
            lo,hi,s0,s1=extrema([l,r])
            assert abs(trace_direct([l,r],s0,4*s0/3,4*s0/3)-lo)<1e-14
            assert abs(trace_direct([l,r],s1,1,1)-hi)<1e-14
            details.append(dict(left=l,right=r,tree_GeV=[lo*c0,hi*c1],at_s=[s0,s1]))
    for n in range(7):
        subset=[s for s in details if s["left"]+s["right"]==n]
        lo=min(s["tree_GeV"][0] for s in subset)
        hi=max(s["tree_GeV"][1] for s in subset)
        groups[str(n)]=dict(tree_GeV=[lo,hi],pole_envelope_GeV=[MATCH[0]*lo,MATCH[1]*hi],
                            matching_overlap_outer=[OBS[0]/hi,OBS[1]/lo],
                            surviving_splits=[[s["left"],s["right"]] for s in subset
                                     if s["tree_GeV"][0]*MATCH[0]<=OBS[1]
                                     and s["tree_GeV"][1]*MATCH[1]>=OBS[0]])
    assert [i for i in groups if groups[i]["surviving_splits"]]==["3"]
    s=.23153
    # Split the same generation indices into complementary ordinary and mirror blocks.
    at_common=22*s*s/9+3*(.25-s/3)+2*(.25-2*s/3)
    universal=scale(MZ,GF)*at_common
    small_max=max(groups[str(i)]["tree_GeV"][1] for i in range(3))
    large_min=min(groups[str(i)]["tree_GeV"][0] for i in range(4,7))
    return dict(hadronic_observed_work_box_GeV=OBS,total_pole_transfer=MATCH,
                maximum_GF_MZ_scale_range_GeV=[c0,c1],rank_pairs=details,ranks=groups,
                strictly_sufficient_total_multiplier_range=[OBS[1]/large_min,OBS[0]/small_max],
                common_three_splits={"tree_GeV":universal,"multiplier_to_data":1.7458/universal})


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    payload=json.dumps(compute(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if check:
        assert json.loads(dest.read_bytes())==json.loads(payload)
    else:
        dest.write_text(payload,encoding="utf-8")
    print("1785 inclusive hadronic Z rank",hashlib.sha256(payload.encode()).hexdigest()[:12])
