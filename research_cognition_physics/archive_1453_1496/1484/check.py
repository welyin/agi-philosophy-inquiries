"""Recompute rank-two singular invariants and enclose every integer/source branch."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SPEC = json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))


def window(key):
    center, minus, plus = SPEC[key]
    n = SPEC["widths"]
    return center-n*minus, center+n*plus


def source_ratio(q, k, cosine):
    length = 3+2*k*k
    s = 6*(q+1/q)-8*k*k*cosine
    if s < 4*length-1.e-11:
        return None
    z = s+math.sqrt(max(0,(s-4*length)*(s+4*length)))
    return 2*length**2/z, z/8


def points(q, k, cosine):
    roots=source_ratio(q,k,cosine)
    if roots is None:
        return []
    length=3+2*k*k
    return [(r, (length/(6*r)-1/(3*q))/(q-1/q)) for r in roots]


def enclosing(q_low,q_high,k,cosine,branch):
    """On each interval the two branches of the positive roots are monotone."""
    length=3+2*k*k
    cutoff=2+4*k*k*(1+cosine)/3
    q_min=(cutoff+math.sqrt(cutoff**2-4))/2
    if q_high < q_min:
        return None
    low=max(q_low,q_min)
    edges=np.linspace(low,q_high,SPEC["subintervals"]+1)
    minimum,maximum=1.,0.
    for a,b in zip(edges[:-1],edges[1:]):
        ra,rb=source_ratio(a,k,cosine)[branch],source_ratio(b,k,cosine)[branch]
        num_low=max(0.,length/(6*max(ra,rb))-1/(3*a))
        num_high=max(0.,length/(6*min(ra,rb))-1/(3*b))
        minimum=min(minimum,num_low/(b-1/b))
        maximum=max(maximum,num_high/(a-1/a))
    # The electron weight on the two-dimensional massive plane is one third.
    pad=1.e-10
    p_low=max(0.,minimum-pad)
    p_high=min(1/3,maximum+pad)
    return [4*p_low*(1-p_low),4*p_high*(1-p_high)]


def branch_table(q_low,q_high,measured):
    table={}
    for label,c in [("trivial",1.),("nontrivial",-.5)]:
        k_max=math.floor(math.sqrt(3*(q_high+1/q_high-2)/(4*(1+c))))
        for k in range(k_max+1):
            for branch,label_r in enumerate(("small","large")):
                bounds=enclosing(q_low,q_high,k,c,branch)
                if bounds is None:
                    continue
                key=f"{label}_k{k}_{label_r}"
                table[key]={"amplitude_enclosure":bounds,
                    "overlap_possible":bool(max(measured[0],bounds[0])<=min(measured[1],bounds[1]))}
    return table


def realisation(q,k,cosine,branch,phase_sign=1):
    r,p=points(q,k,cosine)[branch]
    a=np.array([0.,1,1]); b0=np.array([1.,1,-1]); b=b0+k*a
    routes=[(0,1),(1,1),(2,-1)]+[(j,1) for _ in range(k) for j in (1,2)]
    transfer=np.zeros((len(routes),3))
    for path,(j,sign) in enumerate(routes):
        transfer[path,j]=sign/math.sqrt(len(routes))
    erased=transfer.sum(axis=0)/math.sqrt(len(routes))
    np.testing.assert_allclose(erased,b/len(routes))
    enforce=np.exp(1j*phase_sign*math.acos(cosine))
    mass=r*np.outer(a,a)+enforce*np.outer(b,b)
    values,vectors=np.linalg.eigh(mass.conj().T@mass)
    np.testing.assert_allclose(np.sqrt(values[-1]/values[-2]),q)
    np.testing.assert_allclose(abs(vectors[0,2])**2,p)
    jcp=(vectors[0,0]*vectors[1,1]*vectors[0,1].conjugate()*vectors[1,0].conjugate()).imag
    return {"q":q,"k":k,"r":r,"reactor_amplitude":4*p*(1-p),"solar_angle_square":1-2/(3*(1-p)),
            "jarlskog":float(jcp),"muon_third_probability":float(abs(vectors[1,2])**2),
            "b_source_success":float(erased@erased)}


def run():
    solar=window("KamLAND_dm21_eV2")
    atm=window("T2K_dm32_NO_eV2")
    minos=window("MINOS_effective_atm_eV2")
    measured=window("DayaBay_sin2_2theta13")
    q_range=(math.sqrt(1+atm[0]/solar[1]),math.sqrt(1+atm[1]/solar[0]))
    independent_q=(math.sqrt(minos[0]/solar[1]),math.sqrt(1+minos[1]/solar[0]))
    center=math.sqrt(1+SPEC["T2K_dm32_NO_eV2"][0]/SPEC["KamLAND_dm21_eV2"][0])
    first=branch_table(*q_range,measured)
    second=branch_table(*independent_q,measured)
    assert [k for k,v in first.items() if v["overlap_possible"]]==["nontrivial_k2_large"]
    two_solutions=sorted(k for k,v in second.items() if v["overlap_possible"])
    assert two_solutions==["nontrivial_k2_large","trivial_k1_large"]
    proof1=realisation(center,2,-.5,1,1)
    proof2=realisation(5.23,1,1.,1)
    assert independent_q[0]<=5.23<=independent_q[1]
    assert all(measured[0]<=x["reactor_amplitude"]<=measured[1] for x in (proof1,proof2))
    assert independent_q[0]<=center<=independent_q[1]
    assert abs(proof1["jarlskog"]+realisation(center,2,-.5,1,-1)["jarlskog"])<1.e-12
    x3=[(1-math.sqrt(1-t))/2 for t in measured]
    sol=window("JUNO_sin2_theta12_NO")
    min_p1=(1-x3[1])*(1-sol[1])
    assert min_p1>2/3
    angle=math.acos(math.sqrt(2/3))-math.acos(math.sqrt(min_p1))
    q0=SPEC["DayaBay_sin2_2theta13"][0]
    x0=(1-math.sqrt(1-q0))/2
    p10=(1-x0)*(1-SPEC["JUNO_sin2_theta12_NO"][0])
    return {"q_T2K_KamLAND":q_range,"q_MINOS_KamLAND":independent_q,"T2K_branches":first,
            "MINOS_branches":second,"T2K_central_solution":proof1,"MINOS_other_solution":proof2,
            "solar_TM1_given_reactor":[1-2/(3*(1-t)) for t in x3[::-1]],"JUNO_box":sol,
            "minimum_electron_or_null_rotation_deg":math.degrees(angle),
            "central_electron_or_null_rotation_deg":math.degrees(math.acos(math.sqrt(2/3))-math.acos(math.sqrt(p10)))}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=json.loads(json.dumps(run()))
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    elif args.check:
        def same(x,y):
            if isinstance(x,dict):
                return x.keys()==y.keys() and all(same(x[k],y[k]) for k in x)
            if isinstance(x,list):
                return len(x)==len(y) and all(same(a,b) for a,b in zip(x,y))
            return np.isclose(x,y,atol=2e-12,rtol=2e-10)
        assert same(result,json.loads(path.read_text(encoding="utf-8-sig")))
        print("1484 source branches and mass-only checks verified")
    else:
        print(json.dumps(result,indent=2))
