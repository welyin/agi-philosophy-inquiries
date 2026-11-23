"""Independent inversion of pole and gauge inputs for the exotic partner."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
def read(*p):
    return json.loads(HERE.joinpath(*p).read_bytes())

def zero(f, lo, hi):
    fl=f(lo)
    if fl*f(hi)>0:
        raise ValueError("root not bracketed")
    for _ in range(85):
        mid=(lo+hi)/2
        if f(mid)*fl>0:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def mixing():
    mean=np.array([.0962,.23186,-.26939])   # bR, lR, lL, LEP/SLD table 7.9
    sig=np.array([.0063,.00023,.00022])
    corr=np.array([[1,-.35,-.33],[-.35,1,.34],[-.33,.34,1]])
    C=sig[:,None]*sig*corr
    out=[]
    for sign in [-1,1]:
        def support(p):
            d=np.array([1,-1/3-p,p])
            return d@mean+sign*2*math.sqrt(d@C@d)
        root=zero(support,-.20,.25)
        grad=np.array([1,-1/3-root,root])
        point=mean+sign*2*C@grad/math.sqrt(grad@C@grad)
        assert abs((point[0]-point[1]/3)/(point[1]-point[2])-root)<1e-12
        out.append(root)
    return out


def one_loop(g3,singlet1,exotic):
    K=8*math.pi**2
    old=read("..","1608","inputs.json")["reference"]
    inv=np.array([3/(5*old["gY"]**2),1/old["g2"]**2,1/g3**2])
    # Independent direct piecewise running of the three inverse couplings.
    sm=np.array([41/10,-19/6,-7])
    sd=np.array([4/15,0,2/3])
    dd=np.array([10/3,2,4/3])
    def at(logscale):
        def length(m):
            return max(0, logscale-math.log(m/200))
        return inv-(sm*logscale + sd*(length(singlet1)+length(2100)) + dd*length(exotic))/K
    t=zero(lambda logscale:at(logscale)[0]-at(logscale)[2],0,50)
    assert t>max(math.log(m/200) for m in [singlet1,exotic,2100])
    end=at(t)
    W=(end[1]-end[0])*K*3/4
    A=t-W/2
    assert 0<A<t
    residual=end.copy()
    residual[1]-=(4/3)*2*(t-A)/K
    assert max(residual)-min(residual)<2e-14
    return {"meeting_GeV":200*math.exp(t),"minimum_fractional_count":W/t,
            "two_adjoints_GeV":200*math.exp(A),"common_g_squared":1/end[0]}


def certificate(gate):
    strong=read("..","1608","inputs.json")
    g0=strong["reference"]["g3"]
    step=2*strong["matching"]["g3_per_alpha"]*strong["marginal_width"]["alpha_s"]
    pts=[one_loop(g,s,d) for g,s,d in itertools.product([g0*(1-step),g0*(1+step)],
                                                       [1100,1600],gate)]
    limits={k:[min(row[k] for row in pts),max(row[k] for row in pts)] for k in pts[0]}
    lo,hi=limits["minimum_fractional_count"]
    assert 1<lo<hi<2
    proton=read("..","1572","inputs.json")
    W0,*errors=proton["lattice_W0_abs_GeV2"]
    Wmin=W0-2*sum(errors)
    mass=proton["mp_GeV"]
    phase=mass*(1-(proton["mpi_GeV"]/mass)**2)**2/(32*math.pi)
    Bmin=limits["common_g_squared"][0]*Wmin*proton["Vud_matching_floor"]
    Mmax=limits["meeting_GeV"][1]
    life=(proton["hbar_GeV_s"]/proton["seconds_per_year"])*Mmax**4/(phase*Bmin**2)
    limits.update({"same_gauge_baryon_vertex_unit_transfer_lifetime_ceiling_years":life,
                   "required_amplitude_transfer_ceiling_if_X_at_meeting": math.sqrt(life/proton["SK_2025_lifetime90_yr"])})
    return limits


def different_orderings_same_neutral_readout():
    opt=read("..","1608","inputs.json")
    g=opt["reference"]["g3"]
    dx=2*opt["matching"]["g3_per_alpha"]*opt["marginal_width"]["alpha_s"]
    def triplets(q):
        return [one_loop(s,d,q)["two_adjoints_GeV"] for d,s in
                itertools.product([1100,1600],[g*(1-dx),g*(1+dx)])]
    pn=min(triplets(1701.))
    pi=max(triplets(5460.))
    ma,mb=math.sqrt(.0029*.036)*1e-9,math.sqrt(.0029)*1e-9
    mc,md=math.sqrt(.0022)*1e-9,math.sqrt(.0022*(1+.024))*1e-9
    # For a real 2x2 seesaw choose B^2=m(Mbare+m), giving pole eigenvalues -m,Mbare+m.
    bare_no=pn/1.10
    target=ma*mb/((bare_no+ma)*(bare_no+mb))
    def neutral(bare,a,b):
        return a*b/((bare+a)*(bare+b))
    critical=zero(lambda f:neutral(pi*f,mc,md)/target-1,1/1.10,1.10)
    assert abs(neutral(critical*pi,mc,md)/target-1)<1.e-14
    p=mixing()[0]
    yuk_at_high=math.sqrt(2)/246.22*(5460-4.18**2/5460)*math.sqrt(p/(1-p))
    assert yuk_at_high<math.sqrt(4*math.pi)
    return {"exotic_masses_GeV":[1701.,5460.],"gauge_means_GeV":[pn,pi],
            "neutral_to_gauge_mass_scales":[1/1.10,critical],
            "Y_row_norm_floor_at_high_mass_lowest_p":yuk_at_high,
            "common_q_Z":target}


def run():
    pl,ph=mixing()
    vev=246.22
    mb=5.0  # conservative upper envelope on the light mass, in GeV
    def max_mass(cap):
        a=vev/math.sqrt(2)*cap*math.sqrt((1-pl)/pl)
        return (a+math.hypot(a,2*mb))/2
    unit=max_mass(1)
    wide=max_mass(math.sqrt(4*math.pi))
    required=math.sqrt(2)/vev*(1700-mb**2/1700)*math.sqrt(pl/(1-pl))
    assert unit<1700<wide
    return {"mixing_ellipsoid": [pl,ph], "Y_mass_max_for_y_one_GeV":unit,
            "Y_mass_max_for_y_sqrt4pi_GeV":wide, "Yukawa_floor_at_1700":required,
            "neutral_collision":different_orderings_same_neutral_readout(), "inclusive_gate":certificate([200,wide]),
            "gate_with_Wb_limit":certificate([1700,wide])}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=run()
    f=HERE/"independent_results.json"
    if args.check:
        if result != json.loads(f.read_bytes()):
            raise SystemExit("1634 recalculation differs")
        print("1634 independent PASS: chiral ellipse, beta-function roots, adjoint count and proton-decay transfer")
    else:
        with f.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")
