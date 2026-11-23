"""Gauge matching at a fixed physical scale; gY is the ordinary SM hypercharge."""
import argparse
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
K = 8 * math.pi**2
BY, BW = 41/6, -19/6


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def anomaly(n, a):
    """All left. x=3X, twice the fundamental quadratic index is one."""
    b = n-a
    # (signed colour dimension, signed L dimension, x, number of multiplets)
    states = [(1,3,0,n), (3,3,2,a), (3,-3,-1,b), (-3,1,-2,n),
              (-3,1,1,n), (-3,1,-5,a), (-3,1,4,b)]
    coefficients = [0]*6
    for c,l,x,count in states:
        cd,ld = abs(c),abs(l)
        terms=[(1 if c>0 else -1)*ld if cd==3 else 0,
               (1 if l>0 else -1)*cd if ld==3 else 0,
               ld*x if cd==3 else 0, cd*x if ld==3 else 0,
               cd*ld*x, cd*ld*x**3]
        coefficients=[old+count*term for old,term in zip(coefficients,terms)]
    assert coefficients==[0,2*(3*a-n),0,3*(3*a-n),0,162*(n-3*a)]
    return coefficients


def matched_inverse(g1, g2, beta_squared):
    return 1/(.6*g1**2) - beta_squared/g2**2


def edge(pair, width):
    return [pair[0]-width*pair[1],pair[0]+width*pair[1]]


def calc():
    cfg=read(HERE/"inputs.json")
    low=read(HERE/cfg["weak_matching"])
    target=sum(read(HERE/cfg["target_from_Higgs"])["central_lower_scale_crossing_GeV"])/2
    mz=low["mu_ref_GeV"]
    time=math.log(target/mz)
    mass=cfg["threshold_witness_GeV"]
    width=cfg["error_widths"]
    def corners(g1,g2):
        return product(edge(g1,width),edge(g2,width))
    def branch(beta_sq):
        slope=BY-beta_sq*BW
        f=matched_inverse(low["g1_GUT"][0],low["g2"][0],beta_sq)
        fbox=sorted(matched_inverse(g1,g2,beta_sq) for g1,g2 in corners(low["g1_GUT"],low["g2"]))
        values=[x-slope*time/K for x in (fbox[0],fbox[-1])]
        out={"F_at_target_central":f-slope*time/K,"F_at_target_box":values}
        if beta_sq==3:
            out["zero_scale_GeV_central"]=mz*math.exp(K*f/slope)
            out["zero_scale_GeV_box"]=[mz*math.exp(K*x/slope) for x in (fbox[0],fbox[-1])]
        return out
    close=branch(3)
    weak=branch(1/3)
    f=close["F_at_target_central"]
    lo,hi=close["F_at_target_box"]
    assert hi<0<weak["F_at_target_box"][0]
    def second_read(number,mu):
        extra_b2=number*8/3
        g1,g2=low["g1_GUT"][0],low["g2"][0]
        invY=1/(.6*g1**2)-BY*math.log(mu/mz)/K
        invW=1/g2**2-BW*math.log(mu/mz)/K-extra_b2*max(0,math.log(mu/mass))/K
        return invW/(invW+invY)
    witness={}
    for number in (1,2):
        wanted=8*number*math.log(target/mass)
        inv=f+wanted/K
        witness[str(number)]={"inverse_gX2":inv,"alpha_X":1/(4*math.pi*inv),
                "alpha_X_below_one_all_box":lo+wanted/K>1/(4*math.pi),
                "sine2_at_second_scale":second_read(number,cfg["second_scale_GeV"]),
                "max_mass_all_box_for_positive_GeV":target*math.exp(K*lo/(8*number)),
                "max_mass_all_box_for_alpha_below_one_GeV":target*math.exp(K*(lo-1/(4*math.pi))/(8*number))}
    coeff_real=F(1,6)*2*3
    coeff_dirac=F(4,3)*2*3
    assert (coeff_real,coeff_dirac)==(1,8)
    no_threshold_cross={}
    for key,v in cfg["table2_SM"].items():
        no_threshold_cross[key]=[min(z:= [matched_inverse(a,b,3) for a,b in corners(v["g1"],v["g2"])]),max(z)]
    zero_fields=anomaly(3,1)
    families=[n for n in range(1,6) if any(not any(anomaly(n,a)) for a in range(n+1))]
    assert families==[3] and not any(zero_fields)
    # Positive F not sufficient for small corrections; alpha_X < 1 is only a working cut.
    return {"scale_GeV":target,"beta_minus_sqrt3":close,"beta_minus_inverse_sqrt3":weak,
            "integrated_3dB2_minus_dBY_needed_at_zero":-K*f,
            "optimistic_vs_conservative_integrated_box":[-K*hi,-K*lo],
            "integral_needed_for_alpha_X_below_one":K*(1/(4*math.pi)-f),
            "finite_match_alone_needed_box":[-hi,-lo],
            "integral_weak_charge_for_1TeV_to_target":-K*f/math.log(target/mass),
            "real_triplets_for_positive_entire_box":math.floor(-K*lo/math.log(target/mass))+1,
            "real_triplets_for_alpha_below_one_entire_box":math.floor(K*(1/(4*math.pi)-lo)/math.log(target/mass))+1,
            "two_loop_table_F_box":no_threshold_cross,"dirac_triplet_witness":witness,
            "SM_sine2_at_second_scale":second_read(0,cfg["second_scale_GeV"]),
            "allowed_families_with_full_color_AF":families,"anomalies_3_1":zero_fields}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    opt=parser.add_mutually_exclusive_group(required=True)
    opt.add_argument("--write-new",action="store_true")
    opt.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=calc()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(results,out,indent=2)
            out.write("\n")
    else:
        assert results==read(path)
        print("1497: matter anomalies, transport, thresholds and second gauge reading verified.")
