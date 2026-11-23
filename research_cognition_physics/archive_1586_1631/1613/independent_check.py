"""Round 1613 independent check. Free-neutron data are separate experiments."""
import argparse
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UD = math.cos(math.pi / 14)
# Hill & Plestid Eq. (28), including the illustrative short-distance input
# following that equation.  The axial measurement is PERKEO III.
BASE = (5263.284, .04537, .02704, 1.27641)
ERROR = (.017, .00027, .00007, .00056)
OBS = {
    "material_storage_UCNtau_published": (877.83, .22+.17, .22+.20),
    "NIST_proton_count_2013": (887.7, 1.2+1.9, 1.2+1.9),
    "JPARC_electron_count": (877.2, 1.7+3.6, 1.7+4.0),
}


def lifetime(normalization, inner, outer, axial):
    return normalization / ((1+inner)*(1+outer)*(1+3*axial**2)*UD**2)


def band(width):
    values=[lifetime(*[v+w*width*e for v,e,w in zip(BASE,ERROR,sgn)])
            for sgn in itertools.product((-1,1),repeat=4)]
    return min(values),max(values)


def study():
    checks={}
    for widths in (1,2):
        t0,t1=band(widths)
        observed={}
        for name,(v,d,e) in OBS.items():
            lo,hi=v-d*widths,v+e*widths
            observed[name]=dict(seconds=[lo,hi],gap_seconds=max(0,lo-t1,t0-hi),
                                allowed_rate_scale=[t0/hi,t1/lo],
                                neutron_inferred_ud=[UD*math.sqrt(t0/hi),UD*math.sqrt(t1/lo)])
        checks[str(widths)]=dict(predicted_seconds=[t0,t1],holdouts=observed)
    assert checks["2"]["holdouts"]["NIST_proton_count_2013"]["gap_seconds"]>3.62
    assert checks["2"]["holdouts"]["material_storage_UCNtau_published"]["gap_seconds"]==0
    assert checks["2"]["holdouts"]["JPARC_electron_count"]["gap_seconds"]==0
    old=json.loads((ROOT.parent/"1381"/"inputs.json").read_bytes())
    new=json.loads((ROOT/"inputs.json").read_bytes())
    theta=math.pi/14
    comparisons={}
    for n in [1,2]:
        f=old["K_form"]
        f_low,f_high=f[0]-n*f[1],f[0]+n*f[1]
        ub_upper=(old["F_times_cb"][0]+n*sum(old["F_times_cb"][1:]))/old["F_lower"]*(old["ub_over_cb"][0]+n*sum(old["ub_over_cb"][1:]))
        group_ks=[f_low*math.sqrt(math.sin(theta)**2-ub_upper**2),f_high*math.sin(theta)]
        ks_ratio=new["kse3"]["branching_in_1e4"]
        inferred_err=new["kse3"]["product"]*math.hypot(*ks_ratio[1:])/(2*ks_ratio[0])
        nu=[(tag,data[0]+n*data[1]) for tag,data in [("2020",old["ud_beta"]),("2025",new["new_nuclear"])]]
        beta={tag:dict(max_ud=h, gap_degrees=math.degrees(math.acos(h)-theta),
                       maximum_allowed_width_response=(h/UD)**2) for tag,h in nu}
        axial=new["leptonic"]
        ratio_low=(axial["combination"]-n*(axial["branching_error"]+axial["radiative_error"]))/(axial["decay_constant_ratio"]+n*axial["upper"])
        comparisons[str(n)]=dict(beta=beta,group_semileptonic_product=group_ks,
                                 KSe3_branching_alone_implied_error=inferred_err,
                                 axial_ratio_low=ratio_low, minimum_axial_over_common_scale=ratio_low/math.tan(theta))
    return dict(exact_ud=UD, central_lifetime=lifetime(*BASE), neutron=checks, decay_checks=comparisons)


def round_numbers(x):
    if isinstance(x,dict):
        return {k:round_numbers(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [round_numbers(v) for v in x]
    return round(x,12)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    flag=parser.add_mutually_exclusive_group(required=True)
    flag.add_argument("--write-new",action="store_true")
    flag.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=round_numbers(study())
    path=ROOT/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        assert result==json.loads(path.read_bytes())
    print("1613 independent neutron lifetime, nuclear transitions and K ratios verified")
