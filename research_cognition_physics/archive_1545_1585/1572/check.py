"""Adjoint triangle gauge matching and the unseen charged-pion nucleon decay."""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
TAU=2*math.pi


def inputs(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def decreasing_zero(fun, lower=0., upper=45.):
    assert fun(lower)>0>fun(upper)
    for _ in range(110):
        middle=(lower+upper)/2
        if fun(middle)>0: lower=middle
        else: upper=middle
    return (lower+upper)/2


def inclusive_bound(cfg, gauge, delta):
    n=cfg["width_count"]
    abox=[]
    for label in ("g1_GUT","g2","g3"):
        mean,err=gauge[label]
        abox.append([4*math.pi/(mean+n*err)**2,4*math.pi/(mean-n*err)**2])
    a1,a2,a3=abox
    b1,b2,b3=cfg["SM_b"]
    ds=[[TAU*(a1[0]-a[1]-2*delta),TAU*(a1[1]-a[0]+2*delta)] for a in (a2,a3)]
    x=lambda L,d:math.exp(.75*(d-(b1-b2)*L))
    y=lambda L,d:math.exp(.5*(d-(b1-b3)*L))
    xmax=lambda L:x(L,ds[0][1])
    ymin=lambda L:y(L,ds[1][0])
    hi=decreasing_zero(lambda L:xmax(L)+y(L,ds[1][1])-2)
    # Ratio monotone on this interval; thus ymin-xmax is decreasing.
    assert max(xmax(t)/ymin(t) for t in (0,hi))<1
    lo=decreasing_zero(lambda L:ymin(L)-xmax(L)-2,0,hi)
    mz=gauge["mu_ref_GeV"]
    gx_floor_sq=4*math.pi/(a1[1]+delta-b1*lo/TAU)
    MXmax=mz*math.exp(hi)
    W=cfg["lattice_W0_abs_GeV2"]
    Wfloor=W[0]-2*sum(W[1:])
    mass=cfg["mp_GeV"]
    kin=mass/(32*math.pi)*(1-(cfg["mpi_GeV"]/mass)**2)**2
    unit=cfg["hbar_GeV_s"]/cfg["seconds_per_year"]
    life=unit*MXmax**4/(kin*(gx_floor_sq*Wfloor*cfg["Vud_matching_floor"])**2)
    SK=cfg["SK_2025_lifetime90_yr"]
    ref=cfg["operator_transfer_test_floor"]
    return {"L_envelope":[lo,hi],"MX_envelope_GeV":[mz*math.exp(lo),MXmax],
            "M3_over_M32_envelope":[x(hi,ds[0][0]),x(lo,ds[0][1])],
            "g5sq_floor":gx_floor_sq,"W0_two_width_floor_GeV2":Wfloor,
            "partial_lifetime_max_yr_unit_transfer":life,
            "partial_lifetime_max_yr_named_transfer":life/ref**2,
            "full_amplitude_transfer_necessary_max":math.sqrt(life/SK),
            "MX_min_if_named_amplitude_floor_GeV":MXmax*(SK*ref**2/life)**.25}


def center(cfg,gauge):
    aa=[4*math.pi/gauge[k][0]**2 for k in ("g1_GUT","g2","g3")]
    b=cfg["SM_b"]
    x=lambda l:math.exp(.75*(TAU*(aa[0]-aa[1])-(b[0]-b[1])*l))
    y=lambda l:math.exp(.5*(TAU*(aa[0]-aa[2])-(b[0]-b[2])*l))
    edges=[decreasing_zero(lambda l:y(l)-x(l)-2),decreasing_zero(lambda l:y(l)+x(l)-2)]
    L=edges[1]
    M=gauge["mu_ref_GeV"]*math.exp(L)
    m32=.25*M
    s3,s8=m32*x(L),m32*y(L)
    # Construct the mass term and the sole invariant with a 24_H insertion.
    vev_coefficient=(s8-s3)/5
    common_mass=s3+3*vev_coefficient
    assert math.isclose(common_mass-.5*vev_coefficient,m32,rel_tol=2e-12)
    masses={"weak3":s3,"color8":s8,"dirac32":m32}
    aG=[aa[i]-b[i]*L/TAU-sum(cfg["adjoint_beta"][k][i]*math.log(M/m)/TAU for k,m in masses.items()) for i in range(3)]
    assert max(aG)-min(aG)<2e-12
    assert all(gauge["mu_ref_GeV"]<m<M for m in masses.values())
    return {"central_two_ends_L":edges,"both_real_positive_mass_witness_GeV":masses,
            "MX_GeV":M,"physical_mass_singlet_GeV":abs(common_mass-vev_coefficient),
            "matched_inverse_strengths":aG,"M3_over_M32":x(L),"M8_over_M32":y(L)}


def extra_insertion(cfg, gauge):
    """An independently adjustable Tr(rho Sigma rho Sigma) term."""
    scale = 1.e15
    log_scale = math.log(scale / gauge["mu_ref_GeV"])
    a = [4*math.pi / gauge[label][0]**2 for label in ("g1_GUT", "g2", "g3")]
    b = cfg["SM_b"]
    ratios = [math.exp(p*(TAU*(a[0]-a[i])-(b[0]-b[i])*log_scale))
              for i,p in ((1,.75),(2,.5))]
    q = scale/4
    s3, s8 = [q*ratio for ratio in ratios]
    extra = (s3+s8-2*q)/25
    k = (s8-s3+5*extra)/5
    mass = s3+3*k-9*extra
    assert math.isclose(mass-.5*k-6*extra, q, rel_tol=1e-13)
    spectra = {"weak3":s3,"color8":s8,"dirac32":q}
    assert all(gauge["mu_ref_GeV"] < m < scale for m in spectra.values())
    inverse = [a[i] - (b[i]*log_scale + sum(cfg["adjoint_beta"][key][i] *
               math.log(scale/value) for key,value in spectra.items()))/TAU
               for i in range(3)]
    assert max(inverse)-min(inverse)<1e-12 and sum(ratios)<2
    proton = cfg["mp_GeV"]
    phase = proton*(1-(cfg["mpi_GeV"]/proton)**2)**2/(32*math.pi)
    coefficient = (4*math.pi/inverse[0])*.974*cfg["lattice_W0_abs_GeV2"][0]/scale**2
    life = cfg["hbar_GeV_s"]/cfg["seconds_per_year"]/(phase*coefficient**2)
    return {"MX_GeV":scale, "mass_witness_GeV":spectra, "inverse_strengths":inverse,
            "extra_coefficient_over_M32":extra/q,"x_plus_y":sum(ratios),
            "central_partial_lifetime_yr_transfer_one":life}


def compute():
    c=inputs("inputs.json")
    g=inputs(c["gauge_input"])
    b=c["adjoint_beta"]
    assert max(abs(sum(v[i] for v in b.values())-10/3) for i in range(3))<1e-14
    out={str(d):inclusive_bound(c,g,d) for d in c["matching_cases_abs_each_inverse_coupling"]}
    th={}
    for rate in (c["operator_transfer_test_floor"],1.):
        low,high=0.,0.85
        for _ in range(75):
            mid=(low+high)/2
            if inclusive_bound(c,g,mid)["full_amplitude_transfer_necessary_max"]<rate: low=mid
            else:high=mid
        th[str(rate)]=high
    return {"center":center(c,g),"extra_insertion":extra_insertion(c,g),"outer_envelopes":out,
            "necessary_matching_Linf_at_given_amplitude_floor":th}


def rounded(o):
    if isinstance(o,dict):return {k:rounded(v) for k,v in o.items()}
    if isinstance(o,list):return [rounded(v) for v in o]
    if isinstance(o,float):return float(format(o,".12g"))
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser()
    gr=p.add_mutually_exclusive_group(required=True)
    gr.add_argument("--check",action="store_true")
    gr.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    result=rounded(compute())
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    elif result!=inputs("results.json"):
        raise SystemExit("Saved joint certificate differs")
    print(json.dumps(result))
