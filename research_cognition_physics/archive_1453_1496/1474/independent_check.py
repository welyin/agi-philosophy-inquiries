"""Bound scalar pair production for a septet using positive partonic DY rates."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def scaled(k,interval):
    a,b=interval
    return min(k*a,k*b),max(k*a,k*b)


def square(interval):
    a,b=interval
    return (0. if a<=0<=b else min(a*a,b*b),max(a*a,b*b))


def amplitude_bounds(flavor,weight,w,x):
    # Reduced amplitudes multiplied by the common Z denominator and w(1-w).
    charge,iso = (2/3,1/2) if flavor=="up" else (-1/3,-1/2)
    b = scaled(2*charge, (w[0]*(1-w[0]),w[1]*(1-w[1])))
    bx = min(a*c for a in b for c in x),max(a*c for a in b for c in x)
    imag2 = (0., (.028*max(map(abs,bx)))**2)
    low = high = 0.
    for t in (iso,0.):
        lin=scaled(charge*(2-weight)-2*t,w)
        r=square((lin[0]+t*weight-bx[1],lin[1]+t*weight-bx[0]))
        low += r[0]+imag2[0]
        high += r[1]+imag2[1]
    return low,high


def run():
    # ATLAS table 8, 200 GeV, same-sign dimuon signal region.
    template=11.3
    low_template=11.3-1.3
    limit=3.78
    bounds={}
    for flavor in ("up","down"):
        low=float("inf")
        high=0.
        for i in range(80):
            w=(.22+.02*i/80,.22+.02*(i+1)/80)
            for j in range(80):
                xmax=(91.20/(2*200))**2
                x=(xmax*j/80,xmax*(j+1)/80)
                n=amplitude_bounds(flavor,0,w,x)
                d=amplitude_bounds(flavor,1,w,x)
                assert d[0]>0
                low=min(low,n[0]/d[1])
                high=max(high,n[1]/d[0])
        bounds[flavor]=(low-1.e-10,high+1.e-10)
    u,d=bounds["up"][0],bounds["down"][0]
    def up_fraction(n,br,transfer):
        return (limit/(n*br**2*transfer)-d)/(u-d)
    M,mh,other=300.,125.,1500.
    v,mw=246.22,80.379
    ratio=(mw/M)**2
    c=M**2*(1-4*ratio+12*ratio**2)*math.sqrt(1-4*ratio)/(32*math.pi*v*v)
    scalar_points={}
    for y in (.31,.32,.35):
        x=math.sqrt(1-y*y)
        h=np.array([x,y/4,math.sqrt(15)*y/4])
        l=np.array([-1.5*y/x,1.,0.])
        assert abs(h@l+5*y/4)<1.e-13 and abs(h@h-1)<1.e-14
        K=other**2*np.eye(3)+(mh**2-other**2)*np.outer(h,h)
        radial=K-M**2*np.outer(l,l)
        test=M**2*((h@l)**2/mh**2+(l@l-(h@l)**2)/other**2)
        minimum=np.linalg.eigvalsh(radial).min()
        assert (minimum>0)==(test<1)
        z=15*y*y
        scalar_points[str(y)]={"rank_one_positivity_test":test,"radial_min_GeV2":minimum,
                "h_curvature_must_exceed_GeV":5*y*M/4,
                "partial_width_fraction":z*c,"min_zBR_in_narrow_window":c*z*z/.05,
                "single_production_response_needed":2*.425**2*.05/(c*z*z)}
    assert scalar_points["0.32"]["rank_one_positivity_test"]<1
    assert scalar_points["0.35"]["h_curvature_must_exceed_GeV"]>mh
    # High-energy elastic W+W- check with only the specified scalar exchanges.
    elastic = {}
    for dmass,y,mneutral in ((300.,.32,1500.),(200.,.44,475.)):
        z=15*y*y
        h_l_squared=(5*y/4)**2
        ell_perp_squared=1+2.25*y*y/(1-y*y)-h_l_squared
        curvature=dmass**2*(h_l_squared/mh**2+ell_perp_squared/mneutral**2)
        a0=(2*(mh**2+z*mneutral**2)+z*dmass**2)/(16*math.pi*v*v)
        x=(mw/dmass)**2
        c0=dmass**2*(1-4*x+12*x*x)*math.sqrt(1-4*x)/(32*math.pi*v*v)
        min_rate=z*z*c0/.05
        # Schur complement. This lower bound also allows split heavy neutral masses.
        den=1-dmass**2*h_l_squared/mh**2
        assert den>0
        min_neutral_sum=mh**2+dmass**2*(15*y/4)**2/den
        elastic[str(int(dmass))]=dict(
            y=y,other_neutral_GeV=mneutral,positive_potential_test=curvature,
            abs_a0=a0,min_abs_a0_arbitrary_neutral_masses=(
                2*min_neutral_sum+z*dmass**2)/(16*math.pi*v*v),
            fraction_WW_width=z*c0,minimum_z_times_B=min_rate,
            minimum_response_transfer=2*.425**2/min_rate)
        assert curvature<1
    assert elastic["200"]["abs_a0"]<.5<elastic["300"]["min_abs_a0_arbitrary_neutral_masses"]
    return {"single_elastic_channel_check":elastic,"joint_curvature_and_template_points":scalar_points,
            "y_interval_for_WW_test_at_unit_transfer_and_unbounded_other_masses":[
            (2*.425**2*.05/(225*c))**.25,4*mh/(5*M)],
            "m_GeV":200,"unpolarized_pointwise_ratio_bounds":bounds,
            "template_N_central_lower_and_limit":[template,low_template,limit],
            "guaranteed_no_extra_transfer_or_branch_loss_up_fraction_required_central":up_fraction(template,1,1),
            "one_width_up_fraction_required":up_fraction(low_template,1,1),
            "one_width_up_fraction_required_with_BR_095":up_fraction(low_template,.95,1),
            "one_width_all_up_minimum_BR":(limit/(low_template*u))**.5,
            "one_width_all_down_guaranteed_events":d*low_template}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=json.loads(json.dumps(run()),parse_float=lambda x:float(format(float(x),".11g")))
    p=HERE/"independent_results.json"
    if args.write_new:
        with p.open("x",encoding="utf-8") as fp:
            json.dump(data,fp,indent=2)
            fp.write("\n")
    else:
        assert data==json.loads(p.read_text(encoding="utf-8"))
        print("1474: second production channel and table 8 same-sign dimuon test verified")
