"""Gauge thresholds of the actual added colour partners on two registered spectra."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
TAU=2*math.pi
def read(path):
    return json.loads((HERE/path).read_bytes())


def calculate():
    low=read("../1369/inputs.json")
    baseline=read("../1698/results.json")["shared_witness"]
    rich=read("../1720/results.json")["central_conditional_combined_point"]
    single=read("../1706/results.json")["cluster_reference_scales_GeV"][1]
    g=runpy.run_path(str(HERE/"../1696/check.py"))
    MZ,soft=low["mu_ref_GeV"],baseline["soft_GeV"]
    b,_,offset=g["branch"](True,soft,MZ)
    spectrum=[(mass,2/3,g["COLOR"]) for mass in baseline["colored_fermions_GeV"]]
    spectrum.extend((mass,1/3,g["COLOR"]) for mass in baseline["colored_scalars_GeV"])
    spectrum.extend((mass,2/3,g["WEAK"]) for mass in baseline["weak_fermions_GeV"])
    spectrum.extend((mass,1/3,g["WEAK"]) for mass in baseline["weak_scalars_GeV"])
    extra=[(mass,2/3,g["COLOR"]) for mass in rich["five_color_poles_GeV"]]
    extra.extend((mass,1.,g["WEAK"]) for mass in rich["four_weak_poles_GeV"])
    inputs=np.array([low[key] for key in ("g1_GUT","g2","g3")])
    def run(values,signed,rules,number,soft_mates):
        alpha=4*math.pi/values**2
        length=TAU*(g["W"]@(alpha-offset))/(g["W"]@b)
        candidate=MZ*math.exp(length)
        new=[(single,2/3,g["COLOR"])]*number
        if soft_mates:
            new += [(math.hypot(single,soft),1/3,g["COLOR"])]*number
        # At any candidate meeting, the colour inverse coupling is fixed
        # irrespective of how many colourless weak doublets restore the difference.
        def inv_at(masses):
            return alpha-offset-b*length/TAU-sum(
                weight*charge*math.log(candidate/mass)/TAU for mass,weight,charge in masses)
        before=inv_at(rules)
        after=inv_at(rules+new)
        required=TAU*(after[1]-after[2])
        # A positive weak load T shifts A2 by -T/(2*pi) and A1 by -(3/5)*T/(2*pi).
        completed=after-g["WEAK"]*required/TAU
        np.testing.assert_allclose(completed,[after[2]]*3,rtol=0.,atol=2e-11)
        # Cross-check the strong coupling by piecewise integration; locate loss
        # of alpha_3<1 and zero of the formal one-loop inverse before the candidate.
        events=[(soft,g["MSSM"][2]-g["SM"][2],-g["ADJ"][2]/(12*math.pi))]
        events.extend((mass,weight*charge[2],0.) for mass,weight,charge in rules+new)
        events.append((candidate,0.,0.))
        a,prev,slope=alpha[2],MZ,g["SM"][2]
        crossings={}
        for mass,delta,matching in sorted(events):
            final=a-slope*math.log(mass/prev)/TAU
            for mark in (1.,0.):
                if final<=mark<a:
                    crossings[str(int(mark))]=prev*math.exp((a-mark)*TAU/slope)
                elif final+matching<=mark<final:
                    crossings[str(int(mark))]=mass
            a,prev,slope=final+matching,mass,slope+delta
        assert abs(a-after[2])<1e-10
        # Maximum extra weak loading from lowering the two old weak fermions
        # to MZ with their fixed 30 TeV soft prescription.
        old_relief=sum(2/3*math.log(m/MZ)+1/3*math.log(math.hypot(m,soft)/math.hypot(MZ,soft))
                       for m in baseline["weak_fermions_GeV"])
        other_relief=sum(math.log(m/1.e7) for m in rich["four_weak_poles_GeV"]) if len(rules)>len(spectrum) else 0.
        deficit=required-old_relief-other_relief
        dirac=candidate*math.exp(-required/(2/3*number))
        full=candidate*math.exp(-required/number)
        full_three=candidate*math.exp(-required/3)
        if MZ<full<candidate:
            complete_weak=[(full,1.,g["WEAK"])]*number
            np.testing.assert_allclose(inv_at(rules+new+complete_weak),completed,rtol=0.,atol=3e-11)
        if number==2 and soft_mates:
            assert full_three>1.e7
            np.testing.assert_allclose(inv_at(rules+new+[(full_three,1.,g["WEAK"])]*3),
                                       completed,rtol=0.,atol=3e-11)
        return dict(signs=signed,U=candidate,old_inverse_colour=before[2],new_inverse_colour=after[2],
                    required_weak_log_load=required,inverse_coupling_spreads_if_no_new_weak=after-after[2],
                    thresholds_of_formal_alpha3=crossings,minimum_net_new_weak_log_load=deficit,
                    uncovered_with_two_above_10PeV=deficit-2*math.log(candidate/1.e7),
                    three_equal_complete_weak_poles=full_three,
                    single_complete_pair_max_load=length,single_Dirac_pair_max_load=2*length/3,
                    equal_masses_Dirac_only=dirac if MZ<dirac<candidate else None,
                    equal_masses_complete_weak_pairs=full if MZ<full<candidate else None)
    def tidy(o):
        return json.loads(json.dumps(o,default=lambda a:a.tolist(),allow_nan=False),
                          parse_float=lambda s:float(f"{float(s):.11g}"))
    out={}
    for family,rules,original in (("R1698",spectrum,baseline["inverse_alpha_U"]),
                                   ("R1720_R800",spectrum+extra,rich["common_inverse_coupling"])):
        rows={}
        for n,soft_mates in itertools.product((1,2),(False,True)):
            c=run(inputs[:,0],"centre",rules,n,soft_mates)
            assert abs(c["old_inverse_colour"]-original)<1e-7
            corners=[run(inputs[:,0]+2*np.array(s)*inputs[:,1],list(s),rules,n,soft_mates)
                     for s in itertools.product((-1,1),repeat=3)]
            def extrema(key):
                vals=[v[key] for v in corners]
                return [min(vals),max(vals)]
            rows[f"{n}_colour_"+("with_soft_scalars" if soft_mates else "Dirac_only")]=dict(
                central=c, two_width_U=extrema("U"),
                two_width_equal_complete_weak_masses=extrema("equal_masses_complete_weak_pairs"),
                two_width_three_equal_complete_weak_masses=extrema("three_equal_complete_weak_poles"),
                two_width_uncovered_after_two_above_10PeV=extrema("uncovered_with_two_above_10PeV"),
                two_width_new_inverse_colour=extrema("new_inverse_colour"),
                two_width_weak_load=extrema("required_weak_log_load"),
                min_extra_beyond_one_complete_pair=min(r["minimum_net_new_weak_log_load"]-
                                                       r["single_complete_pair_max_load"] for r in corners),
                min_extra_beyond_one_Dirac_pair=min(r["minimum_net_new_weak_log_load"]-
                                                   r["single_Dirac_pair_max_load"] for r in corners))
        out[family]=rows
    fit=read("inputs.json")["nearby_measured_inputs_in_marginal_widths"]["absolute_and_oscillations"]
    shift=np.array([0.,0.,fit["g3"]])
    weak_inputs=inputs[:,0]+shift*inputs[:,1]
    out["R1698_two_new_colours_and_soft_scalars_at_g3_minus_0p67"]=run(
        weak_inputs,shift.tolist(),spectrum,2,True)
    return tidy(out)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    data=calculate()
    path=HERE/"gauge_results.json"
    if args.check:
        assert data==json.loads(path.read_bytes())
        print("1739 two gauge spectra, common-threshold projection and colour crossings verified")
    else:
        path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf8")
