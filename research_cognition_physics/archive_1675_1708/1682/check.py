"""Mass spectrum excludes a single passive shared pole for any charged return phase."""
import argparse
import json
import math
import runpy
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
OLD=runpy.run_path(str(HERE/CFG["charged_intervals"]))
I=OLD["I"]


def read(path):
    return json.loads((HERE/path).read_bytes())


def run():
    q=read(CFG["quark_masses_angles"])["high_scale"]
    le=read(CFG["charged_leptons"])["scales"]["1e16"][3:]
    qintervals=[OLD["iv"].width(row) for row in q["down"]]
    qangles=[OLD["trig"](OLD["iv"].width(row),True)
             for row in q["angles_theta12_theta13_theta23_delta"][:3]]
    full=[qintervals[0]/qintervals[2],qintervals[1]/qintervals[2],*qangles]
    D,P=OLD["spectral_target"](le,0)
    smallmax=I.point(0)
    max_turn=0.
    small_r_in=[1.,0.]
    N=CFG["cos_delta_cells"]
    for j in range(N):
        x=OLD["bounds"](Fraction(2*j,N)-1,Fraction(2*j+2,N)-1)
        co,u0,_,_,_,aux=OLD["coefficients"](full+[x],(1,1,3))
        h,L,base,wf,er,ei,K=aux
        # The endpoint bracketing in the historical level-set solver uses a certified monotone cubic.
        small,large=OLD["level_set"](D,co)
        def scale(w): return (w*K/((1-w)*L)).sqrt().values()
        lower,upper=scale(small)
        small_r_in=[min(small_r_in[0],lower),max(small_r_in[1],upper)]
        rr=D*small.sqrt()/(h*L.sqrt())
        bound=u0*(1-small).square()+base*rr+((wf*rr).sqrt()
                  +(er.square()+ei.square()).sqrt()*small/L).square()
        smallmax=I(min(smallmax.lo,bound.lo),max(smallmax.hi,bound.hi))
        max_turn=max(max_turn,(K/(2*L)).sqrt().values()[1])
        assert small.hi<OLD["iv"].SCALE//3<large.lo
    assert smallmax.hi<P.lo

    # T_d=1+d²+s_q². The entire family of quark phases obeys C² <= H22 <= s23²+s_q².
    d,s,p,u,t=full
    diag22=t.square()+s.square()
    trace_d=1+d.square()+s.square()
    g=I.point(CFG["test_modulus"])
    trace_cap=g.square()*trace_d+(9-g.square())*diag22
    trial=(d*s)*g/(trace_cap*trace_cap.sqrt())
    assert max_turn<g.values()[0]
    assert trial.lo>D.hi
    rn=I.point(7).sqrt()/5
    polemax=(1+rn)/(2*rn.sqrt())
    qc=1/I.point(3).sqrt()
    fixed=qc*(1+rn)/(rn+qc.square())
    gain=g/polemax
    gain_fixed=g/fixed
    xmax=(g-(g.square()-1).sqrt()).square()
    # The neutral transfer is decreasing in x. Endpoint uses outward arithmetic.
    ket=rn*(1+xmax)/(rn.square()+xmax)
    independent_witness=read(CFG["neutrino_witness"])["joint_witness"]["common_abs_t"]
    charged_witness=read(CFG["charged_partner"])["opposite_k8_charged_witness"]["abs_t"]
    assert polemax.hi<g.lo and gain.lo>I.point(1).hi
    # Directly invert randomly generated explicit heavy blocks and check the scalar formula.
    tests=0
    for xx in (.21,float(rn.values()[0]),.9,2.3):
        for qabs in (.11,float(qc.values()[0]),1.):
            for angle in np.linspace(-np.pi,np.pi,9):
                z=qabs*np.exp(1j*angle)
                block=np.array([[np.sqrt(xx),z],[-z.conjugate(),np.sqrt(xx)]])
                measured=-np.linalg.inv(block)[0,1]*(1+xx)
                np.testing.assert_allclose(measured,(1+xx)*z/(xx+qabs*qabs))
                tests+=1
    return {"phase_cells":N, "all_inputs_halfwidths":CFG["widths"],
            "small_s_outer":small_r_in,"small_P_upper":smallmax.values()[1],"data_P_lower":P.values()[0],
            "turn_modulus_upper":max_turn,"trial_modulus":g.values(),
            "trial_D_lower":trial.values()[0],"data_D_upper":D.values()[1],
            "required_large_s_strict_lower":float(CFG["test_modulus"]),
            "known_charged_feasible_s":charged_witness,"known_joint_feasible_s":independent_witness,
            "single_pole_maximum":polemax.values(),"neutral_two_path_charged_fixed_q_response":fixed.values(),
            "minimum_relative_gain":gain.values(),"minimum_relative_gain_for_fixed_q":gain_fixed.values(),
            "necessary_x_upper_if_unchanged_relative_gains":xmax.values(),
            "necessary_neutral_modulus_lower_if_unchanged_relative_gains":ket.values(),
            "direct_inverse_tests":tests}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    target=HERE/"results.json"
    if args.write:target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    else:assert json.loads(target.read_bytes())==result
    print(json.dumps(result,indent=2))
