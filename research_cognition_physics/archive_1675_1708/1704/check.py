"""Neutral flavor contractions as a unitary return with positive or signed weights."""
from argparse import ArgumentParser
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def gain(q,x):
    return q*(1+x)/(x+q*q)


def segment_intersections(vertices,theta,tol=2e-13):
    """Positive ray intersections with each chord, including intervals along the ray."""
    intervals=[]
    for i,v in enumerate(vertices):
        for w in vertices[i+1:]:
            a,b=v*np.exp(-1j*theta),w*np.exp(-1j*theta)
            if max(abs(a.imag),abs(b.imag))<tol:
                lo,hi=sorted((a.real,b.real))
                if hi>=0:
                    intervals.append((max(lo,0.),hi))
            elif abs(a.imag-b.imag)>tol:
                t=a.imag/(a.imag-b.imag)
                z=((1-t)*a+t*b).real
                if -tol<=t<=1+tol and z>=-tol:
                    intervals.append((max(0.,z),max(0.,z)))
    return intervals


def solve():
    cfg=read("inputs.json")
    previous=read(cfg["quark_lepton_matches"])
    phi=math.pi*cfg["flavon_phase_over_pi"]
    omega=np.exp(4j*phi)
    r=abs((3+2*omega)/5)
    vertices=np.array([1.,omega,omega*omega])
    signed=np.r_[vertices,-vertices]
    report={}
    for key,k in zip(cfg["charged_targets"],cfg["charged_phases_in_units_of_phi"]):
        s=previous["down_one_width_lower"][key]["large_modulus"]
        theta=k*phi
        cos_face=math.cos(theta+math.pi/3)
        qpos=1/(2*cos_face)
        smax=gain(min(qpos,math.sqrt(r)),r)
        roots=np.sort(np.roots([s,-1-r,s*r]))
        np.testing.assert_allclose(roots,previous["one_pole_with_unit_neutral"]["witnesses"][key]["charged_bare_moduli"])
        assert s>smax and all(0<v<1 for v in roots)
        cases=[]
        for q in roots:
            z=q*np.exp(1j*theta)
            face=z*np.exp(1j*math.pi/3)
            nmin=2*face.real-1
            diff=2*face.imag/math.sqrt(3)
            p0,p1=(1-nmin+diff)/2,(1-nmin-diff)/2
            weights=np.array([p0,p1,nmin])
            assert min(weights)>0
            represented=np.dot(weights,vertices*np.array([1,1,-1]))
            np.testing.assert_allclose(represented,z)
            np.testing.assert_allclose(gain(abs(represented),r)*np.exp(1j*np.angle(represented)),
                                       s*np.exp(1j*theta))
            chords=segment_intersections(signed,theta)
            delta=min(max(0.,a-q,q-b) for a,b in chords)
            # Which segments of the regular hexagon contain this point?
            if k==1:
                assert delta>0.02
                min_arms=3
                min_two_sign_budget=None
            else:
                assert delta<1e-12
                min_arms=2
                min_two_sign_budget=(1+q)/2
                assert min_two_sign_budget > nmin
                np.testing.assert_allclose(min_two_sign_budget*(-vertices[2])+
                                           (1-min_two_sign_budget)*vertices[2],z)
            cases.append(dict(q=q,minimum_total_negative_weight=nmin,
                              saturating_weights=weights.tolist(),least_number_of_distinct_arms=min_arms,
                              gap_to_two_arm_chord=delta,negative_weight_if_only_two_arms=min_two_sign_budget))
        report[key]=dict(charged_target=s,charged_phase_deg=math.degrees(theta),
                         maximum_same_port_positive_raw=qpos,
                         maximum_same_port_positive_gain=smax,shortfall=s-smax,solutions=cases)
    return dict(fixed_common_pole_ratio=r,matched_points=report)


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,(np.floating,float)) else x


if __name__=="__main__":
    cli=ArgumentParser()
    cli.add_argument("--check",action="store_true")
    result=stable(solve())
    dest=HERE/"results.json"
    if cli.parse_args().check:
        assert result==read("results.json")
        print("PASS: two actual texture phases, exact positive bounds, signed minima and chord alternatives")
    else:
        with dest.open("x",encoding="utf8") as file:
            file.write(json.dumps(result,indent=2)+"\n")
