"""Exact character-lattice certificate, and a three-flavour phase-lock check."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def det(a,b):
    return a[0]*b[1]-a[1]*b[0]


def rank(rows):
    rows=[[Fraction(x) for x in row] for row in rows if any(row)]
    if not rows:
        return 0
    at=0
    for j in range(len(rows[0])):
        p=next((i for i in range(at,len(rows)) if rows[i][j]),None)
        if p is None:
            continue
        rows[at],rows[p]=rows[p],rows[at]
        for i in range(at+1,len(rows)):
            scale=rows[i][j]/rows[at][j]
            rows[i]=[x-scale*y for x,y in zip(rows[i],rows[at])]
        at+=1
    return at


def minimum_harmonic(rows,anomaly):
    """Order in Z² / row lattice, None denotes infinite order."""
    if rank(rows+[anomaly])!=rank(rows):
        return None
    if rank(rows)==1:
        v=next(r for r in rows if any(r))
        divisor=math.gcd(*v)
        v=[x//divisor for x in v]
        j=next(i for i in range(2) if v[i])
        stride=math.gcd(*(r[j]//v[j] for r in rows))
        return stride//math.gcd(stride,anomaly[j]//v[j])
    index=math.gcd(*(det(a,b) for a,b in combinations(rows,2)))
    extended=math.gcd(index,*(det(r,anomaly) for r in rows))
    return index//extended


def compute():
    cfg=read(HERE/"inputs.json")
    anomaly=cfg["topological_step"]
    case={}
    for name,rows in cfg["charge_menus"].items():
        order=minimum_harmonic(rows,anomaly)
        gain=1+rank(rows)-rank(rows+[anomaly])
        assert bool(gain)==(order is not None)
        if order is not None:
            terms=cfg["examples_known_coefficients"][name]
            assert [sum(n*q[j] for n,q in zip(terms,rows)) for j in range(2)]==[order*x for x in anomaly]
        case[name]={"minimum_integer_k":order,"extra_continuous_phase":gain}
    osc=read(HERE/cfg["oscillation"])
    def window(row):
        return row[0]-2*row[1],row[0]+2*row[-1]
    amplitude=window(osc["DayaBay_sin2_2theta13"])
    p=[(1-math.sqrt(1-x))/2 for x in amplitude]
    solar=window(osc["JUNO_sin2_theta12_NO"])
    electron_min=[(1-p[1])*(1-solar[1]),(1-p[1])*solar[0],p[0]]
    assert min(electron_min)>0
    # Historical edge samples, errors added in magnitude, then transfer allowed.
    edges=read(HERE/cfg["quark_edges"])["edges"]
    network={}
    for tag,central,errors in edges:
        if (central-2*sum(errors))/1.02>0:
            a,b=tag
            network.setdefault(a,set()).add(b)
            network.setdefault(b,set()).add(a)
    visited=set()
    todo=["u"]
    while todo:
        a=todo.pop()
        if a not in visited:
            visited.add(a)
            todo.extend(network[a]-visited)
    assert visited==set("udcsbt")
    assert rank([[1,1,0],[3,1,2]])==2
    # Same two histories, different record overlap. CP-even two-phase contrast.
    ratio=.2
    overlap_extrema={}
    for overlap in [0.,.4,1.]:
        vis=2*ratio*overlap/(1+ratio*ratio)
        overlap_extrema[str(overlap)]=[1-vis,1+vis]
    g=read(HERE/cfg["weak_gauge"])["g2"][0]
    return {"angular_character_count":case,"electron_row_two_width_NO_floors":electron_min,
            "connected_quark_ports":sorted(visited),"split_flavour_extra_phase":1+1-2,
            "two_history_rate_relative_incoherent_average":overlap_extrema,
            "log10_single_semiclassical_amplitude_at_MZ_illustration":-8*math.pi**2/(g*g*math.log(10))}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true"); mode.add_argument("--check",action="store_true")
    value=compute()
    target=HERE/"results.json"
    if parser.parse_args().write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(value,f,indent=2,allow_nan=False); f.write("\n")
    else:
        stored=read(target)
        for field in ["electron_row_two_width_NO_floors","log10_single_semiclassical_amplitude_at_MZ_illustration"]:
            a,b=value[field],stored[field]
            if not isinstance(a,list):
                a,b=[a],[b]
            assert all(math.isclose(x,y,rel_tol=1e-12) for x,y in zip(a,b))
        assert all(value[k]==stored[k] for k in value if k not in ("electron_row_two_width_NO_floors","log10_single_semiclassical_amplitude_at_MZ_illustration"))
    print(json.dumps(value,indent=2))
