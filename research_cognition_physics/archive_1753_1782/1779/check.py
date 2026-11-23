"""Exact anomaly and minimum-cost tests for net exotic weak doublets."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import json
import math
from pathlib import Path


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
HALF=F(1,2)


def label(name):
    kind,index=name.split(":")
    index=int(index)
    if kind=="q":
        return 3,F(index)+F(1,6)
    assert kind=="l" and index>=0
    return 1,-F(index)-HALF


def row(d,Y):
    return dict(d=d,Y=Y,coeff=2*d*Y,largest=max(abs(Y-HALF),abs(Y+HALF)),
                bY_doublet=F(4,3)*d*Y**2,bY_all=d*(F(8,3)*Y**2+F(1,3)))


def anomaly_and_singlets(spectrum):
    anomalies=defaultdict(F)
    singlets=defaultdict(lambda:[0,0])
    witten=0
    for d,Y,n in spectrum:
        # n>0 means ordinary LH doublets plus RH singlets, n<0 the mirror.
        qs=[Y-HALF,Y+HALF]
        anomalies["SU3_cubic"]+=(2-len(qs))*n*(d==3)
        anomalies["SU3sq_Y"]+=(2*Y-sum(qs))*n*(d==3)*HALF
        anomalies["SU2sq_Y"]+=d*Y*n*HALF
        anomalies["grav_Y"]+=d*(2*Y-sum(qs))*n
        anomalies["Y_cubic"]+=d*(2*Y**3-sum(q**3 for q in qs))*n
        witten+=d*abs(n)
        for q in qs:
            singlets[d,q][0 if n>0 else 1]+=abs(n)
    assert anomalies["Y_cubic"]==-3*anomalies["SU2sq_Y"]
    return anomalies,witten%2,singlets


def witness(spec):
    ex=[(*label(name),n) for name,n in spec.items()]
    ordinary=[(*label("q:0"),CFG["ordinary_net"]["q0"]),(*label("l:0"),CFG["ordinary_net"]["l0"])]
    anomalies,parity,slots=anomaly_and_singlets(ordinary+ex)
    assert not any(anomalies.values()) and parity==0
    D=sum(d*abs(n) for d,y,n in ex)
    by=sum(abs(n)*row(d,y)["bY_all"] for d,y,n in ex)
    bdouble=sum(abs(n)*row(d,y)["bY_doublet"] for d,y,n in ex)
    all_dbl=sum(row(d,y)["bY_doublet"]*abs(n) for d,y,n in ordinary+ex)
    singpart=sum(F(2,3)*d*q*q*(r+l) for (d,q),(r,l) in slots.items())
    assert all_dbl+singpart==by+sum(abs(n)*row(d,y)["bY_all"] for d,y,n in ordinary)
    return dict(signed_exotics=spec,D=D,max_abs_Q=max(row(d,y)["largest"] for d,y,n in ex),
                exotic_b2=F(D,3),exotic_bY_doublets=bdouble,exotic_bY_full=by,
                inverse_alphaY_drop_above_known_content_one_decade=float(by)*math.log(CFG["common_running_scale_ratio"])/(2*math.pi))


def rows_with_cap(max_q):
    entries=[]
    for k in range(-int(max_q)-2,int(max_q)+3):
        if k:
            v=row(*label(f"q:{k}"))
            if v["largest"]<=max_q:
                entries.append(v)
    for j in range(1,int(max_q)+1):
        v=row(*label(f"l:{j}"))
        if v["largest"]<=max_q:
            entries.append(v)
    return entries


def min_D(rows,cap):
    # Exact discrete reachability up to a known upper witness.
    reachable=[set() for _ in range(cap+1)]
    reachable[0].add(0)
    for budget in range(1,cap+1):
        for v in rows:
            if budget>=v["d"]:
                for previous in reachable[budget-v["d"]]:
                    for sign in [-1,1]:
                        reachable[budget].add(previous+sign*v["coeff"])
        if -1 in reachable[budget]:
            return budget
    return None


def min_by(bound):
    # Costs multiplied by nine are integers. A single row beyond the entire
    # upper bound is unaffordable; no implicit finite hypercharge cutoff.
    entries=[]
    for kind in ("q","l"):
        for sign in ((-1,1) if kind=="q" else (1,)):
            j=1
            while True:
                v=row(*label(f"{kind}:{sign*j}"))
                if v["bY_all"]>bound:
                    break
                entries.append(v)
                j+=1
    budget=9*bound
    assert budget.denominator==1
    costs=[int(v["bY_all"]*9) for v in entries]
    assert all(F(c,9)==v["bY_all"] for c,v in zip(costs,entries))
    paths=[set() for _ in range(int(budget)+1)]
    paths[0].add(0)
    for z in range(1,len(paths)):
        for c,v in zip(costs,entries):
            if c<=z:
                for x in paths[z-c]:
                    paths[z].update([x-v["coeff"],x+v["coeff"]])
        if -1 in paths[z]:
            return F(z,9)
    raise AssertionError("upper witness did not close")


def run():
    witnesses={k:witness(v) for k,v in CFG["explicit_compensators"].items()}
    cost=min_by(witnesses["mixed"]["exotic_bY_full"])
    assert cost==F(173,9)
    frontier={k:min_D(rows_with_cap(F(k)),15) for k in CFG["max_charge_caps"]}
    assert list(frontier.values())==[15,5,3]
    direct_extra=CFG["direct_product_example"]
    alternate=[(*label("q:0"),4),(*label("l:0"),3),(3,F(direct_extra["extra_colour_doublet_Y"]),direct_extra["extra_net"])]
    an,mod,_=anomaly_and_singlets(alternate)
    assert not any(an.values()) and mod==0
    # The direct-product counterexample violates the G6 coloured doublet congruence.
    q6=6*alternate[-1][1]
    assert q6.denominator==1 and q6%6!=1
    checked=0
    for a in range(-2,3):
        for b in range(-2,3):
            for c in range(-2,3):
                quarks=[(*label("q:-1"),a),(*label("q:0"),3),(*label("q:1"),b)]
                exotic=(*label("l:2"),c)
                l0=3-5*a+7*b-5*c
                residual,parity,_=anomaly_and_singlets(quarks+[exotic,(*label("l:0"),l0)])
                assert not any(residual.values()) and not parity
                checked+=1
    universal=row(*label("q:-1"))
    per_unit=universal["bY_all"]/abs(universal["coeff"])
    assert per_unit==F(59,45)
    return dict(witnesses=witnesses,minimum_extra_weighted_doublets_by_charge=frontier,
                minimum_extra_complete_bY=cost,optional_SM3_plus_one_quark_plus_exotic_bY=F(41,6)+F(11,9)+cost,
                lattice_cases_verified=checked,lower_bY_per_unit_general_mismatch=per_unit,
                direct_product=dict(extra_Y=alternate[-1][1],abs_charge_max=F(2,3),
                                    anomaly_free=True,G6_descends=False))


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,F):
        return str(x)
    return round(x,10) if isinstance(x,float) else x


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    verify=cli.parse_args().check
    result=tidy(run())
    path=HERE/"results.json"
    if verify:
        assert json.loads(path.read_bytes())==result
        print("1779 verified: anomaly lattice, three charge fronts, full singlet cost and quotient control")
    else:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
