"""Six published GLS moments, a profiled response bound and a coherent-receiver control."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np


def frac(x):
    return Fraction(str(x))


# Q^2, GLS, statistical error, systematic error. Kim et al., table III.
TABLE = [
    ("1.26","2.39",".08",".14"), ("2.00","2.49",".08",".10"),
    ("3.16","2.55",".06",".10"), ("5.01","2.78",".06",".19"),
    ("7.94","2.82",".07",".19"), ("12.59","2.80",".13",".18"),
]
ROWS = [tuple(map(frac,row)) for row in TABLE]
HMAX, D = frac(".30"), frac(".02")


def endpoints(row):
    _,s,stat,syst=row
    return s-stat-syst, s+stat+syst


def thresholds(indices):
    cuts=[endpoints(ROWS[i])[1]+D+HMAX/ROWS[i][0] for i in indices]
    top=min(cuts)
    return {"numerator":float(top),"N4_floor":float(top/4),"N6_floor":float(top/6),
            "limiting_Q_squared":float(ROWS[indices[cuts.index(top)]][0])}


def feasibility(n,c,indices):
    # For given n, existence of a single H requires intersection of all H intervals.
    lo=frac(0);hi=HMAX
    for i in indices:
        q=ROWS[i][0]
        a,b=endpoints(ROWS[i])
        lo=max(lo,q*(n*c-b-D))
        hi=min(hi,q*(n+D-a))
    return lo<=hi


def counts():
    idx=list(range(6)); others=[i for i in idx if i!=2]
    full=thresholds(idx)
    partial=thresholds(others)
    min_count=max(endpoints(row)[0]-D for row in ROWS)
    table={}
    for c in (".47",".48",".70",".708",".72",".80"):
        table[c]=[n for n in range(1,15) if feasibility(n,frac(c),idx)]
    assert table[".47"]==[3,4,5,6] and table[".48"]==[3,4,5]
    assert table[".70"]==[3,4] and table[".708"]==[3] and table[".80"]==[3]
    assert feasibility(4,frac(".708"),others)
    bound=endpoints(ROWS[2])[1]+D+HMAX/ROWS[2][0]
    assert feasibility(4,bound/4,idx) and not feasibility(4,bound/4+frac(".00000001"),idx)
    assert feasibility(6,bound/6,idx) and not feasibility(6,bound/6+frac(".00000001"),idx)
    return dict(six=full,omit_Q2_3p16=partial,lower_N=float(min_count),compatible_counts=table)


def response_competition():
    moments=[]
    for row in ROWS:
        q=row[0]
        c3=frac(".98")-frac(".2268")/q
        h=frac(".10")
        s=3*c3-h/q
        a,b=endpoints(row)
        assert a<=s<=b and 6*(c3/2)-h/q==s and 4*(3*c3/4)-h/q==s
        moments.append(float(s))
    # Table I at the bin containing Q^2=3.16. The low-x entry includes data
    # and extrapolation; it is not an estimate of the unmeasured tail alone.
    portions=[frac(x) for x in (".558","1.889",".092",".009")]
    return dict(shared_six_moments=moments,table_I_sum=float(sum(portions)),
                low_x_portion_fraction=float(portions[0]/sum(portions)),
                neutron_charge_interval_in_e_1e_minus21=[-1.5,.7],
                count_3_neutral_ud=[1,2],count_6_neutral_ud=[2,4],
                count_4_neutral_ud_with_changed_charges=[1,3])


def memory_control():
    # q is a two-level multiplicity of the same global neutral, fixed-count object.
    # Joint of the would-be receiver q, untouched partner P, environment E.
    ket=np.zeros((2,2,2))
    ket[0,0,0]=ket[1,1,0]=1/np.sqrt(2)
    good=ket.reshape(4,2)@ket.reshape(4,2).T
    leaked=ket.copy()
    leaked[1,:,:]=leaked[1,:,::-1].copy()
    columns=leaked.reshape(4,2)
    bad=columns@columns.T
    bell=ket[:,:,0].flatten()
    assert np.allclose(np.diag(good),np.diag(bad)) and np.allclose([np.trace(good),np.trace(bad)],1)
    assert np.allclose([bell@good@bell,bell@bad@bell],[1,.5])
    return dict(charge_preserving_partner_Bell_pass=[1,.5],charge_preserving_basis_success=[1,1])


def calculate():
    return dict(round=1361,profile=counts(),witnesses=response_competition(),receiver_control=memory_control())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true");mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    path=Path(__file__).with_name("results.json")
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:assert result==json.loads(path.read_bytes())
    print(json.dumps(result))

