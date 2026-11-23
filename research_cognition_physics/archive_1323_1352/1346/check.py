"""Count open weak rows and electromagnetic charges in the two registered windows."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
BOX = {
    "low": (F(208,100) / (3 * F(120,100)), F(228,100) / 3),
    "high": (F(348,100) / (3 * F(112,100)), F(364,100) / 3),
}
W_BRANCH = F(6741,10000), F(27,10000)
C_TAG = F(489,1000), F(20,1000)

def response(q, up, down):
    return up*q*q + down*(1-q)*(1-q)

def extrema(a, b, up, down):
    values = [response(q, up, down) for q in (a,b)]
    if up+down:
        vertex = F(down, up+down)
        if a <= vertex <= b:
            values.append(response(vertex, up, down))
    return min(values),max(values)

def hits(range1, range2):
    return range1[0] <= range2[1] and range2[0] <= range1[1]

def in_box(q, added, down):
    return all(low <= response(q, up, down) <= high
               for up,(low,high) in ((2-added,BOX["low"]),(2,BOX["high"])))

def integers_and_charges():
    # From 2*q^2<=high_max<2*(.78)^2 hence |q|<.78.
    # Down*(1-q)^2<=low_max<16*(.22)^2 hence down<=15.
    assert 2*F(78,100)**2 > BOX["high"][1]
    assert 16*F(22,100)**2 > BOX["low"][1]
    candidates = []
    for added in (1,2):
        for down in range(16):
            kept = []
            witness = None
            for i in range(-780,780):
                a,b = F(i,1000),F(i+1,1000)
                if (hits(extrema(a,b,2-added,down), BOX["low"])
                        and hits(extrema(a,b,2,down), BOX["high"])):
                    kept.append(i)
                    if witness is None:
                        for point in (a,(a+b)/2,b):
                            if in_box(point, added, down):
                                witness = point
                                break
            if kept:
                assert witness is not None
                segments = []
                for i in kept:
                    if segments and segments[-1][1] == i:
                        segments[-1][1] = i+1
                    else:
                        segments.append([i,i+1])
                candidates.append({
                    "new_up":added, "low_up":2-added, "down":down,
                    "q_enclosing_segments": [[a/1000,b/1000] for a,b in segments],
                    "q_exact_witness": str(witness),
                    "q_witness":float(witness),
                })
    return candidates

def main():
    p,e=W_BRANCH
    lo,hi = p-e,p+e
    # Three known normalized leptonic channels, color dimension three.
    row_range = (lo/(1-lo)/F(11,10), hi/(1-hi))
    assert 1 < row_range[0] < 2 < row_range[1] < 3
    cal=integers_and_charges()
    tag0,err=C_TAG
    # Per-row differences in corrected hadronic widths may be 10 percent.
    def tag_range(a):
        return (F(9*a,9*a+11*(2-a)), F(11*a,11*a+9*(2-a)))
    for r in cal:
        tr=tag_range(r["new_up"])
        r["weak_tag_interval"] = [float(t) for t in tr]
        r["tag_overlap"] = hits(tr,(tag0-err,tag0+err))
    eps=F(2,25)
    with_triplet=[]
    for row in cal:
        permitted=[]
        for af,bf in row["q_enclosing_segments"]:
            a,b=3*F(str(af)),3*F(str(bf))
            permitted += [i for i in range(-3,4) if hits((a,b),(i-eps,i+eps))]
        if permitted:
            with_triplet.append([row["new_up"],row["down"],sorted(set(permitted))])
    assert with_triplet==[[1,2,[2]]]
    # Counterexample on two disconnected weak components, q1 != q2.
    x,y=F(1,3),F(2,3)
    other_low=x*x+(x-1)**2+(y-1)**2
    assert (other_low,other_low+y*y,3*x)==(F(2,3),F(10,9),1)
    assert response(F(2,3),1,2)==F(2,3)
    assert response(F(2,3),2,2)==F(10,9)
    return {"round":1346,"status":"PASS","ratios_per_color":{k:[str(i) for i in v] for k,v in BOX.items()},
            "complete_weak_rows_outer":[float(t) for t in row_range],
            "integral_row_count":2,"q_grid_width":.001,"down_enumerated":[0,15],
            "historical_calibration":cal,"tag_data":[float(tag0-err),float(tag0+err)],
            "net_triplet_extra":{"charge_error_at_most":str(eps),"survivors":[{"new_up":a,"down":d,"three_q_nearest":i} for a,d,i in with_triplet]},
            "two_disconnected_components":{"up_charges":[str(x),str(y)],"low":str(other_low),"high":str(other_low+y*y),"tag":"1/2","triplet_charge":str(3*x)},
            "exact_ideal":{"low": "2/3","high":"10/9","q":"2/3","low_up":1,"new_up":1,"down":2,"tag":"1/2"}}

if __name__=="__main__":
    flags=argparse.ArgumentParser()
    commands=flags.add_mutually_exclusive_group(required=True)
    commands.add_argument("--check",action="store_true")
    commands.add_argument("--write-new",action="store_true")
    args=flags.parse_args()
    output=main()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8",newline="\n") as f:
            f.write(json.dumps(output,ensure_ascii=False,indent=2)+"\n")
    else:
        assert json.loads(target.read_text(encoding="utf-8"))==output
    print(json.dumps(output,ensure_ascii=False))

