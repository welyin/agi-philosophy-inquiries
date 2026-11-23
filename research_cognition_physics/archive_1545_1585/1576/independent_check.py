"""SU(15) irreducible-index certificate and the common-threshold consequence."""
from collections import deque
from fractions import Fraction as F
from pathlib import Path
import argparse
import itertools
import json
import math

ROOT=Path(__file__).resolve().parent
N=15

def representation(weight):
    rows=[sum(weight[i:]) for i in range(N)]
    dim=math.prod(F(rows[i]-rows[j]+j-i,j-i) for i in range(N) for j in range(i+1,N))
    boxes=sum(rows)
    c2=(N*boxes+sum(v*(v+1-2*i) for i,v in enumerate(rows,1))-F(boxes**2,N))/2
    return int(dim), dim*c2/(N*N-1)

def small_reps():
    start=(0,)*(N-1)
    todo=deque([start])
    seen={start}
    keep=[]
    frontier=[]
    while todo:
        w=todo.popleft()
        dim, index=representation(w)
        if dim>=455 and index>=39:
            frontier.append([dim,index])
            continue
        keep.append({"weight":[[k+1,a] for k,a in enumerate(w) if a],"dimension":dim,"index":str(index)})
        for i in range(N-1):
            nxt=tuple(x+(k==i) for k,x in enumerate(w))
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    assert len(keep)==8
    return keep, {"minimal_cut_dimension":min(v[0] for v in frontier),
                  "minimal_cut_index":str(min(v[1] for v in frontier)),"visited":len(seen)}

def bounds(g2,g3,mz):
    a2,a3=4*math.pi/g2**2,4*math.pi/g3**2
    b2,b3=-19/6,-7
    cross_log=2*math.pi*(a2-a3)/(b2-b3)
    mu=mz*math.exp(cross_log)
    au=a2-b2*cross_log/(2*math.pi)
    def max_span(b,inv):
        return math.exp(2*math.pi*inv/b)
    log10t=math.log(1e4/mz)
    return {"M23_GeV":mu,"A23_before_threshold":au,
            "A_parent_if_no_light_extras":au/3,
            "M_effective_positive_GeV":mu/max_span(42,au),
            "M_effective_alpha_parent_le_1_GeV":mu/max_span(42,au-3),
            "Q_over_10TeV_positive_weak":max_span(b2+42,a2-b2*log10t/(2*math.pi)),
            "Q_over_10TeV_positive_strong":max_span(b3+42,a3-b3*log10t/(2*math.pi))}

def solve():
    irr,front=small_reps()
    # Solve signed multiplicities in the net basis (5, 10, 15).
    found=[]
    for a,b in itertools.product(range(-3,4),repeat=2):
        if (2*a+b,a+2*b)==(3,0):
            found.append((a,b))
    assert found==[(2,-1)] and 2*(N-4)-(N+4)-3==0
    dims=(2*105+120+3*15)
    index=2*F(N-2,2)+F(N+2,2)+3*F(1,2)
    assert (dims,index)==(375,23)
    delta=F(2,3)*(3*index-3*(F(3,2)+F(1,2)))
    assert delta==42
    data=json.loads((ROOT.parent/"1369"/"inputs.json").read_text(encoding="utf-8-sig"))
    center=bounds(data["g2"][0],data["g3"][0],data["mu_ref_GeV"])
    corners=[bounds(data["g2"][0]+j*2*data["g2"][1],data["g3"][0]+k*2*data["g3"][1],data["mu_ref_GeV"])
             for j,k in itertools.product([-1,1],repeat=2)]
    return {"small_irreps":irr,"dominance_frontier":front,"minimal_SU15_components":dims,
            "minimal_SU15_index":int(index),"surplus_index_each_SU5_normalized_factor":63,
            "minimal_extra_beta_each":int(delta),"one_loop_without_finite_matching_central":center,
            "two_input_width_corner_ranges":{key:[min(a[key] for a in corners),max(a[key] for a in corners)] for key in center}}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    options=ap.parse_args()
    result=solve()
    output=ROOT/"independent_results.json"
    if options.check:
        assert result==json.loads(output.read_text(encoding="utf-8-sig"))
        print("1576 independent index and matter thresholds verified")
    else:
        output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
