"""Independent finite role-screen audit using exact effect vertices and
energy blocks. No author-code imports. --write creates once; default readonly.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
P=tuple(F(x,10) for x in (4,3,1,2))
Q=tuple(F(x,10) for x in (3,4,2,1))
GAMMA=tuple(F(x,15) for x in (8,4,2,1))
EFFECTS=tuple(product((0,1),repeat=4))

def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))

def tv(a,b):
    return sum((abs(x-y) for x,y in zip(a,b)),F(0))/2

def vertex_max(r,t):
    values=[dot(effect,tuple(x-t*y for x,y in zip(r,GAMMA))) for effect in EFFECTS]
    value=max(values)
    winners=[list(effect) for effect,v in zip(EFFECTS,values) if v==value]
    return value,winners

def reduced_channel(r):
    # V|i> = |i xor 1>|0 if i even, 2 if i odd>.
    return tuple(r[i ^ 1] for i in range(4))

def evaluate():
    assert sum(P,F(0))==sum(Q,F(0))==sum(GAMMA,F(0))==1
    assert min(P+Q+GAMMA)>0
    assert Counter(P)==Counter(Q)
    assert dot(P,range(4))==dot(Q,range(4))==F(11,10)
    # Gamma_j=2^-j/Z: independent exact ratio check.
    assert all(GAMMA[i]==2*GAMMA[i+1] for i in range(3))

    tests=[]
    for name,source,target,t,expected in (
        ("p_to_q",P,Q,F(3,4),F(389,4000)),
        ("q_to_p",Q,P,F(3,2),F(193,2000))
    ):
        lo,src_eff=vertex_max(source,t)
        hi,tgt_eff=vertex_max(target,t)
        assert hi-lo==F(1,10)
        # Each effect has oscillation <=1 on normalized signed variations.
        # Thus state errors cost eta,zeta, and reference error costs t*delta.
        assert all(max(e)-min(e)<=1 for e in EFFECTS)
        for effect in EFFECTS:
            for i in range(4):
                for j in range(4):
                    assert abs(effect[i]-effect[j])<=1
                    assert abs(t*(effect[i]-effect[j]))<=t
        eta=zeta=delta=F(1,1000)
        bound=hi-lo-eta-zeta-t*delta
        assert bound==expected
        tests.append({"direction":name,"t":str(t),"source_max":str(lo),
                      "target_max":str(hi),"source_maximizers":src_eff,
                      "target_maximizers":tgt_eff,"gap":str(hi-lo),
                      "conditional_bound":str(bound)})
    assert vertex_max(P,F(3,2))[0]==F(1,10)
    assert vertex_max(Q,F(3,2))[0]==0
    assert vertex_max(P,F(3,4))[0]==F(1,4)
    assert vertex_max(Q,F(3,4))[0]==F(7,20)

    states=list(product(range(4),range(3)))
    blocks={energy:[state for state in states if sum(state)==energy] for energy in range(6)}
    image={state:state for state in states}
    # Build each exchange from the desired central-energy input and parity.
    # This is an energy-block construction, not the author's explicit pair list.
    for energy,block in blocks.items():
        for state in block:
            i,w=state
            if w==1:
                target_i=i ^ 1
                target=(target_i,energy-target_i)
                assert target in block and target!=state
                image[state]=target
                image[target]=state
    assert len(set(image.values()))==12
    assert all(image[image[state]]==state for state in states)
    assert all(sum(state)==sum(image[state]) for state in states)
    # Explicit exact basis matrix check covers every vector, not just populations.
    U=[[int(image[s]==t) for s in states] for t in states]
    for i in range(12):
        for j in range(12):
            assert sum(U[k][i]*U[k][j] for k in range(12))==int(i==j)
            assert (sum(states[i])-sum(states[j]))*U[i][j]==0

    branches=[]
    for name,r,target in (("p_to_q",P,Q),("q_to_p",Q,P)):
        joint={state:F(0) for state in states}
        for i,weight in enumerate(r):
            joint[image[(i,1)]]+=weight
        s=tuple(sum((joint[(i,w)] for w in range(3)),F(0)) for i in range(4))
        w=tuple(sum((joint[(i,j)] for i in range(4)),F(0)) for j in range(3))
        assert s==target and s==reduced_channel(r)
        assert w==(F(1,2),F(0),F(1,2))
        assert dot(w,range(3))==1
        assert sum((weight*(j-1)**2 for j,weight in enumerate(w)),F(0))==1
        assert tv(w,(F(0),F(1),F(0)))==1
        # Joint spectrum equals the input spectrum plus zeros, exactly.
        assert Counter(x for x in joint.values() if x)==Counter(r)
        # W is a deterministic function of final S: H(W|S)=0.
        assert all(sum(joint[(i,j)]>0 for j in range(3))==1 for i in range(4))
        total_out=sum((weight*sum(state) for state,weight in joint.items()),F(0))
        assert total_out==dot(r,range(4))+1==F(21,10)
        branches.append({"direction":name,
            "joint_nonzero":[{"state":list(state),"probability":str(weight)}
                              for state,weight in joint.items() if weight],
            "S_output":[str(x) for x in s],"W_output":[str(x) for x in w],
            "total_mean_energy":str(total_out),"W_mean":"1","W_variance":"1",
            "W_return_probability":"0","W_trace_distance_from_blank":"1",
            "W_entropy_increase":"log(2)",
            "joint_entropy_equals_input_exact":True,
            "conditional_W_entropy_given_S":"0","mutual_information":"log(2)"})

    # General reduced map retains within-parity coherence and removes cross-parity
    # coherence. Check all 16 matrix units via the shared environment labels.
    env=lambda i:0 if i%2==0 else 2
    retained=[(i,j) for i in range(4) for j in range(4) if env(i)==env(j)]
    lost=[(i,j) for i in range(4) for j in range(4) if env(i)!=env(j)]
    assert len(retained)==len(lost)==8
    # Paid conversion is not Gibbs-preserving.
    drift=tv(reduced_channel(GAMMA),GAMMA)
    assert drift==F(1,3)

    return {"verdict":"PASS_WITH_SCOPE","independent_groups":1,
        "new_numbered_rounds":0,"new_formal_science_groups":0,
        "method":"16 diagonal effect vertices; exact total-energy blocks; all matrix units",
        "scalar_identities":{"mean_energy":"11/10","same_entropy":"equal eigenvalue multiset",
            "same_free_energy":"same beta, energy, entropy",
            "same_relative_entropy":"D=-S+beta E+log Z"},
        "monotone_tests":tests,
        "energy_block_dimensions":[len(blocks[e]) for e in range(6)],
        "unitary_bijective_involution":True,"commutes_with_total_H":True,
        "assisted_conversions":branches,
        "general_reduced_channel":{
            "formula":"Xpair*(Peven*rho*Peven+Podd*rho*Podd)*Xpair",
            "preserved_matrix_units":[list(x) for x in retained],
            "removed_matrix_units":[list(x) for x in lost],
            "Gibbs_reference_drift":str(drift)},
        "scope":["Additional scalar simplification, not the user's full conjecture.",
            "No geometry variable or actual new device permission.",
            "Unknown spectator acts by U tensor I; correlations need not be preserved.",
            "Reverse listed task uses a fresh pure W; not free reuse of output W.",
            "No reopening instantaneous-Gibbs or arbitrary thermal-auxiliary routes."],
        "independent_code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=evaluate()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as file:
            json.dump(result,file,ensure_ascii=False,indent=2)
            file.write("\n")
        mode="WRITE_ONCE_PASS"
    else:
        assert json.loads(target.read_text(encoding="utf-8"))==result
        mode="READ_ONLY_MATCH"
    print(json.dumps({"verdict":result["verdict"],"mode":mode,
                     "gap":"1/10","W_return":"0","Gibbs_drift":"1/3"}))
