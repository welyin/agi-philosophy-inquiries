"""Independent exact audit of lottery gains and finite observation scope."""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import json,sys

def v(*x):return tuple(map(F,x))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(t,a):return tuple(F(t)*x for x in a)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def det3(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def kernel(points,weights):
    return {(i,j):mul(weights[i]*weights[j],sub(points[i],points[j])) for i in range(len(points)) for j in range(len(points))}
def mixD(D,P,Q):return tuple(sum(P[i]*Q[j]*D[i,j][k] for i in range(len(P)) for j in range(len(Q))) for k in range(3))
def moments(P,points,weights):
    W=sum(p*w for p,w in zip(P,weights));B=tuple(sum(P[i]*weights[i]*points[i][k] for i in range(len(P))) for k in range(3));return W,B

def run():
    points=[v(0,0,0),v(1,0,0),v(0,1,0),v(0,0,1),v(0,0,0)]
    weights=list(map(F,[1,1,1,1,2]));D=kernel(points,weights)
    distributions=[tuple(F(int(i==j)) for i in range(5)) for j in range(5)]
    distributions += [tuple(F(i==a or i==b,2) for i in range(5)) for a,b in combinations(range(5),2)]
    pairs=0
    for P,Q in product(distributions,repeat=2):
        Wp,Bp=moments(P,points,weights);Wq,Bq=moments(Q,points,weights)
        assert mixD(D,P,Q)==sub(mul(Wq,Bp),mul(Wp,Bq))
        for t in (F(1,3),F(1,2),F(4,5)):
            R=tuple((1-t)*p+t*q for p,q in zip(P,Q));Wr,Br=moments(R,points,weights)
            assert Wr==(1-t)*Wp+t*Wq and Br==add(mul(1-t,Bp),mul(t,Bq))
        pairs+=1
    # A genuine fourth affine moment coordinate survives duplicate position.
    assert D[0,4]==v(0,0,0)
    assert dot(v(1,0,0),D[0,1])/8==-F(1,8)
    assert dot(v(1,0,0),D[4,1])/8==-F(1,4)
    affine_columns=[sub((weights[i],)+mul(weights[i],points[i]),(weights[0],)+mul(weights[0],points[0])) for i in (1,2,3,4)]
    assert affine_columns==[v(0,1,0,0),v(0,0,1,0),v(0,0,0,1),v(1,0,0,0)]
    badpoints=points[:4]
    A={(i,j):F(2 if {i,j}=={1,2} else 1) for i in range(4) for j in range(4)}
    bad={(i,j):mul(A[i,j],sub(badpoints[i],badpoints[j])) for i in range(4) for j in range(4)}
    opposite=(A[0,1]*A[2,3],A[0,2]*A[1,3],A[0,3]*A[1,2])
    assert opposite==(1,1,2)
    r=v(F(1,2),0,F(-1,2));P=v(0,F(1,3),0,F(2,3));Q=v(0,0,1,0);R=v(1,0,0,0)
    tie=[dot(r,mixD(bad,X,Y))/8 for X,Y in ((P,Q),(Q,R),(P,R))]
    assert tie==[0,0,-F(1,48)] and dot(r,r)==F(1,2)
    r=v(F(1,4),F(1,2),F(3,4));P=v(0,F(1,6),F(2,3),F(1,6));Q=v(0,F(1,3),F(1,6),F(1,2));R=v(F(1,6),F(1,6),0,F(2,3))
    cycles=[dot(r,mixD(bad,X,Y))/8 for X,Y in ((P,Q),(Q,R),(R,P))]
    assert cycles==[F(1,1152),F(1,1152),F(1,576)] and dot(r,r)==F(7,8)
    norm2=max(dot(x,x)/64 for x in bad.values());assert norm2==F(1,8)<F(1,4)
    # Same one-shot moments, shared hidden label in two shots: higher joint moment differs.
    p_plus,p_minus=F(5,8),F(3,8)
    shared=(p_plus*p_plus+p_minus*p_minus)/2
    fresh=((p_plus+p_minus)/2)**2
    assert shared==F(17,64) and fresh==F(1,4) and shared-fresh==F(1,64)
    # Constant-gain lottery compatibility does not restrict a separately supplied mass-weighted meeting.
    first=F(1,2);second=F(1,3)
    assert (first-second)/8==F(1,48)
    return dict(round=1082,passed=True,arithmetic="Fraction exact",assertion_groups=6,
       factorized_moment_pairs=pairs,fourth_affine_moment_visible=True,
       duplicate_position_one_shot_probability_gap="1/8",
       nonfactorized_opposite_products=list(map(str,opposite)),binary_tie_probability_contrasts=list(map(str,tie)),
       strict_cycle_probability_contrasts=list(map(str,cycles)),maximum_effect_norm_squared=str(norm2),
       shared_label_two_shot_gap=str(shared-fresh),unconstrained_meeting_gap="1/48",
       scope="Finite exact certificates; general factorization proved separately. Moments cover balanced one-shot reports and fresh independent routing, not arbitrary instruments or shared-label histories.")
if __name__=="__main__":
    output=run();p=Path(__file__).with_name("independent_results.json")
    if "--write" in sys.argv:
        assert not p.exists();p.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    elif "--check" in sys.argv:assert json.loads(p.read_text(encoding="utf-8"))==output
    print(json.dumps(output,ensure_ascii=False,indent=2))
