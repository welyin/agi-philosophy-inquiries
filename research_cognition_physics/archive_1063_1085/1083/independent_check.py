"""1083 independent exact audit.

All scientific arithmetic uses Fraction or Gaussian rational matrices.
Default prints JSON; --check only reads; --write replaces this file's
independent_results.json. No author assets are imported or modified.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import hashlib, json, os, sys

@dataclass(frozen=True)
class G:
    re: F=F(0)
    im: F=F(0)
    def __post_init__(self):
        object.__setattr__(self,"re",F(self.re))
        object.__setattr__(self,"im",F(self.im))
    def __add__(self,b):
        b=asG(b); return G(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self): return G(-self.re,-self.im)
    def __sub__(self,b): return self+-asG(b)
    def __rsub__(self,b): return asG(b)+-self
    def __mul__(self,b):
        b=asG(b);return G(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def conj(self): return G(self.re,-self.im)

def asG(x): return x if isinstance(x,G) else G(x)
def zero(n,m=None): return [[G() for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n):
    z=zero(n)
    for i in range(n):z[i][i]=G(1)
    return z
def ms(t,A):return [[t*x for x in row] for row in A]
def ma(A,B):return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def mm(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),G()) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A):return sum((A[i][i] for i in range(len(A))),G())
def tensor(A,B):
    return [[A[i][j]*B[k][l] for j in range(len(A[0])) for l in range(len(B[0]))] for i in range(len(A)) for k in range(len(B))]
def partialA(S,n=3):
    return [[sum((S[a*n+b][c*n+b] for b in range(n)),G()) for c in range(n)] for a in range(n)]
def partialB(S,n=3):
    return [[sum((S[a*n+b][a*n+c] for a in range(n)),G()) for c in range(n)] for b in range(n)]
def vec(*x):return tuple(map(F,x))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def scale(t,a):return tuple(t*x for x in a)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def textv(a):return list(map(str,a))
Z=vec(0,0,0)
X,Y,ZAX=vec(1,0,0),vec(0,1,0),vec(0,0,1)

def rho(p):
    x,y,z=p
    return [[G((1+z)/2),G(x/2,-y/2)],[G(x/2,y/2),G((1-z)/2)]]
def effect(d):return ma(ms(F(1,2),eye(2)),ms(2,ma(rho(d),ms(-F(1,2),eye(2)))))
def source(q,p):
    R=rho(p);A=zero(3)
    for i in range(2):
        for j in range(2):A[i][j]=q*R[i][j]
    A[2][2]=G(1-q);return A
def ready(p):return source(F(1),p)
VAC=source(F(0),Z)

def joint_flags(qi,qj,h,pi,pj):
    # Fixed local ready states; only the Bernoulli flags can be correlated.
    terms=[(h,tensor(ready(pi),ready(pj))),
           (qi-h,tensor(ready(pi),VAC)),
           (qj-h,tensor(VAC,ready(pj))),
           (1-qi-qj+h,tensor(VAC,VAC))]
    A=zero(9)
    for t,B in terms:
        assert t>=0
        A=ma(A,ms(t,B))
    return A

def comparator():
    # Q (dim2), source A (dim3), source B (dim3).
    # Both-ready: fair SWAP-test against A / flipped SWAP-test against B.
    # Otherwise: fair output. No readiness postselection.
    A=ms(F(1,2),eye(18))
    ix=lambda q,a,b:q*9+a*3+b
    for q,a,b in product(range(2),repeat=3):
        A[ix(a,q,b)][ix(q,a,b)]+=G(F(1,4))
        A[ix(b,a,q)][ix(q,a,b)]-=G(F(1,4))
    return A

def effective(global_effect,S):
    return [[sum((S[s][sp]*global_effect[q*9+sp][qp*9+s]
                  for s in range(9) for sp in range(9)),G())
             for qp in range(2)] for q in range(2)]
def mixD(K,P,Q):
    return tuple(sum((P[i]*Q[j]*K[i,j][k] for i in range(len(P)) for j in range(len(Q))),F(0)) for k in range(3))
def flags_distribution(correlated):
    if correlated:
        return [(F(1,8),(a,b,b,c)) for a,b,c in product((0,1),repeat=3)]
    return [(F(1,16),bits) for bits in product((0,1),repeat=4)]
def flag_moments(dist):
    q=[sum((p*bits[i] for p,bits in dist),F(0)) for i in range(4)]
    h={(i,j):sum((p*bits[i]*bits[j] for p,bits in dist),F(0)) for i in range(4) for j in range(4)}
    return q,h

class Audit:
    def __init__(self):self.n=0;self.groups=[]
    def eq(self,a,b,label):
        self.n+=1
        if a!=b:raise AssertionError(label+": "+str(a)+" != "+str(b))
    def yes(self,a,label):
        self.n+=1
        if not a:raise AssertionError(label)
    def group(self,name,**data):self.groups.append(dict(name=name,passed=True,**data))

def run():
    a=Audit();M=comparator();V=[Z,X,Y,ZAX]
    a.eq(M,[[M[j][i].conj() for j in range(18)] for i in range(18)],"global comparator is Hermitian")
    # Positive implementation is a convex instrument: four swap-test outcome
    # branches on the ready sector and fair branches on its orthogonal sector.
    ia=eye(18);L=zero(18);SA=zero(18);SB=zero(18)
    ix=lambda q,x,y:q*9+x*3+y
    for q,x,y in product(range(2),repeat=3):
        L[ix(q,x,y)][ix(q,x,y)]=G(1)
        SA[ix(x,q,y)][ix(q,x,y)]=G(1)
        SB[ix(y,x,q)][ix(q,x,y)]=G(1)
    for S in (SA,SB):
        a.eq(mm(S,S),L,"restricted SWAP squares to the ready projector")
        for sign in (F(-1),F(1)):
            proj=ms(F(1,2),ma(L,ms(sign,S)))
            a.eq(mm(proj,proj),proj,"swap-test fine effect is an exact projector")
    a.eq(M,ma(ma(ms(F(1,4),ma(L,SA)),ms(F(1,4),ma(L,ms(-1,SB)))),ms(F(1,2),ma(ia,ms(-1,L)))),"all branches summed without postselection")
    source_cases=[(F(1,2),F(3,4),vec(F(1,3),F(1,2),F(-1,4)),vec(F(-1,2),F(1,3),F(1,4))),
                  (F(1),F(1,2),X,Y),(F(0),F(1),Z,ZAX)]
    for qi,qj,pi,pj in source_cases:
        a.yes(dot(pi,pi)<=1 and dot(pj,pj)<=1,"ready states inside Bloch ball")
        S=tensor(source(qi,pi),source(qj,pj))
        E=effective(M,S)
        a.eq(E,effect(scale(qi*qj/8,sub(pi,pj))),"exact q_i q_j comparator effect including imaginary entries")
        a.eq(tr(S),G(1),"joint source normalized")
    a.group("fixed_quantum_comparator",dimension=18,arithmetic="Gaussian rational",source_cases=len(source_cases),coarse_effect="I/2 + q_i*q_j*(p_i-p_j).sigma/8",failure_branch="fair report, retained readiness record; never postselected")

    # Actual whole-source random SWAP, including source vacua.
    qx,qy=F(3,4),F(1,4);px=vec(F(1,3),0,F(1,2));py=vec(0,F(1,2),F(-1,2));lam=F(2,5)
    qnew=(1-lam)*qx+lam*qy
    Bnew=add(scale((1-lam)*qx,px),scale(lam*qy,py));pnew=scale(1/qnew,Bnew)
    swap=zero(9)
    for i,j in product(range(3),repeat=2):swap[j*3+i][i*3+j]=G(1)
    for S in (tensor(source(qx,px),source(qy,py)),joint_flags(qx,qy,F(1,8),px,py)):
        output=partialA(ma(ms(1-lam,S),ms(lam,mm(mm(swap,S),swap))))
        a.eq(output,source(qnew,pnew),"whole-source random SWAP has affine unnormalized moments")
        a.eq(output,ma(ms(1-lam,partialA(S)),ms(lam,partialB(S))),"meeting formula uses actual marginals even with correlated input flags")
    a.group("same_q_in_actual_meeting",lambda_right=str(lam),q_output=str(qnew),B_output=textv(Bnew),p_output=textv(pnew),scope="Selected source marginal; no equality of full instruments or retained environments is asserted.")

    # Two-anchor residual criterion and its four-scalar reconstruction.
    anchor0=vec(F(1,4),0,0);anchor1=vec(0,F(1,2),0);axis=sub(anchor1,anchor0)
    for deltaW,deltaB in [(F(0),Z),(F(1,7),vec(F(2,7),F(-3,7),F(4,7))),(-F(2,5),Y)]:
        R0=sub(deltaB,scale(deltaW,anchor0));R1=sub(deltaB,scale(deltaW,anchor1))
        recoveredW=-dot(axis,sub(R1,R0))/dot(axis,axis)
        recoveredB=add(R0,scale(recoveredW,anchor0))
        a.eq((recoveredW,recoveredB),(deltaW,deltaB),"three components plus one second-anchor direction reconstruct all residual moments")
    # Same V but changed W: signs agree while actual probability amplitude changes.
    qweak=qnew/2
    a.eq(sub(scale(qweak,pnew),scale(qweak,X)),scale(F(1,2),sub(Bnew,scale(qnew,X))),"weak representative differs by a positive scalar for every comparison")
    a.yes(qweak!=qnew,"weak sign representation does not fix the weight")
    a.group("actual_representation_test",exact_residual="R_i/(kappa*w_i)=Delta B-Delta W*v_i",scalar_checks=4,anchor_difference_squared=str(dot(axis,axis)),weak_scope="Only V is forced by all-r sign agreement; representative weight remains free.")

    # Positive attenuation is an erasure channel, not a free inverse.
    for q,p,alpha in [(F(1),X,F(1,2)),(F(1,2),Y,F(1,3)),(F(3,4),px,F(2,5))]:
        A=source(q,p);atten=ma(ms(alpha,A),ms(1-alpha,VAC))
        a.eq(atten,source(alpha*q,p),"identity/vacuum replacement implements source attenuation")
    target=F(1,4)
    def normalize_known(q,p):return source(target,p)
    # Same CPTP map cannot realize exact normalization for unknown q.
    qmix=F(3,4);pmix=scale(F(1,3),X)
    nonlinear_left=normalize_known(qmix,pmix)
    nonlinear_right=ma(ms(F(1,2),normalize_known(F(1),Z)),ms(F(1,2),normalize_known(F(1,2),X)))
    a.yes(nonlinear_left!=nonlinear_right,"normalization is nonaffine on unknown mixed sources")
    # Equalization is an additional known-q protocol; bare meeting is different.
    bare1=(F(1),scale(F(1,2),X));bare2=(F(3,4),scale(F(2,3),X))
    a.eq(scale(bare1[0],bare1[1]),scale(bare2[0],bare2[1]),"bare outputs have equal B despite different V")
    d1=scale(bare1[0]/8,sub(bare1[1],X));d2=scale(bare2[0]/8,sub(bare2[1],X))
    a.eq(dot(X,sub(d2,d1)),F(1,32),"nonzero anchor detects bare-meeting weight difference")
    a.eq(dot(X,sub(scale(bare1[0]/8,bare1[1]),scale(bare2[0]/8,bare2[1]))),F(0),"zero anchor does not detect this particular difference")
    a.group("attenuation_and_balancing_boundary",known_weight_required=True,universal_unknown_normalization_is_affine=False,bare_normalized_outputs=["e1/2","2*e1/3"],nonzero_anchor_probability_gap="1/32",zero_anchor_probability_gap="0",balanced_source="q=c and p=(1-lambda)*p_x+lambda*p_y after explicit known-q equalization")

    # Same local source states, same comparator; only joint readiness changes.
    qi,hi=flag_moments(flags_distribution(False));qc,hc=flag_moments(flags_distribution(True))
    a.eq(qi,[F(1,2)]*4,"independent flags fair")
    a.eq(qc,qi,"correlated flags have identical local readiness probabilities")
    kernels={}
    for name,h in [("independent",hi),("correlated",hc)]:
        K={}
        for i,j in combinations(range(4),2):
            expected_h=F(1,2) if name=="correlated" and {i,j}=={1,2} else F(1,4)
            a.eq(h[i,j],expected_h,"enumerated joint ready probability")
            S=joint_flags(F(1,2),F(1,2),h[i,j],V[i],V[j])
            a.eq(partialA(S),source(F(1,2),V[i]),"left local state unchanged")
            a.eq(partialB(S),source(F(1,2),V[j]),"right local state unchanged")
            d=scale(h[i,j]/8,sub(V[i],V[j]))
            a.eq(effective(M,S),effect(d),"unchanged fixed comparator yields the h_ij gain")
            K[i,j]=d;K[j,i]=scale(-1,d)
        for i in range(4):K[i,i]=Z
        kernels[name]=K
    a.group("correlated_readiness_with_fixed_local_sources",joint_atoms=8,all_local_q="1/2",h12="1/2",all_other_distinct_h="1/4",local_states_equal=True,processor_equal=True,scope="A single joint preparation realizes all pair marginals. Independence is a model input, not a universal cognitive axiom.")

    P=vec(0,F(1,6),F(2,3),F(1,6));Q=vec(0,F(1,3),F(1,6),F(1,2));R=vec(F(1,6),F(1,6),0,F(2,3))
    r=vec(F(1,4),F(1,2),F(3,4));K=kernels["correlated"]
    cycle=[dot(r,mixD(K,A,B)) for A,B in ((P,Q),(Q,R),(R,P))]
    a.eq(cycle,[F(1,4608),F(1,4608),F(1,2304)],"exact three-lottery strict cycle")
    a.yes(dot(r,r)<=1,"cycle probe legal")
    tieP=vec(0,F(1,3),0,F(2,3));pure2=vec(0,0,1,0);pure0=vec(1,0,0,0);rt=vec(F(1,2),0,F(-1,2))
    tie=[dot(rt,mixD(K,A,B)) for A,B in ((tieP,pure2),(pure2,pure0),(tieP,pure0))]
    a.eq(tie,[F(0),F(0),-F(1,192)],"exact binary tie-substitution failure")
    a.eq(max(dot(d,d) for d in K.values()),F(1,128),"every effective binary effect strictly physical")
    # Matrix probability and vector formula agree for an imaginary-entry case.
    E=effect(K[1,2]);prob=tr(mm(rho(r),E))
    a.eq(prob,G(F(1,2)+dot(r,K[1,2])),"Born probability agrees with Pauli normalization")
    a.group("finite_probability_failures",strict_cycle_contrasts=textv(cycle),binary_tie_contrasts=textv(tie),maximum_effect_vector_norm_squared="1/128",counterexample="Only joint readiness was changed; local p_i,q_i and the physical comparator stay fixed.")

    # Two distinct ways actual construction can fail full representation.
    # General comparison weights here are not source readiness probabilities.
    kap=F(1,8);average=scale(F(1,2),X);wrong=scale(F(1,3),X);rt=vec(F(1,2),F(1,4),0)
    mixed_to_y=scale(kap,sub(average,Y));wrong_to_y=scale(kap,sub(wrong,Y))
    a.eq(dot(rt,mixed_to_y),F(0),"random representative tied with e2")
    a.eq(dot(rt,wrong_to_y),-F(1,96),"strict-between wrong representative violates full comparison")
    p_lottery=F(1,2)+dot(X,scale(kap,average))
    p_double_weight=F(1,2)+dot(X,scale(2*kap,average))
    a.eq((p_lottery,p_double_weight,p_double_weight-p_lottery),(F(9,16),F(5,8),F(1,16)),"correct V with wrong weight preserves sign but changes probability")
    finite_partners=[(F(1),Z),(F(1),X),(F(1),Y),(F(1),ZAX),(F(1),wrong),(F(2),average)]
    anchor_det=X[0]*(Y[1]*ZAX[2]-Y[2]*ZAX[1])-X[1]*(Y[0]*ZAX[2]-Y[2]*ZAX[0])+X[2]*(Y[0]*ZAX[1]-Y[1]*ZAX[0])
    a.eq(anchor_det,F(1),"representation counterexamples retain four affine-rank-three anchors")
    maximum=max(dot(scale(kap*w*t,sub(p,q)),scale(kap*w*t,sub(p,q))) for w,p in finite_partners for t,q in finite_partners)
    a.yes(maximum<F(1,4),"both representation counterexamples use legal effects")
    a.group("two_actual_representation_counterexamples",affine_anchor_determinant=str(anchor_det),wrong_position_probability_contrast="-1/96",correct_position_probabilities=["9/16","5/8"],wrong_weight_probability_gap="1/16",maximum_effect_vector_norm_squared=str(maximum),scope="General positive weights, not q>1 readiness. Betweenness and ordinal representation each fall short of full probability representation.")

    return dict(round=1083,passed=True,arithmetic="Fraction and exact Gaussian-rational matrices",group_count=len(a.groups),assertion_count=a.n,groups=a.groups,
        scope=dict(general_representation_theorem_proved_by_finite_examples=False,
            actual_spatial_dimension_derived=False,neutral_attenuation_is_original_contact_equivalence=False,
            independent_source_preparations_are_universal_cognitive_axiom=False,
            passive_reference_full_instrument_equality_claimed=False,
            output_readiness_postselected=False),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=="__main__":
    result=run();path=Path(__file__).with_name("independent_results.json")
    if "--write" in sys.argv:
        tmp=path.with_name(path.name+".tmp-review1068")
        with tmp.open("w",encoding="utf-8",newline="\n") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n");f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    if "--check" in sys.argv:
        assert result==json.loads(path.read_text(encoding="utf-8")),"saved result mismatch"
    print(json.dumps(result,ensure_ascii=False,indent=2))
