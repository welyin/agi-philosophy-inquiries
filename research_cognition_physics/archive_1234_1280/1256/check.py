"""Round1256: exact recorder cotangent, local source and deletion checks.
Default is read-only; --write is allowed only for the initial result creation.
General field identities are proved in proof.md; fixtures do not prove them.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import argparse, json, hashlib
HERE=Path(__file__).resolve().parent

def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def mv(a,b): return [dot(row,b) for row in a]
def tr(a): return [list(x) for x in zip(*a)]
def det(a):
    return sum(((-1)**j*a[0][j]*(a[1][(j+1)%3]*a[2][(j+2)%3]-a[1][(j+2)%3]*a[2][(j+1)%3])*(-1)**j for j in range(3)),F(0))
def inv(a):
    v=[list(row)+[F(i==j) for j in range(3)] for i,row in enumerate(a)]
    for k in range(3):
        z=next(i for i in range(k,3) if v[i][k]); v[k],v[z]=v[z],v[k]
        pivot=v[k][k];v[k]=[x/pivot for x in v[k]]
        for i in range(3):
            if i!=k:
                c=v[i][k];v[i]=[x-c*y for x,y in zip(v[i],v[k])]
    return [row[3:] for row in v]

def sr(x):
    x=F(x);a,b=isqrt(x.numerator),isqrt(x.denominator)
    assert a*a==x.numerator and b*b==x.denominator
    return F(a,b)

class Dual:
    def __init__(self,v,d=0): self.v,self.d=F(v),F(d)
    @staticmethod
    def cast(x): return x if isinstance(x,Dual) else Dual(x)
    def __add__(a,b):
        b=Dual.cast(b);return Dual(a.v+b.v,a.d+b.d)
    __radd__=__add__
    def __neg__(a): return Dual(-a.v,-a.d)
    def __sub__(a,b): return a+-Dual.cast(b)
    def __rsub__(a,b): return Dual.cast(b)+-a
    def __mul__(a,b):
        b=Dual.cast(b);return Dual(a.v*b.v,a.d*b.v+a.v*b.d)
    __rmul__=__mul__
    def __truediv__(a,b):
        b=Dual.cast(b);return Dual(a.v/b.v,(a.d*b.v-a.v*b.d)/(b.v*b.v))
    def __rtruediv__(a,b): return Dual.cast(b)/a
    def root(a):
        v=sr(a.v);return Dual(v,a.d/(2*v))

def pack(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [pack(v) for v in x]
    return x

def canonical():
    fixtures=[]
    for s in (1,2,3):
        B=[[F(1+s,2),F(1,3),F(-1,7)],[F(0),F(3,2),F(1,5)],[F(0),F(0),F(4,3)]]
        E=inv(B); jac=det(B); assert jac>0
        assert [[dot(row,col) for col in tr(E)] for row in B]==[[F(i==j) for j in range(3)] for i in range(3)]
        n=F(1,10+s);eta=n*jac
        P=[F(s,7),F(-1,11),F(2,13)];g=[F(1,3),F(s,5),F(-1,7)];I=F(7,5)
        dz=[F(1,17),F(-2,19),F(s,23)];dphi=F(3,29)
        dx=[-x for x in mv(E,dz)];dphiL=dphi+dot(g,dx)
        Pi=eta*I;pa=[-eta*dot([P[i]+I*g[i] for i in range(3)],col) for col in tr(E)]
        old=eta*(dot(P,dx)+I*dphiL); new=dot(pa,dz)+Pi*dphi
        assert old==new
        D=[dot(pa,col)+Pi*g[i] for i,col in enumerate(tr(B))]
        assert D==[-eta*x for x in P]
        reconstructed=[-d/eta for d in D]; assert reconstructed==P
        assert [-dot(pa,col)/eta-I*g[i] for i,col in enumerate(tr(B))]==P
        fixtures.append(dict(jacobian=jac,label_density=n,eta=eta,Pi=Pi,pi_labels=pa,one_form=old,D=D,mechanical_P=P))
    return dict(exact_fixtures=fixtures,scope='Same original cotangent one-form, invertible map and mechanical momentum; no independent density added.')

def transport_source():
    B=[[F(2),F(1,3),F(0)],[F(0),F(1),F(1,5)],[F(0),F(0),F(1,2)]];E=inv(B)
    eta=F(1,5)*det(B); p=F(1,4);g=[F(1),F(1,3),F(-1,5)]
    h=[[F(2),F(1,5),F(0)],[F(1,5),F(3),F(1,7)],[F(0),F(1,7),F(4)]]
    m=[F(1,10),F(-1,20),F(1,30)];Q=dot(m,mv(h,m));t=F(2)
    Er=(t*t+Q)/(2*t);sigma=(t*t-Q)/(2*t);assert sigma>0 and Er*Er==Q+sigma*sigma
    Pi=(10*sigma-eta)/(1+p);I=Pi/eta;M=sigma/eta;P=[x/eta for x in m]
    assert M==(1+(1+p)*I)/10 and Pi>0
    pa=[-dot([m[i]+Pi*g[i] for i in range(3)],col) for col in tr(E)]
    D=[-x for x in m];otherD=[F(1,17),F(-1,19),F(1,23)];Dt=[a+b for a,b in zip(D,otherD)]
    Qt=dot(Dt,mv(h,Dt));u=F(4);A=(u*u+Qt)/(2*u);root=(u*u-Qt)/(2*u);Cother=-A-Er
    N=A/root;beta=[x/root for x in mv(h,Dt)];v=[N*x/Er-b for x,b in zip(mv(h,m),beta)]
    U=F(1,4); assert root-U>0
    # Differentiate the full reduced square root, not a frozen reconstructed NE+betaD.
    def H(pa0,Pi0,p0,h0):
        mm=[-dot(pa0,col)-Pi0*g[i] for i,col in enumerate(tr(B))]
        ss=(eta+(1+p0)*Pi0)/10
        ee=Dual.cast(dot(mm,mv(h0,mm))+ss*ss).root()
        dd=[otherD[i]-mm[i] for i in range(3)]
        return U-(Dual.cast(Cother+ee)*(Cother+ee)-dot(dd,mv(h0,dd))).root()
    label_rates=[]
    for a in range(3):
        pp=[Dual(x,F(i==a)) for i,x in enumerate(pa)]
        rate=H(pp,Pi,p,h).d;expected=-dot(B[a],v);assert rate==expected;label_rates.append(rate)
    phase=H(pa,Dual(Pi,1),p,h).d
    material_rate=phase+dot(v,g)
    assert material_rate==N*sigma*(1+p)/(10*Er)
    e_particle=Er/eta
    assert material_rate==N*M/e_particle*(1+p)/10
    dHdp=H(pa,Pi,Dual(p,1),h).d
    source=-N*sigma*Pi/(10*Er)
    assert -dHdp==source==-N*eta*M/e_particle*(I/10)
    assert Er==eta*e_particle
    assert [[m[i]*m[j]/Er for j in range(3)] for i in range(3)]==[[eta*P[i]*P[j]/e_particle for j in range(3)] for i in range(3)]
    # Explicit geometry contribution of both fluids, including original clock root variation.
    dh=[[F(1,7),F(1,11),F(0)],[F(1,11),F(-1,13),F(1,17)],[F(0),F(1,17),F(1,19)]]
    hd=[[Dual(h[i][j],dh[i][j]) for j in range(3)] for i in range(3)]
    geom_derivative=H(pa,Pi,p,hd).d
    assert geom_derivative==N*dot(m,mv(dh,m))/(2*Er)+dot(Dt,mv(dh,Dt))/(2*root)
    # Continuity cancellation for arbitrary jets: both densities transport with v.
    geta=[F(1,7),F(-1,11),F(1,13)];gPi=[F(-1,17),F(1,19),F(1,23)];divv=F(1,29)
    detadt=-dot(v,geta)-eta*divv;dPidt=-dot(v,gPi)-Pi*divv
    dIdt=(eta*dPidt-Pi*detadt)/(eta*eta);gradI=[(eta*b-Pi*a)/(eta*eta) for a,b in zip(geta,gPi)]
    assert dIdt+dot(v,gradI)==0
    # Determinant identity for Z advection, with arbitrary derivative of v and n(Z).
    V=[[F(1,31),F(1,37),F(0)],[F(-1,41),F(1,43),F(0)],[F(0),F(1,47),F(-1,53)]]
    LB=[[-dot(row,col) for col in tr(V)] for row in B]
    detmaterial=det(B)*sum((dot(E[i],tr(LB)[i]) for i in range(3)),F(0))
    assert detmaterial==-det(B)*sum((V[i][i] for i in range(3)),F(0))
    return dict(N=N,beta=beta,v=v,eta=eta,Er=Er,Pi=Pi,label_rates=label_rates,phase_material_rate=material_rate,p_source=source,full_metric_directional_derivative=geom_derivative,material_I_derivative=F(0),scope='Exact derivatives of original1255 root and original worldline source; field transport follows the analytic functional proof.')

def deletion():
    eta=F(1,10);I=F(1);Pi=eta*I;g=F(1);pa=-Pi*g
    m=-(pa+Pi*g);bad=-pa;M=F(1,5)
    goodE=eta*M;badE2=bad*bad+(eta*M)**2
    assert m==0 and bad==F(1,10) and goodE==F(1,50) and badE2==F(26,2500)
    bound=F(2,25); assert badE2>(goodE+bound)**2
    # A phase profile on the original period L=10 has derivative1 at x=0:
    # phi(x)=L/(2*pi)*sin(2*pi*x/L). No global nonperiodic preparation needed.
    return dict(correct_m=m,false_m=bad,correct_energy=goodE,false_energy_squared=badE2,strict_energy_error_lower=bound,proof='26>25 implies (sqrt(26)-1)/50>2/25',periodic_profile='L/(2*pi)*sin(2*pi*x1/L), L=10; derivative1 at x1=0',scope='Rejects omitted internal-gradient term in the same canonical variables; not a new finite future prediction.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=pack(dict(round=1256,status='PASS',scientific_groups=3,canonical=canonical(),transport_source=transport_source(),deletion=deletion(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),new_adopted_axioms=0,CT='unadopted_conjecture'))
    out=HERE/'results.json'
    if args.write:
        assert not out.exists();out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');mode='INITIAL_WRITE'
    else:
        assert json.loads(out.read_text('utf-8-sig'))==result;mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status='PASS',mode=mode,scientific_groups=3)))
if __name__=='__main__':main()
