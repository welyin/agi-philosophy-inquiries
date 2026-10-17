"""1255: exact rational material-clock reduction and initial curvature; read-only default."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt, sqrt
import argparse, hashlib, json

HERE=Path(__file__).resolve().parent
RESULT=HERE/'results.json'

def pack(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [pack(v) for v in x]
    return x

def sqrt_interval(x, digits=40):
    q=10**digits
    z=isqrt(x.numerator*q*q//x.denominator)
    lo,hi=F(z,q),F(z+1,q)
    assert lo*lo<=x<hi*hi
    return lo,hi

def dot(u,v): return sum((a*b for a,b in zip(u,v)),F(0))
def mv(h,u): return [dot(row,u) for row in h]

def canonical_checks():
    rows=[]
    for scale in (1,2,3):
        # Non-diagonal positive inverse metric; q = D^T h^-1 D.
        h=[[F(2),F(1,5),F(0)],[F(1,5),F(3),F(1,7)],[F(0),F(1,7),F(4)]]
        D=[F(scale,10),F(-scale,15),F(scale,20)]
        Q=dot(D,mv(h,D)); u=F(2*scale)
        A=(u*u+Q)/(2*u); r=(u*u-Q)/(2*u)
        n=F(scale+1,3); U=r/(3*n); J=r/n-U
        P=[d/n for d in D]; M=J+U; E=A/n
        N=A/r; beta=[v/r for v in mv(h,D)]
        assert A*A-Q==r*r and r>0 and J>0
        assert n*n*(dot(P,mv(h,P))+M*M)==A*A
        assert N*M/E==1
        assert [N*v/E-b for v,b in zip(mv(h,P),beta)]==[0,0,0]
        # h_ij beta^i beta^j = D_i h^ij D_j / r^2.
        assert -N*N+Q/(r*r)==-1
        C0=-A; dC=F(1,11); dD=[F(1,13),F(-1,17),F(1,19)]
        dh=[[F(1,23),F(1,29),F(0)],[F(1,29),F(-1,31),F(0)],[F(0),F(0),F(1,37)]]
        dU=F(1,41)
        dQ=2*dot(dD,mv(h,D))+dot(D,mv(dh,D))
        dr=(2*C0*dC-dQ)/(2*r)
        pullback=n*dU-dr
        full=n*dU+N*dC+dot(beta,dD)+dot(D,mv(dh,D))/(2*r)
        assert pullback==full
        assert n*J==r-n*U
        rows.append(dict(n=n,C0=C0,D=D,root=r,J=J,N=N,beta=beta,variation=full))
    # Independent finite differences of the actual square-root function.
    h=[[2.,.2,0.],[.2,3.,1/7],[0.,1/7,4.]]
    D=[.1,-1/15,.05]; C=-1.; n=2/3; U=.3
    dh=[[.1,.02,0.],[.02,-.03,0.],[0.,0.,.04]]; dD=[.03,-.02,.01]; dC=.04; dU=.05
    q=dot(D,mv(h,D)); r=sqrt(C*C-q)
    predicted=n*dU-C/r*dC+dot(mv(h,D),dD)/r+dot(D,mv(dh,D))/(2*r)
    def value(t):
        ht=[[h[i][j]+t*dh[i][j] for j in range(3)] for i in range(3)]
        dt=[D[i]+t*dD[i] for i in range(3)]
        return n*(U+t*dU)-sqrt((C+t*dC)**2-dot(dt,mv(ht,dt)))
    errors=[]
    for eps in (1e-3,1e-4,1e-5):
        err=abs((value(eps)-value(-eps))/(2*eps)-predicted)
        assert err<1e-8
        errors.append(err)
    return dict(exact_fixtures=rows,finite_difference_diagnostic_errors=errors,
                proof_scope='Finite-difference diagnostics do not replace symplectic pullback proof.')

def tensor_witness(d):
    k2=F(1,10); kap=F(3); H2=F(111,100); U=F(167,200)
    added=d*d/(2*k2); A=F(187,200)-added; root2=A*A-d*d
    assert A>0 and root2>U*U
    rl,rh=sqrt_interval(root2)
    Jlo,Jhi=rl-U,rh-U
    rho_rod=3*k2/2+added; rho_X=F(1,200); rho_rec=F(1,50); rho_clock=A
    rho=rho_rod+rho_X+rho_rec+rho_clock
    assert rho==H2 and -3*H2/kap+rho==0
    assert d-d==0
    S=[[-k2/2-F(1,200)+added if i==j else F(0) for j in range(3)] for i in range(3)]
    S[0][0]+=d*d/A
    eta=[-1,1,1,1]; R={}
    def set_component(a,b,c,e,v):
        for aa,bb,cc,ee,sign in [(a,b,c,e,1),(b,a,c,e,-1),(a,b,e,c,-1),(b,a,e,c,1),(c,e,a,b,1),(e,c,a,b,-1),(c,e,b,a,-1),(e,c,b,a,1)]:
            key=(aa,bb,cc,ee)
            assert key not in R or R[key]==sign*v
            R[key]=sign*v
    trS=sum(S[i][i] for i in range(3))
    # Spatial Gauss data and the Einstein Ricci tensor reconstruct all initial Riemann components.
    for i in range(3):
        for j in range(3):
            Rij=kap*(S[i][j]+(rho-trS)/2*(i==j))
            set_component(0,i+1,0,j+1,2*H2*(i==j)-Rij)
    for i in range(1,4):
        for j in range(1,4):
            for k in range(1,4):
                for l in range(1,4):
                    R[i,j,k,l]=H2*((i==k)*(j==l)-(i==l)*(j==k))
    ric=[[sum(eta[a]*R.get((a,b,a,e),F(0)) for a in range(4)) for e in range(4)] for b in range(4)]
    scalar=sum(eta[a]*ric[a][a] for a in range(4))
    T=[[F(0) for _ in range(4)] for _ in range(4)]; T[0][0]=rho
    for i in range(3):
        for j in range(3): T[i+1][j+1]=S[i][j]
    for b in range(4):
        for e in range(4): assert ric[b][e]-F(1,2)*(eta[b] if b==e else 0)*scalar==kap*T[b][e]
    def metric(a,b): return eta[a] if a==b else 0
    Weyl={}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for e in range(4):
                    Weyl[a,b,c,e]=R.get((a,b,c,e),F(0))-(metric(a,c)*ric[e][b]-metric(a,e)*ric[c][b]-metric(b,c)*ric[e][a]+metric(b,e)*ric[c][a])/2+scalar*(metric(a,c)*metric(e,b)-metric(a,e)*metric(c,b))/6
    invariant=sum(eta[a]*eta[b]*eta[c]*eta[e]*v*v for (a,b,c,e),v in Weyl.items())
    electric=[Weyl[0,i,0,i] for i in range(1,4)]
    assert sum(electric)==0 and invariant==8*sum(x*x for x in electric)
    return dict(d=d,C0=-A,root_squared=root2,J_interval=[Jlo,Jhi],rho=rho,
                complete_spatial_stress=S,electric_Weyl=electric,Weyl_squared=invariant,
                lapse_interval=[A/rh,A/rl],shift1_interval=[d/rh,d/rl],
                recorder_coordinate_speed_is_minus_shift=True,
                reference_clock_rate_interval=[A/(10*rh),A/(10*rl)])

def calculate():
    can=canonical_checks(); base=tensor_witness(F(0)); moving=tensor_witness(F(1,10))
    assert base['Weyl_squared']==0 and moving['Weyl_squared']==F(16,10443)
    assert moving['electric_Weyl']==[F(-2,177),F(1,177),F(1,177)]
    assert moving['J_interval'][0]>F(1,25) and moving['J_interval'][1]<F(9,200)
    A=F(177,200); d=F(1,10)
    # Copying the zero-momentum J=-h formula leaves a strictly positive constraint residual.
    lo,hi=sqrt_interval(A*A+d*d)
    residual=(lo-A,hi-A)
    assert residual[0]>F(1,200)
    true_rlo,true_rhi=sqrt_interval(A*A-d*d)
    h_gap=(A-true_rhi,A-true_rlo)
    assert h_gap[0]>0
    return pack(dict(round=1255,status='PASS',scientific_groups=3,
        canonical_and_source=can,initial_geometry={'baseline':base,'moving_rod':moving},
        deletion={'wrong_homogeneous_clock_constraint_residual':residual,'residual_lower':'1/200',
                  'physical_H_missing_momentum_correction_interval':h_gap},
        scope='Exact classical spatial-parent reduction and initial curvature witness; no finite-time measured quantum geometry.',
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--write',action='store_true'); args=ap.parse_args()
    data=calculate()
    if args.write:
        assert not RESULT.exists(), 'Refuse to overwrite saved scientific result.'
        RESULT.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('WROTE_RESULTS_1255')
    else:
        saved=json.loads(RESULT.read_text(encoding='utf-8-sig'))
        # Floating derivative checks are diagnostics; exact scientific values must agree byte-for-value.
        for obj in (saved,data): obj['canonical_and_source'].pop('finite_difference_diagnostic_errors')
        assert data==saved
        print('READ_ONLY_MATCH_1255: 3 scoped groups; exact constraints, complete source, initial curvature')
if __name__=='__main__': main()
