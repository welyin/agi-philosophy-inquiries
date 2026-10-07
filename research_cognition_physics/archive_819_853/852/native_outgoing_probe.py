"""852: exact causal transfer to a genuinely later probe test.

The rational 1+1 lattice KG operator is a support/adjoint calibration, not
an integration of the original curved spacetime PDE. The CAR algebra is
imported from the original round 829 code without altering its files.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, importlib.util, json, sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
TARGET=HERE/'native_outgoing_probe_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code
import relational_probe_budget_probe as budget

def zeros(nt,nx):return [[F(0) for x in range(nx)] for t in range(nt)]
def dot(a,b):return sum(x*y for ra,rb in zip(a,b) for x,y in zip(ra,rb))
def transpose(a):return [list(c) for c in zip(*a)]
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def inv(a):
    n=len(a);q=[list(a[i])+[F(i==j) for j in range(n)] for i in range(n)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if q[i][k]);q[k],q[pivot]=q[pivot],q[k]
        d=q[k][k];q[k]=[x/d for x in q[k]]
        for i in range(n):
            if i!=k:
                d=q[i][k];q[i]=[x-d*y for x,y in zip(q[i],q[k])]
    return [row[n:] for row in q]

def run():
    nt,nx=12,27;mass2=F(3,17);speed2=F(2,11);lam=F(2,7)
    def spatial(row):return [mass2*row[x]+speed2*(2*row[x]-row[(x-1)%nx]-row[(x+1)%nx]) for x in range(nx)]
    def op(a):
        out=zeros(nt,nx)
        for t in range(1,nt-1):
            c=spatial(a[t]);out[t]=[a[t+1][x]-2*a[t][x]+a[t-1][x]+c[x] for x in range(nx)]
        return out
    def homogeneous(center,second_scale=F(1)):
        v=zeros(nt,nx)
        for x in range(nx):
            d=min((x-center)%nx,(center-x)%nx)
            v[4][x]=F(d<=2);v[5][x]=second_scale*F(d<=2)
        for t in range(5,nt-1):
            c=spatial(v[t]);v[t+1]=[2*v[t][x]-v[t-1][x]-c[x] for x in range(nx)]
        for t in range(4,0,-1):
            c=spatial(v[t]);v[t-1]=[2*v[t][x]-v[t+1][x]-c[x] for x in range(nx)]
        assert not any(x for row in op(v) for x in row)
        return v
    def retarded(s):
        p=zeros(nt,nx)
        for t in range(1,nt-1):
            c=spatial(p[t]);p[t+1]=[2*p[t][x]-p[t-1][x]-c[x]+s[t][x] for x in range(nx)]
        assert op(p)[1:-1]==s[1:-1]
        return p
    centers=(3,12,21);vs=[homogeneous(c) for c in centers];bs=[];os=[];fs=[]
    chi=[F(0) if t<=5 else F(1,3) if t==6 else F(2,3) if t==7 else F(1) for t in range(nt)]
    for k,c in enumerate(centers):
        b=zeros(nt,nx);o=zeros(nt,nx)
        for t,w in ((3,F(2,7)),(4,F(5,7))):
            for dx,a in ((-1,F(1,4)),(0,F(1,2)),(1,F(1,4))):
                x=c+dx;b[t][x]=w*a
                assert vs[k][t][x]!=0
                o[t][x]=b[t][x]/vs[k][t][x]
                assert all(vs[j][t][x]==0 for j in range(3) if j!=k)
        bs.append(b);os.append(o)
        z=[[chi[t]*x for x in vs[k][t]] for t in range(nt)]
        f=[[-x for x in row] for row in op(z)];fs.append(f)
        advanced=[[(1-chi[t])*x for x in vs[k][t]] for t in range(nt)]
        assert op(advanced)==f
    responses=[retarded([[-lam*x for x in row] for row in o]) for o in os]
    transfers=[[dot(f,p) for p in responses] for f in fs]
    assert transfers==[[(-lam if i==j else F(0)) for j in range(3)] for i in range(3)]
    raw=[retarded([[-lam*x for x in row] for row in b]) for b in bs]
    uncalibrated=[[dot(f,p) for p in raw] for f in fs]
    assert uncalibrated!=transfers
    # A future readout is not an equation of motion smear with compact j:
    # it has nonzero pairing with a free homogeneous solution.
    free_pairings=[dot(fs[i],homogeneous(c,F(2))) for i,c in enumerate(centers)]
    assert all(x!=0 for x in free_pairings)
    source_times=sorted({t for o in os for t,row in enumerate(o) if any(row)})
    output_times=sorted({t for f in fs for t,row in enumerate(f) if any(row)})
    assert max(source_times)<min(output_times)
    # Non-orthogonal material gradients: Q^T B Q reproduces b/(kappa v).
    q=[[F(x) for x in row] for row in ((2,1,0,1),(0,3,1,0),(1,0,2,1),(0,1,0,2))]
    qi=inv(q);stress=[[F((a+1)*(b+1),13)+F(a==b,7) for b in range(4)] for a in range(4)]
    tensor_res=[]
    for k in range(3):
        b=[[F((a+b+k)%5-2,11) for b in range(4)] for a in range(4)]
        v=vs[k][3][centers[k]];kappa=F(17)
        prof=[[x/(kappa*v) for x in row] for row in mm(mm(transpose(qi),b),qi)]
        reconstructed=mm(mm(transpose(q),prof),q)
        err=max(abs(reconstructed[a][c]*kappa*v-b[a][c]) for a in range(4) for c in range(4))
        assert err==0
        assert dot(stress,reconstructed)*v/2==dot(stress,b)/(2*kappa)
        tensor_res.append(str(err))
    # Original 829 six-Majorana content witness, including the mixed complement.
    gamma,_,_,_,comp,sx,sz=code.code_data()
    sy=[code.mul((0,0,3),code.mul(sz[b],sx[b])) for b in range(5)]
    a=[code.mul((0,0,2),x) for x in (sz[0],sz[1],sy[3])]
    logical=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    eps=F(1,10);r0=F(0);r1=F(3,5);d=1024
    def phase(p):
        assert p in (0,2)
        return F(1 if p==0 else -1)
    def core(p):
        typ,ph=comp(p)
        return phase(ph) if typ=='scalar' else F(0)
    def ex(p,r):
        full=phase(p[2]) if p[:2]==(0,0) else F(0)
        return ((1-eps)-2*eps/(d-2))*core(p)+(1-eps)*r*core(code.mul(logical,p))+eps*d/(d-2)*full
    def product(indices):return code.product(a[i] for i in indices)
    assert all(ex(product(ix),r1)==ex(product(ix),r0) for ix in [(i,) for i in range(3)]+[(i,j) for i in range(3) for j in range(3)])
    def cum(r):
        return ex(product((0,1,2)),r)-sum(ex(product(ix),r)*ex(a[k],r) for ix,k in (((0,1),2),((0,2),1),((1,2),0)))+2*ex(a[0],r)*ex(a[1],r)*ex(a[2],r)
    source_delta=cum(r1)-cum(r0);assert source_delta==F(27,50)
    output_delta=source_delta*transfers[0][0]*transfers[1][1]*transfers[2][2]
    assert output_delta==-lam**3*F(27,50)
    # At epsilon=sqrt(hbar) order six, content requires all six fermion legs.
    allowed=[]
    for nf in range(7):
        for nb in range(7):
            for loops in range(4):
                for shifts in range(nb+1):
                    if nf>=6 and nf+nb+2*loops+shifts==6:allowed.append([nf,nb,loops,shifts])
    assert allowed==[[6,0,0,0]]
    balance=budget.run();assert balance==json.loads(budget.TARGET.read_text('utf-8'))
    return dict(round=852,all_checks_passed=True,fresh_test_groups=1,
        diagnostic='exact finite KG support/adjoint calibration plus original 829 CAR algebra; not the original continuum PDE',
        lattice=[nt,nx],coupling=str(lam),source_time_sites=source_times,output_time_sites=output_times,
        calibrated_transfer=[[str(x) for x in row] for row in transfers],
        uncalibrated_transfer=[[str(x) for x in row] for row in uncalibrated],
        nonzero_pairing_with_free_probe=[str(x) for x in free_pairings],
        material_tensor_residuals=tensor_res,
        original_source_content_difference=str(source_delta),outgoing_content_difference=str(output_delta),
        epsilon_six_allowed_legs=allowed,
        dynamical_window_budget_reproduced=True,
        omitted_clock_energy_rates=[x['frozen_switch_missing_energy_rate'] for x in balance['rows']],
        full_curved_PDE_or_loop_coefficients_computed=False,
        local_formal_extension_and_continuum_transfer_justification='research_note_852.md',
        finite_coupling_probabilities_or_autonomous_apparatus_claimed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
