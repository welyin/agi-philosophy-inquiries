"""Independent exact 1188 review. No author code or results imported.
Default read-only; --write creates review_checks.json exclusively.
"""
from fractions import Fraction as F
from pathlib import Path
import json,sys
BASE=Path(__file__).resolve().parent

def trn(a):return [list(r) for r in zip(*a)]
def mm(a,b):return [[sum((x*y for x,y in zip(r,c)),F(0)) for c in zip(*b)] for r in a]
def add(a,b):return [[x+y for x,y in zip(r,t)] for r,t in zip(a,b)]
def sc(s,a):return [[s*x for x in r] for r in a]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def inv(a):
    n=len(a);b=[r[:]+t for r,t in zip(a,eye(n))]
    for i in range(n):
        j=next(j for j in range(i,n) if b[j][i]);b[i],b[j]=b[j],b[i]
        p=b[i][i];b[i]=[v/p for v in b[i]]
        for j in range(n):
            if j!=i:
                p=b[j][i];b[j]=[v-p*w for v,w in zip(b[j],b[i])]
    return [r[n:] for r in b]
def trace(a):return sum((a[i][i] for i in range(len(a))),F(0))
def tt(a,hi):return trace(mm(hi,a))
def nn(a,hi):u=mm(hi,a);return trace(mm(u,u))

def geometry():
    L=[[F(1),F(0),F(0)],[F(1,2),F(1),F(0)],[-F(1,3),F(2,5),F(1)]]
    h=mm(L,trn(L));hi=inv(h)
    kh=[[F(2),F(1,3),-F(1,5)],[F(1,3),F(1),F(0)],[-F(1,5),F(0),F(3)]]
    K=mm(mm(L,kh),trn(L))
    ah=[[[F(1),F(0),F(0)],[F(0),F(-1),F(0)],[F(0),F(0),F(0)]],
        [[F(1),F(0),F(0)],[F(0),F(0),F(0)],[F(0),F(0),F(0)]]]
    shapes=[mm(mm(L,a),trn(L)) for a in ah]
    alpha=F(5,4);gamma=F(3,4);lam=F(2,5);s=F(1,2)
    ambient=F(5,7);Lambda=F(2,9)
    def R(k):return ambient+nn(k,hi)-tt(k,hi)**2
    def kinetic(k):
        p=sc(1/alpha,add(k,sc(lam*tt(k,hi)/(1-3*lam),h)))
        assert sc(alpha,add(p,sc(-lam*tt(p,hi),h)))==k
        return alpha*(nn(p,hi)-lam*tt(p,hi)**2)
    rho=-kinetic(K)+gamma*R(K)-2*Lambda
    assert rho>0
    def H(k):return kinetic(k)-gamma*R(k)+2*Lambda+rho
    assert H(K)==0
    rows=[]
    for A in shapes:
        kp=add(K,sc(s,A));km=add(K,sc(-s,A))
        dR=R(kp)+R(km)-2*R(K)
        assert dR==2*s*s*(nn(A,hi)-tt(A,hi)**2)
        actual=H(kp)+H(km)-2*H(K)
        expected=2*s*s*((1/alpha-gamma)*nn(A,hi)+(gamma+lam/(alpha*(1-3*lam)))*tt(A,hi)**2)
        assert actual==expected
        shifted=lambda k:H(k)+F(89,23)-gamma*F(13,19)
        assert shifted(kp)+shifted(km)-2*shifted(K)==actual
        rows.append({'trace_A':str(tt(A,hi)),'norm_A_squared':str(nn(A,hi)),
          'Gauss_center_difference':str(dR),'constraint_center_difference':str(actual)})
    return {'metric':'nonorthogonal rational frame; det h=1','alpha':str(alpha),'gamma':str(gamma),
      'Hamilton_lambda':str(lam),'Lagrange_lambda':str(lam/(3*lam-1)),
      'same_complete_source':str(rho),'base_constraint_residual':'0','s':str(s),'shape_checks':rows,
      'scope':'local tensor identity, not a prepared physical quantum geometry'}

def selection():
    rows=[[F(2),F(0),F(-2)],[F(1),F(-1),F(0)]]
    assert all(sum(r)==0 for r in rows)
    assert rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]!=0
    for alpha in (F(1),F(5,4)):
        lam=F(1,2);gamma=1/alpha
        assert 1/alpha==lam/(alpha*(3*lam-1))==gamma
    p=[[F(2),F(1),F(0)],[F(1),F(-1),F(0)],[F(0),F(0),F(3)]]
    assert trace(add(p,sc(-trace(p)/3,eye(3))))==0
    assert F(1,2)-(-F(1,2))==1
    assert 3*F(3,5)-1==F(4,5)
    assert F(8)/F(5,2)==F(16,5)
    for lam in (F(2,5),F(9,20),F(1,2),F(11,20),F(3,5)):
        C=(1-2*lam)/(1-3*lam)
        assert abs(C)>=F(5,2)*abs(lam-F(1,2))
    return {'velocity_coefficients':'a=b=c','selected_Hamilton_lambda':'1/2',
      'selected_product':'alpha*gamma=1','lambda_one_third':'Legendre image has trace K=0; three nonzero trace offsets cannot all belong',
      'finite_error_domain':'lambda in [2/5,3/5], alpha>0, s=1/2, rank-one trace-one A',
      'per_slice_total_residual_epsilon_bound':'abs(lambda-1/2) <= (16/5)*alpha*epsilon',
      'scope':'conditional propagation of independently certified total residual error, not an instrument certificate'}

def records():
    # Full history is retained. This coherently computes a basis-label function;
    # it does not clone an unknown state or discard correlations.
    images=[]
    for h in range(3):
        for s in range(3):
            for y in range(3):
                f=(s*h+s*s)%3;out=(h,s,(y+f)%3)
                assert (out[2]-f)%3==y
                images.append(out)
    assert len(set(images))==27
    for rdim in (2,3):assert len({(out,r) for out in images for r in range(rdim)})==27*rdim
    lam=F(2,5);s=F(1,2);C=(1-2*lam)/(1-3*lam)
    assert 2*s*s*C==F(-1,2)
    return {'full_relabeling_dimension':27,'unitarity':'exact reversible basis permutation',
      'preserved':'complete history, coherent slice/label registers, arbitrary reference',
      'independent_canonical_countercomparison_lambda':'2/5','canonical_center_residual':'-1/2',
      'scope':'record permission and CT are different; no complete quantum lambda-R countermodel claimed'}

def evaluate():
    return {'round':1188,'review_status':'PASS_FOR_CONDITIONAL_THREE_SLICE_IDENTITY','author_module_imported':False,
      'groups':{'1_same_event_source_and_geometry':geometry(),'2_selection_degenerate_and_error':selection(),'3_coherent_records_not_CT':records()},
      'scope_limits':['Smooth effective Lorentz history and same complete source/geometric dictionary are assumed.',
       'The local pure-metric quadratic momentum and scalar-curvature Hamiltonian class is assumed.',
       'Nonminimal or higher-derivative momentum mixing requires a new Legendre/source audit.',
       'CT is a candidate condition, not an adopted axiom or a consequence of record relabeling alone.',
       'No physical quantum geometry encoding, finite instrument certificate, nonlinear GR construction or empirical result.'],
      'new_axioms':0,'new_empirical_calibration':0,'gr_hard_deliverables':0}
if __name__=='__main__':
    result=evaluate();out=BASE/'review_checks.json'
    if '--write' in sys.argv:
        with out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WROTE independent 1188 rational audit: 3 groups PASS')
    else:
        assert result==json.loads(out.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH independent 1188: full-source slice identity, selection, coherent relabeling')
