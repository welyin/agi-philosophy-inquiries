"""1282 independent audit. Fraction commutants and Wick moments.
No imports from author calculations; no random scans or singular preparations.
Default reads saved results. --write creates independent_results.json exclusively.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]

def eye(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]
def zeros(n,m):
    return [[F(0) for _ in range(m)] for _ in range(n)]
def tr(a):
    return sum(a[i][i] for i in range(len(a)))
def T(a):
    return [list(x) for x in zip(*a)]
def add(a,b):
    return [[x+y for x,y in zip(u,v)] for u,v in zip(a,b)]
def sub(a,b):
    return [[x-y for x,y in zip(u,v)] for u,v in zip(a,b)]
def sc(s,a):
    return [[s*x for x in row] for row in a]
def mm(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in T(b)] for row in a]
def join(a,b):
    return [u+v for u,v in zip(a,b)]
def rank(a):
    a=[row[:] for row in a]; nr=len(a); nc=len(a[0]); r=0
    for col in range(nc):
        k=next((k for k in range(r,nr) if a[k][col]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]; v=a[r][col]
        a[r]=[x/v for x in a[r]]
        for j in range(nr):
            if j!=r and a[j][col]:
                v=a[j][col]; a[j]=[x-v*y for x,y in zip(a[j],a[r])]
        r+=1
        if r==nr:break
    return r

def qm(q,p):
    # Quaternion product from the oriented multiplication table, not matrix formulas.
    out=[F(0)]*4
    cyclic={(1,2):3,(2,3):1,(3,1):2}
    for i in range(4):
        for j in range(4):
            if i==0:k,sign=j,1
            elif j==0:k,sign=i,1
            elif i==j:k,sign=0,-1
            elif (i,j) in cyclic:k,sign=cyclic[(i,j)],1
            else:k,sign=cyclic[(j,i)],-1
            out[k]+=sign*q[i]*p[j]
    return out

BASIS=eye(4); I=eye(4)
def left(q):
    return T([qm(q,b) for b in BASIS])
def right(q):
    return T([qm(b,q) for b in BASIS])
LEFT=[left(b) for b in BASIS[1:]]
RIGHT=[right(b) for b in BASIS[1:]]

def cross(b0,bv):
    b=sc(b0,I)
    for x,j in zip(bv,RIGHT):b=add(b,sc(x,j))
    return b
def edge(a,d,b0,bv,u):
    b=mm(cross(b0,bv),u)
    return join(sc(a,I),sc(-1,b))+join(sc(-1,T(b)),sc(d,I))
def covariance(s,t,u):
    aa=(s+t)/2; cc=(s-t)/2
    return join(sc(aa,I),sc(cc,u))+join(sc(cc,T(u)),sc(aa,I))

def qmean(a,cov):
    return sum(a[i][j]*cov[i][j] for i in range(len(a)) for j in range(len(a)))
def qproduct(a,b,cov):
    # Isserlis pairing, evaluated entry by entry rather than trace identities.
    aa=[(i,j,a[i][j]) for i in range(len(a)) for j in range(len(a)) if a[i][j]]
    bb=[(k,l,b[k][l]) for k in range(len(b)) for l in range(len(b)) if b[k][l]]
    return sum(v*w*(cov[i][j]*cov[k][l]+cov[i][k]*cov[j][l]+cov[i][l]*cov[j][k])
               for i,j,v in aa for k,l,w in bb)
def selector(n,indices):
    a=zeros(n,n)
    for i in indices:a[i][i]=F(1)
    return a
def num(a):
    return np.array([[float(x) for x in row] for row in a])
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def calculate():
    # Exact real commutant dimensions, plus an explicit complete basis.
    e=[]
    for i in range(4):
        for j in range(4):
            a=zeros(4,4);a[i][j]=F(1);e.append(a)
    cols=[[x for l in LEFT for row in sub(mm(l,a),mm(a,l)) for x in row] for a in e]
    comm_dim=16-rank(T(cols))
    sym=[]
    for i in range(4):
        for j in range(i,4):
            a=zeros(4,4);a[i][j]=a[j][i]=F(1);sym.append(a)
    cols=[[x for l in LEFT for row in sub(mm(l,a),mm(a,l)) for x in row] for a in sym]
    sym_dim=10-rank(T(cols))
    assert (comm_dim,sym_dim)==(4,1)
    assert rank([[x for row in a for x in row] for a in [I]+RIGHT])==4
    for i,j in enumerate(RIGHT):
        assert T(j)==sc(-1,j)
        for l in LEFT:assert mm(j,l)==mm(l,j)
        for k,h in enumerate(RIGHT):
            assert add(mm(j,h),mm(h,j))==sc(-2 if i==k else 0,I)

    quats=[[F(1),F(0),F(0),F(0)],
           [F(1,3),F(2,3),F(2,3),F(0)],
           [F(1,2)]*4]
    us=[left(q) for q in quats]
    assert all(mm(T(u),u)==I for u in us)
    a=F(21,20);d=F(19,20);b0=F(99,100);bv=[F(1,12),F(1,20),F(0)]
    beta=b0*b0+sum(x*x for x in bv)
    assert a*d>beta
    q=a+d-2*b0;r=a+d+2*b0;diff=a-d
    assert diff*diff+4*sum(x*x for x in bv)<=q*r
    s1,s2,t=F(2,3),F(7,3),F(1,5)
    eta=F(1,2000);eps=q+2*eta/(s2-s1);rmax=F(4)
    assert q<=eps and r<=rmax
    max_residual=0.;ratios=[];means=[]
    for idx,u in enumerate(us):
        m=edge(a,d,b0,bv,u)
        g1=us[(idx+1)%3];g2=us[(idx+2)%3]
        changed=mm(mm(g1,u),T(g2))
        transform=join(g1,zeros(4,4))+join(zeros(4,4),g2)
        assert mm(mm(T(transform),edge(a,d,b0,bv,changed)),transform)==m
        for s in (s1,s2):
            mean=qmean(m,covariance(s,t,u))/2
            assert mean==s*q+t*r
            means.append(str(mean))
        # Exact y/z change; its only square-root factor cancels in T^T M T.
        unscaled=join(I,I)+join(T(u),sc(-1,T(u)))
        v=sc(diff/2,I)
        for c,j in zip(bv,RIGHT):v=add(v,sc(c,j))
        target_yz=join(sc(q/2,I),v)+join(T(v),sc(r/2,I))
        assert sc(F(1,2),mm(mm(T(unscaled),m),unscaled))==target_yz
        target=edge(r/4,r/4,r/4,[F(0)]*3,u)
        err=float(np.linalg.norm(num(sub(m,target)),2))
        guard=(float(eps)+np.sqrt(float(eps*rmax)))/2
        assert err<=guard
        ratios.append(err/guard)
        actual=np.linalg.eigvalsh(num(m))
        ev=np.sqrt(float(diff*diff+4*beta))
        expected=np.sort(np.repeat([(float(a+d)-ev)/2,(float(a+d)+ev)/2],4))
        max_residual=max(max_residual,float(np.max(abs(actual-expected))))

    # Gaussian moments and complete source: actual 10-variable original W,
    # including a nonzero L12; no author zero-cross-term fixture reused.
    energy=[]
    kin=F(3,7);el=F(2,9);volume=F(7,5)
    l11,l12,l22=F(2),F(1,3),F(3,2)
    uh,uu=F(2,3),F(5,4)
    assert l11*l22>l12*l12
    kappa=F(5,4);u=us[1]
    for s in (s1,s2):
        cov=covariance(s,t,u)
        p1=selector(8,range(4));p2=selector(8,range(4,8))
        cmat=join(zeros(4,4),sc(F(1,2),u))+join(sc(F(1,2),T(u)),zeros(4,4))
        det=qproduct(p1,p2,cov)-qproduct(cmat,cmat,cov)
        assert det==12*s*t
        grad_sum=F(0)
        for l in LEFT:
            ju=mm(sc(F(1,2),l),u)
            ga=join(zeros(4,4),sc(F(1,2),ju))+join(sc(F(1,2),T(ju)),zeros(4,4))
            grad_sum+=qproduct(ga,ga,cov)
        assert grad_sum==3*s*t
        group=((1/t-1/s)/4)**2*grad_sum
        assert group==3*(s-t)**2/(16*s*t)
        invcov=join(sc((1/s+1/t)/2,I),sc((1/s-1/t)/2,u))+join(sc((1/s-1/t)/2,T(u)),sc((1/s+1/t)/2,I))
        assert mm(cov,invcov)==eye(8)
        coord=tr(invcov)/4
        assert coord==1/s+1/t
        cov10=[row+[F(0),F(0)] for row in cov]+[ [F(0)]*8+[F(1),F(0)], [F(0)]*9+[F(1)] ]
        w=F(0)
        for site in range(2):
            rp=selector(10,range(4*site,4*site+4));sp=selector(10,[8+site])
            mr=qmean(rp,cov10);ms=qmean(sp,cov10)
            w+=volume/4*(l11*(qproduct(rp,rp,cov10)-2*uh*mr+uh*uh)
                 +2*l12*(qproduct(rp,sp,cov10)-uu*mr-uh*ms+uh*uu)
                 +l22*(qproduct(sp,sp,cov10)-2*uu*ms+uu*uu))
        h=s+t
        assert w==volume/2*(l11*(6*h*h-4*uh*h+uh*uh)
                 +2*l12*(2*h-uh)*(1-uu)+l22*(3-2*uu+uu*uu))
        selected=qmean(edge(kappa,kappa,kappa,[F(0)]*3,u),cov)/2
        assert selected==4*kappa*t==1
        total=kin*(coord+F(1,2))+el*group+w+selected
        energy.append(dict(s=str(s),t=str(t),coordinate_per_kin=str(coord),
              group_per_el=str(group),gauss_gradient_c_squared=str(grad_sum),
              W_nonzero_L12=str(w),edge=str(selected),total=str(total)))
    assert energy[0]['total']!=energy[1]['total']

    # Each deletion is within its stated remaining assumptions.
    signed=edge(F(1),F(-1),F(0),[F(0)]*3,I)
    assert signed[4][4]<0
    for s in (s1,s2):assert qmean(signed,covariance(s,t,I))==0
    anis=[[F(0)]*4 for _ in range(4)]
    for j,w in enumerate((1,2,4,8)):anis[j][j]=F(w)
    ma=join(anis,sc(-1,anis))+join(sc(-1,anis),anis)
    for s in (s1,s2):assert qmean(ma,covariance(s,t,I))/2==15*t
    assert mm(mm(T(LEFT[0]),anis),LEFT[0])!=anis
    zcov=sc(t,I)
    quartic=F(1)*qmean(I,zcov)+4*F(2)*qproduct(I,I,zcov)
    assert quartic==4*t+192*t*t
    no_cross=edge(F(2),F(2),F(0),[F(0)]*3,I)
    assert qmean(no_cross,covariance(F(1),F(3),I))/2==16
    assert qmean(no_cross,covariance(F(3),F(1),I))/2==16
    assert edge(F(0),F(0),F(0),[F(0)]*3,I)==zeros(8,8)

    # Source assertions are recomputed, not copied; saved author results are provenance only.
    sources=[HERE/'proof.md',HERE/'check.py',HERE/'results.json',HERE/'NEXT.md',
             HERE.parent/'research_note_1282.md',Path(__file__).resolve(),
             RESEARCH/'archive_554_584/research_note_557.md']
    return dict(round=1282,status='PASS_WITH_SCOPE',new_axioms=0,independent_groups=1,
        exact_commutant=dict(real_dimension=comm_dim,symmetric_dimension=sym_dim),
        finite_case=dict(a=str(a),d=str(d),b0=str(b0),b_vector=[str(x) for x in bv],
            q=str(q),r=str(r),s1=str(s1),s2=str(s2),t=str(t),eta_each=str(eta),
            epsilon=str(eps),r_max=str(rmax),means=means,
            coefficient_bound=float(guard),largest_ratio=max(ratios)),
        exact_source_budgets=energy,
        exact_deletions=dict(positivity=True,left_local_covariance=True,quadratic_class=True,
            nonzero_interaction_unselected=True,fixed_difference_width=True),
        numerical_diagnostic=dict(rational_link_cases=3,eigenvalue_residual=max_residual,
            numerical_samples_are_not_proof=True),
        scope=dict(normal_Gauss_states=True,edge_cost_not_total_energy=True,
            unknown_candidate_instrument_is_effective_input=True,
            gauge_group_generated=False,autonomous_preparation_proved=False,
            potential_operator_norm_bound=False,all_time_channel_bound=False),
        source_sha256={str(p.relative_to(RESEARCH)).replace('\\','/'):sha(p) for p in sources})

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();data=calculate();dest=HERE/'independent_results.json'
    if args.write:
        with dest.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        saved=json.loads(dest.read_text('utf-8'))
        assert saved==data,'Independent result/source mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status=data['status'],mode=mode,round=1282,new_axioms=0,
        commutant=data['exact_commutant'],source_cases=len(data['exact_source_budgets'])),ensure_ascii=False))

if __name__=='__main__':main()
