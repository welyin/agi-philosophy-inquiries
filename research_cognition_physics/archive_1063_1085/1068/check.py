"""Finite formula witnesses for round 1068; default is read-only."""
import argparse
import json
from pathlib import Path
import numpy as np

def density(v): return np.outer(v,v.conj())
def omega(b):
    b=np.array(b,dtype=complex);b/=np.linalg.norm(b)
    return b,density(b.reshape(-1))
def channel(u,rho): return u@rho@u.conj().T
def overlap(pure,rho): return float(np.trace(pure@rho).real)
def phase(n,r,t): return np.diag([np.exp(1j*t)]*r+[1]*(n-r))
def swap(n):
    s=np.zeros((n*n,n*n))
    for i in range(n):
        for j in range(n):s[j*n+i,i*n+j]=1
    return s

def run():
    residuals=[]
    def close(a,b):
        e=float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
        residuals.append(e);assert e<2e-12,e
    b2,w2=omega([[0,1],[2,0]])
    b4,w4=omega([[0,0,1,0],[0,0,0,2],[0,3,0,0],[4,0,0,0]])
    positives=[]
    for n,r,b,w in [(2,1,b2,w2),(4,2,b4,w4)]:
        cases=0
        for t in np.linspace(0,np.pi,9):
            u=phase(n,r,t);close(channel(np.kron(u,u),w),w)
            for start,end in [(0,r),(r,n)]:
                for i in range(start,end):
                    for j in range(start,end):
                        e=np.zeros((n,n),complex);e[i,j]=1
                        close(channel(u,e),e);cases+=1
        assert np.linalg.matrix_rank(b)==n
        positives.append(dict(n=n,rank_P=r,schmidt_rank=n,phase_settings=9,branch_units=cases))
    ranks=[];z=(3+4j)/5
    for n in range(2,9):
        for r in range(1,n):
            pp=np.diag([1]*r+[0]*(n-r)).astype(complex);qq=np.eye(n)-pp
            cross=np.zeros((n,n),complex)
            for i in range(min(r,n-r)):cross[i,r+i]=1;cross[r+i,i]=2
            u=np.diag([z]*r+[1]*(n-r));rr=[]
            for m,c in [(qq,1),(cross,z),(pp,z*z)]:
                close(u@m@u.T,c*m);rr.append(int(np.linalg.matrix_rank(m)))
            assert rr==[n-r,2*min(r,n-r),r]
            assert (max(rr)==n)==(n==2*r)
            ranks.append(dict(n=n,r=r,ranks_QQ_cross_PP=rr))
    ux=np.array([[0,1j],[1j,0]])
    partial=overlap(w2,channel(np.kron(ux,ux),w2));close(partial,16/25)
    b3=np.zeros((3,3),complex);b3[0,1]=1;b3[1,0]=-1
    b3,w3=omega(b3);u3=phase(3,1,np.pi/2)
    close(channel(np.kron(u3,u3),w3),w3);assert np.linalg.matrix_rank(b3)==2
    anti=(np.eye(9)-swap(3))/6
    close(channel(np.kron(u3,u3),anti),anti)
    close(np.einsum('ijkj->ik',anti.reshape(3,3,3,3)),np.eye(3)/3)
    assert np.linalg.matrix_rank(anti)==3
    _,max3=omega(np.eye(3));flip=phase(3,1,np.pi)
    close(channel(np.kron(flip,flip),max3),max3)
    midpoint=overlap(max3,channel(np.kron(u3,u3),max3));close(midpoint,1/9)
    j=np.array([[0,1],[-1,0]]);_,sym4=omega(np.kron(np.eye(2),j))
    weak=np.diag([1j,-1j,1,1]);p=np.diag([1,0,0,0])
    close(weak.conj().T@p@weak,p);close(channel(np.kron(weak,weak),sym4),sym4)
    q=density(np.array([0,1,1,0])/np.sqrt(2))
    branch=overlap(q,channel(weak,q));close(branch,.5)
    close(channel(np.kron(u3,u3.conj()),max3),max3)
    _,bell=omega(np.eye(2));zz=np.diag([1,-1,-1,1])
    correlated=(bell+channel(zz,bell))/2;close(correlated,bell)
    independent=np.diag(np.diag(bell))
    independent_overlap=overlap(bell,independent);close(independent_overlap,.5)
    classical=np.diag([.5,0,0,.5]);close(np.diag(np.diag(classical)),classical)
    u=phase(2,1,.37);ks=[u/np.sqrt(3),u*np.sqrt(2/3)]
    close(sum(k.conj().T@k for k in ks),np.eye(2))
    close(sum(channel(np.kron(a,b),w2) for a in ks for b in ks),w2)
    kraus_rank=int(np.linalg.matrix_rank(np.stack([k.reshape(-1) for k in ks])))
    assert kraus_rank==1
    return dict(round=1068,passed=True,
        scope='Finite matrix witnesses; general theorem and physical premises are in proof.md.',
        positive_contracts=positives,charge_rank_cases=ranks,
        partial_relation_common_SU2_overlap=partial,embedded_n3_schmidt_rank=2,
        mixed_n3_relation_rank=3,discrete_n3_midpoint_overlap=midpoint,
        only_P_statistics_n4_complement_overlap=branch,conjugate_role_n3_invariant=True,
        correlated_noise_bell_overlap=overlap(bell,correlated),
        independent_noise_bell_overlap=independent_overlap,kraus_span_rank=kraus_rank,
        max_formula_residual=max(residuals))

def compare(a,b,path='result'):
    if isinstance(a,dict):
        assert set(a)==set(b),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+str(i))
    elif isinstance(a,float):assert np.isclose(a,b,atol=1e-11,rtol=1e-10),(path,a,b)
    else:assert a==b,(path,a,b)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run();p=Path(__file__).with_name('results.json')
    if args.write:
        with p.open('x',encoding='utf8',newline='\n') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(p.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1068,passed=True,write=args.write,max_formula_residual=out['max_formula_residual'])))
if __name__=='__main__':main()
