"""1301 exact role-flow and matrix certificates.
The all-charge and all-N theorem is analytic; samples only cross-check it.
Default/--check read-only. --write exclusively creates results.json.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import math
import numpy as np
HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'


def bareiss(matrix):
    a=[list(map(int,row)) for row in matrix];n=len(a);prev=1;sgn=1
    for k in range(n-1):
        if a[k][k]==0:
            row=next((j for j in range(k+1,n) if a[j][k]),None)
            if row is None:return 0
            a[k],a[row]=a[row],a[k];sgn=-sgn
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                assert numerator%prev==0
                a[i][j]=numerator//prev
        for i in range(k+1,n):a[i][k]=0
        prev=pivot
    return sgn*a[-1][-1]


def role_components(charges,matrix):
    roles=sorted(set(charges));adj={q:set() for q in roles}
    for i,q in enumerate(charges):
        for j,r in enumerate(charges):
            if q!=r and matrix[i][j]:adj[q].add(r)
    remaining=set(roles);components=[]
    while remaining:
        seen={min(remaining)}
        while True:
            new=seen|set().union(*(adj[q] for q in seen))
            if new==seen:break
            seen=new
        components.append(sorted(seen));remaining-=seen
    return components


def n8_flows():
    def all_flows(total,length,prefix=()):
        if length==0:
            if total==0 and prefix[-1]>0:yield prefix
            return
        low=1 if not prefix or prefix[-1]==0 or length==1 else 0
        for value in range(low,total+1):
            yield from all_flows(total-value,length-1,prefix+(value,))
    flows=[]
    for L in range(2,8):
        for x in all_flows(4,L):
            if sum((-1)**i*v for i,v in enumerate(x)):continue
            A=sum((-1)**i*v*i for i,v in enumerate(x))
            B=sum((-1)**i*v*i*(i+1) for i,v in enumerate(x))
            root=None if A==0 else F(-B,2*A)
            assert (A==0 and B!=0) or root<1
            flows.append({'flow':list(x),'coefficient_of_k':2*A,'constant':B,
                          'root':None if root is None else str(root),
                          'no_k_ge_1':True})
    expected={(2,2),(1,2,1),(1,0,1,2),(1,1,1,1),(2,1,0,1),
              (1,1,0,1,1),(1,0,1,1,0,1)}
    assert {tuple(r['flow']) for r in flows}==expected
    return {'all_seven_shapes':flows,'maximum_length_from_total_flow':7,
            'integer_charge_scan_used_for_lower_bound':False}


def path_template(q,n,x):
    charges=[z for z,c in zip(q,n) for _ in range(c)]
    indices=[];start=0
    for count in n:indices.append(list(range(start,start+count)));start+=count
    unused=[a[:] for a in indices];pairs=[]
    for i,count in enumerate(x):
        for _ in range(count):pairs.append((unused[i].pop(0),unused[i+1].pop(0)))
    assert not any(unused)
    return charges,pairs


def matrix_certificate(charges,pairs,h,denominator=16,expect_connected=True):
    N=len(charges);M0=np.zeros((N,N),dtype=int)
    for i,j in pairs:M0[i,j]=M0[j,i]=1
    assert np.array_equal(M0@M0,np.eye(N,dtype=int))
    assert sum(charges)==sum(q**3 for q in charges)==0
    assert all(q!=0 and -q not in charges for q in charges)
    assert math.gcd(abs(h),*map(abs,charges))==1
    A=np.array([[int(abs(q+r)==abs(h)) for r in charges] for q in charges],dtype=int)
    assert np.all((M0==0)|(A==1))
    num=denominator*M0+A
    degree=int(np.max(np.sum(A,axis=1)));lower=1-F(degree,denominator)
    assert lower>0
    det=bareiss(num.tolist());assert det!=0
    components=role_components(charges,num.tolist())
    assert (len(components)==1)==expect_connected
    # Actual coefficients, not just allowed zeros: every permitted entry nonzero.
    assert np.all((A==0)|(num!=0))
    M=num/denominator
    singular=np.linalg.svd(M,compute_uv=False)
    assert min(singular)>float(lower)-1e-12
    # M(S)=S*Y_minus + S.conjugate()*Y_plus.
    sums=np.array(charges)[:,None]+np.array(charges)[None,:]
    minus=M*(sums==-h);plus=M*(sums==h)
    theta=.237;S=np.exp(1j*h*theta);R=np.diag(np.exp(1j*np.array(charges)*theta))
    Mnew=S*minus+S.conjugate()*plus
    assert np.linalg.norm(R.T@Mnew@R-M)<2e-12
    phases=np.exp(1j*np.arange(N)*.173);D=np.diag(phases)
    assert np.max(np.abs(np.linalg.svd(D.T@M@D,compute_uv=False)-singular))<2e-12
    return {'N':N,'h':h,'charges':charges,'linear_anomaly':0,'cubic_anomaly':0,
            'matching_pairs':[list(p) for p in pairs],
            'max_allowed_row_sum':degree,'epsilon':str(F(1,denominator)),
            'mass_singular_lower_bound':str(lower),'integer_numerator_determinant':str(det),
            'role_components':components,'all_permitted_edges_actually_nonzero':True,
            'reference_stabilizer_order':abs(h),
            'S_only_linear_interpolator_available':[q%abs(h)==0 for q in charges],
            'gauge_covariance_and_rephasing_residual_below_2e_minus_12':True,
            'full_flavor_irreducibility_claimed':False}


def templates_and_controls():
    definitions={
      'N12':([-2,3,-4,5,-6,7],[2,4,2,1,2,1],[2,2,0,1,1]),
      'N16':([-1,2,-3,4,-5],[1,4,6,4,1],[1,3,3,1]),
      'N20':([-1,2,-3,4],[3,8,7,2],[3,5,2]),
      'h4_N8':([9,-5,1,3,-7],[1,1,1,3,2],[1,0,1,2])}
    out={}
    for name,(q,n,x) in definitions.items():
        charges,pairs=path_template(q,n,x)
        out[name]=matrix_certificate(charges,pairs,4 if name=='h4_N8' else 1,
                                     denominator=16 if name in ['N12','h4_N8'] else 4*sum(n))
        out[name]['matching_flow']=x
        out[name]['zero_matching_bridges_activated']=[i for i,v in enumerate(x) if not v]
    charges=[2,-1,10,-9,4,-5,8,-9];pairs=[(0,1),(2,3),(4,5),(6,7)]
    out['disconnected_h1_N8']=matrix_certificate(charges,pairs,1,expect_connected=False)
    assert len(out['disconnected_h1_N8']['role_components'])==3
    return out


def all_n_construction():
    # Analytical coverage is p=3a+4b+5c. These are finite implementation checks.
    defs={12:([2,-3,4,-5,6,-7],[2,4,2,1,2,1],[2,2,0,1,1]),
          16:([-1,2,-3,4,-5],[1,4,6,4,1],[1,3,3,1]),
          20:([-1,2,-3,4],[3,8,7,2],[3,5,2])}
    records=[]
    for p in range(3,13):
        b=int(p%3==1);c=int(p%3==2);a=(p-4*b-5*c)//3
        assert min(a,b,c)>=0 and 3*a+4*b+5*c==p
        charges=[];pairs=[]
        for size,count in [(12,a),(16,b),(20,c)]:
            qq,pp=path_template(*defs[size])
            for _ in range(count):
                offset=len(charges);charges.extend(qq);pairs.extend((i+offset,j+offset) for i,j in pp)
        cert=matrix_certificate(charges,pairs,1,denominator=4*len(charges))
        assert F(cert['mass_singular_lower_bound'])>=F(3,4)
        records.append({'N':4*p,'template_counts':[a,b,c],
                        'strict_mass_bound':cert['mass_singular_lower_bound'],
                        'connected_actual_role_graph':True,'anomalies_exactly_zero':True})
    return {'finite_cross_checks':records,
            'all_p_ge_3_proof':'p mod 3: use 3a, 4+3a, or 5+3a',
            'finite_scan_not_the_all_N_proof':True}


def payload():
    sources=['proof.md','../_shared/round1301_admission.md',
             '../../archive_1009_1043/research_note_1022.md',
             '../_admission/positive_task_screen_after1296.md']
    return {'status':'PASS_AUTHOR_STRUCTURE_CERTIFICATES','N8':n8_flows(),
            'templates_and_deletions':templates_and_controls(),'all_N':all_n_construction(),
            'new_axioms':0,'actual_instrument_O_claimed':False,
            'full_SM_or_generations_selected':False,'full_global_anomalies_claimed':False,
            'source_bindings':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
            'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();out=payload()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
        print('CREATED',out['status'])
    else:
        expected=json.loads(OUT.read_text(encoding='utf-8'))
        if out!=expected:raise SystemExit('MISMATCH: author diagnostics or source binding changed')
        print('READ_ONLY_MATCH',out['status'])
if __name__=='__main__':main()
