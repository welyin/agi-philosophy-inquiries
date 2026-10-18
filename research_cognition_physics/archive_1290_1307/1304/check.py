# 1304: native identity under a finite polynomial mass menu.
# Default/--check read only. --write exclusively creates results.json.
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'
PAIRS=((1,0),(3,0),(0,1),(0,3),(2,1),(1,2))

def moments(xs):
    return tuple(sum(x**r for x in xs) for r in (1,2,3))

def coefficients(q,f):
    return {f'g{r}F{s}':sum(x**r*y**s for x,y in zip(q,f)) for r,s in PAIRS}

def components(vertices,edges):
    vertices=set(vertices)
    adj={v:set() for v in vertices}
    for x,y in edges:
        adj[x].add(y); adj[y].add(x)
    groups=[]
    while vertices:
        start=min(vertices)
        seen={start}; todo=[start]
        while todo:
            for nxt in adj[todo.pop()]-seen:
                seen.add(nxt); todo.append(nxt)
        groups.append(sorted(seen))
        vertices-=seen
    return sorted(groups)

def coordinate_graph(a,b,maxpower):
    vs=[('A',x) for x in set(a)]+[('B',x) for x in set(b)]
    es=[(('A',x),('B',y)) for x in set(a) for y in set(b)
        if 0<abs(x-y)<=maxpower]
    return components(vs,es)

def finite_box():
    buckets=defaultdict(list)
    for a in combinations_with_replacement(range(15),4):
        buckets[moments(a)].append(a)
    assert sum(map(len,buckets.values()))==3060
    rows=[]
    for bag in buckets.values():
        for a,b in combinations(bag,2):
            if min(a+b)!=0 or set(a)&set(b): continue
            groups=coordinate_graph(a,b,2)
            cuts=[t for t in range(14) if min(a+b)<=t<max(a+b)
                  and not any((x<=t<y or y<=t<x) and abs(x-y)<=2
                              for x in a for y in b)]
            matching=any(all(0<abs(a[i]-b[j])<=2 for i,j in enumerate(p))
                         for p in permutations(range(4)))
            assert len(groups)>1 and cuts
            rows.append({'A':list(a),'B':list(b),'moments':list(moments(a)),
                         'distinct_role_components':len(groups),
                         'groups':groups,'edge_free_cuts':cuts,
                         'has_full_matching_at_power_2':matching})
    rows.sort(key=lambda z:(z['A'],z['B']))
    expected=[
        ([0,3,4,7],[1,1,6,6],2),
        ([0,4,7,11],[1,2,9,10],2),
        ([0,5,5,10],[1,2,8,9],3),
        ([0,6,7,13],[1,3,10,12],6),
        ([0,6,8,14],[2,2,12,12],4)]
    assert [(r['A'],r['B'],r['distinct_role_components']) for r in rows]==expected
    assert all(t in r['edge_free_cuts'] for r,t in zip(rows,[3,4,2,1,2]))
    return {'max_power':2,'proved_span_bound':14,'multisets_checked':3060,
            'normalized_disjoint_equal_moment_pairs':rows,
            'connected_candidates':0,'box_bound_is_analytic_not_arbitrary':True}

def determinant(a):
    a=[[F(x) for x in row] for row in a]
    ans=F(1)
    for j in range(len(a)):
        k=next((k for k in range(j,len(a)) if a[k][j]),None)
        if k is None: return F(0)
        if k!=j: a[j],a[k]=a[k],a[j]; ans=-ans
        p=a[j][j]; ans*=p
        for k in range(j+1,len(a)):
            ratio=a[k][j]/p
            for col in range(j+1,len(a)): a[k][col]-=ratio*a[j][col]
            a[k][j]=0
    return ans

def role_check(q,f,M):
    assert all(q) and not any(x+y==0 for i,x in enumerate(q) for y in q[i+1:])
    assert all(f[i]==f[j] for i in range(len(q)) for j in range(len(q)) if q[i]==q[j])
    assert not any(coefficients(q,f).values())
    assert all((f[i]+f[j])*M[i][j]==0 for i in range(len(q)) for j in range(len(q)))
    edges=[(i,j) for i in range(len(q)) for j in range(i+1,len(q)) if M[i][j]]
    fields=components(range(len(q)),edges)
    roles=components(q,[(q[i],q[j]) for i,j in edges])
    return {'q':list(q),'F':list(f),'six_coefficients':coefficients(q,f),
            'distinct_role_components':len(roles),'role_groups':roles,
            'field_support_components':len(fields),'field_groups':fields,
            'F_contract_exact':True,'same_charge_same_F':True}

def block_matrix(D):
    n=len(D)
    return [[F(0) if (i<n)==(j<n) else
             D[i][j-n] if i<n else D[j][i-n] for j in range(2*n)]
            for i in range(2*n)]

def main_example():
    a,b=(1,6,6,11),(2,3,9,10)
    q=a+tuple(-x for x in b); f=(1,)*4+(-1,)*4
    assert moments(a)==moments(b)==(24,194,1764)
    K=[[F(int((i,j) in [(0,1),(1,2),(3,2)])) for j in range(4)] for i in range(4)]
    D=[[F(int(i==j))+K[i][j]/4 for j in range(4)] for i in range(4)]
    assert determinant(D)==1
    kt=[[sum(K[k][i]*K[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    assert kt==[[F(v) if i==j else F(0) for j in range(4)] for i,v in enumerate([0,1,2,0])]
    assert F(2)<F(3,2)**2  # sqrt(2)<3/2, so sigma_min>5/8.
    M=block_matrix(D)
    info=role_check(q,f,M)
    assert info['distinct_role_components']==info['field_support_components']==1
    assert determinant(M)==1
    eps=F(1,8); scale=eps**3/2
    terms=[]
    for i in range(4):
        for j in range(4):
            if not D[i][j]: continue
            delta=a[i]-b[j]
            u,w=max(-delta,0),max(delta,0); p=u+w
            c=D[i][j]*eps**(3-p)/2
            assert 1<=p<=3 and 0<c<=F(1,2)
            assert delta+u-w==0 and f[i]+f[4+j]==0
            assert c*eps**p==scale*D[i][j]
            terms.append({'i':i,'j':4+j,'charge_sum':delta,'u':u,'w':w,
                          'operator_dimension':3+p,'c':str(c),'D':str(D[i][j])})
    assert len(terms)==7 and max(t['operator_dimension'] for t in terms)==6
    df=np.array([[float(x) for x in row] for row in D])
    sv=np.linalg.svd(df,compute_uv=False)
    assert sv[-1]>5/8

    def mass(s,ds=False,db=False):
        out=np.zeros((8,8),dtype=complex)
        for t in terms:
            u,w=t['u'],t['w']; c=float(F(t['c']))
            if ds:
                val=0j if u==0 else c*u*s**(u-1)*s.conjugate()**w
            elif db:
                val=0j if w==0 else c*w*s**u*s.conjugate()**(w-1)
            else: val=c*s**u*s.conjugate()**w
            i,j=t['i'],t['j']; out[i,j]=out[j,i]=val
        return out
    s=.13+.07j; theta=.47; eta=.29
    R=np.diag(np.exp(1j*np.array(q)*theta))
    RF=np.diag(np.exp(1j*np.array(f)*eta))
    current=mass(s); ds=mass(s,ds=True); db=mass(s,db=True)
    gauge=np.max(np.abs(R.T@mass(s*np.exp(1j*theta))@R-current))
    fid=np.max(np.abs(RF.T@current@RF-current))
    Q=np.diag(q)
    ward=np.max(np.abs(Q.T@current+current@Q+s*ds-s.conjugate()*db))
    h=1e-6
    dx=(mass(s+h)-mass(s-h))/(2*h)
    dy=(mass(s+1j*h)-mass(s-1j*h))/(2*h)
    derivative=max(np.max(np.abs(dx-ds-db)),np.max(np.abs(dy-1j*(ds-db))))
    atv=mass(complex(float(eps),0))
    mfloat=np.array([[float(x) for x in row] for row in M])*float(scale)
    assert max(gauge,fid,ward)<1e-12 and derivative<1e-8
    assert np.max(np.abs(atv-mfloat))<1e-15
    assert np.all(np.isfinite(mass(0j))) and np.all(np.isfinite(mass(0j,ds=True)))
    info.update({'A':list(a),'B':list(b),'moments':list(moments(a)),
        'D':[[str(x) for x in row] for row in D],'det_D':str(determinant(D)),
        'det_dimensionless_M':str(determinant(M)),'KtK_diagonal':[0,1,2,0],
        'sigma_strict_rational_lower':'5/8','diagnostic_singular_values_D':sv.tolist(),
        'fixed_eps_diagnostic':str(eps),'fixed_Lambda_diagnostic':1,
        'mass_scale_Lambda_eps3_over2':str(scale),'terms':terms,
        'all_monomial_gauge_weights_exact_zero':True,'all_terms_F_neutral':True,
        'regular_at_S_zero':True,'fixed_coefficients_not_retuned_with_S':True,
        'gauge_covariance_residual':float(gauge),'F_covariance_residual':float(fid),
        'polynomial_source_Ward_residual':float(ward),
        'scalar_derivative_diagnostic_residual':float(derivative),
        'mass_at_v_residual':float(np.max(np.abs(atv-mfloat)))})
    return info

def deletion_controls():
    a,b=(1,4,5,8),(2,2,7,7)
    D=[[F(int(i==j)) for j in range(4)] for i in range(4)]
    c5=role_check(a+tuple(-x for x in b),(1,)*4+(-1,)*4,block_matrix(D))
    assert moments(a)==moments(b)==(18,106,702)
    powers=[abs(x-y) for x,y in zip(a,b)]
    assert powers==[1,2,2,1] and c5['distinct_role_components']==2
    c5.update({'max_dimension':5,'matching_powers':powers,
               'extra_uniform_abs_F_retained':True,'only_sufficiency_claimed':True})
    q=(2,-1,10,-9,4,-5,8,-9); f=(1,-1,3,-3,1,-1,3,-3)
    P=[[F(int(i//2==j//2 and i!=j)) for j in range(8)] for i in range(8)]
    c4=role_check(q,f,P)
    assert determinant(P)==1 and c4['distinct_role_components']==3
    assert c4['field_support_components']==4
    sums=[q[i]+q[i+1] for i in range(0,8,2)]
    assert sums==[1,1,-1,-1]
    c4.update({'max_dimension':4,'pair_charge_sums':sums,'det_P':str(determinant(P)),
               'uniform_abs_F_does_not_follow_without_connectivity':True,
               'same_charge_minus9_identity':-3})
    return {'dim5_with_uniform_F':c5,'dim4_after_connectivity_alone_removed':c4}

def payload():
    sources=['proof.md','NEXT.md','../research_note_1304.md',
        '../_shared/round1304_admission.md','../_shared/round1304_admission_clarification.md',
        '../1301/proof.md','../1302/proof.md','../1302/acceptance.json',
        '../../archive_1009_1043/research_note_1022.md',
        '../_shared/comparison_scope.md']
    result={'status':'PASS_AUTHOR_NATIVE_IDENTITY_MENU_THRESHOLD',
        'finite_complete_box_certificate':finite_box(),'connected_dim6_model':main_example(),
        'deletion_controls':deletion_controls(),'new_axioms':0,
        'PTE_mapping_and_N8_lower_bound_are_prior_results':True,
        'quantified_main_claim':'connected role class: N>=8, N8 minimal maximum operator dimension=6',
        'full_quantum_or_operational_model_claimed':False,'global_UV_completion_claimed':False,
        'exactly_real_world_particles_claimed':False,
        'source_bindings':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
        'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    # Ensure the comparison uses precisely the admitted source records.
    assert result['source_bindings']['../_shared/round1304_admission.md']=='b51d300bcd4bf881d447a31e2b2847bcce239ded86318a36bad1a8563937e169'
    assert result['source_bindings']['../_shared/round1304_admission_clarification.md']=='ece63490c57fc1a4fdb417e791485216b2eb57109645cca6b6b32d448a98a7b3'
    return json.loads(json.dumps(result,ensure_ascii=False))

def main():
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',action='store_true')
    modes.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=payload()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2); stream.write('\n')
        print('CREATED',result['status'])
    else:
        expected=json.loads(OUT.read_text(encoding='utf-8'))
        if result!=expected: raise SystemExit('MISMATCH: sources or results changed')
        print('READ_ONLY_MATCH',result['status'])
    print('box_multisets',result['finite_complete_box_certificate']['multisets_checked'])
    print('role_components',[r['distinct_role_components'] for r in result['finite_complete_box_certificate']['normalized_disjoint_equal_moment_pairs']])
    print('dim6_min_singular',min(result['connected_dim6_model']['diagnostic_singular_values_D']))

if __name__=='__main__': main()
