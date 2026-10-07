"""Exact finite certificate: a five-qubit code in ten original CAR modes.

Integer symplectic Pauli arithmetic checks every Majorana word up to degree
four and all 211^2 error products. This checks the encoding and source
compression; the original BFF and continuum state mapping are analytic.
"""
from pathlib import Path
import argparse,itertools,json
HERE=Path(__file__).resolve().parent;TARGET=HERE/'majorana_code_source_bridge_results.json'
I=(0,0,0)
def mul(a,b):
    x,z,p=a;y,w,q=b
    return x^y,z^w,(p+q+2*(z&y).bit_count())%4
def adj(a):
    x,z,p=a;return x,z,(-p+2*(x&z).bit_count())%4
def commute(a,b):return ((a[0]&b[1]).bit_count()+(a[1]&b[0]).bit_count())%2==0
def product(seq):
    a=I
    for b in seq:a=mul(a,b)
    return a
def code_data():
    gamma=[]
    for i in range(10):
        gamma.extend([(1<<i,(1<<i)-1,0),(1<<i,(1<<(i+1))-1,1)])
    sx=[((1<<(2*b))|(1<<(2*b+1)),0,0) for b in range(5)]
    sz=[(0,1<<(2*b),0) for b in range(5)]
    sy=[mul((0,0,3),mul(sz[b],sx[b])) for b in range(5)]
    st=[(0,(1<<(2*b))|(1<<(2*b+1)),2) for b in range(5)]
    for word in ('XZZXI','IXZZX','XIXZZ','ZXIXZ'):
        st.append(product([{'I':I,'X':sx[b],'Y':sy[b],'Z':sz[b]}[s] for b,s in enumerate(word)]))
    group={I}
    for s in st:group|={mul(a,s) for a in list(group)}
    assert len(group)==512 and (0,0,2) not in group
    assert all(commute(a,b) for a in st for b in st)
    lookup={(a[0],a[1]):a[2] for a in group}
    def syndrome(a):return sum((not commute(a,s))<<i for i,s in enumerate(st))
    def compressed(a):
        if syndrome(a):return 'zero',0
        key=a[:2]
        if key in lookup:return 'scalar',(a[2]-lookup[key])%4
        return 'logical',None
    return gamma,st,group,syndrome,compressed,sx,sz
def run():
    gamma,st,group,syndrome,comp,sx,sz=code_data()
    degrees=[]
    for degree in range(5):
        count={'zero':0,'scalar':0,'logical':0}
        for ids in itertools.combinations(range(20),degree):
            kind,_=comp(product([gamma[i] for i in ids]));count[kind]+=1
        assert count['logical']==0
        degrees.append(dict(degree=degree,**count))
    # Hermitian error basis I, gamma_i, i gamma_i gamma_j.
    errors=[I,*gamma]+[mul((0,0,1),mul(gamma[i],gamma[j])) for i in range(20) for j in range(i+1,20)]
    assert len(errors)==211 and all(adj(e)==e for e in errors)
    scalars=zeros=0;representatives={}
    for e in errors:representatives.setdefault(syndrome(e),e)
    for e in errors:
        r=representatives[syndrome(e)]
        assert comp(mul(adj(r),e))[0]=='scalar'
        for f in errors:
            kind,_=comp(mul(adj(e),f));assert kind!='logical'
            scalars+=kind=='scalar';zeros+=kind=='zero'
    # The exact dual-rail code has vanishing one-/two-Majorana moments.
    assert degrees[1]['zero']==20 and degrees[2]['zero']==190
    logical=product(sx);assert comp(logical)[0]=='logical'
    # First nontrivial logical Majorana word: no odd word preserves parity.
    shortest=None
    for degree in (5,6):
        for ids in itertools.combinations(range(20),degree):
            if comp(product([gamma[i] for i in ids]))[0]=='logical':shortest=list(ids);break
        if shortest:break
    assert shortest is not None and len(shortest)==6
    # Fixed block parity is a non-Gaussian four-point moment, not Wick data.
    four=product(gamma[:4]);kind,phase=comp(four)
    assert kind=='scalar' and phase==0
    # Removing the outer code leaves a one-block logical quadratic source.
    bare_z=sz[0];assert all(commute(bare_z,s) for s in st[:5])
    assert bare_z[:2] not in {product(sub)[:2] for r in range(6) for sub in itertools.combinations(st[:5],r)}
    return dict(round=829,all_checks_passed=True,exact_integer_arithmetic=True,
        CAR_modes=10,Majoranas=20,independent_stabilizers=9,code_dimension=2,
        fixed_total_particle_number=5,detected_word_counts=degrees,
        hermitian_error_basis_size=211,ordered_KL_products_checked=211**2,
        KL_zero_products=zeros,KL_scalar_products=scalars,nontrivial_error_syndromes=len(representatives)-1,
        syndrome_recovery_exact_on_entire_error_span=True,
        all_190_quadratic_Majorana_sources_scalar=True,
        self_dual_two_point_on_code_modes='I_20 / 2 for every logical density matrix',
        scalar_four_Majorana_expectation=1,Gaussian_prediction_from_two_point=0,
        shortest_logical_Majorana_word=shortest,Majorana_distance=6,
        dual_rail_only_quadratic_source_counterexample=True,
        original_BFF_and_source_mapping_checked_numerically=False,
        original_autonomous_encoder_decoder_proven=False,
        scope='Exact code, finite error-span recovery, and all quadratic/fourth-degree source compressions. Original continuum mapping and perturbative scope are analytic; this is not a new code discovery.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
