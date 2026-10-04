"""698 entry: exact degree/tail contracts, not evaluation of the full integral."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_full_holonomy_reduction as old
TARGET=HERE/'polynomial_error_entry_results.json'


def degree_bound():
    colors=[(1,0),(0,1),(-1,-1)]
    weights=[(*c,p,1) for c in colors for p in (1,-1)]
    weights += [(-c[0],-c[1],0,q) for q in (-4,2) for c in colors]
    weights += [(0,0,p,-3) for p in (1,-1)]+[(0,0,0,6),(0,0,0,0)]
    assert len(weights)==16
    bare=2*np.abs(np.array(weights)).sum(axis=0)
    weighted=bare+np.array([4,4,2,0])
    assert weighted.tolist()==[20,20,18,72]
    # Only constant term is needed; M+1 roots suffice, not2M+1 interpolation.
    nodes=weighted+1
    return dict(original16_weights=weights,per_channel_exponents=[-2,2],
        source_coordinate_degree_bounds=bare.tolist(),with_Weyl_degree_bounds=weighted.tolist(),
        constant_term_grid_nodes=nodes.tolist(),grid_points=int(np.prod(nodes)),
        sphere_total_degree=32,offdiagonal_integral_field='Q(sqrt(2))',
        actual_constant_term_or_sphere_contraction_evaluated=False)


def tail(n):
    eps=7*F(4,5)**(n+1)
    assert eps<F(1,128)
    # Exact majorants for all factors used in the proof, without floating exp.
    delta=eps+eps*eps/2
    assert 63*eps/2<F(1,2) and 64*eps/2<F(1,2) and 31*delta<F(1,2)
    assert 64*eps+128*delta<256*eps
    return eps,256*eps


def approximate(a,n):
    q=np.eye(32)-a.conj().T@a/5
    power=np.eye(32,dtype=complex);polynomial=power.copy();coef=1.
    for k in range(1,n+1):
        power=power@q;coef*=float(F(2*k-1,2*k));polynomial+=coef*power
    return a@polynomial/np.sqrt(5)


def source(v,t):
    return np.linalg.det((np.eye(32)+v)/2)**2*np.linalg.det((t+v.conj()@t@v.conj().T)/2)


def run():
    degree=degree_bound()
    selected=next(n for n in range(30,400) if tail(n)[1]<=F(1,10**12))
    eps,bound=tail(selected)
    assert bound/10**6<=F(1,10**18)
    rng=np.random.default_rng(69801);angles=rng.uniform(0,2*np.pi,(3,4));hs,_=old.probe.torus(angles)
    es=rng.normal(size=(3,2,10));es/=np.linalg.norm(es,axis=2)[:,:,None]
    rows=[]
    for k,h in enumerate(hs):
        a,v=old.polar(h,True,True)
        spectrum=np.linalg.eigvalsh(a.conj().T@a)
        assert spectrum.min()>=1-1e-12 and spectrum.max()<=9+1e-12
        t=np.zeros((32,32),complex)
        for j in range(2):t[16*j:16*j+16,16*j:16*j+16]=sum(x*y for x,y in zip(es[k,j],old.b.internal.T))
        exact=source(v,t)
        for n in (32,64,selected):
            vn=approximate(a,n);error=float(np.linalg.norm(v-vn,2));source_error=float(abs(source(vn,t)-exact))
            ev,ef=tail(n)
            assert error<float(ev)+3e-14
            assert source_error<float(ef)+1e-25
            rows.append(dict(sample=k,degree=n,polar_error=error,source_error=source_error,
                polar_analytic_bound=float(ev),source_analytic_bound=float(ef)))
    deps=('research_note_697.md','joint_full_holonomy_reduction.py',
          'joint_full_holonomy_reduction_results.json','research_note_657.md','research_note_687.md')
    return dict(date='2026-10-02',entry_round=698,new_formal_round=False,
        exact_offdiagonal_contract=degree,
        static_polynomial=dict(selected_degree=selected,uniform_polar_error_upper=str(eps),
            uniform_source_error_upper=str(bound),source_error_decimal=float(bound),
            test_vector_B11_error_upper=str(bound/10**6),
            exact_coefficient_field='Q(i,sqrt(5))',floating_roundoff_not_covered_by_analytic_truncation_bound=True),
        numerical_crosschecks=rows,
        scope=dict(actual_complete_integral_not_evaluated=True,negative_integral_sign_not_proved=True,
            full_group_and_sphere_measures_retained=True,old_space_interfaces_inherited=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps},all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(entry_round=698,degree=r['static_polynomial']['selected_degree'],
        error=r['static_polynomial']['source_error_decimal'],all_checks_passed=True)))
