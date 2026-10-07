"""994: leading chemical-equilibrium interface of the declared B993 field content.
Exact rational constraints + illustrative positive linear-response relaxation.
No baryogenesis source, physical reaction rates or cosmological abundance is computed.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'charge_history_results.json'
NAMES=['q','u','d','l_e','l_mu','l_tau','e_e','e_mu','e_tau','H']
W=list(map(F,[18,9,9,2,2,2,1,1,1,4]))
Y=[F(1,6),F(2,3),F(-1,3),F(-1,2),F(-1,2),F(-1,2),F(-1),F(-1),F(-1),F(1,2)]
B=[F(1,3)]*3+[F(0)]*7
LF=[[F(int(k in (3+i,6+i))) for k in range(10)] for i in range(3)]
DELTA=[[B[k]/3-LF[i][k] for k in range(10)] for i in range(3)]

def dot(x,y):return sum(a*b for a,b in zip(x,y))
def exact(x):return dict(exact=str(x),value=float(x))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def row(entries):return [F(entries.get(k,0)) for k in range(10)]
def constraints(active=()):
    # Chemical affinities, physical right-handed convention, mu_H has Y=+1/2.
    r=[row({0:-1,1:1,9:-1}),row({0:-1,2:1,9:1})]
    r += [row({3+i:-1,6+i:1,9:1}) for i in range(3)]
    r += [row({0:9,3:1,4:1,5:1})]
    for i,j in active:
        rr=row({9:2});rr[3+i]+=1;rr[3+j]+=1;r.append(rr)
    assert all(dot(rr,Y)==0 for rr in r)
    return r
def rref(matrix):
    a=[[F(x) for x in rr] for rr in matrix];pivots=[];i=0
    if not a:return a,pivots
    for j in range(len(a[0])):
        p=next((k for k in range(i,len(a)) if a[k][j]),None)
        if p is None:continue
        a[i],a[p]=a[p],a[i];scale=a[i][j];a[i]=[x/scale for x in a[i]]
        for k in range(len(a)):
            if k!=i:
                z=a[k][j];a[k]=[x-z*y for x,y in zip(a[k],a[i])]
        pivots.append(j);i+=1
        if i==len(a):break
    return a,pivots
def nullspace(matrix):
    a,piv=rref(matrix);free=[i for i in range(10) if i not in piv];basis=[]
    for f in free:
        v=[F(0)]*10;v[f]=1
        for i,j in enumerate(piv):v[j]=-a[i][f]
        basis.append(v)
    return basis
def solve(matrix,rhs):
    a,piv=rref([rr+[v] for rr,v in zip(matrix,rhs)])
    assert 10 not in piv and len([j for j in piv if j<10])==10
    ans=[F(0)]*10
    for i,j in enumerate(piv):ans[j]=a[i][-1]
    assert all(dot(rr,ans)==v for rr,v in zip(matrix,rhs))
    return ans
def charges(mu):
    dens=[x*w for x,w in zip(mu,W)]
    b=dot(B,dens);leptons=[dot(ll,dens) for ll in LF]
    return dict(B=exact(b),L=exact(sum(leptons)),B_minus_L=exact(b-sum(leptons)),
        hypercharge=exact(dot(Y,dens)),flavor_delta=[exact(dot(d,dens)) for d in DELTA])
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)

def run():
    yn=[w*y for w,y in zip(W,Y)]
    base=constraints();basis=nullspace(base+[yn]);assert len(basis)==3
    for v in basis:
        c=charges(v)
        assert F(c['B']['exact'])==F(28,79)*F(c['B_minus_L']['exact'])
    delta_rows=[[d*w for d,w in zip(di,W)] for di in DELTA]
    full=base+[yn]+delta_rows
    # Chosen charge amplitude normalizes a susceptibility calculation; not observed eta_B.
    mu=solve(full,[F(0)]*7+[F(1),F(0),F(0)])
    partial=constraints([(0,0)])+[yn]
    mu_partial=solve(partial+delta_rows[1:],[F(0)]*8+[F(1),F(0)])
    assert charges(mu_partial)['B']['exact']=='84/257'
    branches=[]
    for name,active in [('slow_DeltaL2',()),('fast_ee_only',((0,0),)),
                         ('fast_ee_mm_tt',((0,0),(1,1),(2,2)))]:
        rr=constraints(active);bs=nullspace(rr+[yn])
        branches.append(dict(name=name,fast_weinberg_pairs=[list(a) for a in active],
            neutral_equilibrium_dimension=len(bs),
            basis_chemical_potentials=[[str(x) for x in v] for v in bs]))
    assert [x['neutral_equilibrium_dimension'] for x in branches]==[3,2,0]

    # Linear-response illustration, positive normalized rates = 1.
    # x=W*mu; dx/ds=-R^T R mu, s is a dimensionless relaxation parameter.
    # The hypercharge neutrality row is NOT a physical reaction.
    wf=np.array(W,dtype=float);sw=np.sqrt(wf)
    r0=np.array(base,dtype=float)
    charge_mu=np.array(mu,dtype=float)*.001
    # Disequilibrium perturbation carries none of Y or the three Delta_i.
    start=charge_mu+.001*r0[0]/wf
    ynf=np.array(Y,dtype=float);bf=np.array(B,dtype=float)
    xminus=np.array([sum(di[k] for di in DELTA) for k in range(10)],dtype=float)
    kinetics=[]
    for name,active in [('preserve_B_minus_L',()),('wash_all',((0,0),(1,1),(2,2)))]:
        rr=np.array(constraints(active),dtype=float)
        A=rr/sw[None,:];K=A.T@A
        eig,u=np.linalg.eigh(K);assert min(eig)>-1e-12
        eig=np.maximum(eig,0);eig[eig<1e-12]=0
        z0=sw*start
        rows=[]
        for s in (0.,1.,10.,100.,1000.):
            z=u@(np.exp(-eig*s)*(u.T@z0));pot=z/sw;dens=wf*pot
            free=float(z@z/2);diss=float(np.sum((rr@pot)**2))
            assert abs(float(ynf@dens))<1e-12
            rows.append(dict(relaxation_parameter=s,B=float(bf@dens),
                B_minus_L=float(xminus@dens),quadratic_free_energy=free,
                dissipation=diss))
        assert all(b['quadratic_free_energy']<=a['quadratic_free_energy']+1e-18
                   for a,b in zip(rows,rows[1:]))
        if not active:
            assert max(abs(z['B_minus_L']-.001) for z in rows)<1e-12
            assert abs(rows[-1]['B']-float(F(28,79))*.001)<1e-12
        else:assert abs(rows[-1]['B'])<1e-12 and abs(rows[-1]['B_minus_L'])<1e-12
        kinetics.append(dict(name=name,rows=rows,smallest_positive_rate=float(min(eig[eig>0]))))
    files=[Path(__file__).resolve(),HERE/'drafts/STATUS.md',HERE/'drafts/selection.md',
        STAGE/'993/common_candidate_v1.md',STAGE/'research_note_993.md',
        STAGE/'981/drafts/common_parent_contract_v1.md',STAGE/'research_note_991.md',
        STAGE/'research_note_992.md',STAGE.parent/'archive_629_652/research_note_629.md',
        STAGE.parent/'archive_702_741/research_note_709.md',STAGE.parent/'archive_702_741/research_note_710.md']
    return dict(round=994,all_scientific_checks_passed=True,
        kind='conditional_common_thermal_history_charge_audit',species_order=NAMES,
        susceptibility_weights=[str(x) for x in W],hypercharges=[str(x) for x in Y],
        reaction_rows=[[str(x) for x in rr] for rr in base],
        neutral_equilibrium_branches=branches,conversion_B_over_B_minus_L='28/79',
        normalized_preserved_charge_example=charges(mu),
        one_fast_flavor_example=charges(mu_partial),relaxation=kinetics,
        actual_cosmological_rates_computed=False,initial_asymmetry_generated=False,
        real_singlet_carries_B_minus_L=False,thermal_flavor_cases_are_actual_history=False,
        full_goal_completed=False,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(TARGET.read_text('utf-8-sig')))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','reaction_rows',
        'neutral_equilibrium_branches','relaxation')},ensure_ascii=False,indent=2))
