"""680: complete original mass reweighting and the independent normalization gate.
The first group keeps the actual original background and all mass coefficients.
The finite Grassmann examples test logical implications, not the physical model.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
import joint_physical_source_reflection as prior
old=prior.old;base=prior.base
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_mass_reflection_congruence_results.json'

def err(a):return float(np.max(abs(a),initial=0))

def pf_small(a):
    return 1.+0j if len(a)==0 else base.pf(a)

def update(w,c,i,j,alpha):
    # N -> N + alpha (Phi_i^T Phi_j - Phi_j^T Phi_i).
    factor=1-alpha*c[i,j]
    assert abs(factor)>1e-8
    u=c[:,i].copy();v=c[:,j].copy()
    cc=c+alpha/factor*(np.outer(u,v)-np.outer(v,u))
    return w*factor,cc,float(abs(factor))

def full_original_mass():
    links,e,phis=old.fixture();mat=base.fixed_matrices(e,phis)
    u,v,_,_,_=base.kernel(links);eps=v@v.conj().T-u@u.conj().T
    q0=old.soft(eps,mat,0);q1=old.soft(eps,mat,.37)
    w0=base.pf(q0['N'])
    c0=q0['Phi']@np.linalg.solve(q0['N'],q0['Phi'].T)
    w=w0;c=c0.copy()
    side=np.r_[np.repeat([0,0,1,1],32),np.repeat([0,0,1,1],32)]
    entries=[(int(i),int(j),.37*mat[-1][i,j]) for i,j in zip(*np.nonzero(np.triu(mat[-1],1)))]
    assert all(side[i]==side[j] for i,j,_ in entries)
    entries.sort(key=lambda row:(-side[row[0]],row[0],row[1]))
    minfactor=1e100;reconstruct=np.zeros_like(mat[-1])
    for i,j,alpha in entries:
        w,c,f=update(w,c,i,j,alpha);minfactor=min(minfactor,f)
        reconstruct[i,j]+=alpha;reconstruct[j,i]-=alpha
    assert len(entries)==264
    assert err(reconstruct-.37*mat[-1])<1e-14
    direct_w=base.pf(q1['N'])
    direct_c=q1['Phi']@np.linalg.solve(q1['N'],q1['Phi'].T)
    relative=float(abs(w/direct_w-1));covariance=err(c-direct_c)
    assert relative<2e-10 and covariance<2e-10
    rng=np.random.default_rng(68031)
    z=rng.normal(size=(256,6))+1j*rng.normal(size=(256,6))
    z/=np.linalg.norm(z,axis=0)
    sources=[]
    for k in (0,2,4,6):
        zz=z[:,:k]
        pred=w*pf_small(zz.T@c@zz)
        direct=base.pf(q1['N']) if k==0 else old.source_coeff(q1,zz)
        rel=float(abs(pred-direct)/max(abs(pred),abs(direct)))
        assert rel<3e-10
        sources.append(dict(count=k,value=base.old.cpair(direct),relative_error=rel))
    for i,j,alpha in reversed(entries):
        w,c,f=update(w,c,i,j,-alpha);minfactor=min(minfactor,f)
    reverse=dict(weight_relative=float(abs(w/w0-1)),covariance_error=err(c-c0))
    assert max(reverse.values())<3e-10
    return dict(original_full_group_nonflat_16_channels=True,
        all_original_mass_terms_included=len(entries),
        positive_half_terms=sum(int(side[i]) for i,j,a in entries),
        negative_half_terms=sum(1-int(side[i]) for i,j,a in entries),
        reconstruction_error=err(reconstruct-.37*mat[-1]),
        full_weight_relative=relative,full_covariance_error=covariance,
        sources=sources,reverse=reverse,minimum_update_factor=minfactor,
        inverse_only_for_this_resolved_numerical_fixture=True,
        general_singular_scope_proved_without_inverse=True,
        all_Haar_and_S9_integrals_not_sampled=True)

# Exact exterior algebra with four generators, order x1,x2,y1,y2.
def add(*polys):
    out={}
    for p in polys:
        for key,val in p.items():out[key]=out.get(key,Fraction(0))+val
    return {key:val for key,val in out.items() if val}

def scale(p,a):return {k:a*v for k,v in p.items() if a*v}

def mul(p,q):
    out={}
    for a,ca in p.items():
        for b,cb in q.items():
            if a&b:continue
            swaps=sum((b&((1<<i)-1)).bit_count() for i in range(4) if a&(1<<i))
            key=a|b;out[key]=out.get(key,Fraction(0))+ca*cb*(-1)**swaps
    return {k:v for k,v in out.items() if v}

def theta(p):
    out={}
    for mask,coeff in p.items():
        term={0:coeff}  # exact real coefficients; complex extension is antilinear.
        for i in reversed([i for i in range(4) if mask&(1<<i)]):
            term=mul(term,{1<<((i+2)%4):Fraction(1)})
        out=add(out,term)
    return out

def integral(p):return -p.get(15,Fraction(0)) # integral(theta(a)*a)=1

def gram(density):
    basis=[{0:Fraction(1)},{1:Fraction(1)},{2:Fraction(1)},{3:Fraction(1)}]
    return [[integral(mul(density,mul(theta(f),g))) for g in basis] for f in basis]

def pack(a):return [[str(x) for x in row] for row in a]

def exact_normalization_boundary():
    one={0:Fraction(1)};a={3:Fraction(1)};b=theta(a)
    base_density=mul(add(one,scale(a,-1)),add(one,scale(b,-1)))
    q0=gram(base_density)
    assert q0==[[Fraction(x) for x in row] for row in
                [[1,0,0,-1],[0,0,0,0],[0,0,0,0],[-1,0,0,1]]]
    rows=[]
    for lam in (Fraction(-1),Fraction(0),Fraction(1,2),Fraction(1),Fraction(2)):
        r=add(one,scale(a,lam));ri=add(one,scale(a,-lam))
        assert mul(r,ri)==one
        density=mul(base_density,mul(theta(r),r))
        q=gram(density);v=[1-lam,Fraction(0),Fraction(0),Fraction(-1)]
        assert q==[[x*y for y in v] for x in v]
        assert q[0][0]==(1-lam)**2
        rows.append(dict(lambda_mass=str(lam),Gram=pack(q),
                         vacuum_norm=str(q[0][0]),positive_rank_one=True))
    assert rows[3]['vacuum_norm']=='0'
    # Complete-polynomial negative direction, despite a positive scalar weight.
    bad=add(one,scale(a,2),scale(b,2),mul(b,a))
    witness=add(one,scale(a,-1))
    badrows=[]
    for lam in (Fraction(0),Fraction(1),Fraction(3)):
        r=add(one,scale(a,lam));ri=add(one,scale(a,-lam))
        density=mul(bad,mul(theta(r),r))
        f=mul(ri,witness)
        value=integral(mul(density,mul(theta(f),f)))
        assert value==-2
        vacuum=integral(density);assert vacuum>0
        badrows.append(dict(lambda_mass=str(lam),vacuum_norm=str(vacuum),
                           transported_negative_norm=str(value)))
    return dict(exact_rational_arithmetic=True,RP_base_Gram=pack(q0),
        normalized_at_zero=True,mass_multipliers_all_invertible=True,
        RP_preserved_at_all_real_lambda_by_outer_product=True,
        zero_normalization_at_lambda_one=True,rows=rows,
        negative_witness_control=badrows,
        logical_counterexample_not_original_model_counterexample=True)

def run():
    return dict(date='2026-10-02',round=680,tests_run=2,failures=0,errors=0,
        full_original_mass=full_original_mass(),
        normalization_boundary=exact_normalization_boundary(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
            for p in ('research_note_668.md','research_note_669.md','research_note_673.md',
                      'research_note_679.md','joint_physical_source_reflection.py',
                      'round680_drafts/mass_congruence_probe_results.json')},
        scope=dict(complete_Gauss_physical_algebra_mass_RP_equivalence=True,
                   fixed_original_Hb_and_other_parameters=True,
                   arbitrary_finite_original_local_mass_coefficients=True,
                   nonzero_normalization_not_implied=True,
                   actual_massless_RP_not_proved=True,original_HF_not_identified=True,
                   continuum_quantum_GR_and_full_unification_not_complete=True,
                   old_space_contracts_inherited=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=680,tests_run=2,all_checks_passed=True)))
