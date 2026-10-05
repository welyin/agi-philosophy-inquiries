"""Working 801: extract record blocks from one fixed scattering process.

Exact finite diagnostic, not the original continuous Einstein-matter action.
The original joint state is retained. Product preparation is a comparison only.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import json
import sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'798'))
import internal_register_probe as a
C=a.C
ORDER=4
def clean(s):return {k:v for k,v in s.items() if any(v.flat)}
def add(*ss):
    out={}
    for s in ss:
        for k,v in s.items():out[k]=out[k]+v if k in out else v.copy()
    return clean(out)
def scale(s,c):return clean({k:a.scale(v,c) for k,v in s.items()})
def mul(s,t):
    out={}
    for i,x in s.items():
        for j,y in t.items():
            if i+j<=ORDER:out[i+j]=out[i+j]+x@y if i+j in out else x@y
    return clean(out)
def adj(s):return clean({k:a.dag(v) for k,v in s.items()})
def equal(s,t):return not add(s,scale(t,-1))
def exp(g):
    n=next(iter(g.values())).shape[0]
    power=out={0:a.eye(n)}
    for k in range(1,ORDER+1):
        power=mul(power,g)
        out=add(out,scale(power,Q(1,factorial(k))))
    return out
def block(s,i,j):return clean({k:v[np.ix_([i,2+i],[j,2+j])] for k,v in s.items()})
def reconstruct(blocks):
    out={}
    for i in range(2):
        for j in range(2):
            e=a.scale(a.I,0);e[i,j]=C(1)
            out=add(out,{k:np.kron(v,e) for k,v in blocks[i,j].items()})
    return out
def expectation(rho,s):
    out={}
    for k,v in s.items():
        z=np.trace(rho@v)
        if z:out[k]=z
    return out
def sumscalar(*ss):
    out={}
    for s in ss:
        for k,v in s.items():out[k]=out.get(k,C())+v
    return {k:v for k,v in out.items() if v}
def difference(s,t):return sumscalar(s,{k:-v for k,v in t.items()})
def dump(s):return {f'e{k}':str(v) for k,v in sorted(s.items())}
def marginals(rho):
    source=a.matrix([[sum((rho[2*i+j,2*k+j] for j in range(2)),C())
                      for k in range(2)] for i in range(2)])
    record=a.matrix([[sum((rho[2*i+j,2*i+k] for i in range(2)),C())
                      for k in range(2)] for j in range(2)])
    return source,record
def run():
    ident=a.eye(4);one={0:ident}
    p0=a.proj(2,0);p1=a.proj(2,1)
    records=[{0:np.kron(a.I,p)} for p in (p0,p1)]
    h1=np.kron(a.X,a.X)+a.scale(np.kron(a.Z,a.I),Q(1,3))+a.scale(np.kron(a.I,a.Y),Q(1,5))
    h2=a.scale(np.kron(a.Y,a.Z),Q(1,7))
    scattering=mul(exp({2:a.scale(h2,C(0,1))}),exp({1:a.scale(h1,C(0,1))}))
    assert equal(mul(adj(scattering),scattering),one)
    assert equal(mul(scattering,adj(scattering)),one)
    blocks={(i,j):block(scattering,i,j) for i in range(2) for j in range(2)}
    assert equal(scattering,reconstruct(blocks))
    for i in range(2):
        for j in range(2):
            assert equal(add(*(mul(adj(blocks[r,i]),blocks[r,j]) for r in range(2))),
                         {0:a.I} if i==j else {})
    rho=a.scale(ident,Q(1,4))+a.scale(np.kron(a.X,a.Y),Q(1,20))
    rho+=a.scale(np.kron(a.Y,a.Z),Q(1,28))+a.scale(np.kron(a.Z,a.I),Q(1,12))
    rho+=a.scale(np.kron(a.I,a.X),Q(1,24))
    lower=Q(1,4)-Q(1,20)-Q(1,28)-Q(1,12)-Q(1,24)
    assert lower>0 and np.trace(rho)==C(1)
    rs,rr=marginals(rho)
    product=np.kron(rs,rr)
    assert all(a.same(x,y) for x,y in zip(marginals(product),(rs,rr)))
    def omegaij(i,j,b):
        e=a.scale(a.I,0);e[i,j]=C(1)
        return expectation(rho,{k:np.kron(v,e) for k,v in b.items()})
    probabilities=[];differences=[]
    for r in range(2):
        effect=mul(mul(adj(scattering),records[r]),scattering)
        assert equal(mul(effect,effect),effect) and equal(adj(effect),effect)
        exact=expectation(rho,effect)
        byblocks=sumscalar(*(omegaij(i,j,mul(adj(blocks[r,i]),blocks[r,j]))
                             for i in range(2) for j in range(2)))
        assert exact==byblocks
        probabilities.append(dump(exact))
        differences.append(dump(difference(exact,expectation(product,effect))))
        # Under an explicitly product input the usual induced CP map is recovered.
        induced=add(*(scale(mul(adj(blocks[r,i]),blocks[r,j]),rr[j,i])
                      for i in range(2) for j in range(2)))
        assert expectation(rs,induced)==expectation(product,effect)
    assert any(differences)
    effects=[mul(mul(adj(scattering),p),scattering) for p in records]
    assert equal(add(*effects),one)
    # A declared original source need not commute with the register.
    source={0:np.kron(a.Z,a.I)+a.scale(np.kron(a.X,a.Y),Q(1,2))}
    jb={(i,j):block(source,i,j) for i in range(2) for j in range(2)}
    assert equal(reconstruct(jb),source) and jb[0,1]
    evolved=mul(mul(adj(scattering),source),scattering)
    total=sumscalar(*(omegaij(i,j,mul(mul(adj(blocks[k,i]),jb[k,l]),blocks[l,j]))
                      for i in range(2) for j in range(2)
                      for k in range(2) for l in range(2)))
    assert total==expectation(rho,evolved)
    diagonal={ij:v for ij,v in jb.items() if ij[0]==ij[1]}
    # Fill missing entries explicitly, keeping reconstruct's fixed matrix shape.
    diagonal.update({(0,1):{},(1,0):{}})
    truncated=mul(mul(adj(scattering),reconstruct(diagonal)),scattering)
    loss=difference(total,expectation(rho,truncated))
    assert loss
    return dict(working_round=801,all_checks_passed=True,epsilon_order=ORDER,
                finite_dimension=4,strict_initial_density_lower_bound=str(lower),
                one_fixed_scattering_process_used=True,
                internal_corner_blocks_reconstruct_full_process=True,
                column_unitarity_identities=4,
                outgoing_record_projections_and_completeness=True,
                same_joint_state_block_probabilities=probabilities,
                silently_factorized_preparation_defects=differences,
                original_noncommuting_source_blocks_preserved=True,
                loss_from_discarding_offdiagonal_source_blocks=dump(loss),
                product_preparation_CP_formula_verified_only_for_product_input=True,
                original_continuum_action_simulated=False,
                autonomous_terminal_or_state_preparation_proven=False,
                finite_probes_are_continuum_proof=False,
                formal_round_completed=False,new_numbered_test_groups=0)
if __name__=='__main__':
    result=run()
    (HERE/'native_scattering_blocks_probe_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
