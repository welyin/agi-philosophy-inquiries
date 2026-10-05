"""Working 798: parity-safe register embedding and correlated-preparation audit.
Exact finite CAR diagnostics only; no original continuum apparatus claim.
"""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
import json
import sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('record797_for798',HERE.parent/'797/source_record_instrument.py')
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
C, matrix, scale, dag, same=old.C, old.matrix, old.scale, old.dag, old.same
I,X,Y,Z=old.I,old.X,old.Y,old.Z
def eye(n):
    return matrix([[i==j for j in range(n)] for i in range(n)])
def car(n):
    out=[]
    for j in range(n):
        a=scale(eye(2**n),0)
        for mask in range(2**n):
            if (mask>>j)&1:
                a[mask^(1<<j),mask]=C((-1)**((mask&((1<<j)-1)).bit_count()))
        out.append(a)
    return out
def proj(n,k):
    out=scale(eye(n),0)
    out[k,k]=C(1)
    return out
def expect(rho,a):
    v=np.trace(rho@a)
    assert not v.im
    return v.re
def run():
    c=car(1)[0]
    parity1=I-scale(dag(c)@c,2)
    # In a single fermion mode, parity-even matrices have only diagonal entries.
    even_dimensions=sum(same(parity1@e,e@parity1) for e in
                        (I,X,Y,Z))
    assert even_dimensions==2
    aa=matrix([[1,0],[0,2]])
    k0,k1=old.exact_instrument(aa,Q(1))
    # Source factor first, CAR mode second.
    e00,e11=proj(2,0),proj(2,1)
    e01=matrix([[0,1],[0,0]])
    e10=dag(e01)
    v1=np.kron(k0,e00)-np.kron(dag(k1),e01)+np.kron(k1,e10)+np.kron(dag(k0),e11)
    p1=np.kron(I,parity1)
    assert same(dag(v1)@v1,eye(4))
    violation=v1@p1-p1@v1
    assert any(violation.flat)
    a,b=car(2)
    id4=eye(4)
    na,nb=dag(a)@a,dag(b)@b
    parity=(id4-scale(na,2))@(id4-scale(nb,2))
    ev=(id4+parity)
    ev=scale(ev,Q(1,2))
    odd=id4-ev
    e0,e1=proj(4,0),proj(4,3)
    f10=dag(a)@dag(b)
    f01=dag(f10)
    assert same(f01@f10,e0) and same(f10@f01,e1)
    assert same(e0+e1,ev)
    v2=np.kron(k0,e0)-np.kron(dag(k1),f01)+np.kron(k1,f10)+np.kron(dag(k0),e1)+np.kron(I,odd)
    parity_full=np.kron(I,parity)
    assert same(dag(v2)@v2,eye(8)) and same(v2@dag(v2),eye(8))
    assert same(v2@parity_full,parity_full@v2)
    rhos=scale(I+scale(X,Q(1,3))+scale(Z,Q(1,5)),Q(1,2))
    # Exactly the pure, uncorrelated blank required by the Kraus formula.
    initial=np.kron(rhos,e0)
    out=v2@initial@dag(v2)
    for r,er,kr in ((0,e0,k0),(1,e1,k1)):
        assert expect(out,np.kron(I,er))==expect(rhos,dag(kr)@kr)
    # Same source and pointer marginals, different initial correlations.
    rho_plus=scale(np.kron(proj(2,0),e0)+np.kron(proj(2,1),e1),Q(1,2))
    rho_minus=scale(np.kron(proj(2,0),e1)+np.kron(proj(2,1),e0),Q(1,2))
    def marginal_s(rho):
        return matrix([[sum((rho[i*4+j,k*4+j] for j in range(4)),C())
                        for k in range(2)] for i in range(2)])
    def marginal_r(rho):
        return matrix([[sum((rho[i*4+j,i*4+k] for i in range(2)),C())
                        for k in range(4)] for j in range(4)])
    assert same(marginal_s(rho_plus),marginal_s(rho_minus))
    assert same(marginal_r(rho_plus),marginal_r(rho_minus))
    output_effect=dag(v2)@np.kron(I,e1)@v2
    probs=[expect(r,output_effect) for r in (rho_plus,rho_minus)]
    assert probs==[Q(7,20),Q(13,20)]
    return dict(working_round=798,all_checks_passed=True,
                single_mode_even_algebra_dimension=even_dimensions,
                single_mode_dilation_unitary_but_parity_violating=True,
                single_mode_parity_defect_entries=sum(bool(x) for x in violation.flat),
                two_mode_even_register_dilation_unitary=True,
                two_mode_encoded_blank_recovers_797_formula=True,
                same_both_marginals_correlated_record_probabilities=[str(x) for x in probs],
                old_743_to_751_native_record_results_reused_not_recounted=True,
                original_native_Hamiltonian_implements_this_dilation=False,
                original_796_state_preparation_matched=False,
                original_continuum_record_projection_constructed=False,
                formal_round_completed=False,new_numbered_test_groups=0)
if __name__=='__main__':
    result=run()
    (HERE/'internal_register_probe_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

