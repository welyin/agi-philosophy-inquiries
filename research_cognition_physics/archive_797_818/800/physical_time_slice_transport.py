"""800: common physical time-slice dictionary, records and state transport.

Exact finite BV diagnostics only. Continuous relative QME repair, support,
and the represented-algebra statement are proved separately in the report.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
import sys
sys.dont_write_bytecode=True
import numpy as np
import causal_cancellation_probe as p
HERE=Path(__file__).resolve().parent

def physical(x):return x[8:10,8:10]
def expect(rho,s):
    out={}
    for key,x in s.items():
        value=np.trace(rho@physical(x))
        if value:out[key]=value
    return out
def dump(s):return {f'e{key[0]}_z{key[1]}':str(v) for key,v in sorted(s.items())}
def subtract(s,t):
    out={k:s.get(k,p.C())-t.get(k,p.C()) for k in s.keys()|t.keys()}
    return {k:v for k,v in out.items() if v}
def lifted(x):
    return p.to_c(np.kron(np.eye(5,dtype=object),x))
def representative(x,c):
    k=p.to_c(np.kron(p.quartet.elementary(5,0,4)+p.quartet.elementary(5,4,3),x))
    ex=p.delta(k,True);ex=ex+p.sharp(ex)
    return lifted(x)+p.scale(ex,c)
def run():
    probe=p.run()
    hp=representative(p.base.X,Q(1,5))
    hm=representative(p.base.Y,Q(1,7))
    past=p.exp(p.gen(hp));mid=p.exp(p.gen(hm))
    total=p.mul(mid,past)
    def beta(x):return p.mul(p.mul(p.adj(past),x),past)
    def invbeta(x):return p.mul(p.mul(past,x),p.adj(past))
    def late(x):return p.mul(p.mul(p.adj(total),{(0,0):lifted(x)}),total)
    id2=p.base.I
    e00=p.scale(id2+p.base.Z,Q(1,2))
    e11=id2-e00
    e01=p.cmat([[0,1],[0,0]]);e10=p.dag(e01)
    units={(0,0):e00,(1,1):e11,(0,1):e01,(1,0):e10}
    xunits={k:late(v) for k,v in units.items()}
    yunits={k:invbeta(v) for k,v in xunits.items()}
    rho=p.scale(id2+p.scale(p.base.X,Q(1,3))+p.scale(p.base.Y,Q(1,7))+
                p.scale(p.base.Z,Q(1,5)),Q(1,2))
    assert Q(1,9)+Q(1,49)+Q(1,25)<1
    for i in range(2):
        for j in range(2):
            assert p.equal(beta(yunits[i,j]),xunits[i,j])
            assert p.equal(p.adj(yunits[i,j]),yunits[j,i])
            for k in range(2):
                for l in range(2):
                    assert p.equal(p.mul(yunits[i,j],yunits[k,l]),
                                   yunits[i,l] if j==k else {})
    original=late(p.base.Z)
    changed=invbeta(original)
    assert expect(rho,original)==expect(rho,beta(changed))
    wrong=subtract(expect(rho,changed),expect(rho,original))
    assert wrong and min(wrong)==(1,0)
    # Same past dictionary for a finite ordered history, including the first read.
    early=[{(0,0):lifted(e00)},{(0,0):lifted(e11)}]
    late_p=[xunits[0,0],xunits[1,1]]
    history={}
    normalization={}
    for r in range(2):
        for s in range(2):
            d=p.mul(late_p[s],early[r])
            yd=invbeta(d)
            assert p.equal(invbeta(p.mul(p.adj(d),d)),p.mul(p.adj(yd),yd))
            actual=expect(rho,p.mul(p.adj(d),d))
            transported=expect(rho,beta(p.mul(p.adj(yd),yd)))
            assert actual==transported
            history[f'{r}{s}']=dump(actual)
            for k,v in actual.items():normalization[k]=normalization.get(k,p.C())+v
    normalization={k:v for k,v in normalization.items() if v}
    assert normalization=={(0,0):p.C(1)}
    # Conjugation must carry nonzero exact inputs as well.
    ex=representative(p.base.Z,Q(1,9))-lifted(p.base.Z)
    tx=beta({(0,0):ex})
    assert any(any(v.flat) for v in tx.values())
    assert all(not any(p.delta(v).flat) and not any(physical(v).flat) for v in tx.values())
    return dict(round=800,all_checks_passed=True,
                finite_causal_BV_probe=probe,
                state_transport=dict(
                    same_original_state_and_one_common_dictionary=True,
                    corrected_matrix_unit_product_identities=16,
                    original_source_recovered_order_by_order=True,
                    untransported_state_expectation_defect=dump(wrong),
                    four_ordered_record_probabilities=history,
                    history_normalization=dump(normalization),
                    nonzero_exact_inputs_preserved_and_annihilated=True),
                combined_test_groups=1,
                continuous_zero_strip_not_proven_by_matrices=True,
                arbitrary_abstract_local_BV_faithfulness_proven=False,
                auxiliary_cancellation_physically_implemented=False,
                finite_coupling_or_autonomous_apparatus_proven=False)
if __name__=='__main__':
    out=run()
    (HERE/'physical_time_slice_transport_results.json').write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))

