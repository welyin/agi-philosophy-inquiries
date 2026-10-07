"""834: syndrome organization, parity and low-degree apparatus sources.

Exact Pauli/combinatorial certificates, plus a differential-Clifford check.
This constructs an endpoint isometry, not an autonomous native device.
"""
from pathlib import Path
import argparse,itertools,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'controller_source_embedding_results.json'
sys.path[:0]=[str(HERE.parent/'829'),str(HERE.parent/'833')]
import majorana_code_source_bridge as code
import joint_short_time_filtration as algebra


def syndrome_completion():
    gamma,st,_,syndrome,comp,_,_=code.code_data()
    representatives={0:code.I}
    for bit in range(10):
        for e in ((1<<bit,0,0),(0,1<<bit,0)):
            se=syndrome(e)
            for s,r in list(representatives.items()):
                representatives.setdefault(s^se,code.mul(e,r))
    assert len(representatives)==512
    errors=[code.I,*gamma]+[code.mul((0,0,1),code.mul(gamma[i],gamma[j]))
                            for i in range(20) for j in range(i+1,20)]
    correctable={}
    for e in errors:correctable.setdefault(syndrome(e),e)
    representatives.update(correctable)
    parity=(0,1023,0)
    parity_counts={-1:0,1:0}
    for s,r in representatives.items():
        assert syndrome(r)==s
        p=-1 if code.commute(parity,r) else 1  # original code is odd
        predicted=-(-1)**((s&31).bit_count())
        assert p==predicted
        parity_counts[p]+=1
    for e in errors:
        r=representatives[syndrome(e)]
        assert comp(code.mul(code.adj(r),e))[0]=='scalar'
    assert parity_counts=={-1:256,1:256}
    return dict(full_syndromes=512,first_order_correctable_syndromes=len(correctable),
                syndrome_qubits=9,logical_and_syndrome_code_blocks=10,
                encoded_CAR_modes=100,parity_tag_CAR_modes=1,
                total_output_CAR_modes=101,additional_pure_CAR_modes=91,
                syndrome_counts_by_input_parity=parity_counts,
                parity_preserving_basis_images=1024,
                parity_tag_depends_only_on_syndrome=True)


def product_compression():
    gamma,_,_,_,comp,_,_=code.code_data()
    local=[]
    for degree in range(5):
        scalar=0
        for ids in itertools.combinations(range(20),degree):
            kind,_=comp(code.product(gamma[i] for i in ids))
            assert kind!='logical'
            scalar+=kind=='scalar'
        local.append(scalar)
    assert local==[1,0,0,0,5]
    # Nonzero tensor compressions require a nonzero local compression in every
    # occupied block. Jordan-Wigner strings only change phases on fixed-parity
    # blocks; they do not change zero/scalar classification.
    nonzero=[1,0,0,0,0]
    for _ in range(10):
        nonzero=[sum(nonzero[d-k]*local[k] for k in range(d+1)) for d in range(5)]
    assert nonzero==[1,0,0,0,50]
    rows=[]
    for degree in range(5):
        total=math.comb(202,degree)
        scalar=nonzero[degree]
        parity=nonzero[degree-2] if degree>=2 else 0
        # A single tag Majorana has zero compression: the syndrome label is
        # unchanged and its tag parity cannot change independently.
        assert (scalar,parity)==[(1,0),(0,0),(0,1),(0,0),(50,0)][degree]
        rows.append(dict(degree=degree,total_words=total,scalar_words=scalar,
                         input_parity_words=parity,zero_words=total-scalar-parity,
                         other_logical_words=0))
    return rows


def repetition_checks():
    x=np.array([[0.,1.],[1.,0.]]);z=np.diag([1.,-1.]);eye=np.eye(2)
    x2=np.kron(eye,x);z2=np.kron(eye,z)
    rows=[]
    for copies in (1,2,3):
        v=np.zeros((4**copies,4))
        for s in range(4):v[sum(s*4**k for k in range(copies)),s]=1
        local_z=np.kron(z2,np.eye(4**(copies-1)))
        assert np.array_equal(v.T@local_z@v,z2)
        for flips in range(1,copies+1):
            obs=np.ones((1,1))
            for k in range(copies):obs=np.kron(obs,x2 if k<flips else np.eye(4))
            compressed=v.T@obs@v
            expected=x2 if flips==copies else np.zeros((4,4))
            assert np.array_equal(compressed,expected)
        rows.append(dict(copies=copies,off_diagonal_dual_rail_degree=2*copies,
                         diagonal_quadratic_label_survives=True))
    return rows


def internal_degree(x,z):
    # Binary inverse of the full 11-mode Jordan-Wigner Majorana basis.
    degree=0
    for j in range(10):
        odd=((z>>j)&1)^((x>>(j+1)).bit_count()&1)
        even=((x>>j)&1)^odd
        degree+=odd+even
    return degree


def twisted_filtration():
    gamma,_,_,_,comp,_,_=code.code_data()
    g20=(1<<10,1023,0);g21=(1<<10,2047,1)
    cross=[code.mul((0,0,1),code.mul(gamma[i],g20)) for i in (0,1)]
    inside=[code.mul((0,0,1),code.mul(gamma[i],gamma[j])) for i,j in ((4,5),(12,14))]
    h=algebra.add(algebra.term(0,2,code.I,-1),algebra.term(4,0,code.I),
                  *[algebra.term(p,0,e) for p,e in zip((1,2,2,3),cross+inside)])
    parity=(0,1023,0)
    obs=algebra.term(0,0,parity)
    rows=[]
    for n in range(7):
        # R_n = ad_H^n(Pi_C) Pi_C. Right multiplication has no Pauli sign.
        weight=max(k[1]+internal_degree(k[2],k[3]^1023) for k in obs)
        assert weight<=n
        logical=sum(comp((k[2]&1023,k[3]&1023,0))[0]=='logical' for k in obs)
        if n<=5:assert logical==0
        assert all(complex(v).real.is_integer() and complex(v).imag.is_integer()
                   and abs(v)<2**40 for v in obs.values())
        rows.append(dict(commutators=n,terms=len(obs),maximum_twisted_weight=weight,
                         nontrivial_logical_terms=logical))
        if n<6:obs=algebra.comm(h,obs)
    # Negative control: full quartic vertices with three internal legs violate
    # the original quadratic-fermion hypothesis and expose parity at order 2.
    f1=code.product([gamma[0],gamma[1],gamma[4],g20])
    f2=code.product([gamma[5],gamma[12],gamma[14],g21])
    assert code.adj(f1)==f1 and code.adj(f2)==f2
    bad=algebra.add(algebra.term(0,0,f1),algebra.term(0,0,f2))
    second=algebra.comm(bad,algebra.comm(bad,algebra.term(0,0,parity)))
    logical=[k for k in second if comp((k[2]&1023,k[3]&1023,0))[0]=='logical']
    assert logical
    return dict(original_structural_example=rows,
                six_is_an_allowed_order_not_a_saturation_claim=True,
                quartic_negative_control_order_two_logical_terms=len(logical))


def run():
    return dict(round=834,all_checks_passed=True,fresh_test_groups=1,
                syndrome_embedding=syndrome_completion(),
                complete_degree_at_most_four_classification=product_compression(),
                repeated_label_limit=repetition_checks(),
                parity_twisted_filtration=twisted_filtration(),
                endpoint_requires_additional_pure_internal_resource=True,
                original_autonomous_preparation_decoder_or_energy_supply_proven=False,
                encoded_source_pair_requires_direct_product_compression=True,
                graph_continuum_or_standard_model_generation_proven=False,
                scope='An exact full-syndrome endpoint isometry embeds logical and syndrome qubits into ten copies of the known code with a parity tag. All degree-at-most-four carrier operators compress to identity or input carrier parity. A twisted differential/Clifford filtration delays parity disclosure through fifth order for the original quadratic class. No autonomous native isometry, ready resource supply, or matching of absolute sources is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    # JSON normalizes integer dictionary keys.
    r=json.loads(json.dumps(r))
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
