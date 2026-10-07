"""837: an even, spatially assigned record distinguishes source-blind states.

Only the 1024 by 2 code is stored. Partial traces are finite-factor diagnostics;
the embedding into the original free Cauchy net is proved in the report.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
TARGET=HERE/'regional_record_entropy_witness_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code

def apply_pauli(p,v):
    x,z,phase=p
    ids=np.arange(1024)
    signs=np.array([(-1)**(int(j)&z).bit_count() for j in ids])
    out=np.empty_like(v)
    out[ids^x]=(1j**phase)*signs[:,None]*v if v.ndim==2 else (1j**phase)*signs*v
    return out

def build_code():
    gamma,st,_,_,_,sx,_=code.code_data()
    seed=np.zeros(1024,complex)
    seed[sum(1<<(2*b+1) for b in range(5))]=1
    for s in st:seed=(seed+apply_pauli(s,seed))/2
    seed/=np.linalg.norm(seed)
    v=np.column_stack([seed,apply_pauli(code.product(sx),seed)])
    assert np.linalg.norm(v.conj().T@v-np.eye(2))<1e-13
    assert max(np.linalg.norm(apply_pauli(s,v)-v) for s in st)<1e-13
    witness=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    z=v.conj().T@apply_pauli(witness,v)
    assert np.linalg.norm(z@z-np.eye(2))<1e-13 and abs(np.trace(z))<1e-13
    support=[b for b in range(5) if (witness[0]|witness[1])&(3<<(2*b))]
    assert support==[0,1,3]
    return v,witness,z

def partial(v,tau,keep):
    keep=list(keep);rest=[b for b in range(5) if b not in keep]
    tensor=v.reshape((4,)*5+(2,),order='F').transpose(keep+rest+[5])
    tensor=tensor.reshape(4**len(keep),4**len(rest),2,order='F')
    return np.einsum('ael,lm,bem->ab',tensor,tau,tensor.conj())

def marginal(v,z,epsilon,r,keep):
    dim=4**len(keep);rest=4**(5-len(keep))
    return ((1-epsilon)*partial(v,(np.eye(2)+r*z)/2,keep)
            +epsilon/1022*(rest*np.eye(dim)-partial(v,np.eye(2),keep)))

def logm(a):
    w,u=np.linalg.eigh(a)
    assert w.min()>0
    return (u*np.log(w))@u.conj().T

def relative(a,b):
    ans=float(np.trace(a@(logm(a)-logm(b))).real)
    assert ans>=-1e-12
    return max(ans,0.)

def kl(p,q):
    p=np.array(p);q=np.array(q)
    return float(np.dot(p,np.log(p/q)))

def run():
    v,witness,z=build_code()
    pair_residual=0.;rows=[]
    for epsilon,r0,r in ((.1,0.,.6),(.4,.3,-.4),(511/512,0.,.6)):
        for keep in itertools.combinations(range(5),2):
            a=marginal(v,z,epsilon,r,keep);b=marginal(v,z,epsilon,r0,keep)
            pair_residual=max(pair_residual,float(np.linalg.norm(a-b)))
            assert np.linalg.norm(a-b)<1e-13
        keep=(0,1,3)
        a=marginal(v,z,epsilon,r,keep);b=marginal(v,z,epsilon,r0,keep)
        assert abs(np.trace(a)-1)<1e-13 and abs(np.trace(b)-1)<1e-13
        d3=relative(a,b)
        mean=(1-epsilon)*r;mean0=(1-epsilon)*r0
        d_read=kl([(1+mean)/2,(1-mean)/2],[(1+mean0)/2,(1-mean0)/2])
        p=[(1-epsilon)*(1+r)/2,(1-epsilon)*(1-r)/2,epsilon]
        q=[(1-epsilon)*(1+r0)/2,(1-epsilon)*(1-r0)/2,epsilon]
        d5=kl(p,q)
        logical=(1-epsilon)*kl([(1+r)/2,(1-r)/2],[(1+r0)/2,(1-r0)/2])
        assert abs(d5-logical)<1e-14
        assert 0<d_read<=d3+1e-12 and d3<=d5+1e-12
        assert abs(np.trace((np.eye(2)+r*z)/2@z).real-r)<1e-13
        rows.append(dict(outside_code_weight=epsilon,reference_bias=r0,state_bias=r,
                         binary_even_readout_relative_entropy=d_read,
                         three_block_finite_marginal_relative_entropy=d3,
                         all_five_block_flag_relative_entropy=d5,
                         three_block_density_difference_norm=float(np.linalg.norm(a-b)),
                         three_block_minimum_eigenvalue=float(np.linalg.eigvalsh(a).min()),
                         three_outcome_state_probabilities=p,
                         three_outcome_reference_probabilities=q))
    # Compression P L on the two logical eigenstates yields positive rank-one
    # P(I +/- L)/2; the remaining outcome has dimension 1022.
    evals=np.linalg.eigvalsh(z)
    assert np.max(abs(evals-[-1,1]))<1e-13
    # Full ten-mode states are even, as P has fixed odd parity and I-P is even.
    parity=code.product((0,1<<i,0) for i in range(10))
    assert code.commute(parity,witness)
    return dict(round=837,all_checks_passed=True,fresh_test_groups=1,
                exact_CAR_modes=10,dual_rail_spatial_blocks=5,
                code_isometry_residual=float(np.linalg.norm(v.conj().T@v-np.eye(2))),
                all_ten_two_block_marginal_difference_maximum=pair_residual,
                local_even_logical_witness_blocks=[0,1,3],
                all_five_block_effect_ranks=[1,1,1022],
                regional_monotonicity_rows=rows,
                reference_and_state_have_identical_degree_at_most_four_moments=True,
                spatial_embedding_proof_is_analytic=True,
                exact_regional_bound_applies_to_original_free_Cauchy_net=True,
                native_instrument_preparation_or_readout_proven=False,
                interacting_finite_coupling_von_Neumann_net_constructed=False,
                full_region_relative_entropy_equals_finite_factor_assumed=False,
                gravitational_area_response_computed=False,
                scope='A finite even measurement in a region containing the chosen compact Cauchy modes gives a strictly positive regional relative-entropy lower bound, while all source moments of fermion degree at most four agree. Pair-block regions are blind in the declared common-rest extension. Regional modular flow and gravity are not inferred from a finite density matrix.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
