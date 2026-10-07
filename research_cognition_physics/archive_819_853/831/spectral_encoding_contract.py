"""831: finite certificates for the spectral encoding and orbit interface.

The examples are spectral-algebra checks, not spectra of the original
infinite-dimensional graph Hamiltonian. Its compact resolvent is inherited
from 598 and used analytically with the native 830 obstruction.
"""
from pathlib import Path
import argparse, json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'spectral_encoding_contract_results.json'


def norm(a): return float(np.linalg.norm(a, 2))


def grams(W, multiplicity):
    return [a.conj().T@a for a in W.reshape(-1, multiplicity, W.shape[1])]


def moment_examples():
    energies=np.array([-2., -1., 0., 1., 2.])
    p=np.array([0., .5, 0., .5, 0.])
    q=np.array([.125, 0., .75, 0., .125])
    a=np.kron(np.sqrt(p)[:, None], np.eye(2))
    b=np.kron(np.sqrt(q)[:, None], np.eye(2))
    assert norm(a.conj().T@a-np.eye(2))<1e-14
    assert norm(b.conj().T@b-np.eye(2))<1e-14
    differences=[float(abs(p@(energies**k)-q@(energies**k))) for k in range(5)]
    assert max(differences[:4])<1e-14 and abs(differences[4]-3)<1e-14
    gap=max(norm(x-y) for x,y in zip(grams(a,2),grams(b,2)))
    assert abs(gap-.75)<1e-14
    # Equal energy probabilities for each basis input omit logical coherence.
    a2=np.vstack([np.eye(2),np.eye(2)])/np.sqrt(2)
    b2=np.array([[1.,1.],[0.,0.],[1.,-1.],[0.,0.]])/np.sqrt(2)
    ga,gb=grams(a2,2),grams(b2,2)
    assert norm(b2.conj().T@b2-np.eye(2))<1e-14
    diagonal=max(float(np.max(abs(np.diag(x-y)))) for x,y in zip(ga,gb))
    full=max(norm(x-y) for x,y in zip(ga,gb))
    assert diagonal==0 and abs(full-.5)<1e-14
    return dict(equal_first_three_energy_moments=True,
                fourth_moment_difference=differences[4], spectral_Gram_gap=gap,
                basis_only_probability_test_gap=diagonal,
                coherent_spectral_Gram_gap=full,
                both_naive_sufficiency_tests_correctly_rejected=True)


def block_unitary_certificate():
    rng=np.random.default_rng(831)
    def unitary(n):
        z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        u,_=np.linalg.qr(z)
        return u
    # Degenerate eigenspaces contain transformations not given by free H time.
    energies=np.array([0.,1.,3.]); multiplicity=4; dimension=12
    w=unitary(dimension)[:,:2]
    u=np.zeros((dimension,dimension),complex)
    for j in range(3):u[4*j:4*j+4,4*j:4*j+4]=unitary(multiplicity)
    h=np.diag(np.repeat(energies,multiplicity)); v=u@w
    gap=max(norm(x-y) for x,y in zip(grams(w,4),grams(v,4)))
    assert gap<1e-14 and norm(h@u-u@h)<1e-14
    # Reconstruct a unitary independently in each block by matching its
    # partial isometry, then extend over the orthogonal complement.
    reconstructed=np.zeros_like(u)
    for j in range(3):
        sl=slice(4*j,4*j+4);x,y=w[sl],v[sl]
        ax,s,bx=np.linalg.svd(x,full_matrices=True)
        rank=int(np.sum(s>1e-12)); assert rank==2
        image=(y@bx[:rank].conj().T)/s[:rank]
        assert norm(image.conj().T@image-np.eye(rank))<1e-13
        # QR supplies only the unused complement. Keep the mapped columns.
        complete=np.linalg.qr(image,mode='complete')[0]
        target_basis=np.column_stack([image,complete[:,rank:]])
        reconstructed[sl,sl]=target_basis@ax.conj().T
    err=norm(reconstructed@w-v)
    assert err<1e-13
    assert norm(reconstructed.conj().T@reconstructed-np.eye(dimension))<1e-13
    # Idle flag example: same spectral Gram and an energy-conserving flag
    # swap exist, but exp(-itH) cannot change that flag at any time.
    probs=np.array([.2,.3,.5]); flag0=np.array([[1.],[0.]])
    flag1=np.array([[0.],[1.]])
    win=np.kron(np.sqrt(probs)[:,None],np.kron(np.eye(2),flag0))
    want=np.kron(np.sqrt(probs)[:,None],np.kron(np.eye(2),flag1))
    swap=np.kron(np.eye(6),np.array([[0.,1.],[1.,0.]]))
    assert norm(swap@win-want)<1e-14 and norm(h@swap-swap@h)==0
    max_overlap=0.
    for t in (0.,.13,1.,np.sqrt(2),10.,1000.):
        evolved=np.exp(-1j*np.diag(h)*t)[:,None]*win
        max_overlap=max(max_overlap,norm(want.conj().T@evolved))
    assert max_overlap==0
    return dict(matched_Gram_residual=gap, reconstructed_unitary_encoding_error=err,
                spectral_unitary_exists_does_not_imply_free_time_realization=True,
                idle_flag_free_evolution_overlap=max_overlap,
                idle_flag_obstruction_analytic_for_all_times=True,
                idle_flag_example_is_not_original_material=True)


def tail_certificate():
    # Common spectral tails control every block unitary, not just one sampled
    # elapsed time. Samples check the bound; the universal argument is analytic.
    rng=np.random.default_rng(83191);count=32
    p=2.**(-np.arange(1,count+1));p/=p.sum()
    input_blocks=np.sqrt(p)[:,None,None]*np.eye(2)[None,:,:]
    output=[]
    for block in input_blocks:
        u,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
        output.append(u@block)
    output=np.array(output)
    rows=[]
    for cutoff in (2,4,8,16):
        tail=norm(input_blocks[cutoff:].reshape(-1,2))
        residual=norm(output[cutoff:].reshape(-1,2))
        assert abs(tail-residual)<1e-14
        rows.append(dict(spectral_blocks=cutoff,uniform_isometry_tail=tail,
                         output_tail=residual,channel_continuity_factor_bound=2*tail))
    return dict(rows=rows,weights_are_declared_diagnostic_not_original_spectrum=True,
                original_encoding_error_gap_not_numerically_computed=True)


def run():
    return dict(round=831,all_checks_passed=True,fresh_test_groups=1,
                spectral_necessity_counterexamples=moment_examples(),
                spectral_sufficiency_and_free_evolution=block_unitary_certificate(),
                compact_orbit_tail_certificate=tail_certificate(),
                scope='Spectral Gram equality is necessary and sufficient for a commuting unitary between specified pure input/output isometries when H has finite-multiplicity pure point spectrum. It is stronger than finitely many moments and weaker than free-H reachability. Compactness applied analytically to the original 598/623 model and 830 input gives a strictly positive fixed-preparation encoding error gap; its numerical value and any bound uniform in resource/model changes are not obtained.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
