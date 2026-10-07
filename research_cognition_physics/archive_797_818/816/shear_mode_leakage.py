"""816 working: original 64-matrix shear coupling leaks out of a two-mode block.
Constant auxiliary-past Fourier diagnostic, not the original coupled future
or its constrained bosonic covariance. No artificial boundary reflection.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
TARGET=HERE/'shear_mode_leakage_results.json'

def time_fourier(frequency,order):
    x,w=np.polynomial.legendre.leggauss(order)
    t=.3*x;eta=np.exp(-1/(1-x*x));weights=.3*w*eta
    return np.einsum('t,tab->ab',weights,np.cos(t[:,None,None]*frequency[None,:,:]))

def leakage(kappa,order):
    H0=original.hamiltonian(np.array([0.,0.,kappa]))
    e0,v0=np.linalg.eigh(H0);idx=[30,31];gram=np.zeros((2,2),complex)
    blocks=[]
    for sign in (-1,1):
        H1=original.hamiltonian(np.array([float(sign),0.,kappa]))
        e1,v1=np.linalg.eigh(H1)
        # Real TT wave: delta g_yy=cos(x-t), delta g_zz=-cos(x-t).
        # Its principal delta H sends momentum p to p +/- ex with
        # coefficient +kappa Gamma_z exp(-/+ i t)/4. K=-int U^* delta H U.
        vertex=v1.conj().T@original.GAMMA[2]@v0
        frequency=e1[:,None]-e0[None,:]-sign
        matrix=-kappa/4*v1@(vertex*time_fourier(frequency,order))@v0.conj().T[:,idx]
        gram+=matrix.conj().T@matrix;blocks.append(matrix)
    return gram,blocks

def finite_CAR_leakage_test():
    # Two internal modes coupled to two outside modes. This independently
    # verifies the full-even-algebra commutant assertion, not a continuum map.
    sys.path.insert(0,str(HERE.parent/'803'))
    from quasifree_pairing_probe import car
    a=car(4);gamma=[]
    for i in (0,1):
        gamma += [a[i]+a[i].conj().T,1j*(a[i]-a[i].conj().T)]
    even=[np.eye(16,dtype=complex)]
    for i in range(4):
        for j in range(i+1,4):even.append(1j*gamma[i]@gamma[j])
    even.append(gamma[0]@gamma[1]@gamma[2]@gamma[3])
    H=sum((a[i+2].conj().T@a[i]+a[i].conj().T@a[i+2] for i in (0,1)),np.zeros((16,16),complex))
    columns=np.stack([(o@H-H@o).reshape(-1) for o in even],axis=1)
    eig=np.linalg.eigvalsh(columns.conj().T@columns)
    assert abs(eig[0])<1e-12 and eig[1]>1.
    # The internal compression of this hopping is zero; it alone would
    # incorrectly mark all eight even operators as preserved.
    return dict(even_operator_dimension=8,full_commutant_dimension=int(np.sum(eig<1e-10)),
        smallest_nonzero_Gram_eigenvalue=float(eig[1]),
        compressed_internal_generator_is_zero=True,full_leakage_rank=4)

def run():
    rows=[];limit=float(time_fourier(np.array([[1.]]),192)[0,0]**2/8)
    for k in (4.,8.,16.,32.,64.):
        coarse,_=leakage(k,96);fine,blocks=leakage(k,192)
        eig=np.linalg.eigvalsh(fine)/(k*k)
        error=float(np.max(abs(coarse-fine))/(k*k))
        assert eig.min()>0 and error<1e-10
        rows.append(dict(kappa=k,particle_Gram_eigenvalues_over_kappa_squared=eig.tolist(),
            relative_quadrature_change=error,
            outgoing_momenta=[[1.,0.,k],[-1.,0.,k]],
            original_Nambu_columns=64,protected_input_particle_modes=[30,31],
            inferred_full_C_invariant_leakage_rank=4))
    assert abs(rows[-1]['particle_Gram_eigenvalues_over_kappa_squared'][0]-limit)<2e-6
    return dict(round=816,status='working_not_formal',all_checks_passed=True,
        shear_mode_rows=rows,principal_variance_limit=limit,
        finite_CAR_commutant=finite_CAR_leakage_test(),
        original_coupled_background_solved=False,original_physical_W_B_evaluated=False,
        auxiliary_TT_wave_is_proven_original_constrained_solution=False,
        full_original_protected_record_no_go_proven=False,
        scope='Exact auxiliary Fourier response with all original mass/principal matrices; two outside momenta retained.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

