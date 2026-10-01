"""629: original quotient topology and mass phases share one rephasing quotient.

The anomalous Jacobian/index and topology theorems are continuum inputs.
Finite CAR conjugation here checks the old mass dictionary, not that anomaly.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_global_anomaly_bundle as poly
import joint_spinor_subgroup_mass as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_topological_mass_phases_results.json'
NAMES=('Q','uc','dc','L','ec','nuc')
SIZES=(6,3,3,2,1,1)
A=np.array([[3,-2,-3,-12],[2,-1,0,-6],[1,-1,0,-6],
            [1,0,-1,-4],[1,0,0,-2],[0,0,0,-2]],dtype=int)
C=np.array([[1,1,0,0,0,0],[1,0,1,0,0,0],[0,0,0,1,1,0],
            [0,0,0,1,0,1],[0,0,0,0,0,2]],dtype=int)
K=np.array([[2,1,1,0,0,1,0,0,0],[-1,-1,0,0,0,0,1,0,0],
            [-6,-6,-2,-2,0,0,0,0,1]],dtype=int)


def exact_rank(matrix):
    m=[[Fraction(int(x)) for x in row] for row in matrix]
    rank=0
    for col in range(len(m[0])):
        pivot=next((i for i in range(rank,len(m)) if m[i][col]),None)
        if pivot is None:continue
        m[rank],m[pivot]=m[pivot],m[rank]
        v=m[rank][col];m[rank]=[x/v for x in m[rank]]
        for i in range(len(m)):
            if i!=rank:
                v=m[i][col];m[i]=[x-v*y for x,y in zip(m[i],m[rank])]
        rank+=1
        if rank==len(m):break
    return rank


def characteristic_check():
    roots=poly.variables(4);roots+= [poly.scale(poly.add(*roots),-1)]
    a=poly.add(*roots[:3]);u=poly.scale(poly.power(a,2),Fraction(1,2))
    v=poly.elementary(roots[:3],2);w=poly.elementary(roots[3:],2)
    J=old.dictionary();module={};start=0
    for name,size,expected in zip(NAMES,SIZES,A):
        weights=[]
        for col in range(start,start+size):
            ix=np.flatnonzero(J[:,col]);assert len(ix)==1
            s=old.STATES[int(ix[0])]
            weights.append(poly.add(*(roots[j] for j in s)))
        start+=size
        ch2=poly.scale(poly.add(*(poly.power(x,2) for x in weights)),Fraction(1,2))
        formula=poly.add(*(poly.scale(p,int(c)) for p,c in zip((u,v,w),expected[:3])))
        assert ch2==formula and expected[3]==-2*size
        module[name]=dict(index_coefficients=expected.tolist(),ch2=poly.rows(ch2))
    assert (A.sum(axis=0)==np.array([8,-4,-4,-32])).all()
    # Geometric existence/integrality is proved/cited in the note, not inferred
    # from an arbitrary numerical matrix. Rows are actual characteristic data.
    witnesses=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]],int)
    assert exact_rank(witnesses)==4
    assert np.array_equal(A@witnesses[0],np.array([3,2,1,1,1,0]))
    # Pullback V imposes a particular one-dimensional phase direction only.
    c2V=np.array([-2,1,1,0],int)
    assert exact_rank(c2V[None,:])==1
    return dict(integer_coordinates=['a_squared_over_2','c2_E3','c2_E2','p1_TM_over_48'],
        module_index_matrix=A.tolist(),module_polynomial_checks=module,
        witness_names=['S2xS2_nonliftable','S4_colour_bundle_c2_1','S4_weak_bundle_c2_1','K3_trivial_gauge'],
        characteristic_witnesses=witnesses.tolist(),integral_basis_determinant=-1,
        SU5_c2_pullback_coefficients=c2V.tolist(),
        SU5_single_gauge_theta_does_not_cover_full_original_gauge_theta_family=True,
        inputs='Wu even intersection pairing and Rokhlin divisibility on smooth closed spin 4-manifolds; unit SU instanton bundles; K3 signature -16; Omega4^Spin(BG)=Z^4 from cited literature.')


def actual_g(phi):
    W=old.ph_matrix();B=W@old.bdg(phi)@W.conj().T
    delta=B[:32,32:]
    # Holomorphic annihilation coefficients are conjugates of creation delta.
    entries=[(0,13),(2,19),(26,29),(24,31),(30,31)]
    F=float(old.old.original.F(phi))
    amplitudes=np.array([phi[1]/np.sqrt(F)]*4+[phi[4]/np.sqrt(F)])
    g=np.array([delta[i,j].conjugate() for i,j in entries])/amplitudes
    expected=-np.array([old.old.Y[x].conjugate() for x in ('u','d','e','nu')]+[old.old.Y['s']])
    assert np.max(abs(g-expected))<1e-14
    return B,g,amplitudes,entries


def phase_quotient_check():
    T=np.vstack((C,-A.T))
    assert exact_rank(C)==5 and exact_rank(T)==6 and exact_rank(K)==3
    assert not (K@T).any()
    assert np.array_equal(K[:,[5,6,8]],np.eye(3,dtype=int))
    baryon=np.array([1,-1,-1,0,0,0])
    assert not (C@baryon).any()
    assert np.array_equal(-A.T@baryon,[0,0,3,0])
    hypercharge=np.array([1,-4,2,-3,6,0])
    assert not (A.T@hypercharge).any()
    assert np.array_equal(C@(-hypercharge),[3,-3,-3,3,0])
    phi=np.array([0.,.4,0.,0.,.7]);B,g,amplitudes,entries=actual_g(phi)
    phases=np.angle(g);rng=np.random.default_rng(629)
    maximum_mass=maximum_spectrum=maximum_invariant=0.
    for _ in range(24):
        alpha=rng.normal(size=6)*5
        alpha16=np.repeat(alpha,SIZES);alpha32=np.repeat(alpha16,2)
        U=np.diag(np.r_[np.exp(-1j*alpha32),np.exp(1j*alpha32)])
        moved=U@B@U.conj().T;delta=moved[:32,32:]
        actual=np.array([delta[i,j].conjugate() for i,j in entries])/amplitudes
        err=float(np.max(abs(actual-g*np.exp(1j*(C@alpha)))))
        maximum_mass=max(maximum_mass,err)
        maximum_spectrum=max(maximum_spectrum,float(np.max(abs(np.linalg.eigvalsh(moved)-np.linalg.eigvalsh(B)))))
        theta=rng.normal(size=4);point=np.r_[phases,theta]
        moved_point=point+T@alpha
        maximum_invariant=max(maximum_invariant,float(np.max(abs(np.exp(1j*(K@moved_point))-np.exp(1j*(K@point))))))
    assert max(maximum_mass,maximum_spectrum,maximum_invariant)<1e-12
    # Eliminate all five coupling phases, retaining the anomalous theta shift.
    alpha=np.zeros(6);alpha[5]=-phases[4]/2
    alpha[3]=-phases[3]-alpha[5];alpha[4]=-phases[2]-alpha[3]
    alpha[1]=-phases[0];alpha[2]=-phases[1]
    zero=phases+C@alpha;theta_shift=-A.T@alpha
    assert np.max(abs(zero))<1e-14
    full=np.r_[zero,theta_shift];base=np.r_[phases,np.zeros(4)]
    proper_error=float(np.max(abs(np.exp(1j*(K@full))-np.exp(1j*(K@base)))))
    erased=np.zeros(9)
    naive_error=float(np.max(abs(np.exp(1j*(K@erased))-np.exp(1j*(K@base)))))
    assert proper_error<1e-13 and naive_error>.1
    # SU5-correlated gauge phases remain gauge-invariant but are a restriction,
    # and their subspace is not preserved by allowed baryon reparameterization.
    c2_direction=np.array([-2,1,1])
    assert exact_rank(np.array([c2_direction,[0,0,3]]))==2
    return dict(coupling_names=['u','d','e','nu','s'],C=C.tolist(),parameter_action=T.tolist(),
        invariant_character_rows=K.tolist(),ranks=dict(C=5,action=6,invariants=3),
        independent_phase_torus_dimension=3,
        baryon_shift=[0,0,3,0],baryon_preserves_all_original_mass_couplings=True,
        original_holomorphic_couplings=[[float(z.real),float(z.imag)] for z in g],
        maximum_CAR_mass_rephasing_error=maximum_mass,maximum_finite_BdG_spectrum_error=maximum_spectrum,
        maximum_invariant_phase_error=maximum_invariant,
        all_real_mass_representative_theta=theta_shift.tolist(),
        correct_joint_rephasing_error=proper_error,erase_mass_phases_without_theta_error=naive_error,
        continuous_field_rephasing_is_basis_change_not_new_gauge_symmetry=True,
        limits='Complete quotient only of the stated nine-angle family by six constant module rephasings, with nonzero one-generation couplings and no extra baryon-violating operator. Not a classification of all QFT parameters, all field dualities, or experimentally measurable CP observables. The continuum Jacobian is an inherited index theorem, not derived from finite CAR.')


def run():
    deps=('research_note_531.md','research_note_581.md','research_note_598.md',
          'research_note_614.md','research_note_616.md','research_note_628.md',
          'joint_global_anomaly_bundle.py','joint_spinor_subgroup_mass.py',
          'joint_fermion_gauss_completion.py')
    return dict(round=629,tests_run=2,failures=0,errors=0,
        topology_and_original_module_indices=characteristic_check(),
        common_mass_and_topological_phases=phase_quotient_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Conditional smooth closed spin four-dimensional original quotient-group continuum target, full topological background menu, five nonzero one-generation mass couplings and four characteristic theta angles. Mature index Jacobians and original finite CAR dictionary give a joint three-torus of phase classes under constant module rephasings. No unique phase values, full regional measure, graph-continuum matching, observational fit, or GR generation is claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(round=629,tests_run=2,phase_quotient=result['common_mass_and_topological_phases']),ensure_ascii=False,indent=2))
