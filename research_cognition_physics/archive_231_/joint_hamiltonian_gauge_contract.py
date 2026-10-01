"""613: spatial CHN Hamiltonian versus the inherited chiral gauge contract.

The CHN identities are established literature. New interface tests use the
original representation/charges, weak algebra and the same geometric source.
Frozen free/vectorlike diagnostics are not the original interacting SM model.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_chiral_fibre_source as chiral
import joint_gapped_link_locality as link
import joint_fermion_gauss_completion as old
import joint_original_mass_spinor_bridge as bridge

HERE=Path(__file__).resolve().parent
TARGET=HERE/"joint_hamiltonian_gauge_contract_results.json"
I4=np.eye(4,dtype=complex)
BETA=chiral.GAMMA[3]
G5=chiral.G5

def op(a):
    return float(np.linalg.norm(a,2))

def symbol(k):
    k=np.asarray(k,float)
    s=np.sin(k);w=float(np.sum(1-np.cos(k))-1)
    r=float(np.sqrt(w*w+s@s))
    X=w*I4+sum(1j*chiral.GAMMA[i]*s[i] for i in range(3))
    D=I4+X/r
    h=BETA@D
    q=G5@(I4-D/2)
    return h,q,float(2*(1+w/r)),r

def exp_i(A,angle):
    e,v=np.linalg.eigh(A)
    return (v*np.exp(1j*angle*e))@v.conj().T

def kinetic_check():
    rng=np.random.default_rng(613)
    maximal=0.
    for k in rng.uniform(-np.pi,np.pi,(35,3)):
        h,q,e2,_=symbol(k)
        maximal=max(maximal,op(h-h.conj().T),op(q-q.conj().T),
                    op(h@q-q@h),op(q@q+h@h/4-I4),op(h@h-e2*I4))
    assert maximal<2e-14
    corners=[]
    for bits in itertools.product((0,1),repeat=3):
        h,q,e2,_=symbol(np.pi*np.array(bits))
        expected=0. if not any(bits) else 2.
        assert abs(np.sqrt(e2)-expected)<1e-13
        corners.append(dict(bits=list(bits),absolute_energy=float(np.sqrt(e2))))
    # Same inherited hopping symbol, restricted to zero temporal momentum.
    H,_,sites=link.kernel(4)
    k=np.array([np.pi/2,0,0,0])
    wave=np.exp(1j*np.asarray(sites)@k)/np.sqrt(len(sites))
    F=np.kron(wave[:,None],I4)
    X=G5@(F.conj().T@H@F)
    w,v=np.linalg.eigh(X.conj().T@X)
    D=I4+X@((v*(1/np.sqrt(w)))@v.conj().T)
    inherited=op(BETA@D-symbol(k[:3])[0])
    assert inherited<1e-12
    h,q,e2,_=symbol(k[:3])
    return dict(maximal_operator_identity_error=maximal,corner_energies=corners,
                inherited_kernel_error=inherited,selected_energy=float(np.sqrt(e2)),
                kinetic_is_vectorlike_not_Weyl=True,continuous_time_is_a_changed_regulator=True)

def gauge_check():
    rows=[]
    for k in ([np.pi/2,0,0],[.3,.2,-.1],[np.pi,0,0]):
        h,q,e2,_=symbol(k)
        P=(I4-q)/2
        T=[np.kron(P,s/2) for s in chiral.S]
        defect=T[0]@T[1]-T[1]@T[0]-1j*T[2]
        predicted=-1j*np.kron(h@h/16,chiral.S[2]/2)
        error=op(defect-predicted)
        assert error<1e-14
        assert abs(op(defect)-e2/32)<1e-14
        period={}
        for name,charge in bridge.CHARGES.items():
            filt=(I4-q)/2 if name in bridge.LEFT else (I4+q)/2
            period[name]=op(exp_i(charge*filt,2*np.pi)-I4)
        if k[0]==np.pi/2:
            assert all(period[n]>.1 for n in bridge.CHARGES if n!="nu")
        rows.append(dict(k=k,energy_squared=e2,
            projector_defect=op(P@P-P),weak_algebra_defect=op(defect),
            exact_defect_formula_error=error,U1_period_defects=period))
    old_period=op(old.representation(np.eye(3),np.eye(2),np.exp(2j*np.pi))-np.eye(32))
    assert old_period<1e-13
    # Normalization works away from q=0 but cannot extend through the cutoff
    # corner with one continuous momentum-space value.
    approaches=[]
    for eps in (.03,.01,.003):
        chis=[]
        for direction in (np.array([1.,0,0]),np.array([0.,1,0])):
            _,q,_,_=symbol(np.array([np.pi,0,0])+eps*direction)
            qe,qv=np.linalg.eigh(q)
            assert np.min(abs(qe))>1e-5
            chi=(qv*np.sign(qe))@qv.conj().T
            assert op(chi@chi-I4)<1e-13
            chis.append(chi)
        distance=op(chis[0]-chis[1])
        assert distance>1.3
        approaches.append(dict(epsilon=eps,chirality_direction_distance=distance))
    return dict(rows=rows,inherited_U1_period_error=old_period,corner_approaches=approaches,
                no_global_gauging_from_nonidempotent_filter=True,
                not_a_no_go_for_all_chiral_gauge_regulators=True)

def grid_spectrum(L=4):
    return np.concatenate([np.linalg.eigvalsh(symbol(2*np.pi*np.array(n)/L)[0])
                           for n in itertools.product(range(L),repeat=3)])

def thermodynamics(e0,sigma,beta,copies=16):
    e=e0*np.exp(-sigma)
    occupation=1/(1+np.exp(beta*e))
    free=-copies*float(np.sum(np.logaddexp(0,-beta*e)))/beta
    U=copies*float(e@occupation)
    sea=copies*float(np.sum(np.minimum(e,0)))
    return free,U,sea

def source_check():
    L=4;e0=grid_spectrum(L)
    # 16 original internal Weyl entries are deliberately extended to Dirac
    # entries here; this doubles the physical spin content. No projection is
    # silently assumed to have removed the mirrors.
    copies=16;beta=.7;sigma=.2;eps=1e-5
    F,U,sea=thermodynamics(e0,sigma,beta,copies)
    fp,_,sp=thermodynamics(e0,sigma+eps,beta,copies)
    fm,_,sm=thermodynamics(e0,sigma-eps,beta,copies)
    source_error=abs((fp-fm)/(2*eps)+U)
    sub_source_error=abs(((fp-sp)-(fm-sm))/(2*eps)+U-sea)
    assert source_error<2e-6 and sub_source_error<2e-6
    V=L**3*np.exp(3*sigma)
    pressure=sea/(3*V)
    cosmological_density=-sea/V
    residual=pressure-cosmological_density
    assert abs(residual-4*sea/(3*V))<1e-12
    # Direct volume differentiation independently checks the fixed-grid pressure.
    vp=L**3*np.exp(3*(sigma+eps));vm=L**3*np.exp(3*(sigma-eps))
    pressure_error=abs(-(sp-sm)/(vp-vm)-pressure)
    assert pressure_error<1e-7
    return dict(L=L,internal_Dirac_copies=copies,one_particle_dimension=copies*len(e0),
        original_target_modes_per_node=32,vectorlike_modes_per_node=64,
        beta=beta,sigma=sigma,free_energy=F,mean_energy=U,sea_energy=sea,
        free_energy_source_error=source_error,subtracted_source_error=sub_source_error,
        vacuum_pressure=pressure,vacuum_pressure_derivative_error=pressure_error,
        cosmological_counterterm_density_at_selected_scale=cosmological_density,
        residual_pressure_after_energy_cancellation=residual,
        fixed_graph_Gibbs_finite=True,
        original_interacting_or_Gauss_projected_partition_not_computed=True,
        regulator_vacuum_not_observed_cosmological_constant=True)

def run():
    deps=("research_note_578.md","research_note_579.md","research_note_598.md",
          "research_note_602.md","research_note_603.md","research_note_611.md",
          "research_note_612.md","joint_chiral_fibre_source.py",
          "joint_gapped_link_locality.py","joint_fermion_gauss_completion.py",
          "joint_original_mass_spinor_bridge.py")
    return dict(round=613,tests_run=3,failures=0,errors=0,
        kinetic=kinetic_check(),gauge=gauge_check(),sources=source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(mature_CHN_identities_not_claimed_as_new=True,
            original_weak_and_U1_contract_tested=True,
            naive_chiral_filter_rejected_not_all_routes=True,
            vectorlike_finite_graph_thermal_and_source_branch=True,
            quantum_SM_GR_and_continuum_still_open=True))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write-results",action="store_true")
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    else:
        assert json.loads(TARGET.read_text("utf8"))==result
    print(json.dumps(result,ensure_ascii=False))
