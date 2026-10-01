"""616: original neutral CAR obstructs gauge replacement of spin/thermal data.

Spin existence uses cited geometric theorems, not numerical verification.
The exact finite example is the actual old singlet Majorana mass, with other
scalar components zero. Its sources retain the old F and complex Y_s.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as old
from joint_chiral_fibre_source import annihilators

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_spin_thermal_source_results.json'
C=annihilators(2)
PARITY=np.diag([1,-1,-1,1]).astype(complex)
ODD=(np.eye(4)-PARITY)/2

def hamiltonian(s):
    _,delta=old.mass_matrices(np.array([0.,0.,0.,0.,s]))
    d=delta[old.SLICES['nu'],old.SLICES['nu']]
    h=sum(.5*d[i,j]*C[i].conj().T@C[j].conj().T+
          .5*d[i,j].conjugate()*C[j]@C[i] for i in range(2) for j in range(2))
    return h,d[0,1]

def heat(s,beta):
    h,m=hamiltonian(s);e,v=np.linalg.eigh(h)
    K=(v*np.exp(-beta*e))@v.conj().T
    return h,m,K

def gauge_spin_check():
    rng=np.random.default_rng(616);rows=[]
    for _ in range(8):
        A=old.gauge.group_exp(rng.normal(size=8),3)
        B=old.gauge.group_exp(rng.normal(size=3),2)
        z=np.exp(1j*rng.normal());R=old.representation(A,B,z)
        neutral=R[old.SLICES['nu'],old.SLICES['nu']]
        assert np.max(abs(neutral-np.eye(2)))==0
        obstruction=float(np.linalg.norm(-neutral-np.eye(2),2))
        assert obstruction==2
        rows.append(obstruction)
    # Eight lifts for a fixed framed flat spatial 3-torus, not a dimension proof.
    structures=[list(b) for b in itertools.product((0,1),repeat=3)]
    holonomies=[[(-1)**b for b in row] for row in structures]
    assert len({tuple(row) for row in holonomies})==8
    # A nontrivial sign cannot be supplied by any original neutral gauge link.
    residuals=[max(abs(sign-1) for sign in row) for row in holonomies]
    assert residuals.count(0)==1 and residuals.count(2)==7
    return dict(neutral_spin_gauge_descent_defect=rows[0],
        subgroup_sample_count=len(rows),analytic_reason='R_nu(g)=I for every original g',
        no_compensating_original_gauge_element=True,flat_T3_spin_lifts=structures,
        fixed_frame_gauge_residuals=residuals,
        diffeomorphism_identifications_not_classified=True,
        spin_existence_from_declared_4D_oriented_global_hyperbolicity_not_numerics=True)

def thermal_check():
    rows=[];worst=0.
    for s,beta in ((.3,2.),(.8,4.),(1.2,5.)):
        h,m,K=heat(s,beta);r=abs(m)
        F=float(old.original.F(np.array([0.,0.,0.,0.,s])))
        dr=abs(old.Y['s'])*old.original.M/F**1.5
        za=float(np.trace(K).real);zp=float(np.trace(PARITY@K).real)
        assert abs(za-(2+2*np.cosh(beta*r)))<1e-12
        assert abs(zp-(2*np.cosh(beta*r)-2))<1e-12
        odd_a=float(np.trace(K@ODD).real/za)
        odd_p=float(np.trace(PARITY@K@ODD).real/zp)
        assert 0<odd_a<1 and odd_p<0
        step=2e-6;dh=(hamiltonian(s+step)[0]-hamiltonian(s-step)[0])/(2*step)
        source_a=float(np.trace(K@dh).real/za)
        source_p=float(np.trace(PARITY@K@dh).real/zp)
        a_exact=-dr*np.tanh(beta*r/2);p_exact=-dr/np.tanh(beta*r/2)
        za_plus=np.trace(heat(s+step,beta)[2]).real
        za_minus=np.trace(heat(s-step,beta)[2]).real
        zp_plus=np.trace(PARITY@heat(s+step,beta)[2]).real
        zp_minus=np.trace(PARITY@heat(s-step,beta)[2]).real
        fd_a=-(np.log(za_plus)-np.log(za_minus))/(2*step*beta)
        fd_p=-(np.log(zp_plus)-np.log(zp_minus))/(2*step*beta)
        err=max(abs(source_a-a_exact),abs(source_p-p_exact),abs(fd_a-a_exact),abs(fd_p-p_exact))
        worst=max(worst,err);assert err<2e-8
        even=(np.eye(4)+PARITY)/2
        z_even=float(np.trace(even@K).real);z_odd=float(np.trace(ODD@K).real)
        assert abs(z_even-(za+zp)/2)<1e-12 and abs(z_odd-(za-zp)/2)<1e-12
        assert np.max(abs(h@PARITY-PARITY@h))==0
        rows.append(dict(s=s,beta=beta,original_F=F,Majorana_mass=[m.real,m.imag],
            mass_magnitude=r,AP_ordinary_trace=za,P_supertrace=zp,
            ordinary_odd_probability=odd_a,supertrace_odd_expectation=odd_p,
            ordinary_source=source_a,supertrace_source=source_p,
            even_sector_partition=z_even,odd_sector_partition=z_odd,
            source_error=err))
    _,_,K0=heat(0.,2.)
    assert np.trace(K0)==4 and np.trace(PARITY@K0)==0
    return dict(rows=rows,max_source_error=worst,zero_mass_ordinary_trace=4,
        zero_mass_supertrace=0,ordinary_state_and_supertrace_not_interchangeable=True,
        parity_projection_requires_combination_of_boundary_sectors=True,
        subsystem_example_not_full_Gauss_spectrum=True)

def run():
    deps=('research_note_375.md','research_note_531.md','research_note_598.md',
          'research_note_603.md','research_note_612.md','research_note_614.md','research_note_615.md',
          'joint_fermion_gauss_completion.py','joint_chiral_fibre_source.py')
    return dict(round=616,tests_run=2,failures=0,errors=0,gauge_spin=gauge_spin_check(),
        thermal=thermal_check(),dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(geometric_spin_existence_is_conditional_mature_theorem=True,
            original_neutral_CAR_requires_uncancelled_spin_lift=True,
            thermal_trace_boundary_and_source_jointly_fixed=True,
            spatial_spin_choice_global_geometry_and_full_reconstruction_open=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
