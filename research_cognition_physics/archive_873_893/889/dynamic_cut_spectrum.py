"""889: actual dynamic cut spectrum in a declared transported, normal-ordered branch.
The full original matter family is retained formally; the numerical sector is
one electron excitation or one neutral even electron-hole pair. No local Gauss
projection, full thermal/source bridge, or derivation of the new kinetic law.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'888'))
import reference_transport_energy as prior
old=prior.old
TARGET=HERE/'dynamic_cut_spectrum_results.json'
ETA=old.ETA
KAPPA=.001
HALF_INTERVAL=.025
CENTER=-1/12

def energy(N,x,m2):
    # Exact principal-angle cut identity, avoiding subtracting two O(N) floats.
    u=6*np.abs(np.asarray(x,float))
    b=(N/2-u)*old.smooth_step(u/ETA)
    return np.sqrt(b*b+m2)

def sturm_count(diagonal,off,value):
    # Number of eigenvalues below value by LDL^T pivots.
    floor=1e-200
    pivot=float(diagonal[0]-value)
    if abs(pivot)<floor:pivot=-floor
    count=int(pivot<0)
    off2=off*off
    for d in diagonal[1:]:
        pivot=float(d-value)-off2/pivot
        if abs(pivot)<floor:pivot=-floor
        count+=int(pivot<0)
    return count

def low_eigenvalues(diagonal,off,count=3):
    low=float(np.min(diagonal)-2*abs(off))
    high=float(np.max(diagonal)+2*abs(off))
    out=[]
    for index in range(count):
        lo=low;hi=high
        for _ in range(100):
            mid=(lo+hi)/2
            if sturm_count(diagonal,off,mid)<=index:lo=mid
            else:hi=mid
            if hi-lo<2e-12*max(1,abs(mid)):break
        out.append((lo+hi)/2)
    return np.array(out)

def spectrum(N,m2,particles=2,points=2047,levels=1):
    h=2*HALF_INTERVAL/(points+1)
    x=np.linspace(-HALF_INTERVAL+h,HALF_INTERVAL-h,points)
    potential=particles*energy(N,x,m2)
    d=2*KAPPA/h**2+potential
    e=low_eigenvalues(d,-KAPPA/h**2,levels)
    return dict(N=N,particles=particles,grid_points=points,
        energies=e.tolist(),ground_over_log_squared=float(e[0]/np.log(N)**2),
        zero_kinetic_threshold=float(particles*np.sqrt(m2)))

def matrix_identity(m,q,m2):
    rows=[]
    for N in (9,33,129):
        for x in (-.013,-.002,.0,.002,.013):
            ids=[26,27,28,29]
            h=old.H(m,q,N,[N//2,0,0],CENTER+x)[np.ix_(ids,ids)]
            eigen=np.linalg.eigvalsh(h)
            expected=float(energy(N,x,m2))
            err=float(np.max(abs(np.abs(eigen)-expected)))
            assert err<1e-10
            rows.append(dict(N=N,offset=x,original_electron_energy=expected,
                             full_block_spectral_error=err))
    return rows

def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.last.em_charges(m,matter)
    assert abs(float(np.sum(q*q))-768)<1e-12
    m2=float(abs(np.exp(old.XI)*m.MASS[26,28])**2)
    identity=matrix_identity(m,q,m2)
    # Independently verify the tridiagonal solver on a constant potential.
    points=63;h=2*HALF_INTERVAL/(points+1);v=2*np.sqrt(m2)
    d=np.full(points,2*KAPPA/h**2+v)
    computed=low_eigenvalues(d,-KAPPA/h**2,3)
    exact=np.array([v+4*KAPPA/h**2*np.sin(j*np.pi/(2*(points+1)))**2 for j in (1,2,3)])
    assert np.max(abs(computed-exact))<1e-8
    rows=[spectrum(N,m2) for N in (17,65,257,4097,1048577,1073741825)]
    refinement=[spectrum(1048577,m2,points=n,levels=3) for n in (511,1023,2047,4095)]
    diffs=[abs(refinement[j+1]['energies'][0]-refinement[j]['energies'][0]) for j in range(3)]
    assert diffs[-1]<diffs[-2]/3 and diffs[-2]<diffs[-3]/3
    assert diffs[-1]/refinement[-1]['energies'][0]<2e-4
    coeff=KAPPA*(6/ETA)**2*np.pi**2/4
    # Analytic rescaled potential demonstrates the hard-wall limit separately
    # from finite-difference eigenvalues; huge integers need no large matrices.
    scales=[]
    for power in (10,30,60,120):
        N=2**power+1;s=np.log(float(N))
        values=[float(2*energy(float(N),ETA*y/(6*s),m2)/s**2) for y in (.75,1.,1.25)]
        scales.append(dict(N=str(N),logN=float(s),scaled_potential_at_y_075_1_125=values))
    assert scales[-1]['scaled_potential_at_y_075_1_125'][0]<.001
    assert scales[-1]['scaled_potential_at_y_075_1_125'][2]>1e3
    return dict(round=889,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3674,
        argument_scope='For the declared finite-cutoff transported kinetic and pointwise vacuum-normal-ordered branch, the original neutral electron-hole cut sector has actual dynamic eigenvalues asymptotic to kappa*(6/eta)^2*(j*pi/2)^2*(log N)^2. This removes its fixed-energy spurious excitations, without assuming pointwise conditional Gibbs. It is a sector theorem, not full Gauss or Q/E equivalence.',
        kappa=KAPPA,eta=ETA,coordinate_half_interval=HALF_INTERVAL,
        boundary_condition='Dirichlet on fixed interval about alpha=-1/12',
        original_electron_mass_squared=m2,asymptotic_ground_coefficient=coeff,
        original_matrix_checks=identity,actual_dynamic_rows=rows,mesh_refinement=refinement,
        refinement_ground_differences=diffs,
        independent_constant_potential_solver_error=float(np.max(abs(computed-exact))),
        rescaled_potential_checks=scales,
        selected_pair_total_EM_charge=0,selected_pair_even_parity=True,
        full_original_species_retained_in_declared_operator=True,
        added_transported_kinetic_rule=True,added_pointwise_vacuum_subtraction=True,
        own_joint_Gibbs_diagonal_assumed_pointwise_free=False,
        local_Gauss_constraint_solved=False,full_thermal_entropy_and_source_convergence=False,
        original_Q_E_equivalence_proved=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
