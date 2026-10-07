"""805: exact Fourier multiplier of the original 730 auxiliary past.

All original 64 Nambu mass/principal matrices retained. The multiplier is
the continuum auxiliary-state operator, not an instantaneous target vacuum.
No numerical future Einstein evolution or bosonic covariance is supplied.
"""
from pathlib import Path
from functools import lru_cache
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'802'))
import original_bff_vertex as vertex
TARGET=HERE/'original_past_covariance_action_results.json'
PHI=np.array([0.,vertex.original.old.CENTER[0],0.,0.,vertex.original.old.CENTER[1]])
MASS=vertex.mass(PHI)
GAMMA=vertex.original.GAMMA
I=np.eye(64,dtype=complex)
CHARGE=np.block([[np.zeros((32,32)),np.eye(32)],[np.eye(32),np.zeros((32,32))]])
EV,VEC=np.linalg.eigh(MASS@MASS)

@lru_cache(maxsize=None)
def root_inverse(radius_squared):
    return (VEC/np.sqrt(EV+radius_squared))@VEC.conj().T

def hamiltonian(k):
    return MASS+sum((GAMMA[a]*k[a] for a in range(3)),np.zeros((64,64),complex))

def occupied(k):
    return .5*(I-hamiltonian(k)@root_inverse(float(np.dot(k,k))))

def apply_column(k,index=30):
    # Positive/greater P is the complement of the original occupied N.
    column=root_inverse(float(np.dot(k,k)))[:,index]
    unit=np.eye(64)[:,index]
    return .5*(unit+hamiltonian(k)@column)

def run():
    errors=dict(clifford_mass_anticommutation=max(float(np.max(abs(g@MASS+MASS@g))) for g in GAMMA))
    rows=[]
    for k in ([0,0,0],[1,0,0],[0,-2,1],[2,3,-1],[7,-4,3],[-.31,.27,-.19]):
        k=np.array(k,dtype=float);h=hamiltonian(k);n=occupied(k)
        val,vec=np.linalg.eigh(h);direct=vec[:,val<0]@vec[:,val<0].conj().T
        square_error=float(np.max(abs(h@h-(np.dot(k,k)*I+MASS@MASS))))
        projection_error=float(np.max(abs(n@n-n)))
        reality_error=float(np.max(abs(CHARGE@n.conj()@CHARGE+occupied(-k)-I)))
        spectral_error=float(np.max(abs(n-direct)))
        rows.append(dict(momentum=k.tolist(),gap=float(min(abs(val))),square_error=square_error,
            projection_error=projection_error,reality_error=reality_error,spectral_formula_error=spectral_error))
        assert max(square_error,projection_error,reality_error,spectral_error)<2e-11
    assert errors['clifford_mass_anticommutation']<1e-12 and EV.min()>0
    # Explicit smooth (not compactly supported in a proper patch) Fourier input.
    # This calibrates the universal L2 contraction/tail estimate only.
    q=.3;s=(1+q*q)/(1-q*q);norm=s**(-1.5)
    energy_by_shell=np.zeros(7)
    for x in range(-6,7):
        for y in range(-6,7):
            for z in range(-6,7):
                k=np.array([x,y,z]);col=apply_column(k)*norm*q**sum(abs(k))
                energy_by_shell[max(abs(k))]+=float(np.vdot(col,col).real)
    tails=[]
    for cut in (0,1,2,4):
        finite_sum=1+2*q*q*(1-q**(2*cut))/(1-q*q)
        exact_input_tail=np.sqrt(max(0,1-(finite_sum/s)**3))
        measured_output_tail=np.sqrt(sum(energy_by_shell[cut+1:]))
        assert measured_output_tail<=exact_input_tail+2e-14
        tails.append(dict(cube_half_width=cut,input_tail_exact=float(exact_input_tail),
            output_tail_through_shell_6=float(measured_output_tail)))
    return dict(round=805,all_checks_passed=True,original_scalar_reference=PHI.tolist(),
        original_Nambu_dimension=64,original_CAR_modes=32,auxiliary_past_gap=float(np.sqrt(EV.min())),
        coefficient_errors=errors,momentum_checks=rows,Fourier_tail_checks=tails,
        continuum_multiplier_formula='N_-(k)=(I-(Gamma.k+M)(|k|^2 I+M^2)^(-1/2))/2; P_-=I-N_-',
        all_original_species_and_masses_retained=True,
        formula_covers_all_allowed_torus_momenta=True,
        test_input_is_original_local_record=False,
        original_future_PDE_solved=False,original_bosonic_covariance_evaluated=False,
        actual_target_covariance_is_instantaneous_projector=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
