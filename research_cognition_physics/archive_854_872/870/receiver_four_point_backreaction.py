"""870: exact finite CAR certificate and ordered receiver backreaction.
The continuous difference-kernel theorem is in the report. Matrix oscillators
only calibrate its polynomial coefficient on a non-boundary vacuum sector.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'869'))
import direct_material_receiver as previous
TARGET=HERE/'receiver_four_point_backreaction_results.json'
I=previous.I;X=previous.X;Y=previous.Y;Z=previous.Z

def annihilators(n):
    a=np.array([[0,1],[0,0]],complex)
    return [previous.kronall(*([Z]*i+[a]+[I]*(n-i-1))) for i in range(n)]

def receiver_state(r,m=.5):
    out=np.zeros((4,4),complex)
    out[0,0]=out[3,3]=(1-r)/2
    out[1:3,1:3]=(r*I+m*X)/2
    assert np.linalg.eigvalsh(out).min()>-1e-14 and abs(np.trace(out)-1)<1e-14
    return out

def comm(a,b):return a@b-b@a

def run():
    modes=annihilators(4)
    def second(A):return sum((modes[i].conj().T@ modes[j])*A[i,j] for i in range(4) for j in range(4))
    # Two retained modes plus two mode-exterior CAR modes: none dropped in O.
    ext=np.diag([.8*.35,.8*.65,.2*.35,.2*.65]).astype(complex)
    contrast=.5;r0=.5;r1=1.
    rho0=np.kron(receiver_state(r0),ext);rho1=np.kron(receiver_state(r1),ext)
    drho=rho1-rho0;delta=r1-r0
    def expect(rho,A):return np.trace(rho@A)
    c0=np.array([[expect(rho0,modes[i].conj().T@modes[j]) for j in range(4)] for i in range(4)])
    c1=np.array([[expect(rho1,modes[i].conj().T@modes[j]) for j in range(4)] for i in range(4)])
    twoerr=float(np.max(abs(c1-c0)));assert twoerr<1e-15
    B=[np.array([[1.,.23+.11j],[.23-.11j,.7]]),
       np.array([[.8,-.14+.09j],[-.14-.09j,1.2]]),
       np.array([[1.1,.08-.13j],[.08+.13j,.9]])]
    currents=[second(np.kron(b,Y)) for b in B]
    densities=np.array([b[0,0].real for b in B])
    record=2*modes[0].conj().T@modes[0]-np.eye(16)
    chi=np.array([1j*expect(rho0,comm(o,record)) for o in currents])
    assert np.max(abs(chi+2*contrast*densities))<1e-15
    assert max(abs(expect(drho,comm(o,record))) for o in currents)<1e-15
    differences=np.array([[expect(drho,oa@ob) for ob in currents] for oa in currents])
    noiseerr=float(np.max(abs(differences-delta*np.outer(densities,densities))))
    commerr=float(max(abs(expect(drho,comm(oa,ob))) for oa in currents for ob in currents))
    assert noiseerr<2e-15 and commerr<2e-15
    # Direct response calibration alone misses this non-Gaussian difference.
    expected_main_cov=np.array([[.5,contrast/2],[contrast/2,.5]])
    assert np.max(abs(c0[:2,:2]-expected_main_cov))<1e-15
    onebody=np.array([[2.,.2j,.1,0],[-.2j,2.,0,.2],[.1,0,3.,.15j],[0,.2,-.15j,3.]],complex)
    energyerr=float(abs(expect(drho,second(onebody))));assert energyerr<2e-15
    # Ordered three-pulse coefficient, retaining full receiver currents.
    n=7;a=np.zeros((n,n),complex)
    for i in range(1,n):a[i-1,i]=np.sqrt(i)
    Q=(a+a.conj().T)/np.sqrt(2);P=(a-a.conj().T)/(1j*np.sqrt(2))
    quadratures=[.7*Q+.2*P,-.3*Q+.6*P,.4*Q-.1*P]
    y=.3*Q+.8*P
    sigmas=np.array([.7*.8-.2*.3,-.3*.8-.6*.3,.4*.8+.1*.3])
    hs=[np.kron(aa,oo) for aa,oo in zip(quadratures,currents)]
    obs=np.kron(y@y,np.eye(16))
    second_coeff=np.zeros_like(obs)
    for t,h in enumerate(hs):
        second_coeff-=.5*comm(h,comm(h,obs))
        for s in range(t):second_coeff-=comm(hs[s],comm(h,obs))
    vacuum=np.zeros((n,n),complex);vacuum[0,0]=1
    got=expect(np.kron(vacuum,drho),second_coeff)
    wanted=delta*float(sigmas@densities)**2
    orderederr=float(abs(got-wanted))
    assert orderederr<3e-14 and wanted>.001
    # Check receiver first response is unchanged for the full pulse word.
    first_record=sum(1j*comm(h,np.kron(np.eye(n),record)) for h in hs)
    firsterr=float(abs(expect(np.kron(vacuum,drho),first_record)))
    assert firsterr<2e-14
    C,p,values,tangent=previous.original_menu()
    alpha=2;c=previous.GAINS[alpha];d=float(tangent[alpha])
    coefficient=delta*c*c*d*d/(4*contrast*contrast)
    assert coefficient>0
    quasifree_r=(1+contrast**2)/2
    rhoq=np.kron(receiver_state(quasifree_r),ext)
    # All three share covariance. Wick would collapse their distinct variance.
    local_y=second(np.kron(np.diag([1.,0]),Y))
    variances=[float(expect(rr,local_y@local_y).real) for rr in (rho0,rhoq,rho1)]
    assert np.allclose(variances,[r0,quasifree_r,r1])
    return dict(round=870,date='2026-10-06',formal_reports=870,
        cumulative_numbered_groups=3655,fresh_numbered_groups=1,all_checks_passed=True,
        contrast=contrast,occupation_correlation_parameters=[r0,r1],
        complete_initial_CAR_two_point_difference=twoerr,
        arbitrary_quadratic_energy_difference=energyerr,
        susceptibilities=chi.real.tolist(),
        full_current_noise_difference_error=noiseerr,
        full_current_commutator_difference=commerr,
        true_and_same_covariance_quasifree_current_variances=variances,
        ordered_three_pulse_variance_difference=float(got.real),
        analytic_three_pulse_variance_difference=wanted,
        ordered_coefficient_error=orderederr,
        unchanged_first_record_response_error=firsterr,
        original863_central_adjoint_tangent=d,
        original_companion_lambda2_hbar_variance_difference_coefficient=coefficient,
        original_smooth_anchor_nonzero='analytic from863 d_A != 0 and shared source calibration',
        continuum_noise_kernel='Delta N_R(x,y)=(r-r0) abar_u(x) abar_u(y); all exterior-mode terms cancel in the difference',
        continuum_outgoing_boson_covariance='[Delta H_out]_(lambda_alpha^2 hbar)=c_alpha^2 (r-r0)/(4 m^2) (E j_alpha) tensor (E j_alpha)',
        same_action_within_compared_family=True,
        same_complete_initial_free_two_point_and_first_menu=True,
        exact869_pure_plus_x_preparation_kept=False,
        new_species=0,all_state_costs_identical_claimed=False,
        finite_coupling_or_total_mean_source_solution=False,
        full_goal_completed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
