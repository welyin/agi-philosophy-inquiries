"""Round1203: fixed autonomous contact + retained physical pointer certificate.

No scattering-grid simulation: exact rational sufficient bounds transport the
previous certified native contact to a normal positive-energy clock. NumPy
checks the finite material/pointer comparison; it is diagnostic, not the
certificate for the infinite-dimensional clock.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import hashlib, json
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=STAGE.parents[1]

def read(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def up(q,places=12):
    k=10**places
    return F((q.numerator*k+q.denominator-1)//q.denominator,k)
def dq(q): return Decimal(q.numerator)/Decimal(q.denominator)
def expm(H,t):
    e,u=np.linalg.eigh(H)
    return (u*np.exp(-1j*t*e))@u.conj().T

def main():
    old=read(STAGE/'1201/results.json'); prev=read(STAGE/'1202/results.json')
    # Exact upper pi bracket; no fitted clock parameters.
    pi_hi=F(22,7)
    T=F(606); mass=F(10**24); sigma=F(1,10**10); g=F(1,10)
    width=F(547,500)  # > max(e)-min(e) from prior root intervals
    spread=sigma*(1+T*T/(8*mass*mass*sigma**4))
    proof={}
    for v in (1,2):
        lipschitz=F(150)+F(3,2)*v*pi_hi
        clock=lipschitz*T*spread
        edge=F(1,10000*v)
        eps=clock+edge
        backreaction=lipschitz*T*T/(2*mass)
        end=lipschitz*(spread+backreaction)
        start=lipschitz*sigma
        residual=start+end
        # All finite read interaction included in residual; initial not zero.
        bare_error=2*width*eps
        clock_work_error=bare_error+residual+F(1,10**10)
        oldtime=str(600//v)
        work=F(str(prev['native_contacts']['energy_exchange'][oldtime]['switch_work_coefficients']['expectation']))
        proof[str(v)]={
            'v':v,'time':606,'contact_time':600//v,
            'lipschitz_upper_exact':str(lipschitz),
            'clock_error_exact_upper':str(clock),'edge_error_exact':str(edge),
            'isometry_error_exact_upper':str(eps),'isometry_error_upper':str(up(eps)),
            'initial_plus_final_interaction_energy_exact_upper':str(residual),
            'initial_plus_final_interaction_energy_upper':str(up(residual)),
            'single_material_energy_error_upper':str(up(width*eps)),
            'clock_energy_prediction_error_upper':str(up(clock_work_error)),
            'clock_delta_energy_interval':[float(-work-clock_work_error),float(-work+clock_work_error)],
            'clock_energy_loss_lower':float(work-clock_work_error),
            'initial_clock_mean_energy_exact':str(mass*v*v/2+1/(8*mass*sigma*sigma)),
            'initial_clock_energy_variance_exact':str(F(v*v)/(4*sigma*sigma)+1/(32*mass**2*sigma**4)),
            'finite_unknown_material_reference_bound':True}
        assert work-clock_work_error>0
    epssum=sum(F(proof[str(v)]['isometry_error_exact_upper']) for v in (1,2))
    eps=up(epssum)
    # Decimal log evaluation with generous explicit final numerical allowance.
    with localcontext() as ctx:
        ctx.prec=80
        x=dq(eps)
        h2=-x*x.ln()-(1-x)*(1-x).ln()
        correction=h2+x*Decimal(255).ln()+Decimal('1.2')*dq(2*width)*x
        loss=correction+Decimal('.000912')
        loose_up=(loss*Decimal(10)**12).to_integral_value(rounding='ROUND_CEILING')/(Decimal(10)**12)+Decimal('1e-12')
        clock_entropies=[]
        for v in (1,2):
            amp=dq(F(proof[str(v)]['isometry_error_exact_upper']))
            q=amp*amp
            s=-q*q.ln()-(1-q)*(1-q).ln()+q*Decimal(255).ln()
            clock_entropies.append(s)
        entropy_loss=sum(clock_entropies)+Decimal('.000912')
        entropy_up=(entropy_loss*Decimal(10)**12).to_integral_value(rounding='ROUND_CEILING')/(Decimal(10)**12)+Decimal('1e-12')
        entropy_strings=[str((s+Decimal('1e-60')).quantize(Decimal('1e-30'),rounding='ROUND_CEILING')) for s in clock_entropies]

    # Elementary polynomial gate identities, independent of numerical quadrature.
    assert 3*F(1,3)-2*F(1,4)==F(1,2)  # integral q on [0,1]
    assert 6*F(1,2)*(1-F(1,2))==F(3,2) # maximum q'
    assert 2*(1-F(1,2))*F(1,1000)==F(1,1000)

    h=np.array([[0,-.2,0,-np.sqrt(3)/80],[-.2,1,0,0],[0,0,1,0],[-np.sqrt(3)/80,0,0,-.025]])
    e,basis=np.linalg.eigh(h)
    Q=np.zeros((4,4));Q[1,2]=Q[2,1]=1
    Q=basis.T@Q@basis
    h=np.diag(e); eye4=np.eye(4)
    H0=np.kron(h,eye4)+np.kron(eye4,h)
    V=np.kron(Q,Q); H=H0+.1*V
    G=np.zeros((64,64)); copy=np.zeros((64,64))
    for j in range(4):
        p=np.zeros((4,4));p[j,j]=1
        swap=np.eye(4)
        if j:
            swap[0,0]=swap[j,j]=0;swap[0,j]=swap[j,0]=1
        G+=np.kron(np.kron(p,eye4),np.pi/2*(eye4-swap))
        copy+=np.kron(np.kron(p,eye4),swap)
    H0K=np.kron(H0,eye4)
    assert np.linalg.norm(G@H0K-H0K@G)<1e-12
    assert np.linalg.eigvalsh(G).min()>-1e-12
    assert np.linalg.norm(expm(G,1)-copy,2)<1e-12
    embed=np.zeros((64,16))
    for i in range(16): embed[4*i,i]=1
    matrices={};finite={}
    for v in (1,2):
        t=600/v
        u=expm(H,t)
        # Free segments are kept as operators, not deleted for arbitrary BR.
        ufull=copy@np.kron(expm(H0,606-601/v)@u@expm(H0,1/v),eye4)@embed
        np.testing.assert_allclose(ufull.conj().T@ufull,np.eye(16),atol=2e-12)
        pointer=np.zeros((4,4))
        for i in range(4):
            vec=ufull[:,2*4+i].reshape(4,4,4)
            pointer[:,i]=np.sum(np.abs(vec)**2,axis=(0,1))
        native=np.array(old['contacts'][str(int(t))]['K'],dtype=float)
        assert np.max(abs(pointer-native))<2e-11
        matrices[str(int(t))]=pointer
        finite[str(v)]={'material_pointer_dimension':64,'copy_norm_error':float(np.linalg.norm(expm(G,1)-copy,2)),
            'pointer_kernel_vs_certified_native_max_error':float(np.max(abs(pointer-native))),
            'unconditional_clock_evolution_simulated':False}

    R=np.array([[float(F(q)) for q in row] for row in prev['two_marginal_rational']])
    joint=matrices['600']@R@matrices['300'].T
    oldjoint=np.array(prev['native_contacts']['joint_outputs'])
    assert np.max(abs(joint-oldjoint))<2e-11
    # Certified table errors add separately to analytic physical-clock bounds.
    joint_tv=up(epssum+8*F(3,10**10))
    corr_lower=F('.00069502')-2*epssum
    inverse_error=F('.015701')+F('3.14')*F(proof['1']['isometry_error_exact_upper'])
    out={'schema':'autonomous_contact_record_1203_v1','status':'PASS_SCOPED_CLOCK_RECORD_CONNECTION',
        'common_terminal_time':606,'one_fixed_H_with_program_Z':[1,2],
        'clock_mass':str(mass),'clock_sigma':str(sigma),'bounds_by_program':proof,
        'combined_isometry_error_exact_upper':str(epssum),
        'combined_isometry_error_upper':str(up(epssum)),
        'pointer_joint_probabilities':joint.tolist(),'pointer_joint_table_TV_error_upper':str(joint_tv),
        'correlation_prediction_TV_lower_after_two_clock_errors':float(corr_lower),
        'correlation_experimental_resolution_claimed':False,
        'population_identification_error_with_original_budget_upper':str(up(inverse_error)),
        'thermodynamics':{'retained_AKC_clausius_lower':'-.000912',
            'AKC_entropy_rank_bound':256,'each_clock_rank_bound':256,'AK_dimension':256,
            'AK_only_clausius_lower':str(-entropy_up),
            'each_clock_entropy_upper_decimal':entropy_strings,
            'clock_pure_comparison_failure_is_amplitude_squared':True,
            'looser_finite_AK_continuity_lower':str(-loose_up),
            'beta':'6/5','initial_clock_pure_and_independent':True,
            'all_clocks_records_source_reference_retained':True,
            'clocks_state_recovery_or_free_reset_claimed':False},
        'finite_matrix_diagnostic':finite,
        'remaining_inputs':['1202 successful reference preparation and full stock/seed resources',
            'normal positive-energy clock degree and specified Gaussian preparation',
            'same fixed smooth gates and energy-projector record coupling',
            'pure degenerate pointer, correlated known Z/mean-velocity preparation',
            'finite terminal window and pointer interrogation; additional physical errors separate',
            'actual spatial region, causal support and geometric identity'],
        'not_claimed':['autonomous full reference/preparation computation',
            'direct physical measurement of the tiny change of the huge clock energy',
            'physical B-energy pointer in this model','attracting thermalization',
            'universal temperature-flow ordering','SR clock mass or 3D area identity','Unruh or GR'],
        'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in
            [STAGE/'1201/results.json',STAGE/'1201/acceptance.json',STAGE/'1202/results.json',STAGE/'1202/acceptance.json']}}
    (HERE/'results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'combined_error':str(up(epssum)),
       'AK_clausius_lower':str(-entropy_up),'clock_losses_positive':True}))

if __name__=='__main__':main()
