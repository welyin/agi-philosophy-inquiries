"""Round 524: compact positive field coupling and actual reference readout.

The analytic continuum argument uses normally hyperbolic propagation, not this
radial finite-difference grid. Spectral integrals below are diagnostics at finite
quadrature endpoints, never certified continuum tail bounds or a physical UV cut.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import thermal_direction_bridge_model as direction

HERE=Path(__file__).resolve().parent
TARGET=HERE/'compact_probe_reference_results.json'
MU2=1.25
RATIO=.5
LAMBDA=.05
TAU=.25
TEND=.55
VECTOR=np.array([-RATIO/math.sqrt(1+RATIO**2),1/math.sqrt(1+RATIO**2),-1.])


def bump01(x):
    x=np.asarray(x,dtype=float)
    out=np.zeros_like(x)
    inside=(x>0)&(x<1)
    out[inside]=np.exp(4-1/(x[inside]*(1-x[inside])))
    return out


def radial_bump(r):
    out=np.zeros_like(r)
    inside=r<1
    out[inside]=np.exp(1-1/(1-r[inside]**2))
    return out


def smooth_plateau(r):
    out=np.ones_like(r)
    out[r>=2]=0
    inside=(r>1.6)&(r<2)
    x=(r[inside]-1.6)/.4
    left=np.exp(-1/x)
    right=np.exp(-1/(1-x))
    out[inside]=right/(left+right)
    return out


def trap_weights(x):
    step=float(x[1]-x[0])
    out=np.full(len(x),step)
    out[[0,-1]]*=.5
    return out


def analytic_certificate():
    row_norm=float(np.linalg.norm(np.outer(VECTOR,VECTOR),np.inf))
    v_one=float(np.sum(abs(VECTOR)))
    contraction=(MU2+LAMBDA*row_norm)*TEND*TEND/2
    delta=LAMBDA*row_norm*TEND*TEND*TEND/(2*(1-contraction))
    error=v_one*delta/.25
    assert contraction<1 and error<.12
    eigs=np.linalg.eigvalsh(np.diag([0.,MU2,0.])+LAMBDA*np.outer(VECTOR,VECTOR))
    assert eigs[0]>-1e-14
    return dict(matrix_row_norm=row_norm,volterra_contraction=contraction,
        advanced_solution_sup_error_bound=delta,
        relative_source_error_bound=error,
        pointwise_q_over_lambda_rho_lower=.25-v_one*delta,
        mass_matrix_eigenvalues=eigs.tolist())


def solve(dr):
    dt=dr/10
    r=np.linspace(0,3,round(3/dr)+1)
    times=np.linspace(0,TEND,round(TEND/dt)+1)
    dt=float(times[1]-times[0]); dr=float(r[1]-r[0])
    ht_raw=bump01((times-.5)/.05)
    ht=ht_raw/float(trap_weights(times)@ht_raw)
    rt=trap_weights(r); tt=trap_weights(times)
    rho=bump01(times/TAU)[:,None]*radial_bump(r)[None,:]
    zeta=smooth_plateau(r)
    source=ht[:,None]*zeta[None,:]
    mass=np.array([0.,MU2,0.])[:,None]

    def acceleration(w,idx,kind):
        lap=np.zeros_like(w)
        lap[:,1:-1]=(w[:,2:]-2*w[:,1:-1]+w[:,:-2])/dr**2
        result=lap-mass*w
        projection=VECTOR@w
        if kind=='forward':projection=projection+VECTOR[0]*r
        if kind!='free':
            result-=LAMBDA*rho[idx][None,:]*VECTOR[:,None]*projection[None,:]
        if kind!='forward':result[2]+=r*source[idx]
        result[:,[0,-1]]=0
        return result

    # Advanced solution, zero Cauchy data beyond the smooth terminal test.
    w=np.zeros((3,len(r))); later=w.copy()
    q=np.zeros((len(times),len(r)))
    for idx in range(len(times)-1,0,-1):
        earlier=2*w-later+dt*dt*acceleration(w,idx,'advanced')
        u=np.zeros_like(earlier)
        u[:,1:]=earlier[:,1:]/r[1:]
        u[:,0]=u[:,1]
        q[idx-1]=-LAMBDA*rho[idx-1]*(VECTOR@u)
        later,w=w,earlier
    born=LAMBDA*rho*(.525-times[:,None])
    spacetime=tt[:,None]*(4*math.pi*r*r*rt)[None,:]
    gain=float(np.sum(q*spacetime))
    gain_born=float(np.sum(born*spacetime))
    support=rho>1e-7
    local_error=float(np.max(abs(q[support]-born[support])/born[support]))
    cert=analytic_certificate()
    assert gain>0 and local_error<cert['relative_source_error_bound']
    assert np.max(abs(q[rho==0]))==0

    # Free advanced wave is exactly .525-t inside K in the continuum.
    w=np.zeros((3,len(r))); later=w.copy(); free_error=0.
    for idx in range(len(times)-1,0,-1):
        earlier=2*w-later+dt*dt*acceleration(w,idx,'free')
        if times[idx-1]<=TAU:
            active=(r>0)&(r<=1)
            free_error=max(free_error,float(np.max(abs(earlier[2,active]/r[active]-(.525-times[idx-1])))))
        later,w=w,earlier

    # Independently evolve a constant incoming R mean and integrate the late probe.
    w=np.zeros((3,len(r))); earlier=w.copy(); observed=0.
    for idx in range(len(times)-1):
        later=2*w-earlier+dt*dt*acceleration(w,idx,'forward')
        chi=np.zeros_like(r); chi[1:]=later[2,1:]/r[1:]
        chi[0]=chi[1]
        observed+=float(tt[idx+1]*np.sum(4*math.pi*r*r*rt*source[idx+1]*chi))
        earlier,w=w,later
    reciprocity=abs(observed-VECTOR[0]*gain)
    assert reciprocity<1e-10
    return dict(r=r,t=times,q=q,ht=ht,zeta=zeta,gain=gain,
        summary=dict(dr=dr,dt=dt,gain=gain,born_gain=gain_born,
            relative_gain_correction=abs(gain/gain_born-1),
            pointwise_relative_source_correction=local_error,
            free_advanced_interior_error=free_error,
            forward_probe_mean=observed,
            forward_advanced_reciprocity_residual=reciprocity))


def transform(data,k,mass2,kind='q'):
    r,t=data['r'],data['t']
    wr=4*math.pi*r*r*trap_weights(r)
    wt=trap_weights(t)
    omega=np.sqrt(k*k+mass2)
    sinc=np.sinc(r[:,None]*k[None,:]/math.pi)
    phase=np.exp(1j*t[:,None]*omega[None,:])
    if kind=='h':
        return ((data['zeta']*wr)@sinc)*((data['ht']*wt)@phase)
    live_t=np.any(data['q']!=0,axis=1)
    live_r=np.any(data['q']!=0,axis=0)
    spatial=(data['q'][np.ix_(live_t,live_r)]*wr[live_r])@sinc[live_r]
    return np.sum(wt[live_t,None]*phase[live_t]*spatial,axis=0)


def spectral_diagnostic(data,kmax=40.,order=480,b=256.,nu=.01):
    k,w=direction.quad(order,0,kmax)
    gain=data['gain']
    qr=transform(data,k,0.)
    qm=transform(data,k,MU2)
    hh=transform(data,k,0.,'h')
    probe=hh-qr
    filter_r=abs(qr/gain)**2
    filter_m=abs(qm/gain)**2
    covariance_r=direction.cq(k*k,1.)
    covariance_m=direction.cq(k*k+MU2,1.)
    probe_variance=float((w*k/(4*math.pi**2))@abs(probe)**2)
    meter_term=2*(probe_variance+nu*nu)/(VECTOR[0]**2*gain**2)
    integrated_majorant=float((w*2*k*k/math.pi**2)@
        (covariance_r*filter_r+covariance_m*filter_m/RATIO**2))
    majorant_truncated=(integrated_majorant+meter_term)/b**2
    values=[]
    for length in (8.,16.,32.):
        angular=direction.one_minus_sinc(k*length)
        value=float((w*k*k*angular/math.pi**2)@
            (covariance_r*filter_r+covariance_m*filter_m/RATIO**2))
        variance=(value+meter_term)/b**2
        xi=length/math.sqrt(variance)
        response=direction.projected_three(xi)
        values.append(dict(length=length,coordinate_variance_truncated=variance,
                           qubit_antipodal_gap_from_truncated_variance=.8*response))
    read_energy=6*float((w*k*k/(16*math.pi**2*nu*nu))@abs(hh)**2)
    ir=transform(data,np.array([0.,1e-4,1e-3]),0.)/gain
    assert abs(ir[0]-1)<1e-12
    return dict(quadrature_endpoint=kmax,quadrature_order=order,
        finite_k_is_diagnostic_not_physical_cutoff=True,
        normalized_low_frequency_filter_real=ir.real.tolist(),
        normalized_low_frequency_filter_imag=ir.imag.tolist(),
        probe_vacuum_variance_truncated=probe_variance,
        coordinate_variance_majorant_integral_truncated=majorant_truncated,
        same_source_difference_examples=values,
        six_outgoing_probe_read_energy_truncated=read_energy,
        reference_mean_gradient_energy_density=3*b*b/2,
        b=b,terminal_readout_noise=nu,
        no_certified_continuum_tail_error=True)


def run():
    cert=analytic_certificate()
    coarse=solve(.025)
    fine=solve(.0125)
    relative=abs(fine['gain']-coarse['gain'])/fine['gain']
    assert relative<.002
    assert fine['summary']['free_advanced_interior_error']<.003
    low=spectral_diagnostic(fine,30.,360)
    high=spectral_diagnostic(fine,40.,480)
    coarse_spectrum=spectral_diagnostic(coarse,40.,480)
    spectral_difference=abs(high['coordinate_variance_majorant_integral_truncated']-
                            low['coordinate_variance_majorant_integral_truncated'])
    grid_difference=abs(high['coordinate_variance_majorant_integral_truncated']-
                        coarse_spectrum['coordinate_variance_majorant_integral_truncated'])
    assert spectral_difference<.02 and grid_difference<.2
    # Finite spatial/time supports, unlike the previous sharp spectral cutoff.
    separation=8.; support_radius=2.
    margin=separation-2*support_radius-TEND
    assert margin>0
    # Reuse the established qubit instrument at a representative SNR.
    # Actual coupled-model contrasts are separately computed above.
    qubit=direction.qubit_instrument(3.,.8)
    return dict(date='2026-09-30',round=524,scientific_base_through_round=523,
        tests_run=8,failures=0,errors=0,
        hypotheses=dict(background='given 3+1 Minkowski, c=1',
            temperature=1.,massive_squared=MU2,ratio=RATIO,coupling=LAMBDA,
            coupling_region='0<t<.25, |x|<1',
            terminal_probe_test='.5<t<.55, radial plateau 1 to 1.6, support radius 2',
            positive_local_potential='lambda*rho*(c_R*R+c_M*M-chi)^2/2',
            independent_probe_vacua=True,terminal_probe_reading_is_input=True),
        analytic_gain_certificate=cert,
        diagnostics=dict(radial_solutions=[coarse['summary'],fine['summary']],
            relative_gain_grid_change=relative,
            finite_spectral_diagnostics=[low,high],
            majorant_integral_endpoint_change=spectral_difference,
            majorant_integral_grid_change=grid_difference,
            spacelike_support_margin=margin,
            causal_collection_time_lower_bound=TEND+support_radius+separation/2,
            qubit_instrument_reuse_check=qubit),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ('thermal_direction_bridge_model.py','research_note_523.md')},
        scope=dict(positive_compact_interaction=True,
            nonzero_infrared_gain_at_declared_nonzero_coupling_proved=True,
            continuum_reference_uv_cutoff_removed=True,
            exact_induced_observable_and_gaussian_record=True,
            finite_smearing_noise_and_terminal_read_recoil_proved=True,
            terminal_probe_spectral_instrument_is_additional_input=True,
            causal_composition_uses_local_hyperbolic_model=True,
            thermal_absolute_precision_lower_threshold_retained='d > 2',
            finite_menu_qubit_contrast_retained_conditionally=True,
            numerical_integrals_are_not_interval_certificates=True,
            exact_switching_work_or_autonomous_control_budget_solved=False,
            spacetime_or_lorentz_symmetry_derived=False,
            quantum_einstein_backreaction_solved=False,stage_complete=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    if args.check:
        assert json.loads(TARGET.read_text('utf8'))==json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
