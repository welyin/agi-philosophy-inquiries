"""585: local-lapse completions of the original fixed-graph Hamiltonian.

Exact finite-graph form identities, an original Gauss packet, and a separate
continuum bracket diagnostic. No quantum Einstein solution is asserted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_quantum_measure_records as packet

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_local_lapse_source_results.json'
ALPHA = .4


def average(a, alpha=ALPHA):
    return (1-alpha)*a+alpha*sum(np.roll(a,k,axis=d)
        for d in range(3) for k in (-1,1))/6


def local_classical_energy(q, psi):
    """Endpoint/face equal allocation of every original Hamiltonian term."""
    lat=original.lattice;eps=q['eps'];phi=q['phi'];P=q['P']
    kinetic=np.einsum('...i,...ij,...j->...',P,original.inverse(phi),P)/(2*eps**3*psi**6)
    onsite=eps**3*psi**6*original.node_potential(phi)
    local=kinetic+onsite
    W=lat.su2(q['a'],eps);z=np.exp(1j*eps*q['a0'])
    bw,b0=lat.PAR['b'][1:]
    components=dict(scalar_kinetic=float(kinetic.sum()),onsite=float(onsite.sum()),
        gradient=0.,electric=0.,magnetic=0.)
    for mu in range(3):
        midpoint=(psi+np.roll(psi,-1,axis=mu))/2
        neighbor=z[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],
            np.roll(q['f']['X'],-1,axis=mu))
        chi=original.real_phi(neighbor,np.roll(q['f']['s'],-1,axis=mu))
        edge_g=eps*midpoint**2*original.distance_squared(phi,chi)/2
        edge_e=midpoint**-2*(bw*np.sum(q['phat'][...,mu,:]**2,axis=-1)
            +b0*q['p0hat'][...,mu]**2)/eps
        edge=edge_g+edge_e
        local+=(edge+np.roll(edge,1,axis=mu))/2
        components['gradient']+=float(edge_g.sum())
        components['electric']+=float(edge_e.sum())
    dg,coef=lat.old.coefficients(lat.PAR)
    for mu in range(3):
        for nu in range(mu+1,3):
            wm,wn=W[...,mu,:,:],W[...,nu,:,:]
            wc=wm@np.roll(wn,-1,axis=mu)@np.swapaxes(np.roll(wm,-1,axis=nu).conj(),-1,-2)@np.swapaxes(wn.conj(),-1,-2)
            zc=z[...,mu]*np.roll(z[...,nu],-1,axis=mu)/np.roll(z[...,mu],-1,axis=nu)/z[...,nu]
            tw=np.trace(wc,axis1=-2,axis2=-1)
            character=lat.PAR['wq']*(3*tw*zc+3*(zc**-4+zc**2))+lat.PAR['wl']*(tw*zc**-3+zc**6+1)
            vf=dg*(12*lat.PAR['wq']+4*lat.PAR['wl']-character.real)+coef[1]*(4-abs(tw)**2)+coef[2]*(1-(zc**6).real)
            center=(psi+np.roll(psi,-1,axis=mu)+np.roll(psi,-1,axis=nu)
                +np.roll(np.roll(psi,-1,axis=mu),-1,axis=nu))/4
            face=center**-2*vf/eps
            local+=(face+np.roll(face,1,axis=mu)+np.roll(face,1,axis=nu)
                +np.roll(np.roll(face,1,axis=mu),1,axis=nu))/4
            components['magnetic']+=float(face.sum())
    components['total']=sum(components.values())
    return local,components


def full_original_allocation_check():
    rows=[]
    for n in (4,6):
        q=original.shared_source(n)
        x,y,z=np.moveaxis(q['grid'],-1,0)
        psi=1.1+.08*np.sin(x)+.06*np.cos(y)+.03*np.sin(z)
        local,parts=local_classical_energy(q,psi)
        reference=original.graph_energy(q,psi)
        component_error=max(abs(parts[k]-reference[k]) for k in reference)
        assert component_error<2e-10
        smoothed=average(local)
        lapse=1+.2*np.cos(x)+.15*np.sin(y)
        dual_error=abs(float(np.sum(lapse*smoothed)-np.sum(average(lapse)*local)))
        total_error=abs(float(local.sum()-smoothed.sum()))
        assert min(local.min(),smoothed.min(),average(lapse).min())>0
        assert max(dual_error,total_error)<2e-10
        direction=.07*np.sin(x+2*y)+.03*np.cos(z)
        def total(t,smear):
            a,_=local_classical_energy(q,psi+t*direction)
            return float((average(a) if smear else a).sum())
        step=1e-4
        derivative_plain=(total(step,False)-total(-step,False))/(2*step)
        derivative_smear=(total(step,True)-total(-step,True))/(2*step)
        derivative_gap=abs(derivative_plain-derivative_smear)
        assert derivative_gap<3e-8
        local_gap=float(np.max(abs(smoothed-local)))
        assert local_gap>1e-4
        rows.append(dict(n=n,total_original=reference['total'],component_error=component_error,
            total_error=total_error,lapse_duality_error=dual_error,
            spatial_weight_derivative_gap=derivative_gap,
            max_local_source_difference=local_gap,
            nonconstant_lapse_energy_difference=float(np.sum(lapse*(smoothed-local)))))
    return dict(rows=rows,original_572_574_classical_source_used=True,
        full_quantum_form_identity_proved_analytically=True)


def packet_phase_cost(nodes):
    z,w=np.polynomial.legendre.leggauss(nodes)
    x,y=np.meshgrid(z,z,indexing='ij')
    h=packet.HCENTER+packet.HRADIUS*x;s=packet.SCENTER+packet.SRADIUS*y
    amplitude=np.exp(-1/(1-x*x)-1/(1-y*y))
    weights=w[:,None]*w[None,:]*packet.HRADIUS*packet.SRADIUS*h**3
    prob=weights*amplitude**2;prob/=prob.sum()
    f=original.M-(h*h+s*s)/6
    coordinates=np.stack((h,s),axis=-1)
    G=f[...,None,None]*(np.eye(2)-coordinates[..., :,None]*coordinates[...,None,:]/(6*original.M))
    grad_log=np.stack((-2*x/(packet.HRADIUS*(1-x*x)**2),
        -2*y/(packet.SRADIUS*(1-y*y)**2)),axis=-1)-coordinates/(2*f[...,None])
    phase_grad=np.stack((h,np.zeros_like(h)),axis=-1)
    theta=.3;hbar=packet.HBAR;weight=.8**3
    dpsi=grad_log+1j*theta*phase_grad/hbar
    before=np.einsum('...a,...ab,...b->...',grad_log,G,grad_log)
    after=np.einsum('...a,...ab,...b->...',dpsi.conj(),G,dpsi).real
    direct=float(hbar**2*np.sum(prob*(after-before))/(2*weight))
    expected=float(theta**2*np.sum(prob*f*(h*h-h**4/(6*original.M)))/(2*weight))
    assert direct>0 and abs(direct-expected)<2e-12
    return dict(nodes=nodes,theta=theta,hbar=hbar,node_volume=weight,
        direct_kinetic_increment=direct,analytic_kinetic_increment=expected,error=abs(direct-expected))


def same_gauss_state_local_source_check():
    rows=[packet_phase_cost(n) for n in (48,80,128)]
    delta=rows[-1]['analytic_kinetic_increment']
    assert abs(delta-rows[-2]['analytic_kinetic_increment'])<1e-12
    increment=np.zeros((5,5,5));increment[0,0,0]=delta
    changed=average(increment)
    assert abs(changed.sum()-delta)<1e-15
    assert abs(changed[0,0,0]-(1-ALPHA)*delta)<1e-15
    assert np.count_nonzero(changed)==7
    assert abs(changed[1,0,0]-ALPHA*delta/6)<1e-15
    x=np.arange(5)*2*np.pi/5
    lapse=np.broadcast_to(1+.2*np.cos(x)[:,None,None],increment.shape)
    difference=float(np.sum(lapse*(changed-increment)))
    prediction=-.2*ALPHA*delta*(1-np.cos(2*np.pi/5))/3
    assert abs(difference-prediction)<1e-15 and difference<0
    return dict(quadrature=rows,source_increment_at_A=float(changed[0,0,0]),
        increment_at_each_of_six_neighbors=float(changed[1,0,0]),
        total_increment=delta,nonconstant_lapse_increment_difference=difference,
        predicted_difference=float(prediction),
        phase_is_original_577_gauge_invariant_preparation=True,
        no_full_graph_propagation_simulated=True)


def smooth_sampling_check():
    rows=[];volume=(2*np.pi)**3
    for n in (8,16,32):
        eps=2*np.pi/n
        x,y,z=np.meshgrid(*(np.arange(n)*eps for _ in range(3)),indexing='ij')
        density=2+.3*np.cos(x)+.2*np.sin(y)
        lapse=1+.2*np.cos(x)+.1*np.sin(y)
        local=eps**3*density
        gap=float(np.sum(lapse*(average(local)-local)))
        expected=-ALPHA*volume*.04*(1-np.cos(eps))/3
        assert abs(gap-expected)<3e-13
        rows.append(dict(n=n,epsilon=eps,lapse_energy_difference=gap,
            exact_fourier_prediction=float(expected),gap_over_eps_squared=gap/eps**2))
    assert abs(rows[-1]['lapse_energy_difference'])<abs(rows[0]['lapse_energy_difference'])/14
    return dict(rows=rows,limit_gap_over_eps_squared=-ALPHA*volume*.04/6,
        only_fixed_smooth_classical_sampling_not_a_quantum_continuum_limit=True)


def derivative(a):
    modes=np.fft.fftfreq(a.size,d=1/a.size)
    return np.fft.ifft(1j*modes*np.fft.fft(a)).real


def normal_bracket_check():
    n=128;x=np.arange(n)*2*np.pi/n;dx=2*np.pi/n
    chi=.3+.15*np.sin(x);momentum=.2*np.sin(2*x)
    lapse=1+.2*np.cos(x);other=np.ones(n)
    s=np.sqrt(6*original.M)*np.tanh(chi/np.sqrt(6))
    f=original.M-s*s/6
    matrix,u,_=original.lattice.scalar.parameters()
    delta=np.stack((-u[0]*np.ones(n),s*s-u[1]),axis=-1)
    V=np.einsum('...i,ij,...j->...',delta,matrix,delta)/4
    dV=s*(delta@matrix)[:,1]
    dU=(dV/f**2+2*s*V/(3*f**3))*f/np.sqrt(original.M)
    grad=derivative(chi)
    def bracket(N,M):
        dchi_N=-derivative(N*grad)+N*dU
        dchi_M=-derivative(M*grad)+M*dU
        return float(dx*np.sum(dchi_N*M*momentum-N*momentum*dchi_M))
    expected=float(dx*np.sum(momentum*grad*(lapse*derivative(other)-other*derivative(lapse))))
    bare=bracket(lapse,other)
    assert abs(bare-expected)<1e-13 and abs(expected)>1e-3
    rows=[]
    for radius in (.8,.4,.2):
        multiplier=1-ALPHA*(1-np.cos(radius))/3
        smeared_lapse=1+.2*multiplier*np.cos(x)
        got=bracket(smeared_lapse,other)
        defect=got-expected;prediction=(multiplier-1)*expected
        assert abs(defect-prediction)<1e-13 and abs(defect)>1e-6
        rows.append(dict(averaging_radius=radius,fourier_multiplier=float(multiplier),
            canonical_bracket=got,defect=defect,predicted_defect=float(prediction)))
    return dict(unsmeared_canonical_bracket=bare,target_momentum_bracket=expected,
        unsmeared_error=abs(bare-expected),rows=rows,
        original_singlet_geodesic_sector_with_full_potential=True,
        transverse_volume_normalized_to_one=True,
        diagnostic_not_a_quantum_graph_constraint_algebra=True)


def run():
    checks={f.__name__:f() for f in (full_original_allocation_check,
        same_gauss_state_local_source_check,smooth_sampling_check,normal_bracket_check)}
    dependencies=('joint_curved_quantum_source.py','joint_quantum_measure_records.py',
        'research_note_350.md','research_note_351.md','research_note_352.md','research_note_366.md',
        'research_note_572.md','research_note_574.md','research_note_577.md','research_note_584.md',
        'research_round_584_checks.json')
    return dict(round=585,tests_run=len(checks),failures=0,errors=0,evidence=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
        scope='positive local-lapse completions of the same original unit-lapse finite-graph Hamiltonian; original Gauss source gives different local constraint energies while all unit-lapse dynamics and spatial-weight responses agree; classical smooth leading limit is shared; restricted continuum normal-bracket test rejects naive smearing at fixed lapse calibration, not a quantum Einstein solution or no-go for cognition/GR')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
