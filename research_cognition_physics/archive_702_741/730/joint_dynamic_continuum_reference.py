"""730: original dynamical initial jets, continuous-state bridge and finite source calibration.

The Hadamard existence argument is analytic and uses cited continuum results.
The time integrations below are local-symbol matrices, not continuum solutions.
"""
import argparse
import hashlib
import json
from pathlib import Path
from functools import lru_cache
import numpy as np
import joint_matter_ground_source as old
import joint_gravity_material_coordinates as geometry
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_dynamic_continuum_reference_results.json'
GAMMA=old.chiral.kinetic_matrices()
SIG=old.chiral.SIG

def smoothstep(s):
    if s<=0:return 0.
    if s>=1:return 1.
    a=np.exp(-1/s);b=np.exp(-1/(1-s))
    return float(a/(a+b))

@lru_cache(maxsize=1)
def initial():
    f=geometry.fields(8);q=f['q'];psi=f['psi']
    metric=psi[...,None,None]**4*np.eye(3)
    K=psi[...,None,None]**-2*f['tensor']+np.sqrt(q['tau2'])/3*metric
    # Explicit ADM convention K=-dot(g)/2 in unit lapse and zero shift.
    dg=-2*K
    bw,b0=old.matter.original.lattice.PAR['b'][1:]
    da=2*bw*psi[...,None,None]**-2*q['f']['E']
    da0=2*b0*psi[...,None]**-2*q['f']['E0']
    return dict(f=f,q=q,g=metric,K=K,dg=dg,v=f['v'],a=q['f']['a'],a0=q['f']['a0'],
                da=da,da0=da0,psi=psi)

def gauge_h(a,a0):
    hs=[]
    for i in range(3):
        gen=np.zeros((32,32),complex)
        for name,sl in old.matter.SLICES.items():
            size=sl.stop-sl.start
            value=old.chiral.CHARGES[name]*a0[i]*np.eye(size,dtype=complex)
            if name in ('Q','L'):
                weak=sum(a[i,j]*SIG[j]/2 for j in range(3))
                value+=np.kron(np.eye(size//4),np.kron(weak,np.eye(2)))
            gen[sl,sl]=value
        hs.append(gen)
    return hs

@lru_cache(maxsize=1)
def point():
    data=initial();idx=(0,2,1)
    return {key:np.array(data[key][idx]) for key in ('g','dg','v','a','a0','da','da0')} | {
        'phi':np.array(data['q']['phi'][idx]),'index':idx}

def collar(t,width):
    p=point();eta=smoothstep((t+2*width)/width)
    phi0=np.array([0.,old.CENTER[0],0.,0.,old.CENTER[1]])
    g=(1-eta)*np.eye(3)+eta*(p['g']+t*p['dg'])
    phi=(1-eta)*phi0+eta*(p['phi']+t*p['v'])
    a=eta*(p['a']+t*p['da']);a0=eta*(p['a0']+t*p['da0'])
    return g,phi,a,a0,eta

def matrices(t,width,k,gamma=0.):
    g,phi,a,a0,eta=collar(t,width)
    ee,v=np.linalg.eigh(g);assert ee.min()>0 and old.matter.original.F(phi)>0
    frame=(v/np.sqrt(ee))@v.T
    hk=sum(GAMMA[j]*sum(frame[j,i]*k[i] for i in range(3)) for j in range(3))
    gen=gauge_h(a,a0)
    gauge=sum(GAMMA[j][:32,:32]@sum(frame[j,i]*gen[i] for i in range(3)) for j in range(3))
    kinetic=hk+old.bdg(gauge,np.zeros_like(gauge))
    mass=old.bdg(*old.matter.mass_matrices(phi))
    scaled=np.exp(-gamma*eta)*kinetic
    return scaled+mass,-eta*scaled

def projector(B):
    e,v=np.linalg.eigh(B)
    return v[:,e<0]@v[:,e<0].conj().T

def step(B,dB,dt):
    e,v=np.linalg.eigh(B);phase=np.exp(-1j*dt*e)
    U=(v*phase)@v.conj().T
    delta=e[:,None]-e[None,:]
    kernel=-1j*dt*np.exp(-1j*dt*(e[:,None]+e[None,:])/2)*np.sinc(dt*delta/(2*np.pi))
    dU=v@(kernel*(v.conj().T@dB@v))@v.conj().T
    return U,dU

def flow(width,steps=96,gamma=0.,derivative=False):
    k=np.array([.31,-.27,.19]);start=-.24;dt=-start/steps
    ps=[];dps=[];ends=[];sources=[]
    for sign in (1,-1):
        B,_=matrices(start,width,sign*k,gamma);P=projector(B);dP=np.zeros_like(P)
        for j in range(steps):
            B,dB=matrices(start+(j+.5)*dt,width,sign*k,gamma)
            U,dU=step(B,dB,dt);before=P
            P=U@before@U.conj().T
            if derivative:dP=dU@before@U.conj().T+U@dP@U.conj().T+U@before@dU.conj().T
        B,dB=matrices(0,width,sign*k,gamma)
        ps.append(P);dps.append(dP);ends.append(B);sources.append(dB)
    # [particles+k, holes-k, particles-k, holes+k] -> global canonical order.
    order=np.r_[np.arange(32),np.arange(64,96),np.arange(96,128),np.arange(32,64)]
    def assemble(xs):return old.block(xs)[np.ix_(order,order)]
    return tuple(assemble(xs) for xs in (ps,dps,ends,sources))

def original_jet_check():
    d=initial();q=d['q'];metric=d['g'];psi=d['psi']
    recovered=psi[...,None]**6*np.einsum('...ij,...j->...i',old.matter.original.metric(q['phi']),d['v'])
    trace=np.einsum('...ij,...ji->...',np.linalg.inv(metric),d['K'])
    err=max(float(np.max(abs(recovered-q['p']))),float(np.max(abs(trace-np.sqrt(q['tau2'])))))
    eigenfloor=10.;Ffloor=10.;jeterr=0.;h=1e-6
    for width in (.08,.12):
        for t in np.linspace(-.24,.01,31):
            g,phi,a,a0,eta=collar(float(t),width)
            eigenfloor=min(eigenfloor,float(np.linalg.eigvalsh(g).min()))
            Ffloor=min(Ffloor,float(old.matter.original.F(phi)))
        zero=collar(0,width);plus=collar(h,width);minus=collar(-h,width)
        p=point()
        for i,key,dkey in [(0,'g','dg'),(1,'phi','v'),(2,'a','da'),(3,'a0','da0')]:
            jeterr=max(jeterr,float(np.max(abs(zero[i]-p[key]))),
                       float(np.max(abs((plus[i]-minus[i])/(2*h)-p[dkey]))))
    p=point();h0,s0=p['phi'][1],p['phi'][4];vh,vs=p['v'][1],p['v'][4]
    a,b=abs(old.matter.Y['nu']),abs(old.matter.Y['s'])
    theta_dot=a*b*(s0*vh-h0*vs)/(b*b*s0*s0+4*a*a*h0*h0)
    assert err<1e-12 and jeterr<1e-9 and eigenfloor>0 and Ffloor>1.8 and abs(theta_dot)>1e-4
    return dict(original_grid=8,point=list(p['index']),canonical_momentum_and_CMC_error=err,
        original_initial_jet_error=jeterr,minimum_sampled_metric_eigenvalue=eigenfloor,
        minimum_sampled_F=Ffloor,extrinsic_curvature_norm=float(np.linalg.norm(d['K'])),
        actual_point_mass_mixing_time_derivative=float(theta_dot),
        diagnostic_is_first_jet_collar_not_solved_future_Einstein_evolution=True)

def transport_check():
    c0=old.bdg(*old.matter.mass_matrices(
        np.array([0.,old.CENTER[0],0.,0.,old.CENTER[1]])))
    e,v=np.linalg.eigh(c0);beta=(v*np.sign(e))@v.conj().T
    errors=[float(np.max(abs(beta@beta-np.eye(64))))]
    for g in GAMMA:errors.append(float(np.max(abs(beta@g+g@beta))))
    for t in (-.2,-.12,-.03,0.):
        B,_=matrices(t,.08,np.array([.31,-.27,.19]))
        m=beta@B
        errors.append(float(np.max(abs(m.conj().T@beta-beta@m))))
    P1,_,B,_=flow(.08);P2,_,B2,_=flow(.12);Pfine,_,_,_=flow(.08,192)
    C=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    errors.extend([float(np.max(abs(B-B2))),
        float(np.max(abs(C@B.conj()@C+B))),
        float(np.max(abs(C@P1.conj()@C+P1-np.eye(128)))),
        float(np.max(abs(P1@P1-P1)))])
    discrepancy=float(np.linalg.norm(P1-P2,'fro'))
    instantaneous=float(np.linalg.norm(P1-projector(B),'fro'))
    integration_error=float(np.linalg.norm(P1-Pfine,'fro'))
    assert max(errors)<3e-11 and discrepancy>1e-4 and instantaneous>.01 and integration_error<5e-5
    return dict(all_original_modes_in_two_momentum_calibration=64,Nambu_dimension=128,
        max_Clifford_reality_and_transport_error=max(errors),
        auxiliary_past_covariance_difference=discrepancy,
        transported_vs_instantaneous_ground_difference=instantaneous,
        midpoint_96_vs_192_difference=integration_error,
        local_symbol_calibration_only_not_continuum_Hadamard_test=True)

def record_source_check():
    P,dP,B,G=flow(.08,96,derivative=True)
    e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2)
    C=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    ce=C@e.conj();Q=np.outer(e,e.conj())+np.outer(ce,ce.conj());R=np.eye(128)-2*Q
    change=lambda X:(R@X@R-X)/2
    J=change(P);dJ=change(dP)
    direct=float(np.trace(G@J).real/2);state=float(np.trace(B@dJ).real/2);total=direct+state
    dx=2e-5;ep=[];perturbed=[]
    for gamma in (dx,-dx):
        PP,_,BB,_=flow(.08,96,gamma=gamma)
        ep.append(float(np.trace(BB@change(PP)).real/2))
        perturbed.append(PP)
    finite=(ep[0]-ep[1])/(2*dx)
    fdP=(perturbed[0]-perturbed[1])/(2*dx)
    derivative_error=float(np.max(abs(fdP-dP)))
    Rg=old.matter.representation(old.matter.gauge.group_exp(np.array([.1]*8),3),
         old.matter.gauge.group_exp(np.array([.13,-.2,.09]),2),np.exp(.17j))
    physical=old.block([Rg,Rg]);full=old.block([physical,physical.conj()])
    gauge_error=float(np.max(abs(full@R-R@full)))
    assert abs(finite-total)<1e-7 and derivative_error<2e-8 and gauge_error<1e-12 and abs(state)>1e-6
    return dict(actual_sterile_mode_occupation=float(np.vdot(e,P@e).real),
        nonselective_record_energy_change=float(np.trace(B@J).real/2),
        direct_fixed_state_source=direct,preparation_response=state,
        total_record_energy_response=total,finite_difference=finite,
        response_error=abs(finite-total),covariance_derivative_error=derivative_error,
        gauge_singlet_record_error=gauge_error,
        preparation_response_cannot_be_dropped=True,
        no_instantaneous_extended_detector_or_continuum_limit_claim=True)

def run():
    results=dict(original_jets_and_auxiliary_collar=original_jet_check(),
                 complete_matter_transport=transport_check(),
                 actual_record_and_total_source=record_source_check())
    deps=('research_note_573.md','research_note_598.md','research_note_630.md','research_note_633.md',
          'research_note_635.md','research_note_663.md','research_note_667.md','research_note_704.md',
          'research_note_729.md','joint_gravity_material_coordinates.py','joint_fermion_gauss_completion.py',
          'round730_drafts/varying_reference_entry.py')
    return dict(round=730,tests_run=3,failures=0,errors=0,results=results,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Conditional continuum Hadamard reference on the original given smooth positive-F dynamical background, via an auxiliary-past construction that is identical on the target collar; all original matrix masses, connection and self-dual CAR retained. Smooth record differences share the same renormalization. Numeric checks validate original initial jets, local complete-matter time transport and total source bookkeeping only. No lattice-to-continuum theorem, unique reference selection, dynamical quantum Gauss state, autonomous detector or self-consistent semiclassical Einstein solution.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results'],ensure_ascii=False,indent=2))
