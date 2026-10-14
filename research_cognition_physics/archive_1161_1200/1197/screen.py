"""Fixed finite-window numerical screen. Default stdout only; --write is exclusive.
Not an interval certificate, instrument implementation, modular-flow or GR proof.
"""
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
REF=HERE.parent/'1196'/'hartree_reference_result.json'

def calculate():
    reference=json.loads(REF.read_text(encoding='utf-8'))
    rd=reference['fixed_point_floating_diagnostics']
    z=rd['z_s']; b=rd['b']; ell=16/65; n=65; M=32
    eps=.15; nu=.35; Om=.5; mass=2.; N=1+eps*z; v=1+nu*z
    mu=mass*N*ell/v
    X=np.array([[0,1],[1,0]],complex); Z=np.diag([1.,-1.]); I=np.eye(2)
    k0=2*np.pi*np.arange(n)/n
    cp=(I[None]-(mu*Z[None]+np.sin(k0)[:,None,None]*X[None])/
        np.sqrt(mu*mu+np.sin(k0)**2)[:,None,None])/2
    C=np.array([[np.einsum('n,nab->ab',np.exp(1j*k0*(i-j)),cp)/n
        for j in range(M)] for i in range(M)]).transpose(0,2,1,3).reshape(2*M,2*M)
    B=np.zeros_like(C)
    for i in range(M):
        B[2*i:2*i+2,2*i:2*i+2]=mu*(i+.5)*Z
        if i<M-1:
            B[2*i:2*i+2,2*i+2:2*i+4]=-(i+1)*.5j*X
            B[2*i+2:2*i+4,2*i:2*i+2]=(i+1)*.5j*X
    eigen,U=np.linalg.eigh(B)
    envelope=np.sin(np.pi*(np.arange(M)+.5)/M)
    envelope/=np.linalg.norm(envelope)
    W=np.kron(envelope[:,None],I); Cw=W.conj().T@C@W
    models={}; temps={}
    for multiple in (1,2,4):
        # Stable Fermi function, without clipping a covariance eigenvalue.
        x=multiple*np.pi*eigen
        vals=np.empty_like(x); pos=x>=0
        vals[pos]=np.exp(-x[pos])/(1+np.exp(-x[pos]))
        vals[~pos]=1/(1+np.exp(x[~pos]))
        CB=(U*vals)@U.conj().T
        reduced=W.conj().T@CB@W; models[multiple]=reduced
        temps[str(multiple)]={'beta_over_pi':multiple,
          'upper_mode_occupation':float(reduced[0,0].real),
          'compressed_covariance_error':float(np.linalg.norm(Cw-reduced,2))}
    residual=temps['2']['compressed_covariance_error']
    controls={str(j):{'prediction_separation_from_2pi':float(np.linalg.norm(models[j]-models[2],2)),
      'residual_over_separation':float(residual/np.linalg.norm(models[j]-models[2],2))} for j in (1,4)}
    # Complete orbit sum of the original 65^3 momenta; no rounded bins.
    ks=[]; weights=[]
    for i in range(33):
        for j in range(i,33):
            for h in range(j,33):
                labels=[i,j,h]; perm=6
                for q in set(labels): perm//=math.factorial(labels.count(q))
                weights.append(perm*math.prod(1 if q==0 else 2 for q in labels))
                ks.append(math.sqrt(sum((math.sin(2*math.pi*q/65)/ell)**2 for q in labels)))
    assert sum(weights)==65**3
    k=np.asarray(ks); wt=np.asarray(weights)/65**3
    E=np.sqrt((mass*N)**2+(v*k)**2); r=math.sqrt(Om*Om+b*b)
    theta=np.pi/8; T=np.pi/2
    source=np.vstack([-v*k/E*np.cos(2*theta),-v*k/E*np.sin(2*theta),-mass*N/E])
    gs=np.array([Om/r,0.,-b/r]); initial=np.concatenate([source.ravel(),gs,[0.]])
    def rhs(y):
        s=y[:-4].reshape(3,-1); g=y[-4:-1]
        lapse=1+eps*g[2]; speed=1+nu*g[2]
        em=2*mass*np.dot(wt,s[2]); et=2*np.dot(wt,k*s[0]); force=eps*em+nu*et-eps*Om
        ds=np.vstack([-2*mass*lapse*s[1],2*(mass*lapse*s[0]-speed*k*s[2]),2*speed*k*s[1]])
        dg=np.array([-2*force*g[1],2*(force*g[0]+Om*g[2]),-2*Om*g[1]])
        return np.concatenate([ds.ravel(),dg,[lapse]])
    def energy(y):
        s=y[:-4].reshape(3,-1); g=y[-4:-1]
        em=2*mass*np.dot(wt,s[2]); et=2*np.dot(wt,k*s[0])
        return (1+eps*g[2])*em+(1+nu*g[2])*et-Om*g[0]-Om-eps*Om*g[2]
    dynamics=[]
    for steps in (1024,2048):
        y=initial.copy(); dt=T/steps
        for _ in range(steps):
            k1=rhs(y); k2=rhs(y+dt*k1/2); k3=rhs(y+dt*k2/2); k4=rhs(y+dt*k3)
            y+=dt*(k1+2*k2+2*k3+k4)/6
        s=y[:-4].reshape(3,-1)
        dynamics.append({'steps':steps,'delta_tau':float(y[-1]-N*T),
          'probe_PminusY_probability_change':float((np.sin(.2*y[-1])-np.sin(.2*N*T))/2),
          'Gs_z_change':float(y[-2]-z),'complete_Hartree_energy_drift':float(energy(y)-energy(initial)),
          'source_Bloch_length_error':float(np.max(np.abs(np.sum(s*s,axis=0)-1))),
          'Gs_Bloch_length_error':float(abs(np.dot(y[-4:-1],y[-4:-1])-1))})
    beta_rotation=np.diag(np.exp(-1j*theta*np.array([1.,-1.])))
    delta_et=(np.cos(2*theta)-1)*rd['source_hopping_energy_per_site']
    delta_force=nu*delta_et
    closure={'same_initial_macro_covariance_numerical_residual':float(np.linalg.norm(beta_rotation@Cw@beta_rotation.conj().T-Cw,2)),
      'same_initial_macro_state_exact_reason':'real envelope Fourier intensity is even; sin(k) term cancels; Cmacro commutes with beta',
      'source_hopping_energy_change':float(delta_et),'initial_source_force_change':float(delta_force),
      'initial_Gs_z_second_derivative_difference':float(-4*Om*Om*delta_force/r),
      'Hartree_initial_probe_probability_third_derivative_difference':float(-.1*eps*4*Om*Om*delta_force/r),
      'claim_complete_preparation_record_same':False,
      'rejects':'autonomous closure from this two-mode instantaneous state plus identical Gs/Gt; not general coarse-graining'}
    return {'status':'COMPLETED_NUMERICAL_SCREEN_NOT_CERTIFIED',
      'reference_sha256':hashlib.sha256(REF.read_bytes()).hexdigest(),
      'parameters':{'ell':'16/65','Nlat':65,'longitudinal_region_sites':32,'time':'pi/2','source_beta_rotation':'pi/8'},
      'macro_readout':{'envelope':'normalized sin(pi*(j+1/2)/32), j=0..31',
        'transverse':'p_y=p_z=0 as in 1196; full periodic cross-section',
        'source_upper_occupation':float(Cw[0,0].real),
        'source_compressed_covariance_eigenvalues':np.linalg.eigvalsh(Cw).tolist(),
        'temperature_controls':temps,'identifiability_controls':controls,
        'boost_leakage_from_two_mode_span':float(np.linalg.norm(B@W-W@(W.conj().T@B@W),2))},
      'closure_test':closure,
      'source_response':{'cubic_orbits':len(k),'total_momenta':sum(weights),
        'initial_pulse_energy_per_site_copy_excluding_common_baseline':float(energy(initial)-rd['complete_Hartree_energy_per_site_excluding_global_baseline']),
        'runs':dynamics,'probability_step_difference':abs(dynamics[1]['probe_PminusY_probability_change']-dynamics[0]['probe_PminusY_probability_change'])},
      'scope':{'full_modular_flow_identity':False,'whole_parent_region':False,
        'actual_source_preparation_device':False,'new_finite_K_time_certificate':False,
        'RK_step_difference_is_rigorous_error':False,'actual_coarse_source_instrument':False,
        'GR_or_Jacobson_derived':False,'new_round_assigned':False}}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=calculate();out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text=json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:
        with (HERE/'result.json').open('x',encoding='utf-8',newline='\n') as f:f.write(text)
    print(text,end='')
if __name__=='__main__': main()
