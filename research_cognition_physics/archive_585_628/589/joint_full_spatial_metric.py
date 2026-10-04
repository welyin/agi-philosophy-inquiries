"""589: one positive finite-graph matter operator for a full spatial metric.

Geometry is a fixed positive background; no quantum geometry or continuum
limit at fixed hbar is claimed. The magnetic entry is reused unchanged.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_matter_energy_current as current
import joint_record_source_compression as records

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_full_spatial_metric_results.json'
spec = importlib.util.spec_from_file_location('magnetic589', HERE/'round589_drafts/full_metric_magnetic_entry.py')
mag = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mag)
PAR, PAIRS = mag.P, mag.PAIRS


def gram(shape):
    inv = np.linalg.inv(shape)
    return np.stack([np.stack([inv[..., a,c]*inv[..., b,d]-inv[..., a,d]*inv[..., b,c]
                              for c,d in PAIRS], axis=-1) for a,b in PAIRS], axis=-2)


def shape_exp(S, t=1.):
    eig, U = np.linalg.eigh(S)
    return (U*np.exp(t*eig)[..., None, :])@np.swapaxes(U, -1, -2)


def shape_field(grid):
    x,y,z = np.moveaxis(grid, -1, 0)
    S = np.zeros(grid.shape[:-1]+(3,3))
    S[...,0,0] = .16*np.sin(x); S[...,1,1] = .13*np.cos(y)
    S[...,2,2] = -S[...,0,0]-S[...,1,1]
    S[...,0,1] = S[...,1,0] = .18+.07*np.sin(x+y)
    S[...,0,2] = S[...,2,0] = .05*np.cos(z)
    S[...,1,2] = S[...,2,1] = -.09+.04*np.cos(x)
    return shape_exp(S)


def target_log(x, y):
    gx, _ = current.edge_gradients(x,y)
    return -np.einsum('...ab,...b->...a', original.inverse(x), gx)


def adj_batch(W):
    T=mag.old.generators(2)
    return 2*np.einsum('aij,...jk,bkl,...li->...ab', T, W, T, current.dag(W)).real


def ew_magnetic_grams(W, z):
    """All original representations evaluated on the ACTUAL colour-identity source.

    Colour multiplicities are exact; adjoint colour contributes zero here.
    Independent full-colour tests remain in the frozen magnetic entry.
    """
    delta, coeff = mag.old.coefficients(PAR)
    def scalar(p): return (z**p)[...,None,None]
    reps = [(3*delta*PAR['wq'], W*z[...,None,None]),
            (3*delta*PAR['wq'], scalar(-4)), (3*delta*PAR['wq'], scalar(2)),
            (delta*PAR['wl'], W*(z**-3)[...,None,None]),
            (delta*PAR['wl'], scalar(6)), (coeff[1], adj_batch(W)),
            (coeff[2], scalar(6))]
    odd = np.zeros(z.shape[:-1]+(3,3)); even = np.zeros_like(odd)
    for alpha, R in reps:
        O = (R-current.dag(R))/(2j)
        E = np.eye(R.shape[-1])-(R+current.dag(R))/2
        odd += .5*alpha*np.einsum('...fij,...gij->...fg', O.conj(), O).real
        even += .5*alpha*np.einsum('...fij,...gij->...fg', E.conj(), E).real
    return odd, even


def graph_data(q, psi):
    eps=q['eps']; phi=q['phi']; X=q['f']['X']; s=q['f']['s']
    W=original.lattice.su2(q['a'],eps); z=np.exp(1j*eps*q['a0'])
    logs=[]; edges=[]; faces=[]; faceW=[]; facez=[]
    for mu in range(3):
        edges.append((psi+np.roll(psi,-1,axis=mu))/2)
        Xnext=z[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],np.roll(X,-1,axis=mu))
        logs.append(target_log(phi,original.real_phi(Xnext,np.roll(s,-1,axis=mu))))
    for mu,nu in PAIRS:
        wm,wn=W[...,mu,:,:],W[...,nu,:,:]
        faceW.append(wm@np.roll(wn,-1,axis=mu)@current.dag(np.roll(wm,-1,axis=nu))@current.dag(wn))
        facez.append(z[...,mu]*np.roll(z[...,nu],-1,axis=mu)/np.roll(z[...,mu],-1,axis=nu)/z[...,nu])
        faces.append((psi+np.roll(psi,-1,axis=mu)+np.roll(psi,-1,axis=nu)+np.roll(np.roll(psi,-1,axis=mu),-1,axis=nu))/4)
    pe=np.stack(edges,axis=-1); pf=np.stack(faces,axis=-1)
    ell=np.stack(logs,axis=-2); u=ell*pe[...,None]
    grad=np.einsum('...ma,...ab,...nb->...mn',u,original.metric(phi),u)
    electric=PAR['b'][1]*np.einsum('...ma,...na->...mn',q['phat'],q['phat'])
    electric+=PAR['b'][2]*q['p0hat'][..., :,None]*q['p0hat'][...,None,:]
    electric/=pe[..., :,None]*pe[...,None,:]
    odd,even=ew_magnetic_grams(np.stack(faceW,axis=-3),np.stack(facez,axis=-1))
    scale=pf[..., :,None]*pf[...,None,:]
    return dict(eps=eps, grad=grad, electric=electric, odd=odd/scale, even=even/scale,
                scalar_kinetic=float(np.sum(np.einsum('...a,...ab,...b->...',q['P'],original.inverse(phi),q['P'])/(2*eps**3*psi**6))),
                onsite=float(eps**3*np.sum(psi**6*original.node_potential(phi))))


def energy(data,shape):
    inv=np.linalg.inv(shape); B=gram(shape); eps=data['eps']
    gradient=float(eps/2*np.sum(inv*data['grad']))
    electric=float(np.sum(shape*data['electric'])/eps)
    magnetic=float((np.sum(B*data['odd'])+np.sum(np.diagonal(B,axis1=-2,axis2=-1)*np.diagonal(data['even'],axis1=-2,axis2=-1)))/eps)
    out=dict(scalar_kinetic=data['scalar_kinetic'],onsite=data['onsite'],gradient=gradient,electric=electric,magnetic=magnetic)
    out['total']=sum(out.values())
    return out


def continuum_energy(q,psi,shape,with_densities=False):
    c=original.geometry.make_source(q['N']); f=c['f']; inv=np.linalg.inv(shape)
    eps=q['eps']; metric=original.metric(c['phi'])
    grad=.5*psi**2*np.einsum('...mn,...ma,...ab,...nb->...',inv,c['Dphi'],metric,c['Dphi'])
    elec=psi**-2*(PAR['b'][1]*np.einsum('...mn,...ma,...na->...',shape,f['E'],f['E'])+
                         PAR['b'][2]*np.einsum('...mn,...m,...n->...',shape,f['E0'],f['E0']))
    fw=[]; f0=[]
    for mu,nu in PAIRS:
        fw.append(original.geometry.derivative(f['a'][...,nu,:],mu)-original.geometry.derivative(f['a'][...,mu,:],nu)-np.cross(f['a'][...,mu,:],f['a'][...,nu,:]))
        f0.append(original.geometry.derivative(f['a0'][...,nu],mu)-original.geometry.derivative(f['a0'][...,mu],nu))
    fw=np.stack(fw,axis=-2); f0=np.stack(f0,axis=-1); B=gram(shape)
    magnetic=.5*psi**-2*(PAR['K'][1]*np.einsum('...fg,...fa,...ga->...',B,fw,fw)+PAR['K'][2]*np.einsum('...fg,...f,...g->...',B,f0,f0))
    densities=dict(scalar_kinetic=.5*psi**-6*c['pKp'],onsite=psi**6*c['U'],gradient=grad,electric=elec,magnetic=magnetic)
    out={name:float(eps**3*np.sum(value)) for name,value in densities.items()}; out['total']=sum(out.values())
    return (out,densities) if with_densities else out


def magnetic_check():
    got=mag.run()
    assert got==json.loads(mag.TARGET.read_text('utf8'))
    return dict(frozen_entry_reproduced=True,**got)


def scalar_electric_check():
    rng=np.random.default_rng(5892); x=rng.normal(size=5)*.2; y=rng.normal(size=(3,5))*.3
    ell=target_log(x,y); norm=np.einsum('ma,ab,mb->m',ell,original.metric(x),ell)
    log_error=float(np.max(abs(norm-original.distance_squared(x,y))))
    shape=shape_exp(np.array([[.2,.17,-.08],[.17,-.13,.11],[-.08,.11,-.07]]))
    pe=np.array([1.03,.97,1.11]); u=ell*pe[:,None]
    G=u@original.metric(x)@u.T; value=float(np.sum(np.linalg.inv(shape)*G))
    W=mag.old.group_exp(rng.normal(size=3),2); z=np.exp(.31j)
    def rotate(v):
        X=z**3*np.einsum('ij,...j->...i',W,v[...,:2]+1j*v[...,2:4])
        return original.real_phi(X,v[...,4])
    ell2=target_log(rotate(x),rotate(y)); u2=ell2*pe[:,None]
    err=float(np.max(abs(ell2-rotate(ell))))
    value2=float(np.sum(np.linalg.inv(shape)*(u2@original.metric(rotate(x))@u2.T)))
    elec=0.; elec2=0.
    for n,b in ((3,PAR['b'][0]),(2,PAR['b'][1])):
        p=rng.normal(size=(3,n*n-1))/pe[:,None]; R=mag.adj(mag.old.group_exp(rng.normal(size=n*n-1),n))
        elec+=b*float(np.sum(shape*(p@p.T))); pr=p@R.T
        elec2+=b*float(np.sum(shape*(pr@pr.T)))
    p=rng.normal(size=3)/pe; elec+=PAR['b'][2]*p@shape@p; elec2+=PAR['b'][2]*p@shape@p
    assert max(log_error,err,abs(value2-value),abs(elec2-elec))<2e-12
    assert min(value,elec)>0
    return dict(log_squared_distance_error=log_error,log_equivariance_error=err,
                scalar_gauge_error=abs(value2-value),all_group_electric_gauge_error=abs(elec2-elec),
                positive_scalar_value=value,positive_electric_value=float(elec),shape_eigenvalues=np.linalg.eigvalsh(shape).tolist())


def common_source_check():
    rows=[]
    for N in (8,12):
        q=original.shared_source(N); x,y,z=np.moveaxis(q['grid'],-1,0)
        psi=1.1+.05*np.cos(x)+.02*np.sin(y); shape=shape_field(q['grid']); data=graph_data(q,psi)
        new=energy(data,shape); old=original.graph_energy(q,psi); recovered=energy(data,np.eye(3))
        exact=max(abs(old[k]-recovered[k]) for k in old)
        spectra=[np.linalg.eigvalsh(shape),np.linalg.eigvalsh(np.linalg.inv(shape)),np.linalg.eigvalsh(gram(shape))]
        lo=min(1.,*(float(a.min()) for a in spectra)); hi=max(1.,*(float(a.max()) for a in spectra))
        gw,g0=original.lattice.residual(q); gauss=max(float(np.max(abs(gw))),float(np.max(abs(g0))))/q['eps']**3
        assert exact<1e-9 and gauss<1e-11 and lo*old['total']<=new['total']<=hi*old['total']
        assert min(new.values())>0
        rows.append(dict(N=N,exact_conformal_error=exact,Gauss_density_error=gauss,new_energy=new,
                         form_comparison_constants=[lo,hi],conformal_total=old['total'],
                         determinant_shape_error=float(np.max(abs(np.linalg.det(shape)-1)))))
    return dict(rows=rows,same_original_physical_source=True,new_background_not_asserted_Einstein_initial_data=True)


def continuum_check():
    rows=[]
    for N in (8,12,20,32):
        q=original.shared_source(N);x,y,z=np.moveaxis(q['grid'],-1,0)
        psi=1.1+.05*np.cos(x)+.02*np.sin(y);shape=shape_field(q['grid'])
        data=graph_data(q,psi);got=energy(data,shape); expected,densities=continuum_energy(q,psi,shape,True)
        err={k:abs(got[k]-expected[k]) for k in expected}
        B=gram(shape); eps=q['eps']
        local=dict(gradient=.5/eps**2*np.sum(np.linalg.inv(shape)*data['grad'],axis=(-1,-2)),
                   electric=np.sum(shape*data['electric'],axis=(-1,-2))/eps**4,
                   magnetic=(np.sum(B*data['odd'],axis=(-1,-2))+np.sum(np.diagonal(B,axis1=-2,axis2=-1)*np.diagonal(data['even'],axis1=-2,axis2=-1),axis=-1))/eps**4)
        local_errors={k:float(np.max(abs(v-densities[k]))) for k,v in local.items()}
        rows.append(dict(N=N,discrete=got,continuum_quadrature=expected,errors=err,max_local_density_errors=local_errors))
    for k in ('gradient','electric','magnetic'):
        assert rows[-1]['max_local_density_errors'][k]<rows[0]['max_local_density_errors'][k]/2,k
    assert rows[-1]['errors']['total']<rows[0]['errors']['total']/10
    return dict(rows=rows,classical_smooth_symbol_only=True,colour_source_zero=True,
                no_fixed_hbar_quantum_limit=True,no_second_order_rate_claim=True)


def shear_response_check():
    q=original.shared_source(8);x,y,z=np.moveaxis(q['grid'],-1,0)
    psi=1.1+.05*np.cos(x)+.02*np.sin(y);data=graph_data(q,psi)
    S=np.array([[0.,.31,-.12],[.31,0.,.17],[-.12,.17,0.]])
    I=np.eye(3)
    dB=np.array([[-S[a,c]*I[b,d]-I[a,c]*S[b,d]+S[a,d]*I[b,c]+I[a,d]*S[b,c]
                   for c,d in PAIRS] for a,b in PAIRS])
    deriv=dict(scalar_kinetic=0.,onsite=0.,gradient=float(-q['eps']/2*np.sum(S*data['grad'])),
               electric=float(np.sum(S*data['electric'])/q['eps']),
               magnetic=float(np.sum(dB*data['odd'])/q['eps']))
    deriv['total']=sum(deriv.values()); rows=[]
    for t in (.02,.01,.005):
        p,m=energy(data,shape_exp(S,t)),energy(data,shape_exp(S,-t))
        error=max(abs((p[k]-m[k])/(2*t)-deriv[k]) for k in deriv)
        rows.append(dict(step=t,max_derivative_error=error))
    assert rows[-1]['max_derivative_error']<rows[0]['max_derivative_error']/12
    assert abs(deriv['total'])>1e-5 and abs(deriv['magnetic'])>1e-5
    return dict(first_variations=deriv,centered_difference=rows,determinant_preserving=True,
                frozen_conformal_only_response=0.,full_covariant_stress_not_claimed=True)


def instrument_check():
    rng=np.random.default_rng(5896); phi=rng.normal(size=(23,5))*.3
    K=original.inverse(phi); psi=1.1+.1*rng.random(23); w=.6**3*psi**6; hbar=.7
    shape=shape_field(rng.normal(size=(23,3))); pe=1+.1*rng.random((23,3))
    C=shape/(pe[..., :,None]*pe[...,None,:]); grad=rng.normal(size=(23,5))+1j*rng.normal(size=(23,5))
    lg=rng.normal(size=(23,3,3))+1j*rng.normal(size=(23,3,3)); potential=1+rng.random(23)
    lf,df,lc,dc=records.instruments(phi[...,4])
    def form(l,dl):
        v=l[...,None]*grad;v[...,4]+=dl
        kin=hbar*hbar/(2*w)*np.einsum('...a,...ab,...b->...',v.conj(),K,v).real
        link=l[...,None,None]*lg
        return kin+PAR['b'][1]*hbar*hbar/.6*np.einsum('...mn,...ma,...na->...',C,link.conj(),link).real+l*l*potential
    f=sum(form(lf[:,j],df[:,j]) for j in range(4));c=sum(form(lc[:,j],dc[:,j]) for j in range(2))
    predicted=hbar*hbar/(2*w)*K[:,4,4]*records.delta_a(phi[:,4]);error=float(np.max(abs(f-c-predicted)))
    assert error<1e-12 and predicted.min()>0
    # A single original read has A=cos²(s)/(16*(1-sin²(s)/4)) <= 1/16.
    s=np.linspace(-3.46,3.46,501); A=np.cos(s)**2/(16*(1-np.sin(s)**2/4))
    assert A.max()<=1/16+1e-15
    return dict(pointwise_full_form_instrument_gap_error=error,minimum_positive_gap=float(predicted.min()),
                single_read_A_max=float(A.max()),energy_injection_bound_coefficient=original.M/32,
                numerical_check_is_local_differential_identity=True,
                full_Hilbert_statement_proved_in_note=True,off_conformal_histories_not_asserted_equal=True)


def run():
    checks=(magnetic_check,scalar_electric_check,common_source_check,continuum_check,shear_response_check,instrument_check)
    evidence={fn.__name__:fn() for fn in checks}
    names=('joint_curved_quantum_source.py','joint_matter_energy_current.py','joint_record_source_compression.py',
           'joint_quotient_gauge_completion.py','joint_gauss_einstein_initial_data.py',
           'joint_gauss_continuum_sampling.py','round589_drafts/STATUS.md',
           'round589_drafts/full_metric_magnetic_entry.py','round589_drafts/full_metric_magnetic_entry_results.json')
    hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    return dict(round=589,tests_run=len(checks),failures=0,errors=0,evidence=evidence,dependency_hashes=hashes,
                scope='Fixed finite graph and arbitrary fixed positive spatial metric; same quotient matter and instruments. Classical smooth-symbol consistency only; geometry dynamics and quantum refinement remain open.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
