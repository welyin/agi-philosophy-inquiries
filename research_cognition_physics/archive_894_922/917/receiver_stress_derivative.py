"""917 working: derivative of the original receiver stress on the same background.
No replacement by an on-shell spinor; derivative of the full off-shell stress.
"""
from pathlib import Path
import sys,json,argparse,hashlib,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'916'))
import covariant_joint_residual as c
r=c.ex.previous.rebuilt.r;weighted=c.aj.prior.weighted
TARGET=HERE/'receiver_stress_derivative_results.json'

def spinor_derivative(field,x,axes=()):
    if len(axes)<2:return field.evaluate(x,None if not axes else axes[0])
    nt=axes.count(0);theta=x[:,0]/field.T;powers=np.zeros((len(x),4))
    for n in range(nt,4):
        factor=1
        for j in range(nt):factor*=n-j
        powers[:,n]=factor*theta**(n-nt)/field.T**nt
    co=field.co.copy()
    for mu in axes:
        if mu:co*=1j*field.k[:,mu-1,None,None]
    out=[]
    for first in range(0,len(x),128):
        end=first+128;q=np.exp(1j*x[first:end,1:]@field.k.T)@co.reshape((-1,16))
        out.append(np.einsum('bp,bpf->bf',powers[first:end],q.reshape((-1,4,4))))
    return np.concatenate(out)

def extended_frame(p,z):
    # Analytic matrix square root near the actual positive real spatial metric.
    # Fixed real eigenbasis is only the Sylvester preconditioner; residual
    # correction includes imaginary perturbations, including repeated eigenvalues.
    gam=z['gamma'];vals,O=np.linalg.eigh(gam.real);roots=np.sqrt(vals)
    B=np.einsum('bai,bi,bci->bac',O,roots,O).astype(gam.dtype)
    def sylvester0(X):
        local=np.einsum('bai,bac,bcj->bij',O,X,O)
        return np.einsum('bai,bij,bcj->bac',O,local/(roots[:,:,None]+roots[:,None,:]),O)
    for _ in range(3):B+=sylvester0(gam-B@B)
    eye=np.eye(3)
    operator=(np.einsum('ik,blj->bijkl',eye,B)+np.einsum('bik,lj->bijkl',B,eye)).reshape((-1,9,9))
    dg=p['dg'];dB=np.linalg.solve(operator,np.moveaxis(dg[:,:,1:,1:].reshape((-1,4,9)),1,2)).swapaxes(1,2).reshape((-1,4,3,3))
    ig=z['ig'];di=z['di'];alpha=z['alpha'];beta=z['beta']
    da=.5*alpha[:,None]**3*di[:,:,0,0]
    db=-di[:,:,0,1:]/ig[:,0,0,None,None]+ig[:,0,1:][:,None,:]*di[:,:,0,0,None]/ig[:,0,0,None,None]**2
    E=np.zeros_like(p['g']);E[:,0,0]=alpha;E[:,1:,1:]=B;E[:,1:,0]=np.einsum('bai,bi->ba',B,beta)
    dE=np.zeros_like(dg);dE[:,:,0,0]=da;dE[:,:,1:,1:]=dB
    dE[:,:,1:,0]=np.einsum('bmai,bi->bma',dB,beta)+np.einsum('bai,bmi->bma',B,db)
    inv=np.linalg.inv(E)
    omega=np.einsum('bar,brms,bsq->bmaq',E,z['Gamma'],inv)-np.einsum('bmar,brq->bmaq',dE,inv)
    low=np.einsum('aq,bmqc->bmac',r.shared.ETA,omega)
    return dict(E=E,dE=dE,W=alpha[:,None,None]*inv,low=low,inv=inv)

def stress_jet(bg,u,x):
    cache=c.aj.Cache(bg,x);p=cache.jets();z=c.geometry(p);fr=extended_frame(p,z)
    oldfr=r.shared.frame_data(dict(g=p['g']),z)
    frame_error=max(c.maximum(fr[k]-oldfr[k]) for k in ('E','dE','W','low'))
    val=u.evaluate(x);du=np.stack([u.evaluate(x,i) for i in range(4)],axis=1)
    ddu=np.stack([np.stack([spinor_derivative(u,x,(i,j)) for j in range(4)],axis=1) for i in range(4)],axis=1)
    base=r.bilinear_jet(dict(g=p['g']),z,fr,val,du)
    iv=z['invgamma'];dg=p['dg'][:,:,1:,1:];ddg=p['ddg'][:,:,:,1:,1:]
    ell=.5*np.einsum('bij,bmji->bm',iv,dg)
    div=-np.einsum('bij,brjk,bkl->bril',iv,dg,iv)
    dell=.5*(np.einsum('brij,bmji->brm',div,dg)+np.einsum('bij,brmji->brm',iv,ddg))
    root=np.sqrt(z['vol']);psi=val/root[:,None]
    dpsi=(du-.5*ell[:,:,None]*val[:,None,:])/root[:,None,None]
    ddpsi=(ddu-.5*ell[:,:,None,None]*du[:,None,:,:]-.5*ell[:,None,:,None]*du[:,:,None,:]+(.25*ell[:,:,None]*ell[:,None,:]-.5*dell)[:,:,:,None]*val[:,None,None,:])/root[:,None,None,None]
    conn=r.connection(fr);cov=dpsi+np.einsum('bmij,bj->bmi',conn,psi)
    dconn=[];dInv=[];h=1e-24
    for mu in range(4):
        pp={k:v.astype(complex) for k,v in p.items()}
        pp['g']+=1j*h*p['dg'][:,mu];pp['dg']+=1j*h*p['ddg'][:,mu]
        ff=extended_frame(pp,c.geometry(pp))
        dconn.append(r.connection(dict(low=ff['low'].imag/h)))
        dInv.append(ff['inv'].imag/h)
    dconn=np.stack(dconn,axis=1);dInv=np.stack(dInv,axis=1)
    dcov=ddpsi+np.einsum('brmij,bj->brmi',dconn,psi)+np.einsum('bmij,brj->brmi',conn,dpsi)
    Cup=np.einsum('bmA,Aij->bmij',fr['inv'],r.ALPHA);dCup=np.einsum('brmA,Aij->brmij',dInv,r.ALPHA)
    Clow=np.einsum('bmn,bnij->bmij',p['g'],Cup)
    dClow=np.einsum('brmn,bnij->brmij',p['dg'],Cup)+np.einsum('bmn,brnij->brmij',p['g'],dCup)
    dZ=np.einsum('bra,bmab,bnb->brmn',dpsi.conj(),Clow,cov) if False else None
    dZ=np.einsum('bri,bmij,bnj->brmn',dpsi.conj(),Clow,cov)+np.einsum('bi,brmij,bnj->brmn',psi.conj(),dClow,cov)+np.einsum('bi,bmij,brnj->brmn',psi.conj(),Clow,dcov)
    dscalar=2*np.real(np.einsum('bri,ij,bj->br',dpsi.conj(),r.BETA,psi))
    dlag=-np.imag(np.einsum('bri,bmij,bmj->br',dpsi.conj(),Cup,cov)+np.einsum('bi,brmij,bmj->br',psi.conj(),dCup,cov)+np.einsum('bi,bmij,brmj->br',psi.conj(),Cup,dcov))-r.MASS*dscalar
    dT=.5*np.imag(dZ+dZ.swapaxes(-1,-2))+p['dg']*base['lagrangian'][:,None,None,None]+p['g'][:,None,:,:]*dlag[:,:,None,None]
    T=base['stress'];covT=dT-np.einsum('blrm,bln->brmn',z['Gamma'],T)-np.einsum('blrn,bml->brmn',z['Gamma'],T)
    divergence=np.einsum('brm,brmn->bn',z['ig'],covT)
    return dict(stress=T,derivative=dT,divergence=divergence,dirac=base['dirac_residual'],frame_error=frame_error)

def direct_stress(bg,u,x):
    p=bg.jets(x)
    # Only first metric jets are needed by the original frame and stress.
    zero=dict(ddg=np.zeros((len(x),4,4,4,4)))
    z=c.geometry(dict(p,**zero));fr=r.shared.frame_data(dict(g=p['g']),z)
    val=u.evaluate(x);du=np.stack([u.evaluate(x,i) for i in range(4)],axis=1)
    return r.bilinear_jet(dict(g=p['g']),z,fr,val,du)['stress']

def run():
    start=time.time();bg,_,_=c.ex.previous.rebuilt.build_pair();u,_,_=weighted.receiver_modes();src=c.ex.loop.LoopSource(bg,2,8);x=src.points[::64]
    data=stress_jet(bg,u,x);value_error=c.maximum(data['stress']-direct_stress(bg,u,x));rows=[]
    for step in (2e-6,1e-6):
        errs=[]
        for mu in range(4):
            xp=x.copy();xm=x.copy();xp[:,mu]+=step;xm[:,mu]-=step
            fd=(direct_stress(bg,u,xp)-direct_stress(bg,u,xm))/(2*step)
            errs.append(c.maximum(fd-data['derivative'][:,mu]))
        rows.append(dict(step=step,coordinate_derivative_differences=errs))
    assert value_error<1e-12 and data['frame_error']<1e-11
    assert max(rows[-1]['coordinate_derivative_differences'])<1e-6,rows
    return dict(round=917,status='working',date='2026-10-06',samples=len(x),original_receiver_and_full_offshell_stress_retained=True,
      independent_original_stress_error=value_error,original_frame_error=data['frame_error'],
      actual_stress_divergence_sample_max=c.maximum(data['divergence']),actual_stress_divergence_components_max=np.max(abs(data['divergence']),axis=0).tolist(),
      actual_Dirac_residual_sample_max=c.maximum(data['dirac']),actual_stress_first_derivative_sample_max=c.maximum(data['derivative']),
      independent_central_difference_diagnostics=rows,finite_difference_used_only_for_check=True,
      source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),STAGE/'916/covariant_joint_residual.py',STAGE/'906/compact_receiver_propagation.py',STAGE/'914/weighted_receiver_pairing.py')},
      uniform_source_divergence_bound_certified=False,actual_finite_observable_error_certified=False,full_goal_completed=False,elapsed_seconds=round(time.time()-start,3))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))
