"""902 working: original classical background, harmonic Einstein + canonical YM/H5/KG.
Floating spectral/RK diagnostics, not a validated time-dependent error enclosure.
Fermion/receiver backgrounds vanish, as in859/860; no quantum modes are solved here.
"""
from pathlib import Path
import sys,json,argparse,time
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'common_background_evolution_results.json'

def generators():
    T=[]
    for i,j in ((0,1),(0,2),(1,2)):
        a=np.zeros((3,3),complex);a[i,j]=a[j,i]=.5
        b=np.zeros((3,3),complex);b[i,j]=-.5j;b[j,i]=.5j
        T.extend((a,b))
    d3=np.diag([.5,-.5,0]);d8=np.diag([1.,1.,-2.])/np.sqrt(12)
    T=np.array([T[0],T[1],d3,T[2],T[3],T[4],T[5],d8])
    terms=[]
    for a in range(8):
        for b in range(8):
            bracket=1j*(T[a]@T[b]-T[b]@T[a])
            for c in range(8):
                q=float((2*np.trace(T[c]@bracket)).real)
                if abs(q)>1e-13:terms.append((c,a,b,q))
    sig=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                q=-.5j*np.trace(sig[c]@(sig[a]@sig[b]-sig[b]@sig[a]))
                # structure coefficient of i[T_a,T_b] is minus epsilon.
                if abs(q)>1e-13:terms.append((8+c,8+a,8+b,-float(q.real)/2))
    R=np.zeros((12,6,6))
    for a in range(4):
        U=1j*sig[a]/2 if a<3 else 3j*np.eye(2)
        R[8+a,:4,:4]=np.block([[U.real,-U.imag],[U.imag,U.real]])
    return T,terms,R
T,TERMS,REP=generators()
def bracket(a,b):
    out=np.zeros(np.broadcast_shapes(a.shape,b.shape))
    for c,i,j,f in TERMS:out[...,c]+=f*a[...,i]*b[...,j]
    return out

class Model:
    def __init__(self,N,old):
        self.N=N;self.old=old;self.geo=old.geo
        k=np.fft.fftfreq(N,d=1/N);self.k=np.stack(np.meshgrid(k,k,k,indexing='ij'),axis=-1)
        self.L,self.u,_=self.geo.old.scalar.parameters()
        par=self.geo.old.PAR;self.b=np.r_[np.full(8,par['b'][0]),np.full(3,par['b'][1]),par['b'][2]]
        self.K=1/(2*self.b)
    def grad(self,a):
        h=np.fft.fftn(a,axes=(0,1,2));out=[]
        for j in range(3):
            k=self.k[...,j].reshape(self.k.shape[:3]+(1,)*(a.ndim-3))
            out.append(np.fft.ifftn(1j*k*h,axes=(0,1,2)).real)
        return np.stack(out,axis=3)
    def lapweighted(self,a,ig):
        h=np.fft.fftn(a,axes=(0,1,2));out=np.zeros_like(a)
        for i in range(3):
            for j in range(3):
                k=(self.k[...,i]*self.k[...,j]).reshape(self.k.shape[:3]+(1,)*(a.ndim-3))
                dij=np.fft.ifftn(-k*h,axes=(0,1,2)).real
                out+=ig[...,i,j,None,None]*dij
        return out
    def initial(self):
        q,p0,_,_,_=self.old.completed(self.N,1.)
        x=q['grid'][...,0];new=dict(q);new['B']=q['B']+.0004*np.cos(x)**2;new['C']=q['C']-.0004*np.sin(x)**2;new['U']=q['U']+.0002*np.sin(x)**2
        psi,At,stats=self.geo.solve_hamiltonian(new,initial=p0)
        gamma=psi[...,None,None]**4*np.eye(3);K=psi[...,None,None]**-2*At+np.sqrt(q['tau2'])/3*gamma
        g=np.zeros(psi.shape+(4,4));g[...,0,0]=-1;g[...,1:,1:]=gamma
        gd=np.zeros_like(g);gd[...,1:,1:]=-2*K;gd[...,0,0]=2*np.sqrt(q['tau2'])
        gradp=self.grad(psi);gd[...,0,1:]=gd[...,1:,0]=-2*gradp/psi[...,None]
        phi=np.concatenate((q['phi'],(.02*np.sin(x))[...,None]),axis=-1)
        pi=np.concatenate((q['p'],np.zeros(x.shape+(1,))),axis=-1)
        A=np.zeros(psi.shape+(3,12));E=np.zeros_like(A)
        for i,a,aa,ee in ((0,0,.31,.23),(1,1,.27,.19),(2,3,.21,.17)):A[...,i,a]=aa;E[...,i,a]=ee
        A[...,8:11]=q['f']['a'];E[...,8:11]=q['f']['E'];A[...,11]=q['f']['a0'];E[...,11]=q['f']['E0']
        return dict(g=g,gd=gd,phi=phi,pi=pi,A=A,E=E),stats
    def geometry(self,y):
        g,gd=y['g'],y['gd'];ig=np.linalg.inv(g);gamma=g[...,1:,1:];invgamma=np.linalg.inv(gamma)
        det=np.linalg.det(gamma);assert np.min(det)>0 and np.max(ig[...,0,0])<0
        vol=np.sqrt(det);alpha=1/np.sqrt(-ig[...,0,0]);beta=-ig[...,0,1:]/ig[...,0,0,None]
        d=np.concatenate((gd[...,None,:,:],self.grad(g)),axis=3)
        # Z_{mu nu sigma}=d_mu g_{nu sigma}+d_nu g_{mu sigma}-d_sigma g_{mu nu}.
        Z=d+np.swapaxes(d,-3,-2)-np.moveaxis(d,-3,-1)
        Gamma=.5*np.einsum('...rs,...mns->...rmn',ig,Z)
        C=np.einsum('...ab,...abm->...m',ig,d)-.5*np.einsum('...ab,...mab->...m',ig,d)
        # Quadratic part of Ricci - nabla_(mu C_nu), with all second partials zero.
        di=-np.einsum('...ab,...mbc,...cd->...mad',ig,d,ig)
        dGamma=.5*np.einsum('...ars,...mns->...armn',di,Z)
        Rq=np.einsum('...rrmn->...mn',dGamma)-np.einsum('...nrmr->...mn',dGamma)
        Rq+=np.einsum('...rrl,...lmn->...mn',Gamma,Gamma)-np.einsum('...rnl,...lmr->...mn',Gamma,Gamma)
        dC=np.einsum('...mab,...abn->...mn',di,d)-.5*np.einsum('...mab,...nab->...mn',di,d)
        Q=Rq-.5*(dC+np.swapaxes(dC,-1,-2))+np.einsum('...lmn,...l->...mn',Gamma,C)
        db=self.grad(beta);dgamma=d[...,1:,1:,1:]
        lie=np.einsum('...k,...kij->...ij',beta,dgamma)+np.einsum('...ik,...jk->...ij',gamma,db)+np.einsum('...jk,...ik->...ij',gamma,db)
        K=-(gd[...,1:,1:]-lie)/(2*alpha[...,None,None])
        return dict(ig=ig,gamma=gamma,invgamma=invgamma,vol=vol,alpha=alpha,beta=beta,d=d,Gamma=Gamma,C=C,Q=Q,K=K)
    def matter(self,y,z):
        phi,pi,A,E=y['phi'],y['pi'],y['A'],y['E'];vol=z['vol'];iv=z['invgamma'];gam=z['gamma'];alpha=z['alpha'];beta=z['beta']
        ph=phi[...,:5];FF=2-np.sum(ph*ph,axis=-1)/6;assert np.min(FF)>0
        G=np.zeros(phi.shape[:-1]+(6,6));G[...,:5,:5]=np.eye(5)/FF[...,None,None]+ph[..., :,None]*ph[...,None,:]/(6*FF[...,None,None]**2);G[...,5,5]=1
        Gi=np.zeros_like(G);Gi[...,:5,:5]=FF[...,None,None]*(np.eye(5)-ph[..., :,None]*ph[...,None,:]/12);Gi[...,5,5]=1
        v=np.einsum('...ab,...b->...a',Gi,pi)/vol[...,None]
        rp=np.einsum('aAB,...B->...aA',REP,phi)
        D=self.grad(phi)+np.einsum('...ia,...aB->...iB',A,rp)
        da=self.grad(A);curv=da-np.swapaxes(da,-3,-2)+bracket(A[..., :,None,:],A[...,None,:,:])
        up=np.einsum('...ik,...jl,...kla->...ija',iv,iv,curv)
        e=np.einsum('...ij,...ja->...ia',gam,E)/(vol[...,None,None]*self.K)
        delta=np.stack((np.sum(ph[...,:4]**2,axis=-1)-self.u[0],ph[...,4]**2-self.u[1]),axis=-1)
        Ld=np.einsum('ab,...b->...a',self.L,delta);V=.25*np.sum(delta*Ld,axis=-1);U=V/FF**2+.5*phi[...,5]**2
        Vp=np.concatenate((Ld[...,0,None]*ph[...,:4],Ld[...,1,None]*ph[...,4,None]),axis=-1)
        Up=np.concatenate((Vp/FF[...,None]**2+2*V[...,None]*ph/(3*FF[...,None]**3),phi[...,5,None]),axis=-1)
        DG=np.einsum('...AB,...iB->...iA',G,D)
        gradnorm=np.einsum('...ij,...iA,...jA->...',iv,D,DG)
        velnorm=np.einsum('...A,...AB,...B->...',v,G,v)
        e2=np.einsum('...ij,...ia,...ja,a->...',iv,e,e,self.K)
        f2=np.einsum('...ija,...ija,a->...',curv,up,self.K)
        rho=.5*(velnorm+gradnorm)+U+.5*e2+.25*f2
        stress=np.einsum('...iA,...jA->...ij',D,DG)+gam*(.5*(velnorm-gradnorm)-U)[...,None,None]
        stress+=-np.einsum('...ia,...ja,a->...ij',e,e,self.K)+np.einsum('...kl,...ika,...jla,a->...ij',iv,curv,curv,self.K)+gam*(.5*e2-.25*f2)[...,None,None]
        M=np.einsum('...A,...iA->...i',pi,D)+np.einsum('...ja,...ija->...i',E,curv)
        trace=-rho+np.einsum('...ij,...ij->...',iv,stress)
        T=np.zeros(y['g'].shape);T[...,1:,1:]=stress
        T[...,0,1:]=T[...,1:,0]=alpha[...,None]*M/vol[...,None]+np.einsum('...j,...ji->...i',beta,stress)
        T[...,0,0]=alpha**2*rho+2*alpha*np.einsum('...i,...i->...',beta,M)/vol+np.einsum('...i,...j,...ij->...',beta,beta,stress)
        return dict(F=FF,G=G,Gi=Gi,v=v,rp=rp,D=D,curv=curv,up=up,e=e,U=U,Up=Up,DG=DG,rho=rho,stress=stress,M=M,T=T,trace=trace)
    def rhs(self,y,aux=False):
        z=self.geometry(y);m=self.matter(y,z);alpha=z['alpha'];beta=z['beta'];vol=z['vol'];A,E=y['A'],y['E'];phi=y['phi'];iv=z['invgamma']
        dgd=self.grad(y['gd']);wave=2*z['Q']-2*m['T']+y['g']*m['trace'][...,None,None]
        wave-=2*np.einsum('...i,...imn->...mn',z['ig'][...,0,1:],dgd)
        wave-=self.lapweighted(y['g'],z['ig'][...,1:,1:])
        gdd=wave/z['ig'][...,0,0,None,None]
        phidot=alpha[...,None]*m['v']+np.einsum('...i,...iA->...A',beta,m['D'])
        flux=(alpha*vol)[...,None,None]*np.einsum('...ij,...jA->...iA',iv,m['DG'])+beta[..., :,None]*y['pi'][...,None,:]
        dflux=self.grad(flux);pidot=np.einsum('...iiA->...A',dflux)+np.einsum('...ia,aAB,...iB->...A',A,REP,flux)
        # 1/2 d_A G_BC (v^B v^C - gamma^ij D_i phi^B D_j phi^C).
        ph=phi[...,:5];VV=m['v'][...,:5];DD=m['D'][...,:5];FF=m['F']
        tensor=VV[..., :,None]*VV[...,None,:]-np.einsum('...ij,...iA,...jB->...AB',iv,DD,DD)
        tr=np.trace(tensor,axis1=-2,axis2=-1);tp=np.einsum('...AB,...B->...A',tensor,ph);ptp=np.sum(ph*tp,axis=-1)
        force=ph*tr[...,None]/(6*FF[...,None]**2)+tp/(6*FF[...,None]**2)+ph*ptp[...,None]/(18*FF[...,None]**3)
        pidot[...,:5]+=(alpha*vol)[...,None]*force;pidot-=(alpha*vol)[...,None]*m['Up']
        Adot=alpha[...,None,None]*m['e']+np.einsum('...j,...jia->...ia',beta,m['curv'])
        fluxE=(alpha*vol)[...,None,None,None]*self.K*m['up']
        fluxE+=beta[..., :,None,None]*E[...,None,:,:]-beta[...,None,:,None]*E[..., :,None,:]
        dFE=self.grad(fluxE);Edot=np.einsum('...jjia->...ia',dFE)
        Edot+=sum(bracket(A[...,j,None,:],fluxE[...,j,:,:]) for j in range(3))
        Edot-=np.einsum('...iA,...aA->...ia',flux,m['rp'])
        out=dict(g=y['gd'],gd=gdd,phi=phidot,pi=pidot,A=Adot,E=Edot)
        return (out,z,m) if aux else out
    def constraints(self,y):
        z=self.geometry(y);m=self.matter(y,z);gam=z['gamma'];iv=z['invgamma'];K=z['K']
        dgam=self.grad(gam);ZZ=dgam+np.swapaxes(dgam,-3,-2)-np.moveaxis(dgam,-3,-1)
        C=.5*np.einsum('...rs,...mns->...rmn',iv,ZZ);dc=self.grad(C)
        Ric=np.einsum('...kkij->...ij',dc)-np.einsum('...jkik->...ij',dc)+np.einsum('...kkl,...lij->...ij',C,C)-np.einsum('...kjl,...lik->...ij',C,C)
        R=np.einsum('...ij,...ij->...',iv,Ric);Km=np.einsum('...ij,...jk->...ik',iv,K);tau=np.trace(Km,axis1=-2,axis2=-1)
        H=R+tau**2-np.einsum('...ij,...ji->...',Km,Km)-2*m['rho']
        dkm=self.grad(Km);momentum=np.einsum('...jji->...i',dkm)+np.einsum('...jjk,...ki->...i',C,Km)-np.einsum('...kji,...jk->...i',C,Km)-self.grad(tau)+m['M']/z['vol'][...,None]
        gauss=np.einsum('...iia->...a',self.grad(y['E']))+sum(bracket(y['A'][...,i,:],y['E'][...,i,:]) for i in range(3))+np.einsum('...A,...aA->...a',y['pi'],m['rp'])
        return dict(H_max=float(np.max(abs(H))),H_rms=float(np.sqrt(np.mean(H*H))),M_max=float(np.max(abs(momentum))),Gauss_max=float(np.max(abs(gauss))),harmonic_max=float(np.max(abs(z['C']))),minimum_metric_eigenvalue=float(np.min(np.linalg.eigvalsh(gam))),minimum_F=float(np.min(m['F'])))

def add(y,k,f):return {key:y[key]+f*k[key] for key in y}
def rk4(model,y,dt):
    a=model.rhs(y);b=model.rhs(add(y,a,dt/2));c=model.rhs(add(y,b,dt/2));d=model.rhs(add(y,c,dt))
    return {key:y[key]+dt/6*(a[key]+2*b[key]+2*c[key]+d[key]) for key in y}
def run(N=12,steps=4,Tend=.002):
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        model=Model(N,old);y,initial_stats=model.initial();first={k:v.copy() for k,v in y.items()};before=model.constraints(y);start=time.time()
        for n in range(steps):
            y=rk4(model,y,Tend/steps)
        after=model.constraints(y)
        return dict(N=N,steps=steps,Tend=Tend,initial=before,final=after,changes={k:float(np.max(abs(y[k]-first[k]))) for k in y},elapsed_seconds=time.time()-start),model,first,y
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=8);ap.add_argument('--steps',type=int,default=2);ap.add_argument('--time',type=float,default=.002);ap.add_argument('--write',action='store_true');a=ap.parse_args()
    r,_,_,_=run(a.N,a.steps,a.time)
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
