"""904: exact-direction complex-step differentiation of the full902 evolution.
The source902 is read, never changed. Only default zeros dtype and FFT extension
are adapted to retain infinitesimal imaginary parts. No claim of a physical Green
inverse, quantization, or rigorous continuum error is made by this numerical JVP.
"""
from pathlib import Path
import sys,ast,types,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'902'))
import common_background_evolution as base
SOURCE=Path(base.__file__)
class ComplexZeros(ast.NodeTransformer):
    def visit_Call(self,node):
        node=self.generic_visit(node)
        if isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name) and node.func.value.id=='np' and node.func.attr=='zeros' and len(node.args)==1 and not any(k.arg=='dtype' for k in node.keywords):
            node.keywords.append(ast.keyword(arg='dtype',value=ast.Name(id='complex',ctx=ast.Load())))
        return node
module=types.ModuleType('_round904_analytic902');module.__file__=str(SOURCE)
code=ast.fix_missing_locations(ComplexZeros().visit(ast.parse(SOURCE.read_text('utf-8'))))
exec(compile(code,str(SOURCE),'exec'),module.__dict__)
class AnalyticModel(module.Model):
    def grad(self,a):
        if not np.iscomplexobj(a):return base.Model.grad(self,a)
        return base.Model.grad(self,a.real)+1j*base.Model.grad(self,a.imag)
    def lapweighted(self,a,ig):
        ar=np.fft.fftn(a.real,axes=(0,1,2));ai=np.fft.fftn(a.imag,axes=(0,1,2));out=np.zeros_like(a,dtype=complex)
        for i in range(3):
            for j in range(3):
                k=(self.k[...,i]*self.k[...,j]).reshape(self.k.shape[:3]+(1,)*(a.ndim-3))
                dij=np.fft.ifftn(-k*ar,axes=(0,1,2)).real+1j*np.fft.ifftn(-k*ai,axes=(0,1,2)).real
                out+=ig[...,i,j,None,None]*dij
        return out
    def jvp(self,y,u,h=1e-24):
        return {k:v.imag/h for k,v in self.rhs({k:y[k].astype(complex)+1j*h*u[k] for k in y}).items()}

class Family:
    def __init__(self,N,old):
        self.model=base.Model(N,old);self.analytic=AnalyticModel(N,old);self.old=old
        self.q,self.p0,_,_,_=old.completed(N,1.)
        self.base_y,_=self.model.initial()
    def data(self,epsilon=.02,with_tangent=False):
        model=self.model;geo=model.geo;q=self.q;x=q['grid'][...,0]
        new=dict(q);new['B']=q['B']+epsilon**2*np.cos(x)**2;new['C']=q['C']-epsilon**2*np.sin(x)**2;new['U']=q['U']+.5*epsilon**2*np.sin(x)**2
        psi,At,stats=geo.solve_hamiltonian(new,initial=self.p0);tau=np.sqrt(q['tau2'])
        y={k:v.copy() for k,v in self.base_y.items()}
        gam=psi[...,None,None]**4*np.eye(3);K=psi[...,None,None]**-2*At+tau/3*gam
        y['g'][...,1:,1:]=gam;y['gd'][...,1:,1:]=-2*K
        y['gd'][...,0,1:]=y['gd'][...,1:,0]=-2*model.grad(psi)/psi[...,None]
        y['phi'][...,5]=epsilon*np.sin(x)
        if not with_tangent:return y
        AA=np.sum(At*At,axis=(-1,-2))+q['pKp']
        J=5*new['C']*psi**4-new['B']+7*AA*psi**-8+6*new['Y']*psi**-4
        source=2*epsilon*(np.sin(x)**2*psi**5+np.cos(x)**2*psi)
        k2=np.sum(geo.waves(q['N'])**2,axis=-1)
        op=lambda u:-8*geo.laplace(u)+J*u
        pre=lambda u:geo.ifft(geo.fft(u)/(8*k2+np.mean(J)))
        dp,it=geo.cg(op,source,pre,tol=2e-13)
        dgam=4*psi[...,None,None]**3*dp[...,None,None]*np.eye(3)
        dK=-2*psi[...,None,None]**-3*dp[...,None,None]*At+tau/3*dgam
        u={k:np.zeros_like(v) for k,v in y.items()};u['g'][...,1:,1:]=dgam;u['gd'][...,1:,1:]=-2*dK
        u['gd'][...,0,1:]=u['gd'][...,1:,0]=-2*(model.grad(dp)/psi[...,None]-model.grad(psi)*dp[...,None]/psi[...,None]**2)
        u['phi'][...,5]=np.sin(x)
        return y,u,dict(elliptic_tangent_residual=float(np.max(abs(op(dp)-source))),cg_iterations=it,minimum_J=float(J.min()),epsilon=epsilon)

def constraints(model,y):
    # Original902 fields before taking norms; complex extension can differentiate them.
    z=model.geometry(y);m=model.matter(y,z);gam=z['gamma'];iv=z['invgamma'];K=z['K']
    dgam=model.grad(gam);ZZ=dgam+np.swapaxes(dgam,-3,-2)-np.moveaxis(dgam,-3,-1)
    C=.5*np.einsum('...rs,...mns->...rmn',iv,ZZ);dc=model.grad(C)
    Ric=np.einsum('...kkij->...ij',dc)-np.einsum('...jkik->...ij',dc)+np.einsum('...kkl,...lij->...ij',C,C)-np.einsum('...kjl,...lik->...ij',C,C)
    R=np.einsum('...ij,...ij->...',iv,Ric);Km=np.einsum('...ij,...jk->...ik',iv,K);tau=np.trace(Km,axis1=-2,axis2=-1)
    H=R+tau**2-np.einsum('...ij,...ji->...',Km,Km)-2*m['rho']
    dkm=model.grad(Km);momentum=np.einsum('...jji->...i',dkm)+np.einsum('...jjk,...ki->...i',C,Km)-np.einsum('...kji,...jk->...i',C,Km)-model.grad(tau)+m['M']/z['vol'][...,None]
    gauss=np.einsum('...iia->...a',model.grad(y['E']))+sum(module.bracket(y['A'][...,i,:],y['E'][...,i,:]) for i in range(3))+np.einsum('...A,...aA->...a',y['pi'],m['rp'])
    return dict(H=H,M=momentum,Gauss=gauss,harmonic=z['C'])

def tangent_constraints(model,y,u):
    h=1e-24;c=constraints(model,{k:y[k].astype(complex)+1j*h*u[k] for k in y})
    return {k:float(np.max(abs(v.imag/h))) for k,v in c.items()}

def step(model,analytic,y,u,dt):
    a=model.rhs(y);ua=analytic.jvp(y,u)
    yy=base.add(y,a,dt/2);uu=base.add(u,ua,dt/2);b=model.rhs(yy);ub=analytic.jvp(yy,uu)
    yy=base.add(y,b,dt/2);uu=base.add(u,ub,dt/2);c=model.rhs(yy);uc=analytic.jvp(yy,uu)
    yy=base.add(y,c,dt);uu=base.add(u,uc,dt);d=model.rhs(yy);ud=analytic.jvp(yy,uu)
    return ({k:y[k]+dt/6*(a[k]+2*b[k]+2*c[k]+d[k]) for k in y},
            {k:u[k]+dt/6*(ua[k]+2*ub[k]+2*uc[k]+ud[k]) for k in u})

def evolve(family,steps=8,Tend=.01):
    y,u,stats=family.data(with_tangent=True);initial=tangent_constraints(family.analytic,y,u)
    u0={k:v.copy() for k,v in u.items()}
    for _ in range(steps):y,u=step(family.model,family.analytic,y,u,Tend/steps)
    return dict(N=family.model.N,steps=steps,Tend=Tend,initial_solver=stats,initial_linear_constraints=initial,
        final_linear_constraints=tangent_constraints(family.analytic,y,u),
        final_tangent_max={k:float(np.max(abs(v))) for k,v in u.items()},
        induced_H5=float(np.max(abs(u['phi'][...,:5]))),
        tangent_change_max={k:float(np.max(abs(u[k]-u0[k]))) for k in u}),y,u
