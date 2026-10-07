"""909 actual908 weak kernel -> finite Fourier forcing -> full904 mixed response.
This is a numerical weak-source discretization, not a continuous source density,
a constraint-preserving Green construction, or the complete872 quantum feedback.
"""
from pathlib import Path
import sys, time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
for n in ('908','904','907'):sys.path.insert(0,str(STAGE/n))
import coupled_boson_tangent as tangent
import material_background as mb
import relational_loop_source as loop
import material_reference_jets as refs
ev=tangent.base
VOLUME=(2*np.pi)**3

def zeros(N):
    return {k:np.zeros((N,N,N)+shape) for k,shape in
            dict(g=(4,4),gd=(4,4),phi=(6,),pi=(6,),A=(3,12),E=(3,12)).items()}

def maxnorm(u):return {k:float(np.max(abs(a))) for k,a in u.items()}

class FourierMap:
    def __init__(self,N,points):
        self.N=N;self.points=points
        k=np.fft.fftfreq(N,d=1/N);self.k=np.stack(np.meshgrid(k,k,k,indexing='ij'),axis=-1).reshape((-1,3))
        self.keep=np.all(abs(self.k)<N/2,axis=1) if N%2==0 else np.ones(N**3,bool)
        self.phase=np.exp(-1j*self.k@points.T)
        self.phase[~self.keep]=0
    def deposit(self,weights):
        # weights include original anchor/path quadrature weights already.
        shape=weights.shape[1:]
        return (self.phase@weights.reshape((len(weights),-1))/VOLUME).reshape((self.N**3,)+shape)
    def grid(self,coeff):
        a=np.fft.ifftn(coeff.reshape((self.N,)*3+coeff.shape[1:])*self.N**3,axes=(0,1,2))
        assert np.max(abs(a.imag))<2e-10*max(1.,float(np.max(abs(a.real))))
        return a.real
    def evaluate(self,field,axis=None):
        f=np.fft.fftn(field,axes=(0,1,2)).reshape((self.N**3,-1))/self.N**3
        if axis is not None:f=f*(1j*self.k[:,axis,None])
        result=self.phase.conj().T@f
        return result.real.reshape((len(self.points),)+field.shape[3:])

def temporal_rule(times,order,halfwidth=2e-5):
    assert np.max(abs(times))<halfwidth
    nodes=np.sort(halfwidth*np.cos((2*np.arange(order)+1)*np.pi/(2*order)))
    w=np.ones((order,len(times)))
    for a in range(order):
        for b in range(order):
            if a!=b:w[a]*=(times-nodes[b])/(nodes[a]-nodes[b])
    assert np.max(abs(w.sum(axis=0)-1))<2e-14
    return nodes,w

def canonical_metric(points_g,Pg):
    ig=np.linalg.inv(points_g);mu=np.sqrt(-np.linalg.det(points_g))
    s=2*np.einsum('nma,nab,nbc->nmc',points_g,Pg,points_g)
    s=s/mu[:,None,None]
    tr=np.einsum('bij,bij->b',ig,s)
    return (-2*s+points_g*tr[:,None,None])/ig[:,0,0,None,None]

def build_forcing(src,N,channel=2,time_order=3,halfwidth=2e-5):
    if not hasattr(src,'sources'):src.kernels()
    s=src.sources[channel];fm=FourierMap(N,src.points[:,1:]);nodes,w=temporal_rule(src.points[:,0],time_order,halfwidth)
    metric=canonical_metric(src.jets['g'],s['g']);items=[]
    for a in range(time_order):
        wa=w[a];B=zeros(N);S=zeros(N)
        B['gd']=fm.grid(fm.deposit(metric*wa[:,None,None]))
        Pphi=fm.deposit(s['phi']*wa[:,None]);Qphi=fm.deposit(s['dphi'][:,1:,:]*wa[:,None,None])
        B['pi']=fm.grid(Pphi-1j*np.einsum('ki,kia->ka',fm.k,Qphi))
        PA=fm.deposit(s['A'][:,1:,:]*wa[:,None,None]);QA=fm.deposit(s['dA'][:,1:,1:,:]*wa[:,None,None,None])
        B['E']=fm.grid(PA-1j*np.einsum('kj,kjia->kia',fm.k,QA))
        S['pi']=fm.grid(fm.deposit(s['dphi'][:,0,:]*wa[:,None]))
        S['E']=fm.grid(fm.deposit(s['dA'][:,0,1:,:]*wa[:,None,None]))
        items.append((B,S))
    return nodes,items,fm,dict(channel=channel,record=float(src.records[channel]),time_order=time_order,
        time_halfwidth=halfwidth,source_sample_time_range=[float(src.points[:,0].min()),float(src.points[:,0].max())],
        source_points=len(src.points),N=N,even_grid_nyquist_planes_removed=N%2==0)

def projection_check(src,N=16,channel=2):
    nodes,items,fm,stats=build_forcing(src,N,channel,3)
    x=2*np.pi*np.arange(N)/N;X,Y,Z=np.meshgrid(x,x,x,indexing='ij');phase=X+2*Y-Z
    rng=np.random.default_rng(909);v=zeros(N)
    for key,a in v.items():
        c=rng.normal(size=a.shape[3:]);d=rng.normal(size=a.shape[3:])
        if key in ('g','gd'):c=(c+c.T)/2;d=(d+d.T)/2
        a[:]=np.sin(phase).reshape((N,)*3+(1,)*c.ndim)*c+np.cos(X+Z).reshape((N,)*3+(1,)*d.ndim)*d
    B={k:sum(row[0][k] for row in items) for k in v};S={k:sum(row[1][k] for row in items) for k in v}
    gridB=VOLUME/N**3*sum(np.sum(B[k]*v[k]) for k in v)
    gridS=VOLUME/N**3*sum(np.sum(S[k]*v[k]) for k in v)
    s=src.sources[channel];vv={k:fm.evaluate(a) for k,a in v.items()}
    direct=np.sum(canonical_metric(src.jets['g'],s['g'])*vv['gd'])+np.sum(s['phi']*vv['pi'])+np.sum(s['A'][:,1:]*vv['E'])
    for i in range(3):
        direct+=np.sum(s['dphi'][:,i+1]*fm.evaluate(v['pi'],i))+np.sum(s['dA'][:,i+1,1:]*fm.evaluate(v['E'],i))
    ds=np.sum(s['dphi'][:,0]*vv['pi'])+np.sum(s['dA'][:,0,1:]*vv['E'])
    # Independent analytic values for the same trigonometric test, without FFT.
    xp,yp,zp=src.points[:,1:].T;ph=xp+2*yp-zp
    value_error=0.
    rng=np.random.default_rng(909)
    for key,a in v.items():
        c=rng.normal(size=a.shape[3:]);d=rng.normal(size=a.shape[3:])
        if key in ('g','gd'):c=(c+c.T)/2;d=(d+d.T)/2
        q=np.sin(ph).reshape((-1,)+(1,)*c.ndim)*c+np.cos(xp+zp).reshape((-1,)+(1,)*d.ndim)*d
        value_error=max(value_error,float(np.max(abs(q-vv[key]))))
    return dict(spatial_B_pair_grid=float(gridB),spatial_B_pair_direct=float(direct),B_pair_error=float(abs(gridB-direct)),
        temporal_S_pair_grid=float(gridS),temporal_S_pair_direct=float(ds),S_pair_error=float(abs(gridS-ds)),
        independent_trig_value_error=value_error,source_normalization='coordinate density on [0,2pi)^3',
        source_kernel_is_not_a_grid_density=True)

def readonly_scalar_observations(model,analytic,y,u,radius=.025,order=2):
    # A fixed finite spatial smear of linear relational scalar variations.
    # Old reference X=(h,s,dh^2,ds^2), as785; pointwise gauge invariance is analytic.
    nodes,w=loop.bump_rule(order);ids=np.array(list(__import__('itertools').product(range(order),repeat=3)))
    points=mb.CENTER[None,1:]+radius*nodes[ids];weights=np.prod(w[ids],axis=1)
    fm=FourierMap(model.N,points);F=model.rhs(y)
    r=refs.references(analytic,y);dr=refs.derivative(analytic,y,u)
    directions=[F,*[{k:np.take(model.grad(a),i,axis=3) for k,a in y.items()} for i in range(3)]]
    dR=[refs.derivative(analytic,y,v) for v in directions]
    J=fm.evaluate(np.stack([q['old'] for q in dR],axis=-1))
    variation=fm.evaluate(dr['old']);eta=np.linalg.solve(J,variation[...,None])[...,0]
    dn=fm.evaluate(dr['new']);gradnew=fm.evaluate(np.stack([q['new'] for q in dR],axis=-1))
    rel=dn-np.einsum('bam,bm->ba',gradnew,eta)
    # h,s occur in both references and must cancel as a scalar consistency check.
    return dict(relational_probe=float(weights@rel[:,2]),relational_magnetic=float(weights@rel[:,3]),
        h_s_reference_cancellation=float(np.max(abs(rel[:,:2]))),sample_min_abs_old_Jdet=float(np.min(abs(np.linalg.det(J)))),
        radius=radius,quadrature_order=order,uniform_observer_chart_certified=False)

def propagate(model,analytic,y0,nodes,items,Tend=.00015,substeps=1,omit_temporal=False,observer=True):
    start=-2e-5; y=ev.rk4(model,y0,start);u=zeros(model.N);t=start;kicks=[]
    initial_background=model.constraints(y)
    for tn,(B,S) in zip(nodes,items):
        for _ in range(substeps):y,u=tangent.step(model,analytic,y,u,(tn-t)/substeps)
        LS=analytic.jvp(y,S);D={k:B[k]-(0 if omit_temporal else LS[k]) for k in u}
        u={k:u[k]+D[k] for k in u};kicks.append(dict(time=float(tn),B_max=maxnorm(B),S_max=maxnorm(S),LS_max=maxnorm(LS)));t=tn
    for _ in range(substeps):y,u=tangent.step(model,analytic,y,u,(Tend-t)/substeps)
    out=dict(Tend=Tend,substeps_per_interval=substeps,omitted_temporal_potential=omit_temporal,
        final_response_max=maxnorm(u),final_source_free_linear_constraints=tangent.tangent_constraints(analytic,y,u),
        initial_background_constraints=initial_background,final_background_constraints=model.constraints(y),kicks=kicks)
    if observer:out['observations']=readonly_scalar_observations(model,analytic,y,u)
    return out,y,u


def constraint_modes(analytic,y,u,cutoffs=(0,1,2)):
    h=1e-24;fields=tangent.constraints(analytic,{k:y[k].astype(complex)+1j*h*u[k] for k in y})
    out={}
    for key,a in fields.items():
        f=np.fft.fftn(a.imag/h,axes=(0,1,2))/analytic.N**3
        out[key]={str(c):float(np.max(abs(f[np.max(abs(analytic.k),axis=-1)<=c]))) for c in cutoffs}
    return out
