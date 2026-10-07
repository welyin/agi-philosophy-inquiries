"""908: full first-jet weak source of the original three material loop records.
Actual902 Fourier/Hermite background, fixed863-type four-dimensional compact
anchor profile, original (3,1)_{-2} and adjoint records. No thin-loop quantum limit.
"""
from pathlib import Path
import json,itertools
import numpy as np
import magnetic_reference_source as ms
from material_background import CENTER
HERE=Path(__file__).resolve().parent
ref=json.loads((HERE.parent/'907/material_reference_jets_results.json').read_text('utf-8'))['rows'][0]['data']['new']
DESIGN_X=np.array(ref['value']);DESIGN_J=np.array(ref['Jacobian'])
WIDTH=np.array([2e-5,2e-5,2e-5,.005]);L=.08
DIRECTIONS=np.array([[0.,1.,-1.,0.],[0.,0.,0.,1.]])
# Freeze the material curve using901's explicit initial conformal witness.
# A tangent parallelogram in material labels is not the image of a finite
# physical rectangle; the failed initial choice is preserved in drafts.
witness=json.loads((HERE.parent/'901/initial_psi_witness.json').read_text('utf-8'))
WK=np.array([r[:3] for r in witness['coefficients']],float);WK[:,2]*=2
WC=np.array([complex(r[3],r[4])/int(witness['denominator']) for r in witness['coefficients']])
S_COLOR=.008868915

def design_curve(segments):
    physical=L*DIRECTIONS;starts=np.array([np.zeros(4),physical[0],physical.sum(axis=0),physical[1]])
    deltas=np.array([physical[0],physical[1],-physical[0],-physical[1]])
    u=(np.arange(segments)+.5)/segments
    x=(CENTER[None,None,:]+starts[:,None,:]+u[None,:,None]*deltas[:,None,:]).reshape((-1,4))
    velocity=np.repeat(deltas,segments,axis=0)
    phase=x[:,1:]@WK.T;psi=np.cos(phase)@WC.real-np.sin(phase)@WC.imag
    dpsi=(-np.sin(phase)@((WC.real[:,None]*WK)) - np.cos(phase)@((WC.imag[:,None]*WK)))
    base_phase=CENTER[1:]@WK.T;psi0=np.cos(base_phase)@WC.real-np.sin(base_phase)@WC.imag
    a=x[:,1];labels=np.stack((a*0,.06*DESIGN_X[1]*np.sin(a),.02*np.sin(a),S_COLOR*(psi**-8-psi0**-8)),axis=-1)
    tangent=np.stack((a*0,.06*DESIGN_X[1]*np.cos(a)*velocity[:,1],.02*np.cos(a)*velocity[:,1],
        -8*S_COLOR*psi**-9*np.sum(dpsi*velocity[:,1:],axis=-1)),axis=-1)
    return labels,tangent,x

def bump_rule(n,base_order=128):
    x,w=np.polynomial.legendre.leggauss(base_order);w=w*np.exp(-1/(1-x*x));w/=w.sum()
    prev=np.zeros_like(x);p=np.ones_like(x);bet=0.;alpha=[];off=[]
    for i in range(n):
        a=np.sum(w*x*p*p);alpha.append(a);q=(x-a)*p-bet*prev
        if i+1<n:
            b=np.sqrt(np.sum(w*q*q));off.append(b);prev,p=p,q/b;bet=b
    jac=np.diag(alpha)+np.diag(off,1)+np.diag(off,-1);nodes,V=np.linalg.eigh(jac)
    return nodes,V[0]**2

def anchors(order):
    nodes,w=bump_rule(order);ids=np.array(list(itertools.product(range(order),repeat=4)))
    shifts=nodes[ids]*WIDTH;labels=DESIGN_X+shifts@DESIGN_J.T;weights=np.prod(w[ids],axis=1)
    return labels,weights,CENTER+shifts

def exponential(B):
    # B is Hermitian; return exp(-iB) and its differential with respect to B.
    e,V=np.linalg.eigh(B);adj=V.conj().swapaxes(-1,-2)
    E=(V*np.exp(-1j*e)[...,None,:])@adj
    dd=-1j*np.exp(-.5j*(e[..., :,None]+e[...,None,:]))*np.sinc((e[..., :,None]-e[...,None,:])/(2*np.pi))
    return E,(V,adj,dd)

def frechet(data,dB):
    V,adj,dd=data;return V@(dd*(adj@dB@V))@adj

def transports(E):
    n,s=E.shape[:2];I=np.broadcast_to(np.eye(3),(n,3,3)).astype(complex);before=[];U=I.copy()
    for i in range(s):before.append(U);U=E[:,i]@U
    after=[None]*s;R=I.copy()
    for i in range(s-1,-1,-1):after[i]=R;R=R@E[:,i]
    return U,np.stack(before,axis=1),np.stack(after,axis=1)

def connection(A):
    color=np.einsum('...ma,aij->...mij',A[...,:8],ms.ev.T)
    full=color-2*A[...,11,None,None]*np.eye(3)
    return color,full

class LoopSource:
    def __init__(self,bg,anchor_order=2,segments=8,variation=None,epsilon=0.):
        self.bg=bg;self.anchor_order=anchor_order;self.segments=segments
        anchor,weights,base=anchors(anchor_order);self.weights=weights;self.na=len(anchor);self.ns=4*segments
        ell,dell,physical=design_curve(segments)
        labels=(anchor[:,None,:]+ell[None,:,:]).reshape((-1,4))
        guess=(base[:,None,:]+physical[None,:,:]-CENTER).reshape((-1,4))
        self.points,self.inverse_stats=bg.inverse(labels,guess,variation,epsilon)
        J=bg.jacobian(self.points,variation,epsilon);self.Jinv=np.linalg.inv(J)
        self.velocity=np.einsum('bma,ba->bm',self.Jinv,np.tile(dell,(self.na,1)))
        self.jets=bg.jets(self.points,variation,epsilon);p=self.jets
        self.data=ms.geometry(p['g'],p['phi'],p['dphi'],p['A'],p['dA'])
        assert np.max(self.data['q'])<0,np.max(self.data['q'])
        color,full=connection(p['A']);self.channels=[]
        for A in (color,full):
            B=np.einsum('bm,bmij->bij',self.velocity,A).reshape((self.na,self.ns,3,3))/segments
            E,dE=exponential(B);U,before,after=transports(E);self.channels.append((U,before,after,dE))
        trc=np.trace(self.channels[0][0],axis1=-2,axis2=-1);trf=np.trace(self.channels[1][0],axis1=-2,axis2=-1)
        self.records=np.sum(weights[:,None]*np.stack((trf.real,trf.imag,abs(trc)**2-1),axis=-1),axis=0)
        self.stats=dict(**self.inverse_stats,clock_margin=float(np.min(-self.data['q'])),
            maximum_sampled_receiver_radius=float(np.max(np.linalg.norm(self.points[:,1:]-CENTER[None,1:],axis=-1))),
            anchor_count=self.na,path_sample_count=len(self.points))
    def kernels(self):
        color,full=self.channels;trc=np.trace(color[0],axis1=-2,axis2=-1)
        currents=np.zeros((3,len(self.points),4,12))
        for a in (*range(8),11):
            M=ms.ev.T[a] if a<8 else -2*np.eye(3)
            dEf=frechet(full[3],M/self.segments)
            dtf=np.trace(full[2]@dEf@full[1],axis1=-2,axis2=-1)
            if a<8:
                dEc=frechet(color[3],M/self.segments)
                dtc=np.trace(color[2]@dEc@color[1],axis1=-2,axis2=-1)
                da=2*np.real(trc.conj()[:,None]*dtc)
            else:da=np.zeros_like(dtf.real)
            values=np.stack((dtf.real,dtf.imag,da),axis=0)*self.weights[None,:,None]
            currents[...,a]=values.reshape((3,-1,1))*self.velocity[None]
        p=self.jets;self.currents=currents;self.ks=[];self.sources=[]
        for current in currents:
            k=-np.einsum('bma,bvma,bvr->br',current,self.data['F'],self.Jinv)
            self.ks.append(k);self.sources.append(ms.first_jet_source(self.data,p['phi'],p['dphi'],p['A'],current,k))
        return self.sources
    def response(self,variation):
        if not hasattr(self,'sources'):self.kernels()
        v=variation(self.points,self.jets);p=self.jets
        dx=ms.variation(self.data,p['phi'],p['dphi'],p['A'],v['g'],v['phi'],v['dphi'],v['A'],v['dA'])
        source_value=np.array([np.sum(ms.pair(s,v)) for s in self.sources])
        direct=np.array([np.sum(c*v['A']) for c in self.currents]);reference=np.array([np.sum(k*dx) for k in self.ks])
        assert np.max(abs(source_value-direct-reference))<2e-10
        return dict(total=source_value.tolist(),direct=direct.tolist(),reference=reference.tolist(),
            jet_pairing_identity_residual=float(np.max(abs(source_value-direct-reference))))

def zero_variation(points,p):return {k:np.zeros_like(p[k]) for k in ('g','phi','dphi','A','dA')}
def metric_variation(points,p):
    v=zero_variation(points,p);f=.1*(1+.2*np.cos(points[:,1]+points[:,3]))
    v['g'][:,1:,1:]=f[:,None,None]*np.eye(3);return v
def clock_variation(points,p):
    v=zero_variation(points,p);t,x,y,z=points.T
    v['phi'][:,1]=.01*(t+.003*np.sin(x)+.002*np.cos(y))
    v['dphi'][:,0,1]=.01;v['dphi'][:,1,1]=.00003*np.cos(x);v['dphi'][:,2,1]=-.00002*np.sin(y)
    return v
def probe_variation(points,p):
    v=zero_variation(points,p);x,z=points[:,1],points[:,3]
    v['phi'][:,5]=.02*np.cos(x)*np.cos(z)
    v['dphi'][:,1,5]=-.02*np.sin(x)*np.cos(z);v['dphi'][:,3,5]=-.02*np.cos(x)*np.sin(z)
    return v
def color_gauge_variation(points,p):
    v=zero_variation(points,p);phase=points[:,1]+.4*points[:,3];wave=np.array([0.,1.,0.,.4])
    alpha=np.zeros_like(p['A'][:,0]);alpha[:,1]=.1*np.sin(phase)
    da=np.zeros_like(p['A']);da[:,:,1]=.1*np.cos(phase)[:,None]*wave
    dda=np.zeros_like(p['dA']);dda[:,:,:,1]=-.1*np.sin(phase)[:,None,None]*wave[None,:,None]*wave[None,None,:]
    v['A']=ms.bracket(alpha[:,None,:],p['A'])-da
    v['dA']=ms.bracket(da[:,:,None,:],p['A'][:,None,:,:])+ms.bracket(alpha[:,None,None,:],p['dA'])-dda
    return v

def translation_variation(bg,vector):
    # Constant on the whole sampled support; a compact extension exists locally.
    vector=np.asarray(vector)
    def variation(points,p):
        v=zero_variation(points,p);v['g']=np.einsum('m,bmij->bij',vector,p['dg'])
        v['phi']=np.einsum('m,bmA->bA',vector,p['dphi']);v['A']=np.einsum('m,bmna->bna',vector,p['dA'])
        h=1e-24;z=points.astype(complex)+1j*h*vector[None,:];shift=bg.jets(z)
        v['dphi']=shift['dphi'].imag/h;v['dA']=shift['dA'].imag/h
        return v
    return variation

def affine_diffeo_variation(bg):
    # Affine on all loop support, smoothly extendible with a compact cutoff.
    scale=np.array([.01,.02,-.015,.025])
    def variation(points,p):
        xi=(points-CENTER)*scale;v=zero_variation(points,p)
        v['g']=np.einsum('br,brmn->bmn',xi,p['dg'])+(scale[None,:,None]+scale[None,None,:])*p['g']
        v['phi']=np.einsum('br,brA->bA',xi,p['dphi'])
        v['A']=np.einsum('br,brma->bma',xi,p['dA'])+scale[None,:,None]*p['A']
        h=1e-24;v['dphi']=scale[None,:,None]*p['dphi'];v['dA']=(scale[None,:,None,None]+scale[None,None,:,None])*p['dA']
        for r in range(4):
            z=points.astype(complex);z[:,r]+=1j*h;q=bg.jets(z)
            v['dphi']+=xi[:,r,None,None]*(q['dphi'].imag/h)
            v['dA']+=xi[:,r,None,None,None]*(q['dA'].imag/h)
        return v
    return variation
