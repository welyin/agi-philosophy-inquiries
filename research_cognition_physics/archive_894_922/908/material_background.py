"""Actual902 short-time background sampled by a finite Fourier/Hermite interpolant.
The interpolant is a numerical device; no spacetime error certificate is claimed.
All off-shell variations act on spacetime g/phi/A, never on a substituted RHS.
"""
from pathlib import Path
import sys,time
import numpy as np
import magnetic_reference_source as source
ev=source.ev
CENTER=np.array([0.,0.,np.pi/2,np.pi/4])

def unpack(z):
    g=z[...,:16].reshape(z.shape[:-1]+(4,4));phi=z[...,16:22]
    A=np.zeros(z.shape[:-1]+(4,12),dtype=z.dtype);A[...,1:,:]=z[...,22:].reshape(z.shape[:-1]+(3,12))
    return g,phi,A

def pack(y):return np.concatenate((y['g'].reshape(y['g'].shape[:3]+(16,)),y['phi'],y['A'].reshape(y['A'].shape[:3]+(36,))),axis=-1)

class Background:
    def __init__(self,N=16,steps=2,T=.00025,drop=1e-11):
        self.N=N;self.T=T;self.steps=steps;self.drop=drop
        with ev.ResearchRuntime(ev.Layout()).installed():
            import joint_reference_constraint_strata as old
            model=ev.Model(N,old);y0,_=model.initial();ends=[]
            for sign in (-1,1):
                y={k:a.copy() for k,a in y0.items()}
                for _ in range(steps):y=ev.rk4(model,y,sign*T/steps)
                ends.append((pack(y),pack(model.rhs(y))))
            ym,dm=ends[0];yp,dp=ends[1]
            coefficients=np.stack(((yp+ym)/2-T*(dp-dm)/4,3*(yp-ym)/4-T*(dp+dm)/4,T*(dp-dm)/4,(T*(dp+dm)-(yp-ym))/4),axis=-2)
            ft=np.fft.fftn(coefficients,axes=(0,1,2)).reshape((-1,4,58))/N**3
            waves=model.k.reshape((-1,3));mask=np.max(abs(ft),axis=(1,2))>=drop
            self.tail0=float(np.max(np.sum(abs(ft[~mask]),axis=(0,1))))
            self.tail1=[float(np.max(np.sum(abs(ft[~mask])*abs(waves[~mask,i,None,None]),axis=(0,1)))) for i in range(3)]
            self.k=waves[mask];self.c=ft[mask]
            self.initial_constraints=model.constraints(y0)
        self.J0=self.jacobian(CENTER[None])[0];self.X0=self.reference(CENTER[None])[0]
    def evaluate(self,points,axis=None):
        points=np.asarray(points);theta=points[:,0]/self.T
        if axis==0:
            powers=np.stack((theta*0,theta*0+1,2*theta,3*theta*theta),axis=-1)/self.T
            coef=self.c
        else:
            powers=np.stack((theta*0+1,theta,theta*theta,theta**3),axis=-1)
            coef=self.c if axis is None else self.c*(1j*self.k[:,axis-1,None,None])
        result=[];flat=coef.reshape((len(self.k),-1))
        for start in range(0,len(points),128):
            q=points[start:start+128];phase=q[:,1:]@self.k.T
            z=np.cos(phase)@flat.real-np.sin(phase)@flat.imag
            result.append(np.einsum('bp,bpf->bf',powers[start:start+128],z.reshape((-1,4,58))))
        return np.concatenate(result,axis=0)
    def jets(self,points,variation=None,epsilon=0.):
        z=self.evaluate(points);g,phi,A=unpack(z)
        dz=np.stack([self.evaluate(points,i) for i in range(4)],axis=1)
        dg,dphi,dA=unpack(dz)
        out=dict(g=g,phi=phi,A=A,dg=dg,dphi=dphi,dA=dA)
        if variation is not None:
            v=variation(points,out)
            out={k:a+epsilon*v[k] if k in v else a for k,a in out.items()}
        return out
    def reference(self,points,variation=None,epsilon=0.):
        p=self.jets(points,variation,epsilon)
        return source.geometry(p['g'],p['phi'],p['dphi'],p['A'],p['dA'])['X']
    def jacobian(self,points,variation=None,epsilon=0.):
        points=np.asarray(points,dtype=float);cols=[];h=1e-24
        for mu in range(4):
            z=points.astype(complex);z[:,mu]+=1j*h
            cols.append(self.reference(z,variation,epsilon).imag/h)
        return np.stack(cols,axis=-1)
    def inverse(self,labels,guess,variation=None,epsilon=0.,iterations=8):
        x=np.array(guess,copy=True)
        for _ in range(iterations):
            r=self.reference(x,variation,epsilon)-labels;J=self.jacobian(x,variation,epsilon)
            dx=np.linalg.solve(J,r[...,None])[...,0];x-=dx
            if np.max(abs(dx))<2e-11:break
        r=self.reference(x,variation,epsilon)-labels
        final=np.max(abs(np.linalg.solve(self.J0,r.T).T))
        assert final<2e-8,final
        assert np.max(abs(x[:,0]))<self.T,(np.min(x[:,0]),np.max(x[:,0]))
        return x,dict(inverse_coordinate_residual=float(final),iterations=_+1,
            time_min=float(np.min(x[:,0])),time_max=float(np.max(x[:,0])),
            minimum_Jacobian_det=float(np.min(np.linalg.det(J))))
