"""915: algebraic first/second coordinate jets for the old material extractor.
No finite-difference steps occur in xi, dxi or dalpha. Derivatives refer to the
original fixed finite Fourier/Hermite fields, not a certified physical solution.
"""
from pathlib import Path
import sys,itertools
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'914'))
import weighted_ward_pairing as prior
ex=prior.ex;previous=ex.previous;mb=ex.mb

class Jet:
    __array_priority__=1000
    def __init__(self,v,d=None,dd=None):
        self.v=np.asarray(v)
        self.d=np.zeros((4,)+self.v.shape,dtype=self.v.dtype) if d is None else np.asarray(d)
        self.dd=np.zeros((4,4)+self.v.shape,dtype=self.v.dtype) if dd is None else np.asarray(dd)
    @staticmethod
    def cast(x):return x if isinstance(x,Jet) else Jet(x)
    def __getitem__(self,key):
        key=key if isinstance(key,tuple) else (key,)
        return Jet(self.v[key],self.d[(slice(None),)+key],self.dd[(slice(None),slice(None))+key])
    def reshape(self,shape):return Jet(self.v.reshape(shape),self.d.reshape((4,)+tuple(shape)),self.dd.reshape((4,4)+tuple(shape)))
    def __neg__(self):return Jet(-self.v,-self.d,-self.dd)
    def __add__(self,other):
        b=Jet.cast(other);return Jet(self.v+b.v,np.stack([self.d[i]+b.d[i] for i in range(4)]),np.stack([np.stack([self.dd[i,j]+b.dd[i,j] for j in range(4)]) for i in range(4)]))
    __radd__=__add__
    def __sub__(self,other):return self+-Jet.cast(other)
    def __rsub__(self,other):return Jet.cast(other)+-self
    def __mul__(self,other):
        b=Jet.cast(other)
        return Jet(self.v*b.v,np.stack([self.d[i]*b.v+self.v*b.d[i] for i in range(4)]),np.stack([np.stack([self.dd[i,j]*b.v+self.d[i]*b.d[j]+self.d[j]*b.d[i]+self.v*b.dd[i,j] for j in range(4)]) for i in range(4)]))
    __rmul__=__mul__
    def __pow__(self,a):
        assert np.isscalar(a)
        v=self.v**a;fp=a*self.v**(a-1);fpp=a*(a-1)*self.v**(a-2)
        return Jet(v,fp*self.d,np.stack([np.stack([fp*self.dd[i,j]+fpp*self.d[i]*self.d[j] for j in range(4)]) for i in range(4)]))
    def __truediv__(self,other):return self*Jet.cast(other)**-1
    def __rtruediv__(self,other):return Jet.cast(other)*self**-1
    def __array_ufunc__(self,ufunc,method,*args,**kwargs):
        assert method=='__call__' and not kwargs
        if ufunc is np.sqrt:return Jet.cast(args[0])**.5
        if ufunc is np.negative:return -Jet.cast(args[0])
        a=Jet.cast(args[0]);b=Jet.cast(args[1])
        if ufunc is np.add:return a+b
        if ufunc is np.subtract:return a-b
        if ufunc is np.multiply:return a*b
        if ufunc is np.true_divide:return a/b
        raise TypeError(ufunc)
    def __array_function__(self,func,types,args,kwargs):
        if func is np.sum:
            a=Jet.cast(args[0]);axis=kwargs.get('axis',None);keep=kwargs.get('keepdims',False)
            f=lambda x:np.sum(x,axis=axis,keepdims=keep)
            return Jet(f(a.v),np.stack([f(x) for x in a.d]),np.stack([np.stack([f(x) for x in row]) for row in a.dd]))
        if func is np.stack:
            a=[Jet.cast(x) for x in args[0]];axis=kwargs.get('axis',0)
            f=lambda xs:np.stack(xs,axis=axis)
            return Jet(f([x.v for x in a]),np.stack([f([x.d[i] for x in a]) for i in range(4)]),np.stack([np.stack([f([x.dd[i,j] for x in a]) for j in range(4)]) for i in range(4)]))
        if func is np.einsum:
            spec=args[0];a=[Jet.cast(x) for x in args[1:]];values=[x.v for x in a]
            def term(replacements):
                xs=values.copy()
                for i,x in replacements.items():xs[i]=x
                return np.einsum(spec,*xs,**kwargs)
            v=term({});d=[];dd=[]
            for i in range(4):
                d.append(sum(term({n:x.d[i]}) for n,x in enumerate(a)))
                row=[]
                for j in range(4):
                    q=sum(term({n:x.dd[i,j]}) for n,x in enumerate(a))
                    for n,x in enumerate(a):
                        for m,y in enumerate(a):
                            if n!=m:q=q+term({n:x.d[i],m:y.d[j]})
                    row.append(q)
                dd.append(np.stack(row))
            return Jet(v,np.stack(d),np.stack(dd))
        if func is np.linalg.inv:
            a=Jet.cast(args[0]);v=np.linalg.inv(a.v)
            d=np.stack([-v@a.d[i]@v for i in range(4)])
            dd=np.stack([np.stack([v@a.d[j]@v@a.d[i]@v+v@a.d[i]@v@a.d[j]@v-v@a.dd[i,j]@v for j in range(4)]) for i in range(4)])
            return Jet(v,d,dd)
        if func is np.linalg.solve:
            a,b=map(Jet.cast,args);return np.einsum('...ij,...jk->...ik',np.linalg.inv(a),b)
        raise TypeError(func)

class Cache:
    def __init__(self,field,points):self.field=field;self.points=points;self.data={}
    def get(self,axes=()):
        key=tuple(sorted(axes))
        if key not in self.data:self.data[key]=previous.derivative(self.field,self.points,key)
        return self.data[key]
    def dual(self,axes=()):
        return Jet(self.get(axes),np.stack([self.get(axes+(i,)) for i in range(4)]),np.stack([np.stack([self.get(axes+(i,j)) for j in range(4)]) for i in range(4)]))
    def components(self,axes=()):return mb.unpack(self.get(axes))
    def jets(self,extra=()):
        g,phi,A=self.components(extra);dg,dphi,dA=mb.unpack(np.stack([self.get(extra+(i,)) for i in range(4)],axis=1))
        ddg,ddphi,ddA=mb.unpack(np.stack([np.stack([self.get(extra+(i,j)) for j in range(4)],axis=1) for i in range(4)],axis=1))
        return dict(g=g,phi=phi,A=A,dg=dg,dphi=dphi,dA=dA,ddg=ddg,ddphi=ddphi,ddA=ddA)
    def reference_inputs(self,extra=()):
        z=self.dual(extra);dz=np.stack([self.dual(extra+(i,)) for i in range(4)],axis=1)
        return dict(g=z[...,:16].reshape((len(self.points),4,4)),phi=z[...,16:22],dphi=dz[...,16:22])

def extractor(bg,field,points):
    bc=Cache(bg,points);vc=Cache(field,points);bp=bc.reference_inputs();vp=vc.reference_inputs()
    J=np.stack([ex.oldref.reference_delta(bp,bc.reference_inputs((mu,))) for mu in range(4)],axis=-1)
    dx=ex.oldref.reference_delta(bp,vp)
    xi=np.linalg.solve(J,dx[...,None])[...,0]
    p=bc.jets();v=vc.jets();x=xi.v;dxi=np.moveaxis(xi.d,0,1)
    _,rho,gram,Q,extract=ex.internal_data(p)
    dz=ex.feature_delta(p,v);lie=ex.lie_feature(p,x,dxi)
    residual={k:dz[k]-lie[k] for k in dz};alpha=extract(residual)
    # The rational/natural feature formulas are real-analytic. Complex-step
    # only differentiates these finite algebraic formulas; xi second jets
    # above supply the analytic derivative of dxi, avoiding nested imag calls.
    h=1e-24;da=[]
    for mu in range(4):
        dp=bc.jets((mu,));dv=vc.jets((mu,))
        pp={k:p[k].astype(complex)+1j*h*dp[k] for k in p}
        vv={k:v[k].astype(complex)+1j*h*dv[k] for k in v}
        xx=x.astype(complex)+1j*h*xi.d[mu]
        ddx=dxi.astype(complex)+1j*h*np.moveaxis(xi.dd[mu],0,1)
        _,_,_,_,ext=ex.internal_data(pp)
        dzz=ex.feature_delta(pp,vv);lzz=ex.lie_feature(pp,xx,ddx)
        da.append(ext({k:dzz[k]-lzz[k] for k in dzz}).imag/h)
    return (x,dxi,J.v,p,v,alpha,rho,residual),np.stack(da,axis=1),dict(background_multiindices=len(bc.data),response_multiindices=len(vc.data),xi_second_jet_symmetry_error=ex.maximum(xi.dd-xi.dd.swapaxes(0,1)))
