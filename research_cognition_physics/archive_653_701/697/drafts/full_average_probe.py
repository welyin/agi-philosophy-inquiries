"""Exploratory complete S9/full-group integration; not a strict sign proof.
Original centre histories, exact2x2 polar algebra, all16 channels retained.
"""
import importlib.util
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import joint_gauss_boundary_functional as b
spec=importlib.util.spec_from_file_location('entry697',HERE/'center_holonomy_entry.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)

def torus(angles):
    a,c,p,x=angles.T
    color=np.stack((np.exp(1j*a),np.exp(1j*c),np.exp(-1j*(a+c))),axis=-1)
    weak=np.stack((np.exp(1j*p),np.exp(-1j*p)),axis=-1);z=np.exp(1j*x)
    left=np.concatenate(((color[:,:,None]*weak[:,None,:]*z[:,None,None]).reshape(-1,6),
        color.conj()*z[:,None]**-4,color.conj()*z[:,None]**2,weak*z[:,None]**-3,
        z[:,None]**6,np.ones((len(z),1))),axis=1)
    j=b.mass.dictionary.dictionary();permutation=np.argmax(abs(j),axis=1)
    h=left[:,permutation]
    density=np.ones(len(z))
    for i in range(3):
        for k in range(i):density*=abs(color[:,i]-color[:,k])**2
    density*=abs(weak[:,0]-weak[:,1])**2/12
    return h,density

def values(h,fields,first,second):
    n=len(h);r=b.prior.rep(np.eye(3),-np.eye(2),1);odd=(np.diag(r).real<0)
    w0=2*first*odd;w1=2*second*odd
    det=w0*w1+h;phase=det/abs(det)
    norm=np.sqrt(w0*w0+w1*w1+2+2*abs(det))
    v=np.zeros((n,2,2,16),complex)
    v[:,0,0]=(w0+phase*w1)/norm;v[:,0,1]=(-1-phase*h.conj())/norm
    v[:,1,0]=(h+phase)/norm;v[:,1,1]=(w1+phase*w0)/norm
    physical=np.prod(((1+v[:,0,0])*(1+v[:,1,1])-v[:,0,1]*v[:,1,0])/4,axis=1)**2
    ts=np.einsum('nta,aij->ntij',fields,np.array(b.internal.T))
    aux=np.zeros((n,32,32),complex)
    for t in range(2):
        for s in range(2):
            block=np.zeros((n,16,16),complex)
            for k in range(2):block+=ts[:,k]*v[:,t,k].conj()[:,:,None]*v[:,s,k].conj()[:,None,:]
            if t==s:block+=ts[:,t]
            aux[:,16*t:16*(t+1),16*s:16*(s+1)]=block/2
    return physical*np.linalg.det(aux)

def crosscheck():
    rng=np.random.default_rng(69741);angles=rng.uniform(0,2*np.pi,(3,4));h,_=torus(angles)
    e=rng.normal(size=(3,2,10));e/=np.linalg.norm(e,axis=2)[:,:,None];errors=[]
    for first,second in ((False,False),(True,False),(True,True)):
        val=values(h,e,first,second)
        for k in range(3):
            original,gap=entry.weight(first,second,np.eye(16),np.diag(h[k]),e[k])
            err=float(abs(original-val[k])/max(abs(original),abs(val[k]),1e-30));errors.append(err)
            assert err<1e-9,(first,second,k,original,val[k],err)
    return max(errors)

def sample(count,seed=69751,batch=256):
    rng=np.random.default_rng(seed);sums=np.zeros(3);sq=np.zeros((3,3));mi=np.inf;ma=-np.inf;imag=0.
    for start in range(0,count,batch):
        n=min(batch,count-start);angles=rng.uniform(0,2*np.pi,(n,4));h,density=torus(angles)
        e=rng.normal(size=(n,2,10));e/=np.linalg.norm(e,axis=2)[:,:,None]
        z=np.stack([values(h,e,*flags) for flags in ((False,False),(True,False),(True,True))],axis=1)*density[:,None]
        imag=max(imag,float(abs(z.imag).max()));z=z.real
        sums+=z.sum(axis=0);sq+=z.T@z;mi=min(mi,float(z.min()));ma=max(ma,float(z.max()))
    mean=sums/count;cov=(sq/count-np.outer(mean,mean))/(count-1)
    return dict(n=count,mean=mean.tolist(),standard_error=np.sqrt(np.diag(cov)).tolist(),
        ratio=float(mean[1]/np.sqrt(mean[0]*mean[2])),covariance_of_mean=cov.tolist(),
        signed_range=[mi,ma],max_imaginary=imag,not_a_certified_integral=True)

if __name__=='__main__':
    print('crosscheck',crosscheck(),flush=True)
    for n in (4096,65536):print(sample(n),flush=True)
