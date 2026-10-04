"""Complete-sphere determinant majorant using actual five complex planes."""
import json
import math
import sys
import time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import joint_full_holonomy_reduction as old
import sphere_majorant_probe as reference


def polar_batch(h):
    odd=np.diag(old.b.prior.rep(np.eye(3),-np.eye(2),1)).real<0
    w=2*odd;delta=w*w+h;eta=delta/abs(delta);den=np.sqrt(2*w*w+2+2*abs(delta))
    v=np.empty((len(h),2,2,16),complex)
    v[:,0,0]=w*(1+eta)/den;v[:,1,1]=v[:,0,0]
    v[:,0,1]=(-1-eta*h.conj())/den;v[:,1,0]=(h+eta)/den
    physical=abs(np.prod(((1+v[:,0,0])*(1+v[:,1,1])-v[:,0,1]*v[:,1,0])/4,axis=1))**2
    return v,physical


def plane_coefficients(v):
    output=[]
    for j in range(5):
        t0,t1=old.b.internal.T[2*j:2*j+2]
        ii,jj=np.nonzero(t0);phase=t0[ii,jj].conj()*t1[ii,jj]
        assert len(ii)==16 and np.array_equal(abs(t0),abs(t1))
        aa=np.zeros(len(v));bb=np.zeros(len(v));cc=np.zeros(len(v));dd=np.zeros(len(v))
        for t in range(2):
            for s in range(2):
                x=((t==0 and s==0)+v[:,t,0,ii].conj()*v[:,s,0,jj].conj())/2
                y=((t==1 and s==1)+v[:,t,1,ii].conj()*v[:,s,1,jj].conj())/2
                aa+=np.sum(abs(x)**2,axis=1)/32;bb+=np.sum(abs(y)**2,axis=1)/32
                cc+=np.sum((x.conj()*y).real,axis=1)/32
                dd+=np.sum((x.conj()*y*phase).real,axis=1)/32
        output.append((aa,bb,cc,dd))
    return output


def sphere_moment(planes,n=16):
    size=len(planes[0][0]);poly=np.zeros((n+1,n+1,size));poly[0,0]=1
    for aa,bb,cc,dd in planes:
        det=aa*bb-cc*cc-dd*dd
        assert det.min()>-2e-15
        new=np.zeros_like(poly)
        for total in range(n+1):
            for a in range(total+1):
                b=total-a
                value=poly[a,b].copy()
                if a:value+=2*aa*new[a-1,b]
                if b:value+=2*bb*new[a,b-1]
                if a and b:value-=4*det*new[a-1,b-1]
                new[a,b]=value
        poly=new
    poch=[1]
    for k in range(n):poch.append(poch[-1]*(5+k))
    return math.factorial(n)/2**n*sum(poly[a,n-a]/(poch[a]*poch[n-a]) for a in range(n+1))


def function(angles):
    h,haar=old.probe.torus(angles);v,physical=polar_batch(h)
    return physical*sphere_moment(plane_coefficients(v))*haar


def check():
    rng=np.random.default_rng(69831);angles=rng.uniform(0,2*np.pi,(6,4));h,_=old.probe.torus(angles)
    vs,_=polar_batch(h);planes=plane_coefficients(vs);values=sphere_moment(planes);errors=[]
    for k,h0 in enumerate(h):
        _,v=old.polar(h0,True,True);_,q=reference.coefficients(v)
        expected=np.zeros((20,20))
        for j,(aa,bb,cc,dd) in enumerate(planes):
            indices=[2*j,2*j+1,10+2*j,10+2*j+1]
            expected[np.ix_(indices,indices)]=[[aa[k],0,cc[k],dd[k]],[0,aa[k],-dd[k],cc[k]],
                [cc[k],-dd[k],bb[k],0],[dd[k],cc[k],0,bb[k]]]
        err=float(np.max(abs(q-expected)));assert err<2e-15
        err2=float(abs(values[k]/reference.moment(q)-1));assert err2<2e-12
        errors.append([err,err2])
    return dict(cases=6,max_Q_error=max(x[0] for x in errors),max_moment_relative=max(x[1] for x in errors))


def integrate(nodes,batch=256):
    total=math.prod(nodes);acc=0.;maximum=0.
    for start in range(0,total,batch):
        index=np.arange(start,min(total,start+batch));angles=[]
        for n in nodes:
            angles.append((index%n)*2*np.pi/n);index//=n
        val=function(np.stack(angles,axis=1));acc+=val.sum();maximum=max(maximum,float(val.max()))
    return dict(nodes=nodes,points=total,value=acc/total,largest_node_value=maximum,
        quadrature_error_not_certified=True)


if __name__=='__main__':
    result=dict(check=check(),integrals=[])
    for nodes in ([12,12,12,48],[16,16,16,64],[20,20,20,80]):
        r=integrate(nodes);result['integrals'].append(r);print(json.dumps(r),flush=True)
    with (HERE/'block_majorant_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
