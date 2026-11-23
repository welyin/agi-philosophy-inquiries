"""Independent audit of round 1339: Choi readout, coherent SWAP, Gram masses and weak currents."""
import argparse
import json
from pathlib import Path
import numpy as np


def unit(vector):
    return vector/np.linalg.norm(vector)


def kron(items):
    v=np.ones(1)
    for x in items: v=np.kron(v,x)
    return v


def j12(v):
    return float((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)


def check_choi():
    rng=np.random.default_rng(1339)
    u,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    partner=unit(rng.normal(size=(3,5))+1j*rng.normal(size=(3,5)))
    # Axes x,p,o,i: x is the unknown port with partner p, o/i the Choi resource.
    state=np.einsum('xp,oi->xpoi',partner,u/np.sqrt(3))
    bell=np.eye(3)/np.sqrt(3)
    success=np.einsum('xi,xpoi->op',bell,state)
    assert np.allclose(success,u@partner/3)
    p=np.vdot(success,success).real
    wedge=np.zeros((9,3))
    for k,(i,j) in enumerate([(1,2),(2,0),(0,1)]):
        wedge[3*i+j,k]=1/np.sqrt(2)
        wedge[3*j+i,k]=-1/np.sqrt(2)
    dual=wedge.T@np.kron(u,u)@wedge
    assert np.allclose(dual,np.linalg.det(u)*u.conj())
    return dict(success=float(p),failure=float(1-p),dual_error=float(np.linalg.norm(dual-np.linalg.det(u)*u.conj())))


def formed(v):
    # Output axes (branch, data, fundamental 0, dual 0, fundamental 1, dual 1, fundamental 2, dual 2).
    left=np.empty((4,3,729,3),complex)
    refs=kron([vec for w in v for vec in [w,w.conj()]])
    for i in range(3):
        x=np.kron(np.eye(3)[:,i],refs).reshape((3,)*7)
        for branch in range(4):
            y=x if branch==0 else x.swapaxes(0,2*branch-1)
            left[branch,:,:,i]=y.reshape(3,729)/2
    assert np.allclose(left.reshape(-1,3).conj().T@left.reshape(-1,3),np.eye(3))
    return left


def predicted_current(x, y, tu, td, c):
    d,a=np.linalg.eigh(x)
    e,b=np.linalg.eigh(y)
    def coeff(s,t):
        total=1+t*t
        z=2*(t*s)**2/(total+np.sqrt(total*total-4*(t*s)**2))
        k=t*t*s/(total-z)
        return k,np.sqrt(1+k*k-2*k*s),np.sqrt(z)
    ku,nu,mu=coeff(d,tu)
    kd,nd,md=coeff(e,td)
    factor=(1-(ku*d)[:,None]-(kd*e)[None,:]+c*ku[:,None]*kd)/(nu[:,None]*nd)
    return (a.conj().T@b)*factor,mu,md


def direct_current(x,y,tu,td,c):
    # The nine-dimensional span covers W and the two port isometries.
    # All other complete left modes have the common bare mass; they are orthogonal spectators.
    eye=np.eye(3)
    null=np.zeros((3,3))
    pu=np.vstack([eye,null,null])
    s=np.sqrt(1-c*c)
    pd=np.vstack([c*eye,s*eye,null])
    residual=(y-c*x)/s
    v,q=np.linalg.eigh(eye-x.conj().T@x-residual.conj().T@residual)
    assert min(v)>0
    w=np.vstack([x,residual,(q*np.sqrt(v))@q.conj().T])
    complement=np.linalg.qr(w,mode='complete')[0][:,3:]
    def light(p,t):
        left,sv,right=np.linalg.svd(np.column_stack([t*p,complement]))
        return left[:,::-1][:,:3],sv[::-1][:3]
    lu,mu=light(pu,tu)
    ld,md=light(pd,td)
    return lu.conj().T@ld,mu,md


def readouts(weights,v,w,mats):
    chi=unit(np.asarray(weights,float))
    bells=kron([np.eye(3).reshape(-1)/np.sqrt(3)]*3)
    result=np.einsum('a,r,axri->xi',chi.conj(),bells.conj(),w)
    target=sum(weights[i]*mats[i] for i in range(4))
    scale=1/(2*np.sqrt(27)*np.linalg.norm(weights))
    assert np.allclose(result,scale*target)
    return chi,result,target


def calculate():
    u=np.diag([1,1j,-1])
    r=np.sqrt([1/2,1/3,1/6])
    v=np.column_stack([r,u@r,u@u@r])
    g=v.conj().T@v
    cyc=g[0,1]*g[1,2]*g[2,0]
    a,b=abs(g[0,1])**2,abs(g[0,2])**2
    detg=np.linalg.det(g).real
    assert np.allclose([a,b,cyc.real,cyc.imag,detg],[2/9,1/9,0,2/27,4/9])
    w=formed(v.T)
    mats=[np.eye(3)]+[np.outer(vec,vec.conj()) for vec in v.T]
    cu,x,A=readouts([1,1,2,0],v,w,mats)
    comparisons=[]
    for name,weights,coef in [('new_direction',[1,0,1,3],-36),('shared_direction',[1,0,3,1],-12)]:
        cd,y,B=readouts(weights,v,w,mats)
        c=cu@cd
        ev_a,va=np.linalg.eigh(A)
        ev_b,vb=np.linalg.eigh(B)
        vm=va.conj().T@vb
        diff=np.prod([e[j]-e[i] for e in [ev_a,ev_b] for j in range(3) for i in range(j)])
        detc=np.linalg.det(A@B-B@A)
        assert np.isclose(detc,1j*coef*detg*cyc.imag)
        assert np.isclose(j12(vm),detc.imag/(2*diff))
        hc_a=A@A; hc_b=B@B
        detc2=np.linalg.det(hc_a@hc_b-hc_b@hc_a)
        squared_diff=np.prod([e[j]**2-e[i]**2 for e in [ev_a,ev_b] for j in range(3) for i in range(j)])
        assert np.isclose(j12(vm),detc2.imag/(2*squared_diff))
        exact,mup,mdp=predicted_current(x,y,1/3,.25,c)
        direct,mud,mdd=direct_current(x,y,1/3,.25,c)
        assert np.allclose(abs(direct),abs(exact))
        assert np.allclose([mup,mdp],[mud,mdd])
        assert abs(j12(direct)-j12(exact))<1e-10
        comparisons.append(dict(name=name,portal_overlap=c,up_spectrum=ev_a.tolist(),down_spectrum=ev_b.tolist(),
                    leading_joint_trace=float(np.trace(A@B).real),comm_det_imag=detc.imag,
                    squared_comm_det_imag=detc2.imag,leading_J=j12(vm),finite_J12=j12(direct),
                    finite_up_masses=mud.tolist(),finite_down_masses=mdd.tolist(),finite_left_abs=abs(direct).tolist(),
                    finite_left_singulars=np.linalg.svd(direct,compute_uv=False).tolist(),
                    max_finite_abs_error=float(np.max(abs(abs(direct)-abs(exact))))))
    assert np.isclose(comparisons[0]['leading_J'],3*comparisons[1]['leading_J'])
    return dict(choi=check_choi(),gram=dict(a=a,b=b,t=cyc.real,h=cyc.imag,determinant=detg),comparisons=comparisons)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a: compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,str): assert a==b
    else: assert np.isclose(a,b,atol=1e-9,rtol=2e-8),(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    res=calculate()
    dest=Path(__file__).with_name('independent_results.json')
    if args.check:compare(res,json.loads(dest.read_text(encoding='utf-8')))
    else:
        with dest.open('x',encoding='utf-8') as f:
            json.dump(res,f,indent=2,allow_nan=False); f.write('\n')
    print([(s['name'],s['leading_J'],s['finite_J12']) for s in res['comparisons']])
