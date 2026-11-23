"""Test neutral single-flavon monomials with jointly movable charged inputs."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
OLD=runpy.run_path(str(HERE/CFG["parent"]))
NU=OLD["NU"]
Q0=OLD["E0"]
L0=OLD["L0"]
QW=OLD["OLD"]["qrows"][:,1]
LW=OLD["OLD"]["lrows"][:,1]
PHI=-np.pi/6
W=np.exp(4j*PHI)
PORTS={"third_3_2":([2*PHI,3*PHI+np.angle(3+2*W)],1),
       "second_2_minus1":([10*PHI+np.angle(2-W),7*PHI],-1),
       "second_original_11_5deg":([10*PHI+np.radians(11.5),7*PHI],-1)}
WINDOWS=np.asarray(NU["read"](NU["INPUT"]["angle_windows"])["joint_windows"])


def edges(q,g):
    # Batched Cholesky in the real convention of the archived texture.
    d,sm,b,t12,t13,t23=q.T
    k=np.cos(t12)*np.tan(t23)*np.sin(t13)/np.sin(t12)
    gamma=np.radians(g)
    delta=gamma+np.arcsin(k*np.sin(gamma))
    p,u,t=np.sin(t12),np.sin(t13),np.sin(t23)
    cp,cu,ct=np.cos(t12),np.cos(t13),np.cos(t23)
    w=np.exp(1j*delta)
    v=np.empty((len(q),3,3),complex)
    v[:,0,:]=np.column_stack([cp*cu,p*cu,u/w])
    v[:,1,:]=np.column_stack([-p*ct-cp*t*u*w,cp*ct-p*t*u*w,t*cu])
    v[:,2,:]=np.column_stack([p*t-cp*ct*u*w,-cp*t-p*ct*u*w,ct*cu])
    h=(v*(np.column_stack([d/b,sm/b,np.ones(len(b))])**2)[:,None,:])@v.conj().transpose(0,2,1)
    left=np.exp(1j*np.column_stack([np.zeros(len(b)),np.angle(h[:,0,1]),np.angle(h[:,0,2])]))
    h=left[:,:,None]*h*left[:,None,:].conj()
    a=np.sqrt(h[:,0,0].real)
    c,e=h[:,0,1].real/a,h[:,0,2].real/a
    f=np.sqrt(h[:,2,2].real-e*e)
    return np.array([a,-d*sm/(b*b*a*f),c,e,f,(h[:,1,2]-c*e)/f])


def charged(edges,k,s):
    a,b,c,e,f,z=edges
    y=np.zeros((len(s),3,3),complex)
    t=s*np.exp(1j*k*PHI)
    y[:,1,0]=a
    y[:,0,1]=b
    y[:,1,1]=3*c
    y[:,2,1]=t*z
    y[:,1,2]=t*e
    y[:,2,2]=t*f
    return y


AD=runpy.run_path(str(HERE.parent/"1675/phase_box_check.py"))
Jet=AD["Jet"]


def box_domain(radius,orientation,transport_margin):
    I=OLD["OV"]["I"]
    def bounds(m,e):
        w=I.point(radius)*I.point(e)
        return I.point(m)+I(-w.hi,w.hi)
    raw=[bounds(*row) for row in OLD["OLD"]["qrows"]]
    raw=raw[:3]+[OLD["OV"]["trig"](r,True) for r in raw[3:]]
    raw.extend(bounds(*row) for row in OLD["OLD"]["lrows"])
    # Old 1675 encloses gamma plus a one-degree allowance including sign reversal.
    if transport_margin:
        angles=tuple(AD["domain"]()[0 if orientation==1 else 1][-1])
    else:
        angles=(.413,.486) if orientation==1 else (-.487,-.414)
        pi=OLD["OV"]["bounds"]("3.14159265358979323846","3.14159265358979323847")
        endlo=OLD["OV"]["trig"](pi*I.point(131)/360,False).values()[0]
        endhi=OLD["OV"]["trig"](pi*I.point(61)/180,False).values()[1]
        assert .413<endlo-.001001 and endhi<.486 and -.487 < -endhi-.001001 and -endlo < -.414
    return np.vstack([[[np.nextafter(x.values()[0],-np.inf),np.nextafter(x.values()[1],np.inf)]
                       for x in raw],angles])


def residual_jet(box,rising,k,orientation):
    variables=[]
    for j in range(10):
        der=np.zeros((len(box),10))
        der[:,j]=1.
        variables.append(Jet((box[:,j,:1],box[:,j,1:]),(der,der)))
    down,strange,bottom,p,u,t,el,mu,tau,x=variables
    d2,s2=(down/bottom).sq(),(strange/bottom).sq()
    cp,cu,ct=[(1-a.sq()).sqrt() for a in (p,u,t)]
    base=d2*cp.sq()+s2*p.sq()
    m=(s2-d2)*cp*p
    a2=cu.sq()*base+u.sq()
    z0=u*(1-base)
    c2=cu.sq()*((ct*m).sq()+(t*z0).sq()+2*ct*m*t*z0*x)/a2
    e2=cu.sq()*((t*m).sq()+(ct*z0).sq()-2*t*m*ct*z0*x)/a2
    ud=(p*ct).sq()+(cp*t*u).sq()+2*p*ct*cp*t*u*x
    us=(cp*ct).sq()+(p*t*u).sq()-2*cp*ct*p*t*u*x
    ub=(t*cu).sq()
    fminor=s2*ud+d2*us+d2*s2*ub
    f2,b2=fminor/a2,d2*s2/fminor
    h22=d2*ud+s2*us+ub
    h33=1+d2+s2-a2-h22
    z2=h22-b2-c2
    K,L=a2+b2+9*c2,z2+h33
    h=(down/bottom)*(strange/bottom)
    tr=el.sq()+mu.sq()+tau.sq()
    D=el*mu*tau/(tr*tr.sqrt())
    P=(el.sq()*mu.sq()+el.sq()*tau.sq()+mu.sq()*tau.sq())/tr.sq()
    w=AD["roots"]((D*K*L.sqrt()/h).sq(),rising)
    v=w/L
    rr=D*w.sqrt()/(h*L.sqrt())
    unit=a2*b2*((1-w)/K).sq()
    mid=a2*(z2+f2)+b2*h33
    if k%12 in (0,6):
        ck,sk=(1 if k==0 else -1),0
    else:
        ck=.5 if k%12 in (2,10) else -.5
        # Even k have sin(k phi)= ±sqrt(3)/2 exactly.
        sk=(x*0+3).sqrt()*(.5 if k%12 in (8,10) else -.5)
    Udown=d2+s2+d2*s2
    twoReal=c2*f2+e2*z2+a2*b2+mid-Udown
    imag=(orientation*cp*p*ct*t*cu.sq()*u*(1-x.sq()).sqrt()
          *(s2-d2)*(1-s2)*(1-d2)/a2)
    return unit+(mid+9*c2*f2)*rr+e2*z2*v.sq()-3*v*rr.sqrt()*(ck*twoReal-2*sk*imag)-P


def certify(dom,rising,k,orient,cap=19):
    def evaluate(box):
        return residual_jet(box,rising,k,orient)
    todo=dom[None].copy()
    center=todo.mean(axis=2)
    pos=evaluate(np.repeat(center[:,:,None],2,axis=2)).v[0].item()>0
    leaves=0
    lower=1.
    for depth in range(cap):
        if not len(todo):return (leaves,depth,lower)
        batch=[]
        for cell in np.array_split(todo,max(1,math.ceil(len(todo)/128))):
            val=evaluate(cell)
            for _ in range(3):
                active=cell[:,:,0]!=cell[:,:,1]
                inc=(val.d[0]>=0)&active
                dec=(val.d[1]<=0)&active&~inc
                if not (inc|dec).any():break
                lo=inc if pos else dec
                hi=dec if pos else inc
                cell[:,:,1]=np.where(lo,cell[:,:,0],cell[:,:,1])
                cell[:,:,0]=np.where(hi,cell[:,:,1],cell[:,:,0])
                val=evaluate(cell)
            mid=cell.mean(axis=2)
            vm=evaluate(np.repeat(mid[:,:,None],2,axis=2)).v
            deriv=AD["prod"](val.d,(AD["down"](cell[:,:,0]-mid),AD["up"](cell[:,:,1]-mid)))
            for j in range(10):
                vm=AD["add"](vm,(deriv[0][:,j:j+1],deriv[1][:,j:j+1]))
            gap=(np.maximum(vm[0],val.v[0]) if pos else -np.minimum(vm[1],val.v[1])).ravel()
            passed=gap>0
            leaves+=int(sum(passed))
            if passed.any():lower=min(lower,float(gap[passed].min()))
            if passed.all():continue
            stale=cell[~passed]
            axes=np.argmax((deriv[1]-deriv[0])[~passed],axis=1)
            idx=np.arange(len(stale))
            left,right=stale.copy(),stale.copy()
            left[idx,axes,1]=right[idx,axes,0]=mid[~passed][idx,axes]
            batch.extend((left,right))
        todo=np.concatenate(batch) if batch else []
        if len(todo)>5000:break
    return ("open",leaves,len(todo))


def check_jets(dom,orient):
    old=OLD["OLD"]
    count=0
    for frac in (.35,.77):
        xx=dom[:,0]+frac*(dom[:,1]-dom[:,0])
        q=np.r_[xx[:3],np.arcsin(xx[3:6])]
        delta=orient*np.arccos(xx[-1])
        phase=orient*old["gamma_deg"](xx[-1],q)
        ed,_,pars=OLD["parameters"](phase,q,xx[6:9])
        np.testing.assert_allclose(edges(q[None,:],np.array([phase]))[:,0],ed[:5]+[complex(*ed[-2:])],rtol=1e-8)
        tr=sum(xx[6:9]**2)
        target=sum((xx[i]*xx[j])**2 for i,j in ((6,7),(7,8),(6,8)))/tr**2
        for k in CFG["neutral_powers"]:
            for b,rising in enumerate((True,False)):
                y=charged(np.array([*ed[:5],complex(*ed[-2:])])[:,None],k,np.array([pars[b][0]]))[0]
                vals=np.linalg.svd(y,compute_uv=False)
                shape=sum((vals[i]*vals[j])**2 for i,j in ((0,1),(1,2),(0,2)))/sum(vals**2)**2
                z=residual_jet(np.repeat(xx[None,:,None],2,axis=2),rising,k,orient)
                assert z.v[0].item()-2e-13<=shape-target<=z.v[1].item()+2e-13
                count+=1
    return count


def witness(k,phase,pulls,port):
    q=Q0+pulls[:6]*QW
    leptons=L0+pulls[6:]*LW
    ed,delta,pars=OLD["parameters"](phase,q,leptons)
    mag,scale,_=pars[1]
    y=scale*OLD["OLD"]["LOW"]["mtx"](ed,[1,1,3],mag*np.exp(1j*k*PHI)).T
    eig,U=np.linalg.eigh(y@y.conj().T)
    np.testing.assert_allclose(np.sqrt(eig),leptons,rtol=1e-9)
    t12,t13,t23=q[3:]
    vud=np.cos(t12)*np.cos(t13)
    vcd=-np.sin(t12)*np.cos(t23)-np.cos(t12)*np.sin(t23)*np.sin(t13)*np.exp(1j*delta)
    vcb=np.sin(t23)*np.cos(t13)
    vub=np.sin(t13)*np.exp(-1j*delta)
    apex=-vud*vub.conjugate()/(vcd*vcb)
    ck=runpy.run_path(str(HERE.parent/"1675/independent_check.py"))
    phi_mz=ck["low_energy"](q,leptons,delta,192)
    resp={"quark_six_pulls_and_lepton_three_pulls":pulls.tolist(),"high_gamma_deg":phase,
          "one_loop_gamma_mod180_at_MZ_deg":phi_mz,"sin2beta_matching":np.sin(-2*np.angle(1-apex)),
          "abs_t":mag,"scale":scale,"physical_charged_yukawas":np.sqrt(eig).tolist()}
    R=5.83
    ds=7.388e-5
    resp["R"]=R
    resp["d21"]=ds
    resp["d31"]=R**2*ds
    for name,(ph,sg) in PORTS.items():
        v=NU["neutrino_mix"](U,ph,sg,R)
        angles=NU["OLD"]["angles"](v)
        resp[name]=angles.tolist()
        if name==port:
            assert np.all((angles>WINDOWS[:,0]) & (angles<WINDOWS[:,1]))
    return resp


def main():
    proofs={}
    checks=0
    for radius,margin in [(1,True),(2,False)]:
        section={}
        for sg in (1,-1):
            dom=box_domain(radius,sg,margin)
            checks+=check_jets(dom,sg)
            for k in CFG["neutral_powers"]:
                for rising in (True,False):
                    label=f"{sg:+d}_k{k}_{'small' if rising else 'large'}"
                    if not rising and ((sg,k)==(1,2) or (radius==2 and (sg,k)==(-1,8))):
                        section[label]="charged solutions verified by explicit witnesses"
                        continue
                    ans=certify(dom,rising,k,sg)
                    assert isinstance(ans[0],int),label
                    section[label]=list(ans)
        proofs[str(radius)]=section
    signs=np.array([1,-1,-1,-1,1,-1,-1,1,1])
    shifts=1.90*signs
    q=Q0+shifts[:6]*QW
    le=L0+shifts[6:]*LW
    g=OLD["OLD"]["zero"](lambda z:OLD["residual"](z,0,1,q,le),58.,80.)
    real=witness(0,g,shifts,"second_2_minus1")
    assert 63.2-1.9*2.2<real["one_loop_gamma_mod180_at_MZ_deg"]<63.2+1.9*2.3
    sig8=np.array([-1,-1,1,-1,1,1,1,-1,1])
    def residual(v):
        return OLD["residual"](-116.8,8,1,Q0+v*sig8[:6]*QW,L0+v*sig8[6:]*LW)
    r8=OLD["OLD"]["zero"](residual,0,2)
    other=witness(8,-116.8,r8*sig8,None)
    assert other["sin2beta_matching"]<0<real["sin2beta_matching"]
    return {"direct_matrix_cross_checks":checks,"charged_exclusion_certificates":proofs,
            "positive_real_transfer_full_witness":real,"opposite_k8_charged_witness":other}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=OLD["OLD"]["tidy"](main())
    path=HERE/"results.json"
    if args.write:path.write_text(json.dumps(results,indent=2)+"\n",encoding="utf8")
    else:assert json.loads(path.read_bytes())==results
    print("PASS: full charged covers, direct SVDs, common witnesses and one-loop transfer")
