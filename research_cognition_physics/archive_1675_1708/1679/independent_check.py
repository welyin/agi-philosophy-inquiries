"""Even charged phase, direct charged SVD and frozen neutrino endpoints."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ARCH=HERE.parent
OLD=runpy.run_path(str(ARCH/'1678/independent_check.py'))
PIV=runpy.run_path(str(ARCH/'1661/independent_check.py'))
Q=np.array(OLD['read']('../1397/inputs.json')['high_scale']['down']+
           OLD['read']('../1397/inputs.json')['high_scale']['angles_theta12_theta13_theta23_delta'][:3])
L=np.array(OLD['read']('../1402/inputs.json')['scales']['1e16'][3:])
ANGLES=np.array(OLD['read']('../1659/results.json')['joint_windows'])
PHI=-np.pi/6
BOX=runpy.run_path(str(ARCH/'1675/phase_box_check.py'))
IV=runpy.run_path(str(ARCH/'1678/check.py'))


def angles_frame(x,leg=0,power=2,gamma=63.2):
    raw=Q[:,0]+x[:,:6]*Q[:,1]
    d=raw[:,:3]/raw[:,2,None]
    a,b,c=raw[:,3:].T
    sa,sb,sc=np.sin([a,b,c])
    ca,cb,cc=np.cos([a,b,c])
    angle=np.radians(gamma)
    delta=angle+np.arcsin(ca*np.tan(c)*sb/sa*np.sin(angle))
    phase=np.exp(1j*delta)
    v=np.stack([ca*cb,sa*cb,sb*phase.conj(),-sa*cc-ca*sc*sb*phase,ca*cc-sa*sc*sb*phase,sc*cb,
                sa*sc-ca*cc*sb*phase,-ca*sc-sa*cc*sb*phase,cc*cb],axis=1).reshape(-1,3,3)
    herm=(v*(d*d)[:,None,:])@v.conj().transpose(0,2,1)
    perm=[0,2,1]
    ch=np.linalg.cholesky(herm[:,perm,:][:,:,perm])
    m=ch[:,perm,:][:,:,[2,0,1]]
    m[:,1,:]*=np.exp(-1j*np.angle(m[:,1,1]))[:,None]
    m[:,2,:]*=np.exp(-1j*np.angle(m[:,2,1]))[:,None]
    m[:,:,0]*=-np.exp(-1j*np.angle(m[:,1,0]))[:,None]
    m[:,:,2]*=np.exp(-1j*np.angle(m[:,2,2]))[:,None]
    t=(1.17+.10*x[:,6])*np.exp(1j*power*PHI)
    m[:,1,1]*=3
    m[:,2,:]*=t[:,None]
    m[:,:2,2]*=t[:,None]
    basis,yy,_=np.linalg.svd(m.transpose(0,2,1))
    e=basis[:,:,::-1]
    yy=yy[:,::-1]
    # Minimax overall scale for the three lepton marginal pulls.
    ia,ib=np.array([0,0,1]),np.array([1,2,2])
    candidates=((L[ia,0]*yy[:,ib]-L[ib,0]*yy[:,ia])/
                (L[ia,1]*yy[:,ib]+L[ib,1]*yy[:,ia]))
    j=abs(candidates).argmax(axis=1)
    scales=(L[ia[j],0]*L[ib[j],1]+L[ib[j],0]*L[ia[j],1])/(
        yy[np.arange(len(x)),ia[j]]*L[ib[j],1]+yy[np.arange(len(x)),ib[j]]*L[ia[j],1])
    pulls=(scales[:,None]*yy-L[:,0])/L[:,1]
    R=5.83+.12*x[:,7]
    rr=PIV['roots'](R)[0]
    aa=np.outer(PIV['A'],PIV['A'])
    bb=np.outer(PIV['B'],PIV['B'])
    cfg=[(2*PHI,3*PHI+np.angle(3+2*np.exp(4j*PHI)),1),
         (10*PHI+np.angle(2-np.exp(4j*PHI)),7*PHI,-1),
         (10*PHI+np.radians(11.5),7*PHI,-1)][leg]
    p=np.exp(1j*np.array([0,*cfg[:2]]))
    n=(aa+rr[:,None,None]*np.exp(2j*np.pi*cfg[2]/3)*bb)*p[:,None]*p
    _,mode=np.linalg.eigh(n.conj().transpose(0,2,1)@n)
    u=e.conj().transpose(0,2,1)@mode
    pm=abs(u)**2
    y=pm[:,0,2]
    obs=np.stack([4*y*(1-y),pm[:,0,1]/(1-y),pm[:,1,2]/(1-y)],axis=1)
    return dict(obs=obs,raw=raw,pulls=pulls,scale=scales,R=R,pmns=u,weak=v,t=t)


def cost(x,leg=0,power=2,gamma=63.2):
    m=angles_frame(x,leg,power,gamma)
    cent=ANGLES.mean(axis=1)
    half=ANGLES[:,1]-cent
    return np.maximum(np.abs(m['pulls']).max(axis=1)/2,(abs(m['obs']-cent)/half).max(axis=1))


def search(leg=0,power=2,gamma=63.2):
    rg=np.random.default_rng(34267)
    lo=np.r_[[-2.]*6,-1.8,-1.]
    hi=np.r_[[2.]*6,1.8,1.]
    n=240
    x=rg.uniform(lo,hi,(n,8))
    v=cost(x,leg,power,gamma)
    for i in range(260):
        a,b,c=x[rg.integers(0,n,size=(3,n))]
        cand=np.clip(a+rg.uniform(.45,.9)*(b-c),lo,hi)
        y=np.where(rg.random((n,8))<.9,cand,x)
        vy=cost(y,leg,power,gamma)
        acc=vy<v
        x[acc],v[acc]=y[acc],vy[acc]
    idx=np.argmin(v)
    print(leg,power,gamma,v[idx], x[idx].tolist(),flush=True)
    return x[idx]


def inverse_bound(cells,small,sign,leg=0,cut=None,above=False,power=2,orientation=1):
    # For the small s branch test the charged shape directly. For large s test
    # the electron weight on the exact null ray, as a positive sum over adjugate columns.
    J=BOX['Jet']
    coords=[]
    for i in range(cells.shape[1]):
        der=np.zeros((len(cells),cells.shape[1]))
        der[:,i]=1
        coords.append(J((cells[:,i,0,None],cells[:,i,1,None]),(der,der)))
    d,sm,bottom,p,u,t,ele,mu,tau,cd=coords[:10]
    d,sm=d/bottom,sm/bottom
    aa,ss=d.sq(),sm.sq()
    cp,cu,ct,sd=[(1-y.sq()).sqrt() for y in (p,u,t,cd)]
    base=aa*cp.sq()+ss*p.sq()
    mix=(ss-aa)*p*cp
    a2=base*cu.sq()+u.sq()
    h12=(cu*(ct*mix+t*u*(1-base)*cd),-cu*t*u*(1-base)*sd)
    h13=(cu*(-t*mix+ct*u*(1-base)*cd),-cu*ct*u*(1-base)*sd)
    def mag(z):
        return z[0].sq()+z[1].sq()
    def product(z,w):
        return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
    def dot(z,w):
        return z[0]*w[0]+z[1]*w[1]
    def grouped(z,w):
        return z[0]*w[0]+z[1]*w[1],z[1]*w[0]-z[0]*w[1]
    c2,e2=mag(h12)/a2,mag(h13)/a2
    v2=[(-p*ct-cp*t*u*cd,-cp*t*u*sd),(cp*ct-p*t*u*cd,-p*t*u*sd),(t*cu,t*0)]
    v3=[(p*t-cp*ct*u*cd,-cp*ct*u*sd),(-cp*t-p*ct*u*cd,-p*ct*u*sd),(ct*cu,t*0)]
    weights=(ss,aa,aa*ss)
    minor=sum(w*mag(z) for w,z in zip(weights,v2))
    b2,f2=aa*ss/minor,minor/a2
    ap=[grouped(v,k) for v,k in zip(v2,v3)]
    adj=[sum(w*z[j] for w,z in zip(weights,ap)) for j in range(2)]
    zz=mag(adj)/(a2*minor)
    raw=grouped(h13,h12)
    S=-6*dot(raw,adj)/a2.sq()
    K=a2+b2+9*c2
    h33=aa*mag(v3[0])+ss*mag(v3[1])+mag(v3[2])
    L=zz+h33
    h=d*sm
    tr=ele.sq()+mu.sq()+tau.sq()
    D=ele*mu*tau/(tr*tr.sqrt())
    P=(ele.sq()*mu.sq()+ele.sq()*tau.sq()+mu.sq()*tau.sq())/tr.sq()
    ww=BOX['roots']((D*K*L.sqrt()/h).sq(),small)
    cf=(c2*f2).sqrt()
    jh=cp*p*ct*t*cu.sq()*u*sd*(ss-aa)*(1-ss)*(1-aa)
    er,ei=S/(6*cf),orientation*jh/(a2*cf)
    # k=2 hence cos=1/2, sin=-sqrt(3)/2. Enclose the irrational phase.
    phase=(J((.5,.5),(0.,0.)),J((BOX['down'](-math.sqrt(3)/2),BOX['up'](-math.sqrt(3)/2)),(0.,0.)))
    if power==8:
        phase=(-phase[0],-phase[1])
    else:
        assert power==2
    trr,tii=product(phase,(er,ei))
    if small:
        r=D*ww.sqrt()/(h*L.sqrt())
        v=ww/L
        sh=a2*b2/K.sq()*(1-ww).sq()+(a2*(zz+f2)+b2*h33)*r+(tii*v).sq()
        return sh+((9*c2*f2*r).sqrt()-trr*v).sq()-P
    s=(ww*K/((1-ww)*L)).sqrt()
    a,b,c,e,f=a2.sqrt(),-b2.sqrt(),c2.sqrt(),e2.sqrt(),f2.sqrt()
    # Neutral third leg uses e^(i alpha)=e^(-i pi/3), e^(i beta)=(-sqrt(3)-2i)/sqrt(7).
    # Both components below follow exp(3i phi)(3+2 exp(4i phi))/sqrt(7).
    def cn(x):
        return J((x-2.e-15,x+2.e-15),(0.,0.))
    alpha=(cn(.5),cn(-math.sqrt(3)/2))
    beta=(cn(-math.sqrt(3/7)),cn(-2/math.sqrt(7)))
    if leg:
        beta=(cn(-math.sqrt(3)/2),cn(.5))
        if leg==1:
            alpha=(cn(1/(2*math.sqrt(7))),cn(3*math.sqrt(3)/(2*math.sqrt(7))))
        else:
            from fractions import Fraction
            iv=IV['OV']
            pi=iv['bounds']('3.14159265358979323846','3.14159265358979323847')
            ga=pi*iv['I'].point(Fraction(143,360))
            def certified(sine):
                low,high=iv['trig'](ga,sine).values()
                return J((BOX['down'](low),BOX['up'](high)),(0.,0.))
            alpha=(certified(False),certified(True))
    cc=e/f
    v=(alpha[0]+cc*beta[0],alpha[1]+cc*beta[1])
    j0=(3*c-s*trr/f,s*tii/f)
    # 4||h1||² - 6|<n,h1>|², h1=(3C-conj(t Z)*E/F,-B,B*E/F).
    first=b.sq()*(4*(1+cc.sq())-mag(v))-4*b*dot(j0,v)
    # The second adjugate column has exactly pure first weak flavour. For
    # the third column, 4||h3||²-6|<n,h3>|² >= -2||h3||² by Cauchy.
    third=a2*(zz/f2+b2/(s.sq()*f2))
    if cut is not None:
        # Exact heavy Takagi ray of the 2-by-2 Gram matrix, in real
        # orthogonal directions (0,1,1) and (1,1,-1). No small cancellation
        # from H/m minus its light eigenvalue projector.
        rr=coords[10]
        real=-2+22*rr
        im=cn((1 if leg else -1)*2*math.sqrt(3))
        gap=4-16*rr+55*rr.sq()
        factor=4*rr/(gap+(gap.sq()+24*rr.sq()*(real.sq()+im.sq())).sqrt())
        ccx=(factor*real,factor*im)
        norm=2+3*mag(ccx)
        # Complex coefficients of the bra w^dagger times the two flavour phases.
        bx,by=ccx[0],-ccx[1]
        v1=(bx,by)
        v2=product((1+bx,by),alpha)
        v3=product((1-bx,-by),beta)
        a1=product(v1,j0)
        amp1=(a1[0]-b*v2[0]+b*cc*v3[0],a1[1]-b*v2[1]+b*cc*v3[1])
        # Third column after dividing the adjugate by conj(t)F and by A.
        first3=(er/(e*f),-ei/(e*f))
        last3=(-b*phase[0]/(s*f),-b*phase[1]/(s*f))
        amp3a,amp3b=product(v1,first3),product(v3,last3)
        amp3=(amp3a[0]+amp3b[0],amp3a[1]+amp3b[1])
        den=mag(j0)+b.sq()*(1+cc.sq())+a2*(1+mag(first3)+mag(last3))
        num=mag(amp1)+a2*(mag(v1)+mag(amp3))
        out=num-cut*norm*den
        return out if above else -out
    if sign<0:
        # Complementary lower enclosure 6|<n,h1>|² - 2||h1||².
        first=2*mag(j0)+4*b*dot(j0,v)+b.sq()*(mag(v)-2*(1+cc.sq()))
    return first-2*third


def box_proof(power=2):
    original=BOX['classify']
    # Bind the existing interval certifier to this round's residual without modifying the old module.
    from types import FunctionType
    direction=1 if power==2 else -1
    dom=BOX['domain']()[0 if direction==1 else 1]
    results={}
    for leg in ([0] if power==2 else [1,2]):
        for small,sg,label in ((True,1,'charged_small'),(False,1,'null_upper_large'),(False,-1,'null_lower_large')):
            def eval_box(cells,s,ss):
                return inverse_bound(cells,s,ss,leg=leg,power=power,orientation=direction)
            classify=FunctionType(original.__code__,{**original.__globals__,'residual':eval_box})
            cert=classify(dom,small,sg)
            mid=dom.mean(axis=1)
            val=eval_box(np.repeat(mid[None,:,None],2,axis=2),small,sg).v
            assert (val[1][0,0]<0 if small else val[0][0,0]>0)
            results[f'leg{leg}_{label}']=cert
    return results


def reactor_proof(width=1,power=2):
    from types import FunctionType
    iv=IV['OV']
    def margin(row):
        return iv['I'].point(row[0])+iv['bounds'](-width,width)*iv['I'].point(row[1])
    quarks=[margin(row) for row in Q]
    ls=[margin(row) for row in L]
    vals=quarks[:3]+[iv['trig'](x,True) for x in quarks[3:]]+ls
    direction=1 if power==2 else -1
    sub=np.array([[BOX['down'](w.values()[0]),BOX['up'](w.values()[1])] for w in vals]+
                 [([.397,.501] if direction==1 else [-.502,-.398])])
    original=BOX['classify']
    ratio=np.array([5.7092,5.9532])
    assert np.all(ratio*np.array([1,-1]) < np.array(IV['NU']['gap_box'](2))*[1,-1])
    # r_small decreases, r_large increases with R. Extra 1e-7 protects the
    # printed endpoints in addition to roundoff in the analytic quadratic.
    roots=PIV['roots'](ratio)
    result={}
    for leg in ([1,2] if power==2 else [0]):
        for branch,rr in enumerate(roots):
            dom=np.vstack([sub,[min(rr)-1.e-7,max(rr)+1.e-7]])
            # Compare the proxy to conservative cutoffs. Error < .000025 in
            # going to the true electron ray. Large roots must miss BOTH
            # reactor branches.
            for cut,above in ([(.02045,False)] if branch==0 else [(.0231,True),(.975,False)]):
                def target(cells,small,sign):
                    return inverse_bound(cells,False,sign,leg,cut,above,power,direction)
                classify=FunctionType(original.__code__,{**original.__globals__,'residual':target,'N':11})
                mid=dom.mean(axis=1)
                assert target(np.repeat(mid[None,:,None],2,axis=2),False,1).v[0][0,0]>0
                result[f'leg{leg}_root{branch}_{cut}']=classify(dom,False,1)
    return result


def direct_tests():
    # Independent SVD check of interval amplitudes at both roots, all ports,
    # both CP directions and two charged powers. Only trace/determinant are
    # imposed on these test points; the inverse-density proxy is exact.
    ntest=0
    for power,direction,ga in ((2,1,63.2),(8,-1,-116.8)):
        raw=Q[:,0]
        sh,long,_,delta=OLD['get_matrices'](raw,ga)
        vals=IV['parameters'](ga,raw,L[:,0])[2]
        mag,scale,_=vals[1]
        y=(sh+mag*np.exp(1j*power*PHI)*long).T
        inv=np.linalg.inv(y@y.conj().T)
        inv/=np.trace(inv)
        coord=np.r_[raw[:3],np.sin(raw[3:]),L[:,0],np.cos(delta)]
        for leg in range(3):
            phases=[(2*PHI,3*PHI+np.angle(3+2*np.exp(4j*PHI)),1),
                    (10*PHI+np.angle(2-np.exp(4j*PHI)),7*PHI,-1),
                    (10*PHI+np.radians(11.5),7*PHI,-1)][leg]
            p=np.exp(1j*np.array([0,*phases[:2]]))
            for rr in PIV['roots'](5.83):
                n=p[:,None]*(np.outer(PIV['A'],PIV['A'])+rr*np.exp(phases[2]*2j*np.pi/3)*
                              np.outer(PIV['B'],PIV['B']))*p
                w=np.linalg.svd(n)[2].conj().T[:,0]
                prob=float(np.vdot(w,inv@w).real)
                x=np.r_[coord,rr]
                cell=np.repeat(x[None,:,None],2,axis=2)
                for offset in (-1.e-6,1.e-6):
                    z=inverse_bound(cell,False,1,leg,prob+offset,False,power,direction)
                    assert (z.v[0][0,0]>0 if offset>0 else z.v[1][0,0]<0)
                    ntest+=1
    return ntest


def compute():
    lep=(L[0,0]+2*L[0,1])/(L[1:,0]-2*L[1:,1])
    err=float(lep@lep)
    assert err<.000025
    low,high=(1-np.sqrt(1-ANGLES[0]))/2
    assert .02045+err<low and .0231-err>high and .975+err<1-high
    assert 2/3+err<(1-high)*(1-ANGLES[1,1])
    assert 1/3-err>high*(1-ANGLES[1,0])
    return dict(proxy_operator_error_upper=err, two_reactor_branches=[[low,high],[1-high,1-low]],
                direct_amplitude_checks=direct_tests(), third_port_k2_positive_two_widths=box_proof(2),
                second_ports_k8_negative_two_widths=box_proof(8), second_ports_k2_positive_one_width=reactor_proof(1,2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--leg',type=int,default=0)
    parser.add_argument('--gamma',type=float,default=63.2)
    parser.add_argument('--power',type=int,default=2)
    parser.add_argument('--certificate',action='store_true')
    parser.add_argument('--reactor',type=int)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.check or args.write:
        out=compute()
        target=HERE/'independent_results.json'
        if args.check:
            ref=json.loads(target.read_bytes())
            def compare(x,y):
                if isinstance(x,dict):
                    for k in x: compare(x[k],y[k])
                else:
                    np.testing.assert_allclose(x,y,rtol=1.e-10,atol=0)
            compare(out, ref)
            print("PASS: inverse-density proxy, six port/root cases, two CP orientations")
        else:
            with target.open('x',encoding='utf8') as f:
                f.write(json.dumps(out,indent=2)+'\n')
    elif args.certificate:
        print(box_proof(args.power))
    elif args.reactor:
        print(reactor_proof(args.reactor,args.power))
    else:
        search(args.leg,args.power,args.gamma)
