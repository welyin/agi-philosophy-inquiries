"""R1692. Enclosures of the charged electron and of both heavy Takagi rays."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OLD = runpy.run_path(str(HERE.parent/"1675"/"phase_box_check.py"))
PORT = runpy.run_path(str(HERE.parent/"1682"/"independent_check.py"))
CFG = json.loads((HERE/"inputs.json").read_bytes())
J, down, up = OLD["Jet"], OLD["down"], OLD["up"]


def constant(value):
    return J((down(value),up(value)), (0.,0.))


def box(x, y):
    return J((down(x),up(y)), (0.,0.))


def conjprod(a,b):
    return a[0]*b[0]+a[1]*b[1], a[1]*b[0]-a[0]*b[1]


def mul(a,b):
    return conjprod(a,(b[0],-b[1]))


def norm(a):
    return a[0].sq()+a[1].sq()


def quotient(a,b):
    den = norm(b)
    return [z/den for z in conjprod(a,b)]


def charged_geometry(cells):
    total = cells.shape[1]
    v = []
    for j in range(total):
        der = np.zeros((len(cells),total))
        der[:,j] = 1.
        v.append(J((cells[:,j,0,None],cells[:,j,1,None]),(der,der)))
    d,s,p,u,t,co,mag = v
    a,b = d.sq(),s.sq()
    cp,cu,ct,si = [(1-e.sq()).sqrt() for e in (p,u,t,co)]
    base = a*cp.sq()+b*p.sq()
    mix = (b-a)*cp*p
    a2 = base*cu.sq()+u.sq()
    h12 = (cu*(ct*mix+t*u*(1-base)*co),-cu*t*u*(1-base)*si)
    h13 = (cu*(-t*mix+ct*u*(1-base)*co),-cu*ct*u*(1-base)*si)
    v2 = [(-p*ct-cp*t*u*co,-cp*t*u*si),(cp*ct-p*t*u*co,-p*t*u*si),(t*cu,t*0)]
    v3 = [(p*t-cp*ct*u*co,-cp*ct*u*si),(-cp*t-p*ct*u*co,-p*ct*u*si),(ct*cu,t*0)]
    weights = [b,a,a*b]
    minor = sum(ww*norm(y) for ww,y in zip(weights,v2))
    babs, c, e, f = [(z).sqrt() for z in (a*b/minor,norm(h12)/a2,norm(h13)/a2,minor/a2)]
    gg = [conjprod(l,r) for l,r in zip(v2,v3)]
    adj = [sum(ww*v[j] for ww,v in zip(weights,gg)) for j in range(2)]
    cross = conjprod(h13,h12)
    er = -(cross[0]*adj[0]+cross[1]*adj[1])/(a2.sq()*c*f)
    ei = cp*p*ct*t*cu.sq()*u*si*(b-a)*(1-b)*(1-a)/(a2*c*f)
    tr,ti = mul((constant(math.sqrt(3)/2),constant(-.5)),(er,ei))
    jr, ji = 3*c-mag*tr/f, mag*ti/f
    slope = ji/jr
    re = babs / (jr*(1+slope.sq()))
    return [re,-re*slope,e/f]


def domain():
    data = OLD["domain"]()[0]
    def divide(i,j):
        return [down(data[i,0]/data[j,1]),up(data[i,1]/data[j,0])]
    return np.array([divide(0,2),divide(1,2),*data[3:6],data[9],[1.09,CFG["upper_s"]]])


def sumc(a,b):
    return a[0]+b[0], a[1]+b[1]


def scale(a,t):
    return a[0]*t,a[1]*t


def sr_model(mag):
    q = constant(3.).sqrt().reciprocal()
    n = constant(7.).sqrt()/constant(5.)
    x = q*(1-mag*q)/(mag-q)
    return n*(1+x)/(n.sq()+x)


def heavy_readout(h,R,e2,k,large,schur=1.):
    a,b,g = 1+h.sq(), 10+h.sq(), 3+h.sq()
    G = 1+5*h.sq()
    shift=(schur if isinstance(schur,J) else constant(schur))-2
    length=(shift.sq()+12).sqrt()
    cos=shift/length
    sin=constant(12.).sqrt()/length
    T=R+R.reciprocal()
    coeff=T*G-2*g.sq()*cos
    radical=(((T-2)*G-2*g.sq()*(1+cos))*((T+2)*G+2*g.sq()*(1-cos))).sqrt()
    r = (coeff+radical)/(2*b.sq()) if large else 2*a.sq()/(coeff+radical)
    # Heavy eigenvector in span(a*,b*) is proportional to a*+z b*.
    z=quotient((g*(a*cos+r*b),-g*a*sin),(G*R-cos*g.sq()-r*b.sq(),sin*g.sq()))
    zc=(z[0],-z[1])
    alpha=(constant(.5),-constant(3.).sqrt()/2)
    beta=(-constant(3./7).sqrt(),-constant(7.).sqrt().reciprocal()*2)
    unit=(1,0)
    second=mul(sumc(unit,scale(zc,3)),alpha)
    third=mul(sumc(unit,zc),beta)
    amp=sumc(zc,mul(e2,sumc(second,scale(third,-k*h))))
    den=(a+2*g*z[0]+b*norm(z))*(1+norm(e2)*(1+k.sq()))
    return norm(amp)/den


def run_threshold(test, boxes):
    """Move unambiguous extrema to faces; bound the remainder by the mean-value theorem."""
    leaves = 0
    mingap = 10.
    for depth in range(28):
        if not len(boxes):
            return dict(certified_boxes=leaves,depth=depth,margin=mingap)
        pending=[]
        for cells in np.array_split(boxes,math.ceil(len(boxes)/256)):
            val=test(cells)
            for _ in range(3):
                fixed=cells[:,:,0]==cells[:,:,1]
                centre_low=(val.d[0]>=0)&~fixed
                centre_high=(val.d[1]<=0)&~fixed&~centre_low
                if not (centre_low|centre_high).any():
                    break
                cells[:,:,1]=np.where(centre_low,cells[:,:,0],cells[:,:,1])
                cells[:,:,0]=np.where(centre_high,cells[:,:,1],cells[:,:,0])
                val=test(cells)
            mid=cells.mean(axis=2)
            pc=test(np.repeat(mid[:,:,None],2,axis=2))
            delta=OLD["prod"](val.d,(down(cells[:,:,0]-mid),up(cells[:,:,1]-mid)))
            mv=pc.v
            for j in range(cells.shape[1]):
                mv=OLD["add"](mv,(delta[0][:,j,None],delta[1][:,j,None]))
            low=np.maximum(val.v[0],mv[0])[:,0]
            good=low>0.
            if good.any():
                leaves+=int(good.sum())
                mingap=min(mingap,float(low[good].min()))
            if not good.all():
                left=cells[~good].copy(); right=left.copy()
                k=np.argmax((delta[1]-delta[0])[~good],axis=1); j=np.arange(len(k))
                left[j,k,1]=right[j,k,0]=mid[~good][j,k]
                pending.extend([left,right])
        boxes=np.concatenate(pending) if pending else []
        if len(boxes)>32768:
            raise ValueError("increase resolution or revise bound")
    raise ValueError("threshold unproved")


def neutrino_jet(cells,large):
    axis=[]
    for i in range(cells.shape[1]):
        der=np.zeros((len(cells),cells.shape[1]));der[:,i]=1
        axis.append(J((cells[:,i,0,None],cells[:,i,1,None]),(der,der)))
    h,R,re,im,k=axis[:5]
    return heavy_readout(h,R,(re,im),k,large,axis[5] if len(axis)>5 else 1.)


def initial_domains(ranges,parts):
    span=[np.linspace(lo,hi,n+1) for (lo,hi),n in zip(ranges,parts)]
    ix=np.indices(parts).reshape(len(parts),-1).T
    return np.stack([np.stack((s[ix[:,j]],s[ix[:,j]+1]),axis=1)
                     for j,s in enumerate(span)],axis=1)


def direct_comparisons():
    rng=np.random.default_rng(1692)
    archive=PORT["BASE"]
    largest=0.
    for _ in range(32):
        raw=archive["Q"][:,0]+rng.uniform(-2,2,6)*archive["Q"][:,1]
        gamma=rng.uniform(*CFG["quark_arc_degrees"])
        sh,long,_,delta=PORT["OLD"]["get_matrices"](raw,gamma)
        mag=rng.uniform(1.091,1.27)
        y=(sh+mag*np.exp(-1j*np.pi/6)*long).T
        # The rank-two charged matrix supplies the proxy using its actual columns.
        proxy=np.linalg.svd(y[:,1:],full_matrices=True)[0][:,-1]
        proxy=proxy/proxy[0]
        d,s=raw[:2]/raw[2]
        x=np.array([d,s,*np.sin(raw[3:]),math.cos(delta),mag])
        values=[(j.v[0].item()+j.v[1].item())/2 for j in
                charged_geometry(np.repeat(x[None,:,None],2,axis=2))]
        k=values[2]
        np.testing.assert_allclose(proxy,[1, complex(*values[:2]),-k*complex(*values[:2])],atol=1e-11)
        R=rng.uniform(*CFG["R_outward"])
        z=rng.uniform(*CFG["excluded_schur_interval"])
        h=PORT["neutral_ratio"](mag)[1]
        legs=PORT["UNIT"].copy();legs[2]*=h
        eta=(z-2+2j*math.sqrt(3))/math.hypot(z-2,2*math.sqrt(3))
        aa,bb,g=1+h*h,10+h*h,3+h*h
        co=(1+5*h*h)*(R+1/R)-2*g*g*eta.real
        for large,r in enumerate(sorted(np.roots([bb**2,-co,aa**2]))):
            mat=legs[:,None]*(PORT["AA"]+r*eta*PORT["BB"])*legs
            sv,vec=np.linalg.eigh(mat.conj().T@mat)
            np.testing.assert_allclose(math.sqrt(sv[-1]/sv[-2]),R,rtol=5e-11)
            p=abs(np.vdot(proxy,vec[:,-1]))**2/np.vdot(proxy,proxy).real
            check=heavy_readout(constant(h),constant(R),tuple(map(constant,values[:2])),
                                constant(k),large,constant(z)).v
            largest=max(largest,abs(float((check[0]+check[1])/2)-p))
    assert largest<1e-10
    return dict(complex_charged_matrices=32,neutral_matrices=64,max_projector_error=round(largest,13))

def compute():
    dom = domain()
    d,s,p,u,t = [box(*v) for v in dom[:5]]
    e,m,tau = [box(*v) for v in OLD["domain"]()[0][6:9]]
    trace=e.sq()+m.sq()+tau.sq()
    determinant=e*m*tau/(trace*trace.sqrt())
    h33=(1-t.sq())*(1-u.sq())
    smax=(d*s/(determinant*h33*h33.sqrt())).sqrt()
    prior=json.loads((HERE.parent/"1682"/"results.json").read_bytes())
    assert prior["required_large_s_strict_lower"]>dom[6,0] and smax.v[1]<dom[6,1]
    hbounds=[1.075,1.29]
    assert sr_model(constant(dom[6,0])).v[0]>hbounds[0]
    assert sr_model(constant(smax.v[1])).v[1]<hbounds[1]
    Rbounds=CFG["R_outward"]
    Rlo,Rhi = PORT["BASE"]["IV"]["NU"]["gap_box"](2)
    assert Rbounds[0]<Rlo<Rhi<Rbounds[1]
    limits=[[.043,.063],[-.0075,-.001],[.74,.97]]
    charged_certs={}
    for k,(lo,hi) in enumerate(limits):
        for side,cut in enumerate([lo,hi]):
            def test(cells):
                v=charged_geometry(cells)[k]
                return constant(cut)-v if side else v-constant(cut)
            charged_certs[f"{k}_{side}"]=run_threshold(test,initial_domains(dom,[1,1,1,1,1,1,1]))
    neut_certs={}
    ranges=[hbounds,Rbounds,*limits,CFG["excluded_schur_interval"]]
    for large in (False,True):
        for side,cut in enumerate([.02365,.2]):
            def test(cells):
                v=neutrino_jet(cells,large)
                return constant(cut)-v if side else v-constant(cut)
            neut_certs[f"{int(large)}_{side}"]=run_threshold(test,initial_domains(ranges,[12,4,1,1,1,4]))
    A2=d.sq()+s.sq()*p.sq()+u.sq()
    mu2=(m/tau).sq()*constant(dom[6,0]).sq()*h33
    overlap=(box(*limits[0]).sq()+box(*limits[1]).sq()).sqrt()
    leakage=A2*overlap/mu2
    eps=constant(leakage.v[1])
    # e = cos(theta)*e0 + sin(theta)*e_perp. Enclose the entire projected ray.
    pmin=(constant(.02365).sqrt()*(1-eps.sq()).sqrt()-eps).sq()
    pmax=(constant(.2).sqrt()+eps).sq()
    prob=box(pmin.v[0],pmax.v[1])
    rea_low=4*constant(prob.v[0])*(1-constant(prob.v[0]))
    assert prob.v[1]<.5 and rea_low.v[0]>.0899
    return dict(s_strict_upper=float(smax.v[1]),s_domain=dom[6].tolist(),h_domain=hbounds,
                R_domain=Rbounds,charged_proxy_geometry=limits,charged_certificates=charged_certs,
                excluded_schur_interval=CFG["excluded_schur_interval"],direct_checks=direct_comparisons(),
                two_neutral_root_certificates=neut_certs,electron_ray_sine_error=float(eps.v[1]),
                reactor_probability_enclosure=[float(prob.v[0]),float(prob.v[1])],
                reactor_amplitude_strict_lower=float(rea_low.v[0]),
                registered_reactor_amplitude_upper=float(PORT["BASE"]["ANGLES"][0,1]))


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    data=compute()
    p=HERE/"results.json"
    if args.check:
        assert data==json.loads(p.read_bytes())
        print("R1692 OK")
    else:
        p.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
