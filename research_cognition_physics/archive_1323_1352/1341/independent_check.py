"""Inverse adjacent-ray test using the ordinary support-plane geometry and exact historical boxes."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np


def projector(x):
    return np.outer(x,x.conj())


def pdg(us,cb,ub,phase):
    c13=np.sqrt(1-ub**2)
    s12,s23=us/c13,cb/c13
    c12,c23=np.sqrt(1-s12*s12),np.sqrt(1-s23*s23)
    e=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,ub/e],
                     [-s12*c23-c12*s23*ub*e,c12*c23-s12*s23*ub*e,s23*c13],
                     [s12*s23-c12*c23*ub*e,-c12*s23-s12*c23*ub*e,c23*c13]])


def from_planes(a,b):
    spec_u,u=np.linalg.eigh(a)
    spec_d,d=np.linalg.eigh(b)
    normals=np.column_stack([u[:,0],d[:,0]])
    if np.linalg.svd(normals,compute_uv=False)[-1]<1e-9:
        raise ValueError('Common plane: the middle ray must be varied.')
    middle=np.linalg.qr(normals,mode='complete')[0][:,2]
    def resolve(y,masses,normal):
        other=np.linalg.qr(np.column_stack([normal,middle]),mode='complete')[0][:,2]
        q=np.column_stack([middle,other])
        t=q.conj().T@(y-masses[0]*np.eye(3))@q
        off=t[0,1]; low=t[1,1].real
        common=t[0,0].real-abs(off)**2/low
        individual=low+abs(off)**2/low
        ray=q@(np.array([off,low])/np.hypot(abs(off),low))
        h=abs(np.vdot(middle,ray))
        recon=masses[0]*np.eye(3)+common*projector(middle)+individual*projector(ray)
        return h,np.linalg.norm(recon-y),common,individual
    hu,eu,cu,pu=resolve(a,spec_u,u[:,0]);hd,ed,cd,pd=resolve(b,spec_d,d[:,0])
    assert min(cu,pu,cd,pd)>0
    return hu,hd,max(eu,ed)


def exact_boxes():
    us0,us1=F('.2236'),F('.2278')
    ub0,ub1=F('.00401'),F('.00461')
    cb0,cb1=F('.041'),F('.0422')
    c13lo=F('.99998');c23lo=F('.999');c12hi=F('.976')
    s12lo,s12hi=us0,us1/c13lo; s23lo,s23hi=cb0,cb1/c13lo
    assert 1-ub1**2>=c13lo**2 and 1-s23hi**2>=c23lo**2 and 1-s12lo**2<=c12hi**2
    cdlo=s12lo*c23lo-c12hi*s23hi*ub1
    cdhi=s12hi+c12hi*s23hi*ub1
    tdlo=s12lo*s23lo-c12hi*ub1
    tdhi=s12hi*s23hi+c12hi*ub1
    assert min(cdlo,tdlo)>0
    klo,khi=(tdlo/cdhi)**2,(tdhi/cdlo)**2
    ud=F('.00085')
    rho_u=(F('.703')-ud)/(F('168.7')-ud)
    rho_d=(F('.040')-F('.00414'))/(F('2.98')-F('.00414'))
    kdlo,kdhi=(ub0/us1)**2,(ub1/us0)**2
    def ulo(r):return (1-r)**2/((1+khi)*(1+r*r/klo))
    dlo=ulo(rho_u)
    dhi=(1-rho_d)**2/((1+kdlo)*(1+rho_d**2/kdhi))
    threshold=F('.0118')
    implied=(1-threshold)*ud+threshold*F('168.7')
    assert dlo>F('.9739')**2 and dhi<F('.85282')**2
    assert ulo(threshold)>dhi and implied>F('1.99149')
    return dict(k_up_lower=float(klo),k_up_upper=float(khi),k_down_lower=float(kdlo),k_down_upper=float(kdhi),
                up_overlap_lower=np.sqrt(float(dlo)),down_overlap_upper=np.sqrt(float(dhi)),
                assumed_rho_upper=float(threshold),hypothetical_up_lower=np.sqrt(float(ulo(threshold))),
                heldout_charm_strict_lower_GeV=float(implied),exact_charm_lower_fraction=str(implied),
                phase_all=True)


def compute():
    rng=np.random.default_rng(134122)
    emax,hmax=0.,0.
    for _ in range(32):
        u=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
        ray=rng.normal(size=3)+1j*rng.normal(size=3)
        ray=ray/np.linalg.norm(ray)
        v=[ray,u@ray,u@u@ray]
        b1,b2,w1,w2,w3,w4=rng.uniform(.2,3,6)
        A=b1*np.eye(3)+w1*projector(v[0])+w2*projector(v[1])
        B=b2*np.eye(3)+w3*projector(v[1])+w4*projector(v[2])
        hu,hd,e=from_planes(A,B)
        h=abs(np.vdot(v[0],v[1]))
        emax=max(emax,e);hmax=max(hmax,abs(hu-h),abs(hd-h))
    assert max(emax,hmax)<1e-8
    a=np.diag([.00127,.619,171.7])
    d=np.diag([.0029,.055,2.89])
    def phase_residual(delta):
        v=pdg(.2257,.0416,.00431,delta)
        beta=np.angle(-v[1,0]*v[1,2].conjugate()/(v[2,0]*v[2,2].conjugate()))
        return np.sin(2*beta)-.687
    central=[]
    for lo,hi in [(0.,1.),(1.,3.)]:
        assert phase_residual(lo)*phase_residual(hi)<0
        for _ in range(64):
            mid=(lo+hi)/2
            if phase_residual(lo)*phase_residual(mid)<=0:hi=mid
            else:lo=mid
        delta=(lo+hi)/2
        v=pdg(.2257,.0416,.00431,delta)
        hu,hd,e=from_planes(a,v@d@v.conj().T)
        k=abs(v[2,0]/v[1,0])**2
        def predict_h(r):
            return (1-r)/np.sqrt((1+k)*(1+r*r/k))
        left,right=0.,1.
        for _ in range(65):
            mid=(left+right)/2
            if predict_h(mid)>hd:left=mid
            else:right=mid
        m=.00127+(171.7-.00127)*(left+right)/2
        central.append(dict(delta=delta,up_overlap=hu,down_overlap=hd,predicted_charm_GeV=m))
    return dict(constructed_orbits=32,max_support_reconstruction_error=float(emax),
                max_derived_overlap_error=float(hmax),exact_boxes=exact_boxes(),central=central)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,(bool,str)):assert a==b
    else:assert np.isclose(a,b,rtol=3e-8,atol=1e-9),(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    res=compute()
    out=Path(__file__).with_name('independent_results.json')
    if args.check:compare(res,json.loads(out.read_text(encoding='utf8')))
    else:
        with out.open('x',encoding='utf8') as f:
            json.dump(res,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(res,allow_nan=False))
