"""Certified spectral exclusion on the full two-width rectangle. No scan of fitted poles."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("interval1486",HERE.parent/"1486/check.py")
iv=importlib.util.module_from_spec(spec)
spec.loader.exec_module(iv)
iv.SCALE=10**26
I=iv.I
one=I.point(1)
zero=I.point(0)


def read(name):
    return json.loads((HERE/name).read_bytes())


def bounds(a,b):
    return I(I.point(a).lo,I.point(b).hi)


def trig(x, sine):
    # All input radians lie in (0,pi/2). Evaluate endpoint alternating Taylor bounds using fractions.
    def endpoint(n):
        s=Fraction(n,iv.SCALE)
        k=1 if sine else 0
        total=sum((-1)**j*s**(k+2*j)/math.factorial(k+2*j) for j in range(12))
        next_term=s**(k+24)/math.factorial(k+24)
        return I.point(total), I.point(total+next_term)
    left=endpoint(x.lo)
    right=endpoint(x.hi)
    return I(left[0].lo,right[1].hi) if sine else I(right[0].lo,left[1].hi)


def cmul(z,w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def star(z):
    return z[0],-z[1]


def abs2(z):
    return z[0].square()+z[1].square()


def weighted_product(v,w,weights):
    vals=[cmul(vi,star(wi)) for vi,wi in zip(v,w)]
    return tuple(sum((weight*row[j] for weight,row in zip(weights,vals)),zero) for j in (0,1))


def coefficients(box, cs):
    d,k,p,u,t,delta=box
    # Sines are the subdivided variables, cosines stay attached through the exact positive root.
    cp,cu,ct=[(1-x.square()).sqrt() for x in [p,u,t]]
    dz=delta
    cd=dz
    sd=(1-cd.square()).sqrt()
    ca,cb,cc=cs
    aa,ss=d.square(),k.square()
    base=aa*cp.square()+ss*p.square()
    mix=(ss-aa)*p*cp
    A=base*cu.square()+u.square()  # squared A, not its square root
    # H12 and H13 (V diag(lambda) V†), split before the real theta23 rotation.
    h12=(cu*(ct*mix+t*u*(1-base)*cd), -cu*t*u*(1-base)*sd)
    h13=(cu*(-t*mix+ct*u*(1-base)*cd), -cu*ct*u*(1-base)*sd)
    C,E=abs2(h12)/A,abs2(h13)/A
    v2=((-p*ct-cp*t*u*cd, -cp*t*u*sd),(cp*ct-p*t*u*cd,-p*t*u*sd),(t*cu,zero))
    v3=((p*t-cp*ct*u*cd,-cp*ct*u*sd),(-cp*t-p*ct*u*cd,-p*ct*u*sd),(ct*cu,zero))
    weights=(ss,aa,aa*ss)
    minor=sum((w*abs2(z) for w,z in zip(weights,v2)),zero)
    F=minor/A
    B=aa*ss/minor
    # H11 H23-H21 H13 = -adj(H)23; its normalized norm is F² |Z|².
    adj=weighted_product(v2,v3,weights)
    Z=abs2(adj)/(A*minor)
    y=cmul(star(h12),h13)
    signed=-2*cc*(y[0]*adj[0]+y[1]*adj[1])/A.square()
    X,Y,W=cb**2*A,ca**2*B,cc**2*C
    K=X+Y+W
    # Sum(E²+F²) = H33 yields a considerably tighter identity than the separated terms.
    h33=aa*abs2(v3[0])+ss*abs2(v3[1])+abs2(v3[2])
    L=Z+h33
    h=(ca*cb)*d*k
    # Use D to eliminate K exactly from (1-w)w/(KL); split the complex
    # 23 minor into real and imaginary squares to avoid interval cancellation.
    base=X*(Z+F)+Y*h33
    jH=cp*p*ct*t*cu.square()*u*sd*(ss-aa)*(1-ss)*(1-aa)
    cf=(C*F).sqrt()
    e_im=jH/(A*cf)
    e_re=signed/(2*cc*cf)
    aux=(h,L,base,W*F,e_re,e_im,K)
    return [(ca*cb)*d*k/(K*L.sqrt()), X*Y/K.square(),
            (X*(Z+F)+Y*h33+W*F)/(K*L),
            signed/(K.sqrt()*L*L.sqrt()), E*Z/L.square(),aux]


def level_set(D,coef):
    # Endpoints of both monotone branches of w(1-w)^2. Dyadics with denominator a fixed multiple of 3.
    a=(D/coef).square()
    assert 0<a.lo<=a.hi<iv.SCALE*4//27
    unit=3*(1<<32)
    def isolate(value,rising):
        lo,hi=(0,unit//3) if rising else (unit//3,unit)
        while hi-lo>1:
            mid=(lo+hi)//2
            v=mid*(unit-mid)**2*iv.SCALE
            if (v<value*unit**3)==rising:lo=mid
            else:hi=mid
        return lo,hi
    for rising in (True,False):
        low=isolate(a.lo if rising else a.hi,rising)[0]
        high=isolate(a.hi if rising else a.lo,rising)[1]
        yield bounds(Fraction(low,unit),Fraction(high,unit))


def classifications(box,cs,D,P):
    co,u0,u2,u3,u4,aux=coefficients(box,cs)
    h,L,base2,wf,er,ei,K=aux
    result=[]
    for w in level_set(D,co):
        b=1-w
        base=u0*b.square()+u2*b*w+u4*w.square()
        off=u3*(b*w).sqrt()*w
        r=D*w.sqrt()/(h*L.sqrt())
        v=w/L
        unit=u0*b.square()+base2*r+(ei*v).square()
        linear=(wf*r).sqrt()
        for shape, positive in zip((base+off,base-off),
                    (unit+(linear+er*v).square(),unit+(linear-er*v).square())):
            shape=I(max(shape.lo,positive.lo),min(shape.hi,positive.hi))
            if shape.hi < P.lo:
                result.append((P.lo-shape.hi)/iv.SCALE)
            elif shape.lo > P.hi:
                result.append((shape.lo-P.hi)/iv.SCALE)
            else:
                result.append(0.)
    return result


def spectral_target(leptons, scale_width):
    ls=[iv.width(row) for row in leptons]
    # A positive electron-only threshold alpha brings model e to observed e.
    # Hence the uncorrected model has to land in the observed interval divided by alpha.
    err=I.point(scale_width)
    ls[0]=ls[0]/I(one.lo-err.hi,one.hi+err.hi)
    e,l=ls[0]/ls[2],ls[1]/ls[2]
    tr=1+e.square()+l.square()
    D=e*l/(tr*tr.sqrt())
    P=(e.square()*l.square()+e.square()+l.square())/tr.square()
    return D,P


def certificate(dom,cs,D,P):
    # Bisect the largest unrefined relative box, weighted by sensitivity of 1489's Schur reconstruction.
    weights=[.6,1.,.18,.75,.5,.55]
    full=[x.hi-x.lo for x in dom]
    todo=[(dom,0)]
    leaves=0
    max_depth=0
    visited=0
    gaps=[1.]*4
    while todo:
        box,depth=todo.pop()
        visited+=1
        out=classifications(box,cs,D,P)
        if min(out)>0:
            gaps=list(map(min,zip(gaps,out)))
            leaves+=1
            max_depth=max(depth,max_depth)
            continue
        if depth>=34:
            raise ArithmeticError(f"unresolved box at {depth}: {[b.values() for b in box]} {out}")
        coordinate=max(range(6),key=lambda i:weights[i]*(box[i].hi-box[i].lo)/full[i])
        for half in box[coordinate].parts(2):
            split=list(box)
            split[coordinate]=half
            todo.append((split,depth+1))
    return {"certified_leaves":leaves,"nodes":visited,"maximum_bisection_depth":max_depth,
            "minimum_positive_shape_gap_by_small_large_w_and_both_signs":gaps}


def gj_witness(source,leptons):
    # Choose actual inputs inside the original unexpanded rectangles, not arbitrary separate ratios.
    raw=source["down"]+source["angles_theta12_theta13_theta23_delta"]
    offsets=[1,-1,-1,-1,1,-1,1]
    point=[I.point(v)+I.point("1.999")*sign*I.point(err) for (v,err),sign in zip(raw,offsets)]
    box=[point[0]/point[2],point[1]/point[2]]+[trig(i,True) for i in point[3:6]]+[trig(point[6],False)]
    co,u0,u2,u3,u4,aux=coefficients(box,(1,1,3))
    h,L,_,_,_,_,K=aux
    t=I.point("1.19656937264")
    scale=I.point(".00802157218263")
    trace=K+L*t.square()
    w=L*t.square()/trace
    b=1-w
    shape=u0*b.square()+u2*b*w-u3*(b*w).sqrt()*w+u4*w.square()
    det2=(h*t).square()/(trace.square()*trace)
    results=[]
    for i,row in enumerate(leptons):
        region=iv.width(row)
        intervals=[]
        for endpoint in (region.lo,region.hi):
            lam=(I(endpoint,endpoint)/scale).square()/trace
            poly=((lam-1)*lam+shape)*lam-det2
            intervals.append(poly.values())
        if i==1: assert intervals[0][0]>0>intervals[1][1]
        else: assert intervals[0][1]<0<intervals[1][0]
        results.append(intervals)
    return {"number_of_marginal_widths_per_quark_specification":"1.999 times signs",
            "signs_d_s_b_theta12_theta13_theta23_delta":offsets,
            "t":t.values(),"k_times_yb":scale.values(),
            "normalized_characteristic_polynomial_at_each_lepton_window_end":results}


def compare_full_svd(source,models):
    import numpy as np
    loader=importlib.util.spec_from_file_location("r1657",HERE.parent/"1657/check.py")
    prev=importlib.util.module_from_spec(loader)
    loader.loader.exec_module(prev)
    edges=prev.make_texture(source)
    d,k,bottom=[I.point(row[0]) for row in source["down"]]
    rad=[I.point(row[0]) for row in source["angles_theta12_theta13_theta23_delta"]]
    box=[d/bottom,k/bottom]+[trig(x,True) for x in rad[:3]]+[trig(rad[3],False)]
    count=0
    for cs in models.values():
        _,u0,u2,u3,u4,(h,L,_,_,_,_,K)=coefficients(box,cs)
        for t in [-9,-4,-1,-.018,-.002,.003,.018,1.3,4,9]:
            x=I.point(t)
            trace=K+L*x.square()
            w=L*x.square()/trace
            b=1-w
            P=u0*b.square()+u2*b*w-u3*(1 if t>0 else -1)*(b*w).sqrt()*w+u4*w.square()
            D=h*abs(t)/(trace*trace.sqrt())
            singular=np.linalg.svd(prev.matrix(edges,cs,t),compute_uv=False)
            sq=singular**2/np.sum(singular**2)
            ref=np.array([sq[0]*sq[1]+sq[1]*sq[2]+sq[0]*sq[2],np.sqrt(sq.prod())])
            np.testing.assert_allclose([P.values()[0],D.values()[0]],ref,rtol=1e-10,atol=1e-13)
            count+=1
    return count


def compute():
    cfg=read("inputs.json")
    q=read(cfg["quark"])["high_scale"]
    l=read(cfg["leptons"])["scales"][cfg["scale"]]
    flavors=read(cfg["clebsch"])["clebsch"]
    down=list(map(iv.width,q["down"]))
    ang=list(map(iv.width,q["angles_theta12_theta13_theta23_delta"]))
    dom=[down[0]/down[2],down[1]/down[2]]+[trig(i,True) for i in ang[:3]]+[trig(ang[3],False)]
    D,P=spectral_target(l[3:],cfg["extra_electron_relative_change"])
    cert={name:certificate(dom,flavors[name],D,P) for name in cfg["certificate_models"]}
    return {"expanded_dimensionless_quark_domains":[x.values() for x in dom],"target_D":D.values(),
            "target_P":P.values(),"certificates":cert,"GJ_verified_witness":gj_witness(q,l[3:]),
            "central_full_matrix_SVD_tests":compare_full_svd(q,flavors)}


if __name__=="__main__":
    cli=argparse.ArgumentParser();cli.add_argument("--check",action="store_true")
    mode=cli.parse_args()
    out=compute()
    if mode.check:
        assert out==read("results.json")
    else:
        (HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out["certificates"],indent=1))
