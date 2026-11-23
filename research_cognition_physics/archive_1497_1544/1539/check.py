"""Single-edge cofactor, finite forbidden-edge tolerance, and common-vertex realization."""
import argparse
import cmath
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def bisect(fun,left=0.,right=1.):
    assert fun(left)<=0<=fun(right)
    for _ in range(90):
        mid=(left+right)/2
        if fun(mid)>0:
            right=mid
        else:
            left=mid
    return (left+right)/2


def interval(pair,n):
    v,d=pair
    return v-n*d,v+n*d


def masses(a,d,x):
    return np.sqrt([a+x*x,a+d+x*x,x*x])


def antisym(n):
    e,m,t=n
    return np.array([[0,t,-m],[-t,0,e],[m,-e,0]],dtype=complex)


def reconstruct(a,d,s,t,atm,x,charged):
    p=np.array([(1-t)*(1-s),(1-t)*s,t])
    ms=masses(a,d,x)
    lens=p/ms
    cos_phase=(lens[2]**2-lens[0]**2-lens[1]**2)/(2*lens[0]*lens[1])
    b=cmath.exp(1j*math.acos(cos_phase))
    phases=np.array([1,b,-(lens[0]+b*lens[1])/lens[2]])
    np.testing.assert_allclose(abs(phases),1,atol=2e-14)
    def plane(i,j,sine):
        v=np.eye(3)
        v[i,i]=v[j,j]=math.sqrt(1-sine**2)
        v[i,j],v[j,i]=sine,-sine
        return v
    u=plane(1,2,math.sqrt(atm))@plane(0,2,math.sqrt(t))@plane(0,1,math.sqrt(s))
    u=u@np.diag(np.sqrt(phases))
    mass=u.conj()@np.diag(ms)@u.conj().T
    ei=np.linalg.inv(mass)[:,0]
    assert abs(ei[0])<2e-12
    n=ei/np.linalg.norm(ei)
    n[0]=0
    f=antisym(n)
    dl=np.diag(np.square(charged)/charged[-1]**2)
    one=f@dl-dl@f
    alpha=(mass@n)[0]/(one@n)[0]
    rem=mass-alpha*one
    fi=np.linalg.pinv(f)
    k=fi@rem@fi.T
    np.testing.assert_allclose(alpha*one+f@k@f.T,mass,atol=4e-15)
    block=np.linalg.svd(mass[1:,1:],compute_uv=False)
    np.testing.assert_allclose(block[1],0,atol=1e-15)
    np.testing.assert_allclose(np.prod(ms),block[0]/np.linalg.norm(ei)**2,rtol=2e-13)
    mbb=abs(mass[0,0])
    formula=a*(a+d)*lens[2]**2-a*d*lens[0]**2+(a+d)*d*lens[1]**2
    np.testing.assert_allclose(mbb**2,formula,rtol=2e-12)
    return dict(masses_eV=ms.tolist(),mbb_eV=mbb,mbeta_eV=math.sqrt(sum(p*ms**2)),
                charged_left_fractions=((1-abs(n)**2)/2).tolist(),
                one_loop_response_eV=float(np.linalg.norm(mass@n)),two_loop_block_eV=block[0],
                reconstruction_error_eV=float(np.linalg.norm(mass-alpha*one-f@k@f.T)))


def compute():
    cfg=read("inputs.json")
    prev=read(cfg["oscillation"])
    n=cfg["widths"]
    reactor=lambda z:(1-math.sqrt(1-z))/2
    t1,t2=map(reactor,interval(prev["DayaBay_amplitude"],n))
    case=prev["IO"]
    s1,s2=interval(case["solar_sine_squared"],n)
    a1,a2=interval(case["atmospheric_gap"],n)
    d1,d2=interval(case["solar_gap"],n)
    left=math.sqrt(a1)/math.sqrt(((1-t1)/t1)**2-1)
    right=math.sqrt(a2)/math.sqrt(((1-t2)*(1-2*s2)/t2)**2-1)
    beta1=math.sqrt(left**2+(1-t2)*a1+(1-t2)*s1*d1)
    beta2=math.sqrt(right**2+(1-t1)*a2+(1-t1)*s2*d2)
    # M_ee lower bound from the first triangle term; upper bound by rms.
    lo=(1-t2)*(1-s2)*math.sqrt(left**2+a1)-(1-t1)*s2*math.sqrt(right**2+a2+d2)-t2*right
    s,t=case["solar_sine_squared"][0],reactor(prev["DayaBay_amplitude"][0])
    a,d=case["atmospheric_gap"][0],case["solar_gap"][0]
    def ratio(x,sign):
        return (1-t)*x*((1-s)/math.sqrt(x*x+a)+sign*s/math.sqrt(x*x+a+d))-t
    central=[bisect(lambda x:ratio(x,k)) for k in (1,-1)]
    y=np.array([v[0] for v in prev["yukawa"]["charged_leptons"]])
    err=np.array([v[1] for v in prev["yukawa"]["charged_leptons"]])
    l,u=(y-n*err)**2,(y+n*err)**2
    # Lower bound on d_tau-d_mu; charged Yukawa normalizations cancel.
    gain=math.hypot(u[2]-l[0],u[1]-l[0])/(l[2]-u[1])/math.sqrt(cfg["minimum_other_null_flavour_shares"])
    s_no=interval(prev["NO"]["solar_sine_squared"],n)[1]
    gap=2*(1-t2)*(1-s_no)-1
    edge=bisect(lambda x:x*(gain+math.sqrt(1+(gain*x)**2))-gap)
    witness=reconstruct(a,d,s,t,case["atmospheric_sine_squared"][0],cfg["m3_witness_eV"],y)
    return dict(IO_m3_central_interval_eV=central,IO_m3_outer_eV=[left,right],
                beta_outer_eV=[beta1,beta2],mbb_outer_eV=[lo,beta2],
                NO_inverse_or_adj_gap=gap,NO_edge_amplitude_min=edge,
                NO_edge_fraction_min=edge**2,NO_e_leg_max=(1-edge**2)/2,
                holdout_mbb_upper_eV=cfg["kamland_zen_mbb_upper_eV"],
                holdout_beta_upper_eV=cfg["katrin_beta_upper_eV"],witness=witness)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    opt=p.parse_args()
    result=compute()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        def compare(x,y):
            if isinstance(x,dict):
                assert x.keys()==y.keys()
                for k in x:
                    compare(x[k],y[k])
            else:
                np.testing.assert_allclose(x,y,rtol=3e-10,atol=2e-14)
        compare(result,read("results.json"))
        print("1539: cofactor, ordering, finite missing-edge bound and two-loop reconstruction verified.")
