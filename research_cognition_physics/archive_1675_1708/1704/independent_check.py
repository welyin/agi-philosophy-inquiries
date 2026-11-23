"""Optimal two-port overlap and finite-momentum root identification."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
CHI=np.exp(-2j*np.pi*np.arange(3)/3)


def best_ratio(z):
    # For the outside edge between the first two vertices, the extremum
    # uses its two complex affine coordinates.
    a,b=CHI[:2]
    w=np.array([(z-b)/(a-b),(a-z)/(a-b),0.])
    c=1/np.sum(abs(w))
    i=np.sqrt(c*abs(w))
    f=i*np.exp(-1j*np.angle(w))
    np.testing.assert_allclose([np.vdot(f,i),np.vdot(f,CHI*i)],[c,c*z],atol=5e-15)
    p=CHI.conj()+z*CHI
    far=np.argmax([abs(p[1]-p[2]),abs(p[2]-p[0]),abs(p[0]-p[1])])
    ray=[p[j]-p[far] for j in range(3) if j!=far]
    cosine=(ray[0]*ray[1].conjugate()).real/abs(ray[0]*ray[1])
    assert cosine<-.5
    np.testing.assert_allclose(3/(abs(ray[0])+abs(ray[1])),c,atol=2e-14)
    return float(c)


def best_absolute(a):
    z=np.exp(1j*math.pi/3)*a
    x,y=z.real,z.imag
    c=(x+math.sqrt(x*x-4*(1-y*y)*(abs(z)**2-.75)))/(2*(1-y*y))
    assert abs(a)<c<1
    np.testing.assert_allclose(best_ratio(a/c),c,atol=2e-14)
    return c


def shape(t,g2,q,r,ref):
    return (q/ref)*(t+g2*(r+ref*ref))/(t+g2*(r+q*q))


def two_probe(q,r,ref,g2):
    s=shape(0.,g2,q,r,ref)
    t=np.array([.4,1.7])
    F=shape(t,g2,q,r,ref)
    normalized_rates=(F/s)**2
    aux=t/(1-np.sqrt(normalized_rates))
    slope=(aux[1]-aux[0])/(t[1]-t[0])
    q_fit=ref*s*(1-1/slope)
    B=(aux[0]-slope*t[0])/slope
    g2_fit=B/(r+q_fit*q_fit)
    np.testing.assert_allclose([q_fit,g2_fit],[q,g2],rtol=2e-13)
    return {"normalized_rates":normalized_rates.tolist(),"q_fit":q_fit,"g2_fit":g2_fit}


def main():
    src=json.loads((HERE/"../1702/results.json").read_bytes())["one_pole_with_unit_neutral"]
    out={}
    for key,row in src["witnesses"].items():
        k=int(key[-1])
        r=src["r"]
        qlo,qhi=row["charged_bare_moduli"]
        roots=[]
        for q,g2 in zip((qlo,qhi),(.8,1.9)):
            z=q*np.exp(-1j*k*np.pi/6)
            cr,ca=best_ratio(z),best_absolute(z)
            roots.append({"q":q,
                          "same_calibrated_ratio_reference_amplitude_max":cr,
                          "same_calibrated_ratio_reference_probability_loss_min":1-cr*cr,
                          "same_absolute_raw_amplitude_reference_amplitude_max":ca,
                          "same_absolute_raw_amplitude_reference_probability_loss_min":1-ca*ca,
                          "unit_reference_two_probes":two_probe(q,r,1.,g2),
                          "neutral_reference_two_probes":two_probe(q,r,r,g2)})
        single={}
        for name,ref in [("ordinary",1.),("neutral",r)]:
            s=shape(0.,1,qlo,r,ref)
            a,b=qlo/ref,qhi/ref
            w=(max(a,b)+s)/2 if max(a,b)<s else (min(a,b)+s)/2
            fits=[(q/ref-w)/((w-s)*(r+q*q)) for q in (qlo,qhi)]
            assert min(fits)>0
            np.testing.assert_allclose([shape(1.,g,q,r,ref) for q,g in zip((qlo,qhi),fits)],w)
            single[name]={"same_F_at_one_transfer":w,"positive_g2_for_two_roots":fits}
        out[key]={"roots":roots,"one_point_degeneracy":single}
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    val=main()
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==val
        print("1704 independent verified")
    else:
        path.write_text(json.dumps(val,indent=2)+"\n",encoding="utf8")
