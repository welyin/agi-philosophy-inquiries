"""An electron-only oscillation readout of a tagged Pati-Salam source row."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
BASE=json.loads((HERE/"../1386/inputs.json").read_text(encoding="utf-8"))
PAIRS=((0,1),(0,2),(1,2))


def guarantees():
    v12,v13,v23=[v+2*err for v,err in BASE["quark_sines_and_phase"][:3]]
    c12,c13,c23=[math.sqrt(1-s*s) for s in (v12,v13,v23)]
    v=min(c12*c13,c12*c23-v12*v23*v13,c23*c13)
    amp,err=BASE["DayaBay_amplitude"]
    p3min=(1-math.sqrt(1-amp+2*err))/2
    s12min=min(BASE[k]["solar_sine_squared"][0]-2*BASE[k]["solar_sine_squared"][1]
               for k in ("NO","IO"))
    pmax=(1-p3min)*(1-s12min)
    def near(q):
        return max(0.,v*math.sqrt(q)-math.sqrt(1-v*v)*math.sqrt(1-q))**2
    lo,hi=0.,1.
    for _ in range(90):
        u=(lo+hi)/2
        if (near(1-u)-pmax)/math.sqrt(6)>u:
            lo=u
        else:
            hi=u
    return {"minimum_diagonal_V":v,"maximum_electron_column_share":pmax,
            "near_min_at_q_99":near(.99),"near_minus_average_min_at_q_99":near(.99)-pmax,
            "oscillation_std_min_at_q_99":(near(.99)-pmax)/math.sqrt(6),
            "uniform_down_impurity_or_up_oscillation_std_floor":lo}


def design(times,squares):
    phases=2.533865*np.asarray(times)[:,None]*np.array([squares[j]-squares[i] for i,j in PAIRS])
    return np.column_stack([np.ones(len(times)),np.cos(phases),np.sin(phases)])


def electron_spectrum(times,squares,p,q,angles):
    amplitudes=np.sqrt(p*q)*np.exp(1j*np.asarray(angles))
    signal=np.abs(np.exp(-2.533865j*np.outer(times,squares))@amplitudes)**2
    return signal


def recover(beta,p):
    amp2=(beta[1:4]**2+beta[4:]**2)/4
    products=np.array([amp2[k]/(p[i]*p[j]) for k,(i,j) in enumerate(PAIRS)])
    nonzero=np.flatnonzero(products>1e-16)
    if len(nonzero)==3:
        ab,ac,bc=products
        return np.sqrt([ab*ac/bc,ab*bc/ac,ac*bc/ab])
    if len(nonzero)==1:
        i,j=PAIRS[nonzero[0]]
        q=np.zeros(3)
        q[i]=(beta[0]-p[j])/(p[i]-p[j])
        q[j]=1-q[i]
        return q
    assert len(nonzero)==0
    return np.eye(3)[np.argmin(np.abs(p-beta[0]))]


def calculate():
    raw=BASE["NO"]
    amp=BASE["DayaBay_amplitude"][0]
    p3=(1-math.sqrt(1-amp))/2
    solar=raw["solar_sine_squared"][0]
    p=np.array([(1-p3)*(1-solar),(1-p3)*solar,p3])
    delta=raw["solar_gap"][0]
    masses=np.array([0,delta,delta+raw["atmospheric_gap"][0]])
    period=2*math.pi/(2.533865*delta)
    times=np.linspace(0,period,257)
    matrix=design(times,masses)
    sigma=np.linalg.svd(matrix,compute_uv=False)
    assert sigma[-1]>10 and sigma[0]/sigma[-1]<1.5
    examples=[]
    for q in (np.array([.18,.51,.31]),np.array([.3,0,.7]),np.array([0,1.,0])):
        prob=electron_spectrum(times,masses,p,q,[.31,-.86,1.57])
        coeff=np.linalg.lstsq(matrix,prob,rcond=None)[0]
        found=recover(coeff,p)
        np.testing.assert_allclose(found,q,atol=5e-13)
        variance=2*sum(p[i]*q[i]*p[j]*q[j] for i,j in PAIRS)
        np.testing.assert_allclose(variance,.5*np.dot(coeff[1:],coeff[1:]),atol=1e-14)
        examples.append({"q":q.tolist(),"reconstructed_q":found.tolist(),"oscillation_std":math.sqrt(variance)})
    # Same fully averaged electron fraction, distinct source rows, distinct coherent curves.
    direction=np.array([p[1]-p[2],p[2]-p[0],p[0]-p[1]])
    qa=np.ones(3)/3
    qb=qa+.35*direction
    np.testing.assert_allclose(np.dot(p,qa),np.dot(p,qb))
    assert min(qb)>0 and np.linalg.norm(qa-qb)>.1
    a=electron_spectrum(times,masses,p,qa,[0,0,0])
    b=electron_spectrum(times,masses,p,qb,[0,0,0])
    short=design(np.linspace(0,100.,257),masses)
    return {"finite_box":guarantees(),"central_solar_period_km_per_GeV":period,
            "design_condition_full_period":float(sigma[0]/sigma[-1]),
            "one_bin_uniform_error_to_coeff_norm_upper":float(math.sqrt(len(times))/sigma[-1]),
            "condition_100_km_per_GeV":float(np.linalg.cond(short)),
            "source_row_examples":examples,"equal_decohered_electron_fraction":float(p@qa),
            "largest_bin_separation_for_same_average":float(np.max(np.abs(a-b)))}


def compare(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for k in x:
            compare(x[k],y[k])
    elif isinstance(x,list):
        for u,v in zip(x,y,strict=True):
            compare(u,v)
    else:
        np.testing.assert_allclose(x,y,rtol=.005 if abs(x)>1e13 else 2e-7,atol=1e-11)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    args=parser.add_mutually_exclusive_group(required=True)
    args.add_argument("--write-new",action="store_true")
    args.add_argument("--check",action="store_true")
    opt=parser.parse_args()
    values=calculate()
    name=HERE/"independent_results.json"
    if opt.write_new:
        with name.open("x",encoding="utf-8") as stream:
            json.dump(values,stream,indent=2)
            stream.write("\n")
    else:
        compare(values,json.loads(name.read_text(encoding="utf-8")))
    print("1546 independent: CKM/oscillation tradeoff and coherent electron-only inversion pass")
