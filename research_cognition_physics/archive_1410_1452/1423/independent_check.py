"""Separate an unknown coherent matter edge from the vacuum-oriented loop using two energies."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def rotate(i,j,sine,delta=0.):
    matrix=np.eye(3,dtype=complex)
    matrix[i,i]=matrix[j,j]=np.sqrt(1-sine*sine)
    matrix[i,j]=sine*np.exp(-1j*delta)
    matrix[j,i]=-matrix[i,j].conjugate()
    return matrix


def loop(m):
    return float((m[0,1]*m[1,2]*m[2,0]).imag)


def measure(m, energy, baseline, vacuum_discriminant):
    vals, u=np.linalg.eigh(m)
    diff=np.array([vals[1]-vals[0],vals[2]-vals[0],vals[2]-vals[1]])
    k=1.26693268*baseline/energy
    propagator=(u*np.exp(-2j*k*vals))@u.conj().T
    delta_p=float(abs(propagator[1,0])**2-abs(propagator[0,1])**2)
    factor=float(16*vacuum_discriminant*np.prod(np.sin(k*diff)/diff))
    assert abs(factor)>1e-5
    extracted=delta_p/factor
    assert abs(extracted-loop(m)/vacuum_discriminant)<2e-12
    return extracted,delta_p,factor


def compute():
    angles=json.loads((HERE.parent/"1386"/"inputs.json").read_text(encoding="utf-8"))
    splittings=json.loads((HERE.parent/"1370"/"inputs.json").read_text(encoding="utf-8"))
    phase=-angles["quark_sines_and_phase"][3][0]
    y=(1-np.sqrt(1-angles["DayaBay_amplitude"][0]))/2
    e0,e1=1.7,2.5
    rho,ye=2.84,.5
    potential=1.526e-4*rho*ye
    outputs={}
    for order in ("NO","IO"):
        xx=angles[order]["solar_sine_squared"][0]
        zz=angles[order]["atmospheric_sine_squared"][0]
        mix=rotate(1,2,np.sqrt(zz))@rotate(0,2,np.sqrt(y),phase)@rotate(0,1,np.sqrt(xx))
        sp=splittings[order]
        ds,da=sp["solar"][0],sp["atm_abs"][0]
        lam=np.array([0,ds,da if order=="NO" else ds-da])
        disc=float(abs(ds*lam[2]*(lam[2]-ds)))
        m=(mix*lam)@mix.conj().T
        # Choose row phases so both e-mu and mu-tau edges are real.
        diag=np.diag(np.exp(1j*np.array([-np.angle(m[0,1]),0.,np.angle(m[1,2])])))
        m=diag@m@diag.conj().T
        assert abs(m[0,1].imag)+abs(m[1,2].imag)<1e-17
        epsilon=np.zeros((3,3),complex)
        epsilon[0,2]=(m[0,2].conjugate()-m[0,2])/(potential*e0)
        epsilon[2,0]=epsilon[0,2].conjugate()
        electron=np.diag([1.,0.,0.])
        assert np.allclose(m+potential*e0*(electron+epsilon),m.conj()+potential*e0*electron, atol=1e-17)
        readings=[measure(m+potential*e*(electron+epsilon),e,1284.852,disc) for e in (e0,e1)]
        fitted=(e1*readings[0][0]-e0*readings[1][0])/(e1-e0)
        vacuum=loop(m)/disc
        assert np.isclose(fitted,vacuum) and np.isclose(readings[0][0],-vacuum)
        outputs[order]={"vacuum_normalized_loop":vacuum,"single_energy_exact_mirror_abs_epsilon":float(abs(epsilon[0,2])),
                        "two_energy_reconstruction":fitted,"T_difference_at_E0_E1":[w[1] for w in readings],
                        "probability_per_normalized_loop_at_E0_E1":[w[2] for w in readings],
                        "uniform_extracted_loop_error_below_for_sign":abs(vacuum)*(e1-e0)/(e1+e0)}
    return {"round":1423,"constant_density_g_cm3":rho,"electron_fraction":ye,"energies_GeV":[e0,e1],
            "error_amplification":(e1+e0)/(e1-e0),"branches":outputs}


def match(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(match(a[k],b[k]) for k in a)
    return bool(np.allclose(a,b,atol=1e-13,rtol=1e-8))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=compute()
    if args.check:
        assert match(result,json.loads((HERE/"independent_results.json").read_text(encoding="utf-8")))
        print("1423 homogeneous mirror and two-energy signed loop inversion verified")
    else:
        print(json.dumps(result,indent=2))
