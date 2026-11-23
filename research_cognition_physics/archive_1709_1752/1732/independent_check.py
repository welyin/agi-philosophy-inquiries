"""Exact closed-port moment resources and same-first-step counterexamples."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def moments(unitary,n):
    np.testing.assert_allclose(unitary.conj().T@unitary,np.eye(len(unitary)),atol=2e-14)
    v=np.eye(len(unitary),dtype=complex)[:,0]
    out=[v[0]]
    for _ in range(n):
        v=unitary@v
        out.append(v[0])
    return np.array(out)


def shift(c,k,phase=1):
    """A lossless k+1-state dilation; the last resource sets the first return."""
    s=math.sqrt(1-c*c)
    U=np.zeros((k+1,k+1),complex)
    U[0,0],U[1,0]=c,s
    for j in range(1,k):
        U[j+1,j]=1
    U[0,k],U[1,k]=-s*phase,c*phase
    return U


def complex_prefix(c,h,z):
    s=math.sqrt(1-c*c)
    b=math.sqrt(max(0,1-abs(h)**2))
    w=np.array([-s,c,0.],complex)
    ancilla=np.array([0,0,1.])
    return np.column_stack([np.array([c,s,0.]),-h*w+b*ancilla,z*(-b*w-h.conjugate()*ancilla)])


def check(c):
    a,r=c*c,1-c*c
    powers=[c**j for j in range(7)]
    # Rank-one alternating legs, with a phase convention deliberately not real.
    u=np.array([1.,0.]); d=np.array([c*np.exp(.37j),math.sqrt(r)])
    v=d.copy()
    for j in range(1,7):
        ray=u if j%2 else d
        v=ray*np.vdot(ray,v)
        np.testing.assert_allclose(np.vdot(ray,v),c**j*np.exp(.37j*(j%2)),atol=3e-15)
    dimensions={}
    for k in range(1,7):
        idx=np.arange(k+1)
        gram=c**np.abs(idx[:,None]-idx[None,:])
        np.testing.assert_allclose(np.linalg.det(gram),r**k,rtol=2e-13)
        eig=np.linalg.eigvalsh(gram)
        phase=np.exp(.48j)
        seq=moments(shift(c,k,phase),k+1)
        np.testing.assert_allclose(seq[:k+1],powers[:k+1],atol=2e-14)
        np.testing.assert_allclose(seq[-1],c**(k+1)-r*phase,atol=2e-14)
        dimensions[str(k)]={"dimension":k+1, "gram_determinant":r**k,
                            "complex_moment_error_floor_if_one_dimension_short":eig[0]/k}
    # Keep the first column exactly fixed. Pick a two-state phase for which U^3 is scalar.
    turning=2*math.acos(1/(2*c))
    low=moments(shift(c,1,1),3)[3]
    high=moments(shift(c,1,np.exp(1j*turning)),3)[3]
    np.testing.assert_allclose([abs(low),abs(high)],[abs(4*c**3-3*c),1.],atol=1e-14)
    # Impose |A2|=c^2 without specifying its phase. In dimension two |h|=1.
    h2=np.exp(1j*math.acos(-r/(2*a)))
    seq2=moments(shift(c,1,-h2),3)
    np.testing.assert_allclose(abs(seq2[1:3]),powers[1:3],atol=2e-14)
    np.testing.assert_allclose(abs(seq2[3])**2,c**6+1-c**4,atol=2e-14)
    # Solve for a genuine three-state completion with A3=c^3 and |A2|=c^2.
    coeff=4*a+2
    t=2/(coeff+math.sqrt(coeff*coeff+4*r))
    re=-r*t/(2*a)
    h3=complex(re,math.sqrt(t-re*re))
    z=-c*(2*h3-h3**2)/(1-t)
    np.testing.assert_allclose(abs(z),1.,atol=3e-15)
    seq3=moments(complex_prefix(c,h3,z),4)
    np.testing.assert_allclose(abs(seq3[:4]),powers[:4],atol=2e-14)
    np.testing.assert_allclose(seq3[3],powers[3],atol=3e-15)
    return {"c":c,"c_cubed":powers[3],"c_seventh":c**7,
            "success_probability_if_three_arms_are_filtered":c**6,
            "two_dimensional_same_full_first_column_third_modulus_range":[abs(low),abs(high)],
            "mass_product_range_if_x_is_unfiltered_third_modulus":[c*abs(low)**2,c],
            "minimum_dimension_per_real_prefix":dimensions,"return_circle_radius_for_exact_prefix_0_1_2":r,
            "two_state_third_modulus_with_first_two_moduli_fixed":abs(seq2[3]),
            "qutrit_with_all_three_moduli_and_real_last":{"phase_of_second_prefix_deg":np.angle(seq3[2],deg=True),
                    "relative_phase_h":[h3.real,h3.imag],"unused_branch_closure_z":[z.real,z.imag],
                    "fourth_moment":[seq3[4].real,seq3[4].imag]}}


def run():
    live=json.loads((HERE/"results.json").read_bytes())["solutions"]
    return {name:check(next(row["Vud_10PeV"] for row in rows if row["k"]==3))
            for name,rows in live.items()}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    results=run()
    target=HERE/"independent_results.json"
    if args.check:
        assert json.loads(target.read_bytes())==results
        print("1732 independent open and closed transfer certificates verified")
    else:
        target.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf-8")
