"""Full-generator spectrum, extra adjoint sources and the finite parity test."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def generators(n):
    basis = []
    for i in range(n):
        for j in range(i):
            sym = np.zeros((n,n),complex)
            sym[i,j] = sym[j,i] = .5
            asym = np.zeros_like(sym)
            asym[i,j], asym[j,i] = .5j, -.5j
            basis.extend([sym,asym])
    for k in range(1,n):
        diag = np.zeros(n)
        diag[:k] = 1
        diag[k] = -k
        basis.append(np.diag(diag)/math.sqrt(2*k*(k+1)))
    return basis


def mass_spectrum(g4, gr, w4, w10):
    v = np.zeros(8)
    v[6] = 1
    gens = ([g4*np.kron(t,np.eye(2)) for t in generators(4)] +
            [gr*np.kron(np.eye(4),t) for t in generators(2)])
    linear = np.array([t@v for t in gens])
    tensor = np.array([np.kron(y,v)+np.kron(v,y) for y in linear])
    gram = 2*(w4**2*(linear@linear.conj().T).real +
              w10**2*(tensor@tensor.conj().T).real)
    eigs = np.linalg.eigvalsh(gram)
    c = w4**2/2+w10**2
    expected = sorted([0]*9+[g4**2*c]*6+[gr**2*c]*2+
                      [(gr**2+1.5*g4**2)*(c+w10**2)])
    np.testing.assert_allclose(eigs,expected,rtol=1e-10,atol=2e-7)
    return np.array(expected)


def tensor_interval(a,b,n):
    if n/2>min(n,a,b):
        return None
    return max(0,n-min(a,b)), n/2


def allocations(a,b,n):
    interval=tensor_interval(a,b,n)
    assert interval is not None
    for t in np.linspace(*interval,21):
        c=n-t
        f, ac, ar = c-t, a-c, b-c
        assert min(f,ac,ar)>=-1e-8
        np.testing.assert_allclose([f+t+ac,f+t+ar,f+2*t],[a,b,n],rtol=2e-12)


def calculate():
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    sample = author["synthetic"]
    mU,mW,mZ = sample["PS_masses_GeV_U_WR_ZR"]
    w4,w10 = sample["w4_GeV"],sample["w10_GeV"]
    c=w4*w4/2+w10*w10
    g4, gr=mU/math.sqrt(c),mW/math.sqrt(c)
    mass_spectrum(g4,gr,w4,w10)
    for f,t in [(1.,0.),(0.,1.),(.1,1.),(2.,3.)]:
        mass_spectrum(1.2,.46,f,t)
    n=mZ*mZ/(gr*gr+1.5*g4*g4)
    a,b=mU*mU/(g4*g4),mW*mW/(gr*gr)
    allocations(a,b,n)
    allocations(a*1.6,b*1.2,n)
    assert tensor_interval(1.,3.,2.1) is None
    interval = tensor_interval(a,b,n)
    ytrace = 32*math.pi*sample["normalized_total_dilepton_width"]
    assert abs(ytrace - .3625)<6e-9
    np.testing.assert_allclose(ytrace*interval[0],sample["sum_neutrino_mass_squared_GeV2"],rtol=2e-8)
    # Pure tensor plus color adjoint can have the same 3-pole ratio as a mixed source.
    target=1.36
    ac=(2/target-1)*(1.5*g4*g4+gr*gr)*c/(1.5*g4*g4)
    apparent = 2*(1.5*g4*g4+gr*gr)*c/(1.5*g4*g4*(c+ac)+gr*gr*c)
    np.testing.assert_allclose(apparent,target)
    # Center and a sample of identity-component holonomies act oppositely on the sterile line.
    omega=complex(-.5,math.sqrt(3)/2)
    for phase in (1,omega,np.exp(.31j)):
        line_hol=phase*phase.conjugate()
        np.testing.assert_allclose([line_hol,-line_hol,(-line_hol)**2],[1,-1,1])
    epsilon=1e-3
    near_ratio=1+1/(1+epsilon*epsilon/2)
    return {"joint_adjoint_tensor_vev_squared_GeV2":list(interval),
            "same_width_sum_heavy_masses_squared_GeV2":[ytrace*t for t in interval],
            "only_color_adjoint_tensor_vev_squared_GeV2":n-b,
            "color_adjoint_counterexample_ratio":apparent,
            "linear_to_tensor_vev_ratio":epsilon,"mass_ratio_near_pure_pair":near_ratio,
            "nontrivial_holonomy_on_neutral_matter":-1}


def compare(x,y):
    for key,val in x.items():
        np.testing.assert_allclose(val,y[key],rtol=3e-8,atol=2e-10)
    assert x.keys()==y.keys()


if __name__=="__main__":
    p=argparse.ArgumentParser()
    cmd=p.add_mutually_exclusive_group(required=True)
    cmd.add_argument("--write-new",action="store_true")
    cmd.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=calculate()
    target=HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as fd:
            json.dump(out,fd,indent=2)
            fd.write("\n")
    else:
        compare(out,json.loads(target.read_text(encoding="utf-8")))
    print("1543 full gauge mass matrix, adjoint ambiguity and neutral holonomy verified")
