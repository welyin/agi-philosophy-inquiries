"""Primal and single-ray lower certificates for the absolute adjoint Wilson coefficient."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def golden(f,lo,hi):
    k=(math.sqrt(5)-1)/2
    left,right=hi-k*(hi-lo),lo+k*(hi-lo)
    fl,fr=f(left),f(right)
    for _ in range(70):
        if fl<fr:
            hi,right,fr=right,left,fl
            left=hi-k*(hi-lo)
            fl=f(left)
        else:
            lo,left,fl=left,right,fr
            right=lo+k*(hi-lo)
            fr=f(right)
    return (lo+hi)/2


def optimal(d,e,real_only):
    def norm(t):
        return np.linalg.norm(d-t*e,2)
    # These searches only find upper-bound candidates; the global proof below
    # has no restriction on either component of t.
    if real_only:
        t=complex(golden(norm,-3.,3.))
    else:
        imag=lambda x:golden(lambda y:norm(x+1j*y),-3.,3.)
        x=golden(lambda x:norm(x+1j*imag(x)),-3.,3.)
        t=x+1j*imag(x)
    def gradient(tt):
        u,_,vh=np.linalg.svd(d-tt*e)
        q=np.vdot(u[:,0],e@vh[0].conj())
        return np.array([q.real]) if real_only else np.array([q.real,q.imag])
    directions=[1.] if real_only else [1.,1j]
    for _ in range(3):
        h=1e-6
        jac=np.column_stack([(gradient(t+h*p)-gradient(t-h*p))/(2*h) for p in directions])
        t-=sum(a*b for a,b in zip(directions,np.linalg.solve(jac,gradient(t))))
    _,s,vh=np.linalg.svd(d-t*e)
    v=vh[0].conj()
    v*=np.exp(-1j*np.angle(v[np.argmax(abs(v))]))
    Ev,Dv=e@v,d@v
    tv=np.vdot(Ev,Dv)/np.vdot(Ev,Ev).real
    if real_only:
        tv=tv.real
    # For *all* t, ||D-tE|| >= ||(D-tE)v|| >= dist(Dv,span(Ev)).
    lower=np.linalg.norm(Dv-tv*Ev)
    upper=s[0]
    assert -2e-16<upper-lower<5e-15
    guard=1e-13
    pair=lambda z:[float(z.real),float(z.imag)]
    return {"t":pair(t),"z":pair(1/t),"singular_values":s.tolist(),
            "unit_5bar_source_real_imag":[pair(x) for x in v],
            "single_ray_minimizing_t":pair(tv),
            "absolute_residual_bracket":[lower-guard,upper+guard]}


def calculate():
    f=runpy.run_path(str(HERE/"../1715/check.py"))
    _,rep=f["charged_background"](.1)
    d=f["unpack"](rep["down_texture_10PeV"])
    e=f["unpack"](rep["charged_texture_10PeV"]).T
    opt={name:optimal(d,e,real) for name,real in (("complex",False),("real",True))}
    ld,ds,rd=np.linalg.svd(d)
    le,es,re=np.linalg.svd(e)
    margins=abs(ds[:,None]*es[None,:]-ds[None,:]*es[:,None])/(es[:,None]+es[None,:])
    i,j=np.unravel_index(np.argmax(margins),margins.shape)
    spec_only=margins[i,j]
    r=(ds[i]+ds[j])/(es[i]+es[j])
    eligned=ld@le.conj().T@e@re.conj().T@rd
    np.testing.assert_allclose(np.linalg.norm(d-r*eligned,2),spec_only,atol=1e-15)
    gauge=read("../1720/results.json")
    U=gauge["baseline_unification_GeV"]
    choice=gauge["central_conditional_combined_point"]
    beta=math.sqrt(1-choice["R"]**-2)
    g=math.sqrt(4*math.pi/choice["common_inverse_coupling"])
    cutoff=read("../1518/inputs.json")["cutoff_GeV"]
    rank2=min(row["lift_singular_values"][0] for row in
              read("../1723/results.json")["exact_and_conditional_implementations"].values() if row["rank"]==2)
    strengths={"unrestricted_complex":opt["complex"]["singular_values"][0]/beta,
               "real":opt["real"]["singular_values"][0]/beta,
               "common_z1":np.linalg.norm(e-d,2)/beta,"exactly_two_paths":rank2,
               "independent_unitary_frames":spec_only/beta}
    wilsons={k:{"minimum_MX_for_unit_Wilson_GeV":v*g*cutoff,
                "minimum_MX_for_sqrt4pi_budget_GeV":v*g*cutoff/math.sqrt(4*math.pi),
                "Wilson_at_MX_equal_U":v*g*cutoff/U}
             for k,v in strengths.items()}
    source_budget=read("../1723/results.json")["five_path_tree_match"]["required_common_vertex_norm_at_this_boundary"]
    def leakage(allowed_vertex):
        pulled=max(0.,spec_only-allowed_vertex)/source_budget
        max_overlap=1-.5*pulled**2
        return 1-max_overlap**2
    return {"tagged":opt,"independently_rotated_targets":{"mass_only_lower":spec_only,"optimal_modulus_t":r,
                "active_singular_pair_heavy_to_light":[int(i),int(j)]},
            "conditions":{"MX_example_GeV":U,"g5_example":g,"cutoff_GeV":cutoff,"right_projection":beta},
            "absolute_Wilson_bounds":wilsons,
            "unmatched_matter_subspaces":{"common_vertex_norm_budget":source_budget,
                "minimum_loss_in_one_conversion_probability_without_adjoint":leakage(0.),
                "minimum_loss_with_unit_Wilson_if_MX_is_U":leakage(U/(g*cutoff))}}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list): return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.15g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    do_check=p.parse_args().check
    results=tidy(calculate())
    output=HERE/"independent_results.json"
    if do_check:
        assert results==json.loads(output.read_bytes())
        print("1726 verified global complex/real Wilson minima, relaxed-frame and overlap bounds")
    else:
        with output.open("x",encoding="utf8") as stream:
            json.dump(results,stream,indent=2,allow_nan=False)
            stream.write("\n")
