"""Strong-CP matching test for three named five-bar sources on the actual 1748 solution."""
import argparse
import json
from pathlib import Path
import runpy
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = runpy.run_path(str(HERE.parent / "1748/check.py"))
EYE = np.eye(3)


def proj(line):
    return np.outer(line, line.conj())


def unitary_step(generator, step):
    values, vec = np.linalg.eigh(generator/1j)
    return (vec*np.exp(1j*step*values)) @ vec.conj().T


def background(step):
    b=PARENT["OLD"]["mother"](step)
    down=b["yW"][1,:,:3] @ b["right_quark_to_GJ"]
    strange=b["tagged_light_columns"][:,1]
    mu=np.linalg.svd(b["yW"][2,:,:3])[0][:,1]
    ps,pm=proj(strange),proj(mu)
    # Eigenspaces of i times the commutator set the sharp orientation-only upper bound.
    probe=(ps@pm-pm@ps)/(2j)
    vals,basis=np.linalg.eigh(probe)
    overlap=abs(np.vdot(strange,mu))
    np.testing.assert_allclose(vals,[-overlap*np.sqrt(1-overlap**2)/2,0,
                                     overlap*np.sqrt(1-overlap**2)/2],atol=1e-14)
    edmlim=json.loads((HERE.parent/"1389/results.json").read_bytes())["conditional_edm"]["theta_limit_rad"]
    rows=[]
    for root in (0,1):
        for signs in PARENT["CFG"]["family_signs"]:
            x,_,_,_=PARENT["terms"](b,root,signs)
            k=x.conj().T@x
            w=ps@pm@k
            z=np.trace(w)
            gram=np.linalg.norm(w,"fro")**2-abs(z)**2
            assert gram>=0
            np.testing.assert_allclose(z.imag,np.trace(probe@k).real,atol=1e-15)
            l=np.linalg.eigvalsh(k)
            envelope=vals[-1]*(l[-1]-l[0])
            best=(basis*l)@basis.conj().T
            np.testing.assert_allclose(np.trace(probe@best),envelope,atol=1e-14)
            for a,c in ((.7,.2),(-1.,.4),(2.,2.),(.03,0.)):
                det=1+(a+c)*z.real-a*c*gram+1j*(a-c)*z.imag
                np.testing.assert_allclose(np.linalg.det(EYE+a*w+c*w.conj().T),det,atol=1e-14)
            tangent=np.tan(edmlim)
            unit=tangent/abs(z.imag)
            center_branch_edges=[-tangent/(abs(z.imag)+tangent*z.real),
                                 tangent/(abs(z.imag)-tangent*z.real)]
            assert center_branch_edges[0]<0<center_branch_edges[1]
            for a in center_branch_edges:
                np.testing.assert_allclose(abs(np.angle(1+a*z)), edmlim,rtol=1e-12)
            # A physical relative unitary across all colored light right slots leaves
            # masses and the left weak current fixed but shifts the full determinant.
            q=unitary_step((w-w.conj().T).conj(),.25)
            light=down@q
            np.testing.assert_allclose(light@light.conj().T,down@down.conj().T,atol=3e-16)
            full=np.block([[down,np.zeros((3,2))],[x.conj(),np.eye(2)]])
            extend=np.eye(5,dtype=complex)
            extend[:3,:3]=q
            transformed=full@extend
            np.testing.assert_allclose(transformed@transformed.conj().T,full@full.conj().T,atol=3e-14)
            # Direct full-rank numerical determinant ratio in the same bare-theta basis.
            actual=np.angle(np.linalg.slogdet(transformed)[0]/np.linalg.slogdet(full)[0])
            np.testing.assert_allclose(actual,-.5*z.imag,atol=4e-13)
            rows.append(dict(root=root+1,signs=signs,ordered_trace=[z.real,z.imag],
                             complementary_area=gram,maximum_CP_trace_for_any_orientation=envelope,
                             fraction_of_sharp_max=abs(z.imag)/envelope,
                             pure_order_tangent_bound=unit,exact_real_epsilon_EDM_interval=center_branch_edges,
                             right_unitary_delta_theta_for_epsilon_quarter=actual))
    for i in (0,4):
        for a,c in ((i,i+1),(i+2,i+3)):
            np.testing.assert_allclose(rows[a]["ordered_trace"][1],-rows[c]["ordered_trace"][1],atol=5e-15)
    return dict(weak_threshold_GeV=b["Mweak"],strange_mu_overlap=overlap,
                conditional_neutron_angle_rad=edmlim,ports=rows)


def native(x):
    if isinstance(x,dict):
        return {k:native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [native(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    cli=ap.parse_args()
    result=native(background(.125))
    fine=native(background(.0625))
    result["step_halved_max_absolute_trace_im_change"]=max(abs(a["ordered_trace"][1]-b["ordered_trace"][1])
        for a,b in zip(result["ports"],fine["ports"]))
    path=HERE/"independent_results.json"
    if cli.check:
        assert result==json.loads(path.read_bytes())
        print("1752 three-source invariant, sharp orientation envelope and determinant phases verified")
    else:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
