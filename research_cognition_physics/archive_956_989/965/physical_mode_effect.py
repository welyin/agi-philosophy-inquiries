"""965: bounded physical-mode effect in the fixed length representation.
The field quadrature is shifted by the material polarization. We evaluate the
displaced effect using exact finite-to-(0,1) oscillator matrix elements, not a
truncated unitary. The model and the original numerical files are unchanged.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/"physical_mode_effect_results.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text("utf-8-sig"))
def displacement_rows(alpha,n):
    p=np.arange(n)
    # <0|D(alpha)|p>, <1|D(alpha)|p>, real alpha.
    fact=np.array([math.sqrt(math.factorial(int(k))) for k in p])
    zero=math.exp(-alpha*alpha/2)*(-alpha)**p/fact
    one=alpha*zero
    for k in range(1,n):
        one[k]+=math.exp(-alpha*alpha/2)*k*(-alpha)**(k-1)/fact[k]
    return zero,one

def run():
    spec=importlib.util.spec_from_file_location("canonical965",HERE/"material_field_window.py")
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    m,Q,Hm,S,W0,local=core.material()
    saved=read(HERE/"material_field_window_results.json")
    T=saved["T"];g=.003;omega=.5
    rows=[]
    for old in saved["fixed_menu_rows"]:
        n=old["cutoff"]+1
        ann=np.diag(np.sqrt(np.arange(1,n)),1);number=np.diag(np.arange(n,dtype=float))
        H=np.kron(Hm,np.eye(n))+np.kron(np.eye(36),omega*number)
        H+=g*np.kron(S,ann+ann.T)+g*g/omega*np.kron(S@S,np.eye(n))
        val,V=np.linalg.eigh(H);U=(V*core.phases(val,T))@V.T
        phi=np.zeros(n,complex);phi[:2]=1/np.sqrt(2)
        ef=np.array([1.,1j*core.phases(np.array([omega]),T)[0]])/np.sqrt(2)
        weights=[]
        for s in np.diag(S):
            r0,r1=displacement_rows(g*s/omega,n)
            weights.append(ef[0].conjugate()*r0+ef[1].conjugate()*r1)
        weights=np.array(weights)
        # Each row is a compression of the same normalized infinite displaced
        # effect; its finite norm <=1. The omitted tail is not discarded in
        # the definition of the observable.
        row_norm_max=float(np.max(np.sum(abs(weights)**2,axis=1)))
        assert row_norm_max<=1+1e-13
        probs=[]
        for label in (0,1):
            psi=np.kron(np.kron(local[:,label],local@np.array([1.,1.])/np.sqrt(2)),phi)
            block=psi.reshape(6,4,6,4,n).transpose(0,2,4,1,3).reshape(36*n,16)
            evolved=(U@block).reshape(36,n,16)
            projected=np.einsum("in,inj->ij",weights,evolved)
            probs.append(float(np.sum(abs(projected)**2)))
        contrast=abs(probs[0]-probs[1])
        eps=old["residual_bound"]["total"]
        rows.append(dict(cutoff=old["cutoff"],field_probabilities=probs,
            contrast=contrast,infinite_model_contrast_lower=contrast-2*eps-1e-10,
            conservative_reported_lower=contrast-4e-6,
            maximum_compressed_effect_norm=row_norm_max,
            difference_from_bare_mode_contrast=contrast-old["probabilities"]["fixed_field_effect_contrast"]))
    assert rows[-1]["conservative_reported_lower"]>.04
    # Explicit independent checks of oscillator displacement coefficients.
    assert np.array_equal(displacement_rows(0.,5)[0],np.eye(5)[0])
    assert np.array_equal(displacement_rows(0.,5)[1],np.eye(5)[1])
    for alpha in (-.012,0.,.012):
        z,o=displacement_rows(alpha,30)
        assert abs(np.dot(z,z)-1)<1e-13
        assert abs(np.dot(o,o)-1)<1e-13 and abs(np.dot(z,o))<1e-13
    return dict(round=965,date="2026-10-07",all_checks_passed=True,
        definition="b=a+(g/omega)S=D((g/omega)S)^dagger a D((g/omega)S)",
        effect="D^dagger (I_material tensor |e_T><e_T|) D; e_T=(|0>+i exp(-i omega T)|1>)/sqrt(2)",
        rows=rows,
        scope=dict(effect_is_bounded_and_same_for_both_preparations=True,
            bare_length_mode_is_not_automatically_outgoing_photon=True,
            photon_detector_or_cavity_boundary_not_constructed=True,
            unitary_representation_change_is_not_full_parent_gauge_matching=True,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),HERE/"material_field_window.py",HERE/"material_field_window_results.json")})
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        assert read(TARGET)["source_hashes"]==out["source_hashes"]
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))

