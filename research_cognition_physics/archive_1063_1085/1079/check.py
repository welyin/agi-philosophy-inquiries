"""1079: a fixed, source-independent physical centering of a binary reader.
Exact source bounds reuse frozen 1077/1078. Matrix checks are diagnostics.
No source-dependent axis selection, new reader Hamiltonian, or image work.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json
import numpy as np
BASE=Path(__file__).resolve().parent
I=np.eye(2,dtype=complex)
P=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def root(a):
    w,v=np.linalg.eigh(a);assert w.min()>-1e-12
    return (v*np.sqrt(np.maximum(w,0)))@v.conj().T
def measure_prepare(e):
    # The original final Z read on the tested qubit has these reduced-Q maps.
    # Remaining original network outputs are traced only for this diagnostic.
    return [[np.outer(I[:,s],root(effect)[a]) for a in range(2)]
            for s,effect in enumerate((e,I-e))]
def partial_m(rho):
    return np.einsum("iarjas->irjs",rho.reshape(2,2,2,2,2,2)).reshape(4,4)
def depolarize_first(rho):
    r=np.trace(rho.reshape(2,2,2,2),axis1=0,axis2=2)
    return np.kron(I/2,r)
def apply_q(ks,rho):
    return sum(np.kron(k,I)@rho@np.kron(k.conj().T,I) for k in ks)
def run():
    groups=[]
    def check(name,condition):
        assert bool(condition),name
        groups.append(name)
    old=load("round1078_certificate",BASE.parent/"1078/check.py").run()
    bound=F(old["uniform_spectral_derivative_boxes"][1][1][1])/4
    check("strict same-source squared-radius derivative bound",bound<F(-1,400))
    check("old Bloch embedding and new half-response remain nondegenerate",
          F(old["reused_bloch_embedding_lower"])/2==F(59,240000))
    raw=json.loads((BASE.parent/"1077/independent_results.json").read_text(encoding="utf-8"))
    c,b=np.array(raw["effect_center"])[0],np.array(raw["effect_center"])[1:]
    e=c*I+sum(x*s for x,s in zip(b,P))
    ks=measure_prepare(e)
    f=.5*e+.5*(I-np.trace(e).real/2*I)
    target=I/2+sum(x*s for x,s in zip(b,P))/2
    check("fixed coarse effect equals actual half traceless response",
          np.max(abs(f-target))<1e-14 and abs(np.trace(f)-1)<1e-14)
    check("raw reduced Z instrument complete",
          np.max(abs(sum(k.conj().T@k for ss in ks for k in ss)-I))<1e-14)
    # Keep replacement memory M and arbitrary entangled reference R explicitly.
    swap=np.zeros((8,8),complex)
    for q in range(2):
        for m in range(2):
            for r in range(2):swap[4*m+2*q+r,4*q+2*m+r]=1
    rng=np.random.default_rng(1079)
    max_coarse=max_memory=max_effect=0.
    fine_values=[]
    for case in range(8):
        v=rng.normal(size=4)+1j*rng.normal(size=4);v/=np.linalg.norm(v)
        rho=np.outer(v,v.conj())
        joint=np.einsum("irjs,ab->iarjbs",rho.reshape(2,2,2,2),I/2).reshape(8,8)
        replaced=swap@joint@swap.conj().T
        # trace Q, preserving M R; output indices a r b s.
        mr=np.trace(replaced.reshape(2,2,2,2,2,2),axis1=0,axis2=3).reshape(4,4)
        max_memory=max(max_memory,float(np.max(abs(mr-rho))))
        full=[[sum(np.kron(k,np.eye(4))@inp@np.kron(k.conj().T,np.eye(4)) for k in kk)/2
               for kk in ks] for inp in (joint,replaced)]
        check("full retained-memory branch probabilities normalized "+str(case),
              abs(sum(np.trace(a).real for rr in full for a in rr)-1)<1e-12
              and min(np.linalg.eigvalsh(a).min() for rr in full for a in rr)>-1e-12)
        coarse=partial_m(full[0][0]+full[1][1])
        expected=(apply_q(ks[0],rho)+apply_q(ks[1],depolarize_first(rho)))/2
        max_coarse=max(max_coarse,float(np.max(abs(coarse-expected))))
        probability=np.trace(np.kron(f,I)@rho).real
        max_effect=max(max_effect,abs(np.trace(coarse).real-probability))
        if case==0:fine_values=[[float(np.trace(a).real) for a in rr] for rr in full]
    check("replacement preserves old unknown Q R in retained memory",max_memory<1e-12)
    check("complete reference coarse-map identity",max_coarse<1e-12 and max_effect<1e-12)
    # Fine records keep c: branch=1, raw+= probability c/2 for any probe.
    check("fine neutral-branch record still carries raw bias",
          abs(fine_values[1][0]-c/2)<1e-12)
    # Uniform probabilistic robustness, independent of unknown source:
    # D_eps=(1-eps)D+eps id; half diamond distance to D=3eps/4.
    worst=0.
    for delta in (-.05,0.,.05):
        for eps in (0.,.04):
            weight=.5+delta
            f_err=weight*e+(1-weight)*(I-((1-eps)*c*I+eps*e))
            err=float(np.linalg.norm(f_err-f,2))
            bound_err=abs(delta)+(1-weight)*3*eps/4
            check("finite coin/depolarization effect bound "+str((delta,eps)),err<=bound_err+1e-13)
            worst=max(worst,err)
    # Bias alone demonstrates centering is not a claim about arbitrary unfair coins.
    bad=.6*e+.4*(I-c*I)
    check("unfair branch retains predicted scalar bias",
          abs(np.trace(bad).real/2-(.5+.1*(2*c-1)))<1e-14)
    delta=F(1,1000)
    factor=2*delta*(F(1,2)+delta)**2
    unfair_minor=[factor*F(v) for v in old["uniform_minor_enclosure"]]
    check("nonzero arbitrarily small bias does not preserve exact spectral rank one",
          unfair_minor[0]<=unfair_minor[1]<0)
    check("Pauli coordinate changes do not change centering rule",
          all(np.max(abs((s@f@s)-(I/2+(s@e@s-c*I)/2)))<1e-14 for s in P))
    return {"round":1079,"status":"passed","assertion_groups":len(groups),"groups":groups,
      "strict_source":{"domain":"original closed 480/1077 cube",
         "squared_radius_u_derivative_upper":str(bound),"strict_upper":"-1/400",
         "bloch_z_strict_lower":"2/25","embedding_lower":"59/240000",
         "Q_response_separation_lower":"59/9120000",
         "ambient_isospectral_dimension":2,"closed_cube_dimension_upper":2,
         "center_angular_rank":2},
      "effect_formula":"F=I/2+tf(E)/2",
      "effect_center_diagnostic":[.5,*list(b/2)],
      "effect_eigenvalues_diagnostic":np.linalg.eigvalsh(f).tolist(),
      "raw_bias_diagnostic":float(c),
      "retained_memory_residual":max_memory,"reference_map_residual":max_coarse,
      "effect_probability_residual":max_effect,"fine_record_probabilities_diagnostic":fine_values,
      "finite_error_example_max":worst,
      "exactness_boundary":{"example_bias":"1/1000",
          "unfair_invariant_minor":list(map(str,unfair_minor)),
          "unfair_exact_spectral_rank":2,
          "small_statistical_error_does_not_preserve_exact_level_dimension":True},
      "robustness":"half diamond instrument error <= eps_I+abs(delta)+(.5-delta)*eps_D, after tracing replacement memory; full lifted error needs its own certificate",
      "source_hashes":{name:hashlib.sha256((BASE.parent/name).read_bytes()).hexdigest() for name in
         ("1078/check.py","1078/results.json","1077/proof.md","1077/independent_results.json")},
      "scope":{"fixed_source_independent_preprocessing":True,
         "new_classical_axis_feedback":False,"source_parameter_estimation":False,
         "full_fine_task_isospectral_or_covariant":False,
         "probe_original_marginal_preserved":False,
         "original_H_autonomously_generates_added_controls":False,
         "actual_reorientation_generated":False,"spatial_dimension_derived":False,
         "new_accepted_cognitive_axioms":0,"empirical_results":0}}
if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
