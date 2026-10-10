"""1080 independent original-H/K diagnostic; default stdout JSON only.
Original 480/474 dense matrices are reused through the frozen migration runtime.
Finite differences and floating rank are diagnostics, not strict certificates.
"""
from pathlib import Path
import hashlib,json,math,sys
import numpy as np
BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from research_layout import Layout,ResearchRuntime
GROUPS=[]
def checked(name,condition):
    if not bool(condition):raise AssertionError(name)
    GROUPS.append(name)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def tensor(values):
    result=np.ones((1,1),complex)
    for a in values:result=np.kron(result,a)
    return result
def graph(columns):
    y=columns.reshape(64,6,28)
    return np.einsum("dga,dha->gh",y,y.conj())/384
def dgraph(columns,derivative):
    y=columns.reshape(64,6,28);dy=derivative.reshape(64,6,28)
    return (np.einsum("dga,dha->gh",dy,y.conj())
            +np.einsum("dga,dha->gh",y,dy.conj()))/384
def gram(a,b):return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))
def run():
    with ResearchRuntime(Layout(),"late").installed():
        import correlated_reference_local_access as source
        import relational_direction_distance_response as reader
        h=source.frame.system()[1];k=source.prior.s12();xi=source.pure_columns()
        ev,vec=np.linalg.eigh(h)
        def uh(t):return (vec*np.exp(-1j*t*ev))@vec.conj().T
        tau=1e-10;pv,ps=np.linalg.eigh(tau*h+(math.pi/2)*k)
        pulse=(ps*np.exp(-1j*pv))@ps.conj().T
        I=np.eye(2,dtype=complex);x,y,z=reader.P
        raw=np.asarray([tensor([I,s,I,I+x,I+y,I+z]) for s in (I,x,y,z)])
        m=np.kron(tensor([I,(I+z)/2,I,I,I,I]),np.eye(6))
        ur=uh(.7);o=ur.conj().T@m@ur
        B=np.einsum("diej,ked->kij",o.reshape(64,6,64,6),raw)/64
        checked("original dimensions, source normalization and reader Hermiticity",
            h.shape==(384,384) and xi.shape==(384,28)
            and abs(np.vdot(xi,xi).real-384)<1e-12
            and max(np.max(abs(a-a.conj().T)) for a in B)<1e-13)
        # Exact factorization is the all-input reason for graph invariance.
        kd=k[::6,::6]
        checked("K exactly factors on data alone and squares to identity",
            np.array_equal(k,np.kron(kd,np.eye(6)))
            and np.array_equal(kd@kd,np.eye(64)))
        def coeff(cols):return np.real(np.einsum("kij,ji->k",B,graph(cols)))
        def derivative(cols,a):
            dg=dgraph(cols,-1j*a@cols)
            return np.real(np.einsum("kij,ji->k",B,dg)),dg
        records=[];columns=[];kappas=(2/3,5/4);kangle=.37
        ku=math.cos(kangle)*np.eye(384)-1j*math.sin(kangle)*k
        worst_k=worst_kfinite=worst_hk=0.
        for t,u in ((.2,.4),(.2,.6)):
            cols=uh(u)@pulse@uh(t)@xi;columns.append(cols)
            e=coeff(cols);vh,dgh=derivative(cols,h);vk,dgk=derivative(cols,k)
            g0=graph(cols);combined=[]
            worst_k=max(worst_k,float(np.max(abs(dgk))))
            worst_kfinite=max(worst_kfinite,float(np.max(abs(graph(ku@cols)-g0))))
            for kap in kappas:
                response,dg=derivative(cols,h+kap*k)
                err=max(float(np.max(abs(response-vh))),float(np.max(abs(dg-dgh))))
                worst_hk=max(worst_hk,err)
                combined.append({"kappa":kap,"effect_derivative":response.tolist(),
                    "max_residual_from_H":err})
            records.append({"source_parameters":[t,u,math.pi/2],
                "effect_coefficients":e.tolist(),"append_H_effect_derivative":vh.tolist(),
                "append_K_effect_derivative":vk.tolist(),"H_plus_kappa_K":combined,
                "graph_trace":float(np.trace(g0).real),
                "graph_coherence_frobenius":float(np.linalg.norm(g0-np.diag(np.diag(g0))))})
        checked("K finite graph invariance and H+kappa K first derivative",
            max(worst_k,worst_kfinite,worst_hk)<2e-13)
        b1,b2=[np.array(a["effect_coefficients"][1:]) for a in records]
        v1,v2=[np.array(a["append_H_effect_derivative"][1:]) for a in records]
        r1,r2=np.dot(b1,b1),np.dot(b2,b2);dot=np.dot(b1,b2)
        N=((np.dot(v1,b2)+np.dot(b1,v2))*r1*r2
            -dot*(np.dot(b1,v1)*r2+np.dot(b2,v2)*r1))
        dgram=float(N/(r1**1.5*r2**1.5));g0=gram(b1,b2)
        matrix=np.vstack([np.dot(b,b)*np.eye(3)-np.outer(b,b) for b in (b1,b2)])
        target=np.r_[np.cross(b1,v1),np.cross(b2,v2)]
        omega,_,rank,singular=np.linalg.lstsq(matrix,target,rcond=None)
        residual=matrix@omega-target
        checked("two-source H response contradicts common rotation diagnostically",
            N < -9e-9 and dgram < -2.9e-4 and np.max(abs(residual))>3e-4 and rank==3)
        # Two positive appended drifts; no inverse-time operation is used.
        finite=[]
        for eps in (1e-3,5e-4):
            ueps=uh(eps);bb=[coeff(ueps@a)[1:] for a in columns]
            after=gram(*bb);slope=(after-g0)/eps
            finite.append({"positive_appended_time":eps,"direction_gram":after,
                "one_sided_slope":slope,"error_from_analytic_slope":abs(slope-dgram)})
        checked("two positive finite drifts converge to the analytic derivative",
            all(a["one_sided_slope"]<0 for a in finite)
            and finite[1]["error_from_analytic_slope"]<finite[0]["error_from_analytic_slope"]
            and max(a["error_from_analytic_slope"] for a in finite)<5e-6)
        # Matrix test fixes the ordering. General commutation follows algebraically.
        max_comm=max(float(np.max(abs(k@np.kron(np.eye(64),a)
                         -np.kron(np.eye(64),a)@k))) for a in B)
        checked("K commutes with all original reader graph coefficients",max_comm==0.)
        hashes={"independent_check.py":sha(Path(__file__)),
            # Runtime intentionally reports historical __file__; hash actual archive assets.
            "480_original_source":sha(ROOT/"research_cognition_physics/archive_467_530/480/correlated_reference_local_access.py"),
            "474_original_reader":sha(ROOT/"research_cognition_physics/archive_467_530/474/relational_direction_distance_response.py"),
            "1077_independent_check.py":sha(BASE.parent/"1077/independent_check.py"),
            "1077_proof.md":sha(BASE.parent/"1077/proof.md")}
        for name in ("check.py","results.json","proof.md"):
            if (BASE/name).exists():hashes["author_"+name]=sha(BASE/name)
        return {"round":1080,"status":"passed","assertion_groups":len(GROUPS),"groups":GROUPS,
            "source_tau":tau,"fixed_reader_time":.7,"sources":records,
            "normalized_direction_gram":g0,"normalized_direction_gram_derivative":dgram,
            "polynomial_N_diagnostic":float(N),
            "two_source_shared_omega_diagnostic":omega.tolist(),
            "stacked_equation_singular_values_diagnostic":singular.tolist(),
            "stacked_equation_max_residual_diagnostic":float(np.max(abs(residual))),
            "stacked_equation_residual_norm_diagnostic":float(np.linalg.norm(residual)),
            "finite_positive_drift_checks":finite,
            "K_invariance":{"finite_K_angle":kangle,"max_graph_derivative_residual":worst_k,
                "max_finite_graph_residual":worst_kfinite,
                "max_H_plus_kappa_K_first_derivative_residual":worst_hk,
                "reader_commutator_residual":max_comm},
            "hashes":hashes,
            "scope":{"same_fixed_H_K_source_and_reader":True,
                "source_dependent_positive_gain_allowed_in_tested_necessity":True,
                "source_dependent_rotation_selected":False,"all_graph_coherences_retained":True,
                "second_source_point_is_analytic_identity_witness_not_new_cube_permission":True,
                "two_fixed_points_not_control_optimization_scan":True,
                "finite_and_floating_checks_are_rigorous_witness_certificate":False,
                "K_first_order_zero_excludes_all_higher_order_control_words":False,
                "known_history_481_strategies_excluded":False,
                "macroscopic_approximate_reorientations_excluded":False,
                "cognitive_axioms_or_spatial_dimension_refuted":False}}
if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
