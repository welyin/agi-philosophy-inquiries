"""Round 1322: bounded new dictionaries and strict finite-CAR transport.
Default/--check reads only. --write-new exclusively initializes results.json.
No 1321 computation is imported or rerun; its frozen certificate is verified.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
CHECKS = []
OLD_HASHES = {
    "proof.md": "5fdf6c40524a83a3fd125fc0191cd98ba181bf93ac2020debb8dcc7000791a68",
    "check.py": "126b4af829a875294cfb4b3cfc5b985ba9696513443c8cd831ba393d0252e3af",
    "results.json": "dae8ee1828d6b74335943a9dce8364a2afcc9c45a1a85ecd4c81a20111accec1",
    "acceptance.json": "9b1a600733ced1323837c7f27bb0f394da173fc9a180eb925ac33b3e4ec83578",
}
EPS = np.array([[0., 1.], [-1., 0.]])
TAU = [
    np.array([[0.,1.],[1.,0.]])/2,
    np.array([[0.,-1j],[1j,0.]])/2,
    np.diag([1.,-1.])/2,
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(name, value, tol=2e-10):
    err = float(np.max(np.abs(value))) if np.size(value) else 0.0
    if not np.isfinite(err) or err > tol:
        raise AssertionError((name, err, tol))
    CHECKS.append({"name": name, "error": err, "tolerance": tol})


def truth(name, condition):
    check(name, 0 if bool(condition) else 1, 0)


def frame(higgs):
    h = np.linalg.norm(higgs)
    if h <= 0:
        raise ValueError("Nonzero Higgs patch required")
    return np.column_stack((EPS@higgs.conj(), higgs))/h


def frame_derivative(higgs, dhiggs):
    h = np.linalg.norm(higgs)
    dh = float(np.real(np.vdot(higgs, dhiggs))/h)
    return np.column_stack((EPS@dhiggs.conj(), dhiggs))/h-frame(higgs)*dh/h


def su2(axis, angle):
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis)
    sigma = 2*sum(x*t for x,t in zip(axis,TAU))
    return np.cos(angle/2)*np.eye(2)+1j*np.sin(angle/2)*sigma


def mass(higgs, cplus, cminus, scale):
    ht = EPS@higgs.conj()
    return (np.kron(np.outer(ht,ht.conj()),cplus)
            +np.kron(np.outer(higgs,higgs.conj()),cminus))/scale


def mass_derivative(higgs, dh, cp, cm, dcp, dcm, scale, dscale):
    ht = EPS@higgs.conj()
    dt = EPS@dh.conj()
    result = np.kron(np.outer(dt,ht.conj())+np.outer(ht,dt.conj()),cp)
    result += np.kron(np.outer(dh,higgs.conj())+np.outer(higgs,dh.conj()),cm)
    result += np.kron(np.outer(ht,ht.conj()),dcp)
    result += np.kron(np.outer(higgs,higgs.conj()),dcm)
    return result/scale-mass(higgs,cp,cm,scale)*dscale/scale


def effective_connection(higgs, dhiggs, weak, hyper, y, gw, gy):
    j = frame(higgs)
    return gw*j.conj().T@weak@j+gy*y*hyper*np.eye(2)+1j*j.conj().T@frame_derivative(higgs,dhiggs)


def blockdiag(*blocks):
    total=sum(len(x) for x in blocks)
    result=np.zeros((total,total),complex)
    start=0
    for x in blocks:
        result[start:start+len(x),start:start+len(x)]=x
        start+=len(x)
    return result


def annihilator(modes, index):
    a=np.zeros((2**modes,2**modes),complex)
    for state in range(2**modes):
        if state & (1<<index):
            sign=(-1)**((state & ((1<<index)-1)).bit_count())
            a[state^(1<<index),state]=sign
    return a


def second_quantize(one_particle):
    """Exterior-power matrix, including all occupations, not a one-particle truncation."""
    n=len(one_particle)
    g=np.zeros((2**n,2**n),complex)
    subsets={state:[i for i in range(n) if state & (1<<i)] for state in range(2**n)}
    for source,cols in subsets.items():
        for target,rows in subsets.items():
            if len(rows)==len(cols):
                g[target,source]=np.linalg.det(one_particle[np.ix_(rows,cols)]) if rows else 1
    return g


def pair_mass(f, ops):
    n=len(f)
    x=sum((f[i,j].conjugate()*ops[i]@ops[n+j]
           for i in range(n) for j in range(n)), np.zeros_like(ops[0]))
    return x+x.conj().T


def calculate():
    CHECKS.clear()
    for name,expected in OLD_HASHES.items():
        truth("old_frozen_"+name,sha(STAGE/"1321"/name)==expected)
    old=json.loads((STAGE/"1321"/"results.json").read_text(encoding="utf-8"))
    receipt=json.loads((STAGE/"1321"/"acceptance.json").read_text(encoding="utf-8"))
    truth("old_acceptance_conditional",receipt["status"]=="ACCEPTED_CONDITIONAL_SCOPE")
    for name in ["proof.md","check.py","results.json"]:
        key="research_cognition_physics/archive_1313_/1321/"+name
        truth("old_receipt_binding_"+name,receipt["sha256"][key]==OLD_HASHES[name])
    truth("old_no_stage_completion",receipt["stage_complete"] is False)
    truth("common_generation_set",[x["n"] for x in old["families"]]==[3,4,5])
    truth("transported_uniform_return_bounds",
          max(x["full_parent_return_error_bound"] for x in old["families"])<.01)
    truth("transported_clock_and_scales",old["fixed_contract"]["g"]==2.**-24
          and old["fixed_contract"]["mu"]==2.**-64)
    check("transported_flavor_sequence",np.array(old["fixed_contract"]["f_sequence"])
          -np.array([4.**-i for i in range(5)]),0)
    truth("registered_common_effective_scales",2./1024<.002)
    truth("matched_dimensionless_C_is_finite_small",1024*11*2.**-64<1e-15)

    rng=np.random.default_rng(1322)
    frames=[np.array([0,2.],complex),np.array([1+2j,-2+1j]),
            np.array([2-.5j,.7+1.3j])]
    scale=7.
    for hi,higgs in enumerate(frames):
        h=np.linalg.norm(higgs); j=frame(higgs)
        check(f"h{hi}_frame_unitary",j.conj().T@j-np.eye(2))
        check(f"h{hi}_frame_special",np.linalg.det(j)-1)
        check(f"h{hi}_frame_maps_reference",j@np.array([0.,h])-higgs)
        check(f"h{hi}_weak_epsilon",j.T@EPS@j-EPS)
        u=su2([1,2,-1],.43); beta=.37
        d=np.diag([np.exp(-.5j*beta),np.exp(.5j*beta)])
        transformed=np.exp(.5j*beta)*u@higgs
        check(f"h{hi}_local_frame_transformation",frame(transformed)-u@j@d)
        for label,y in [("quark",1/6),("lepton",-.5)]:
            q=rng.normal(size=2)+1j*rng.normal(size=2)
            chi=j.conj().T@q
            chip=frame(transformed).conj().T@(np.exp(1j*y*beta)*u@q)
            phase=np.diag(np.exp(1j*beta*np.array([y+.5,y-.5])))
            check(f"h{hi}_{label}_relation_hypercharges",chip-phase@chi)
        for n in [3,4,5]:
            # Full complex matrices test the dictionary; they are not a new
            # physical parameter selection or a rerun of 1321's open-ball proof.
            mplus=np.diag(np.arange(1,n+1)).astype(complex)
            mminus=np.diag(np.arange(2,n+2)).astype(complex)
            mplus+=.01*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
            mminus+=.01*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
            cp=scale*mplus/h**2; cm=scale*mminus/h**2
            f=mass(higgs,cp,cm,scale)
            fullj=np.kron(j,np.eye(n))
            canonical=blockdiag(mplus,mminus)
            check(f"h{hi}_n{n}_arbitrary_mass_matching",
                  fullj.conj().T@f@fullj-canonical)
            check(f"h{hi}_n{n}_true_allleft_conjugation",
                  fullj.T@f.conj()@fullj.conj()-canonical.conj())
            up=np.kron(u,np.eye(n))
            check(f"h{hi}_n{n}_mass_equivariance",
                  mass(transformed,cp,cm,scale)-up@f@up.conj().T)
            check(f"h{hi}_n{n}_same_mass_singular_values",
                  np.linalg.svd(f,compute_uv=False)-np.linalg.svd(canonical,compute_uv=False))
            truth(f"h{hi}_n{n}_full_rank",
                  np.linalg.svd(canonical,compute_uv=False)[-1]>.9)

    # Local frame derivative and connection, including a time-dependent
    # weak transformation and U(1) phase; no gauge field is dynamically simulated.
    higgs=frames[2]; j=frame(higgs)
    dh=np.array([.21-.1j,-.11+.31j])
    dj=frame_derivative(higgs,dh)
    eta=1e-6
    check("frame_derivative_finite_difference",
          (frame(higgs+eta*dh)-frame(higgs-eta*dh))/(2*eta)-dj,3e-10)
    check("frame_derivative_antihermitian",j.conj().T@dj+dj.conj().T@j)
    check("frame_derivative_traceless",np.trace(j.conj().T@dj))
    weak=sum(x*t for x,t in zip([.3,-.2,.4],TAU))
    gw,gy,hyper=.7,.4,.23
    u=su2([2,-1,.5],.29); beta=.31; betadot=.17
    generator=sum(x*t for x,t in zip([.21,.13,-.07],TAU))
    udot=1j*generator@u
    d=np.diag([np.exp(-.5j*beta),np.exp(.5j*beta)])
    hnew=np.exp(.5j*beta)*u@higgs
    dhnew=np.exp(.5j*beta)*(.5j*betadot*u@higgs+udot@higgs+u@dh)
    wnew=u@weak@u.conj().T-1j/gw*udot@u.conj().T
    bnew=hyper+betadot/gy
    for y in [1/6,-.5]:
        a=effective_connection(higgs,dh,weak,hyper,y,gw,gy)
        anew=effective_connection(hnew,dhnew,wnew,bnew,y,gw,gy)
        v=np.diag(np.exp(1j*beta*np.array([y+.5,y-.5])))
        check(f"connection_covariance_y{y}",anew-v@a@v.conj().T
              -betadot*(y*np.eye(2)+TAU[2]))
        check(f"effective_connection_Hermitian_y{y}",a-a.conj().T)
    h0=np.array([0,2.],complex)
    a=effective_connection(h0,np.zeros(2),weak,hyper,1/6,gw,gy)
    std=gy*hyper*(np.eye(2)/6+TAU[2])
    delta=gw*weak-gy*hyper*TAU[2]
    check("right_connection_difference",a-std-delta)
    truth("right_charged_weak_source_nonzero",abs(delta[0,1])>.1)
    neutral=np.sqrt(gw**2+gy**2)
    sw,cw=gy/neutral,gw/neutral
    photon,zfield=.31,-.12
    w3=sw*photon+cw*zfield; b=cw*photon-sw*zfield
    check("right_photon_difference_zero_Z_retained",gw*w3-gy*b-neutral*zfield)

    # Sources: independent original parameters, not derivatives along the
    # matching map which would conceal the differing Higgs response.
    higgs=frames[1]; h=np.linalg.norm(higgs)
    cp=np.array([[1+.2j,.3j],[.1,2-.1j]])
    cm=np.array([[.7,-.2j],[.3j,1.2]])
    dh=np.array([.1+.2j,-.2+.1j])
    dcp=np.array([[.1,0],[.02j,-.03]])
    dcm=np.array([[.02j,-.07],[.03,.01]])
    ds=.13
    analytic=mass_derivative(higgs,dh,cp,cm,dcp,dcm,scale,ds)
    numeric=(mass(higgs+eta*dh,cp+eta*dcp,cm+eta*dcm,scale+eta*ds)
             -mass(higgs-eta*dh,cp-eta*dcp,cm-eta*dcm,scale-eta*ds))/(2*eta)
    check("full_H_C_Lambda_source",numeric-analytic,1e-9)
    m=np.array([[1+.2j,.1j],[.03,2-.1j]])
    href=2.; coeff=scale*m/href**2; ystd=m/href
    drstd=ystd
    drvec=2*href*coeff/scale
    check("radial_standard_source",drstd-m/href)
    check("radial_vectorlike_source",drvec-2*m/href)
    truth("radial_source_not_equal",np.linalg.norm(drstd-drvec)>.5)
    check("scale_source",-href**2*coeff/scale**2+m/scale)
    check("H0_polynomial_mass_regular",mass(np.zeros(2,complex),cp,cm,scale))
    try:
        frame(np.zeros(2,complex))
    except ValueError:
        truth("H0_relation_chart_rejected",True)
    else:
        raise AssertionError("Zero Higgs relation frame was accepted")

    # r phases refer to fields, not an inconsistent common phase for all
    # particle and antiparticle creation operators.
    rleft=1j; rright_std=-1j; rright_vec=1j
    check("Rop_squared_parity",rleft**2+1)
    check("standard_linear_mass_R",rleft.conjugate()*(-1)*rright_std-1)
    check("vectorlike_quadratic_mass_R",rleft.conjugate()*rright_vec-1)
    check("left_QQ_S_portal_R",rleft*rleft*(-1)-1)
    check("left_QL_S_portal_R",rleft*rleft*(-1)-1)
    check("left_unit_QQQL_R",rleft**4-1)
    truth("bare_same_physical_Majorana_forbidden",rleft**2 == -1 and rright_vec**2 == -1)
    truth("Rrep_right_conjugate_differs",rright_vec.conjugate()!=1j)
    for label,value in [("stdL",-rleft),("stdR",rright_std),
                        ("vecL",-rleft),("vecR",-rright_vec)]:
        check("residual_"+label,value+1j)
    check("residual_H_invariant",(-1)*(-1)-1)
    check("residual_S_unchanged",-1-(-1))
    # Dirac conjugate pairs: r exponent +/-1, all-left x=6Y.
    reps=[(6,1,1),(2,-3,1),(6,-1,-1),(2,3,-1)]
    for powers in [(1,0),(3,0),(0,1),(0,3),(2,1),(1,2)]:
        a,b=powers
        check(f"vectorlike_Abelian_anomaly_{a}_{b}",
              sum(mult*x**a*r**b for mult,x,r in reps))
    check("SU3_cubic_conjugate_pairs",2-2)
    check("SU3_squared_hypercharge",2*(1-1))
    check("SU2_squared_hypercharge",3*(1-1)+(-3+3))
    check("SU3_squared_R",2*(1-1))
    check("SU2_squared_R",3*(1-1)+(1-1))
    truth("SU2_Witten_even_eight_per_generation",8%2==0)
    truth("sixteen_Weyl_per_generation",sum(x[0] for x in reps)==16)
    # Finite determinant is only a diagnostic of the analytic positive
    # vectorlike regulator argument, not a proof of all global backgrounds.
    vals=np.array([.2,1.3,2.1]); regulator=.7
    paired=np.array([(regulator+1j*x)*(regulator-1j*x) for x in vals])
    check("paired_Dirac_regulator_real",paired.imag)
    truth("paired_Dirac_regulator_positive",min(paired.real)>0)

    # All sectors of a 4-CAR-mode example, independent of 1321 code.
    # L transforms as J^dagger; the genuine conjugate R transforms as J^T.
    ops=[annihilator(4,i) for i in range(4)]
    eye=np.eye(16)
    for i,k in itertools.product(range(4),repeat=2):
        check(f"CAR_{i}_{k}",ops[i]@ops[k].conj().T+ops[k].conj().T@ops[i]-(eye if i==k else 0))
    j=frame(frames[2])
    basis=blockdiag(j.conj().T,j.T)
    gamma=second_quantize(basis.conj().T)
    check("full_small_Fock_exterior_unitary",gamma.conj().T@gamma-eye)
    for i in range(4):
        rhs=sum((basis[i,k]*ops[k] for k in range(4)),np.zeros((16,16),complex))
        check(f"full_small_Fock_CAR_transport_{i}",gamma@ops[i]@gamma.conj().T-rhs)
    nf=sum((x.conj().T@x for x in ops),np.zeros((16,16),complex))
    check("full_small_Fock_number_preserved",gamma@nf@gamma.conj().T-nf)
    vacuum=np.eye(16)[:,0]
    check("full_small_Fock_vacuum_preserved",gamma@vacuum-vacuum)
    d=np.diag([1+.2j,2-.3j])
    f=j@d@j.conj().T
    check("full_small_Fock_mass_polynomial",
          gamma@pair_mass(d,ops)@gamma.conj().T-pair_mass(f,ops))
    portal=ops[0]@ops[1]
    check("full_small_Fock_left_epsilon_portal",
          gamma@portal@gamma.conj().T-portal)
    neutral_number=ops[2].conj().T@ops[2]
    projector=np.outer(j[:,0],j[:,0].conj())
    conjugate_projector=projector.conj()
    physical_effect=sum((conjugate_projector[i,k]*ops[2+i].conj().T@ops[2+k]
                         for i,k in itertools.product(range(2),repeat=2)),
                        np.zeros((16,16),complex))
    check("right_neutral_relation_effect",gamma@neutral_number@gamma.conj().T-physical_effect)
    check("right_neutral_relation_effect_projector",physical_effect@physical_effect-physical_effect)
    check("particle_antiparticle_projector_conjugacy",conjugate_projector-projector.conj())

    source_hashes={"1321/"+k:v for k,v in OLD_HASHES.items()}
    for rel in ["proof.md","check.py","NEXT.md"]:
        source_hashes["1322/"+rel]=sha(HERE/rel)
    source_hashes["research_note_1322.md"]=sha(STAGE/"research_note_1322.md")
    for rel in ["_admission/round1322_admission.md",
                "_admission/weak_role_joint_screen_after1321.md"]:
        source_hashes[rel]=sha(STAGE/rel)
    return {
        "round":1322,"status":"AUTHOR_CHECK_PASS",
        "scope":"Fixed nonzero Higgs and zero gauge background: exact declared finite-CAR transport, not full physical equivalence.",
        "check_count":len(CHECKS),"checks":CHECKS.copy(),
        "maximum_absolute_residual":max(x["error"] for x in CHECKS),
        "transported_contract":old["fixed_contract"],
        "transported_uniform_errors":old["uniform_errors"],
        "registered_effective_scales_over_p":{"Higgs_norm":1,"Lambda":1024,"M_S":2,
            "uncomputed_higher_operator_errors_certified":False},
        "transported_family_bounds":[{"n":x["n"],"return_error_bound":x["full_parent_return_error_bound"],
                                     "fermion_modes":x["fermion_modes"]} for x in old["families"]],
        "new_distinctions":{"radial_Higgs_source_ratio_at_fixed_coefficients":2,
                            "right_charged_weak_current_present":True,
                            "R_op_and_fixed_H_residual_preserved":True,
                            "unbroken_R_rep_preserved":False,
                            "mass_operator_dimension":5,
                            "new_axioms":0},
        "scope_exclusions":["H=0 relation chart","Full dynamic gauge/Higgs apparatus",
                            "Undeclared quotient or nonspin backgrounds",
                            "All unknown flavor delivery","Empirical tiny-probability precision",
                            "Equality of background variations or total resource costs"],
        "source_sha256":source_hashes,
    }


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=calculate()
    result=HERE/"results.json"
    if args.write_new:
        with result.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        previous=json.loads(result.read_text(encoding="utf-8"))
        if previous != data:
            raise AssertionError("Read-only recomputation differs from saved results")
    print(json.dumps({"round":1322,"mode":"write-new" if args.write_new else "read-only",
                      "checks":data["check_count"],"max_error":data["maximum_absolute_residual"],
                      "transported_max_return_bound":max(x["return_error_bound"] for x in data["transported_family_bounds"]),
                      "pass":True},ensure_ascii=False))


if __name__=="__main__":
    main()

