"""966: finite coherent comparison mechanism; no fitted plaquette term.
Default: recompute and compare saved evidence. --write: exclusive creation.
All matrices use the exact +1 Gauss sector in the X basis (64 states).
"""
from pathlib import Path
import argparse, hashlib, json, math
from fractions import Fraction
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=STAGE.parent.parent
OUT=HERE/"comparison_loop_results.json"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def norm(a):
    return float(np.linalg.norm(a,2))

def permutation(mask):
    a=np.zeros((64,64))
    for k in range(64): a[k^mask,k]=1
    return a

def exact_upper(a):
    # All entries here are dyadic rationals from integer matrices; every
    # arithmetic operation is exact in binary64 at these small dimensions.
    rows=max(sum(abs(Fraction(float(x))) for x in row) for row in a)
    cols=max(sum(abs(Fraction(float(x))) for x in col) for col in a.T)
    # max(||A||_1, ||A||_inf) bounds sqrt(their product) and hence ||A||_2.
    upper=max(rows,cols)
    return upper, {"one":str(cols),"infinity":str(rows),"upper":str(upper)}

def build():
    # k = 8*data + link. data bits are each node's incoming-port X bit.
    # Link bits: 01,12,20. Outgoing-port X bit is fixed by Gauss=+1.
    edges=((0,1),(1,2),(2,0))
    I=np.eye(64)
    Xlink=[np.diag([1-2*((k>>e)&1) for k in range(64)]) for e in range(3)]
    charges=[]
    for i in range(3):
        incident=[e for e,(u,v) in enumerate(edges) if i in (u,v)]
        charges.append(np.diag([
            ((k>>incident[0])^(k>>incident[1]))&1 for k in range(64)]))
    H0=sum(charges)
    T=[permutation((1<<e)|(1<<(3+v))) for e,(u,v) in enumerate(edges)]
    V=-sum(T)
    D=sum((I-x)/4 for x in Xlink)  # mu=epsilon^3/2, Hrel=mu sum n_e
    B=permutation(7)
    record=[permutation(1<<(3+i)) for i in range(3)]
    S=record[0]@record[1]@record[2]
    cols=[8*d+l for d in range(8) for l in (0,7)]
    P=I[:,cols]
    Q=I-P@P.T
    R=np.diag([0 if H0[k,k]==0 else 1/H0[k,k] for k in range(64)])
    K2=P.T@V@(-R@V@P)
    K3=P.T@V@R@V@R@V@P + P.T@D@P
    W1=-R@V@P
    W2=-R@V@W1
    W3=-R@(V@W2-W1@K2)
    A4=V@W3+D@W1-W2@K2-W1@K3
    A5=D@W2-W3@K2-W2@K3
    A6=D@W3-W3@K3
    return locals()

def compute():
    m=build()
    I,H0,V,D,P,R,B,S,Q=[m[k] for k in ("I","H0","V","D","P","R","B","S","Q")]
    K2,K3=m["K2"],m["K3"]
    W1,W2,W3=[m["W"+str(i)] for i in (1,2,3)]
    A4,A5,A6=[m["A"+str(i)] for i in (4,5,6)]
    assert P.shape==(64,16)
    assert np.array_equal(P.T@V@P,np.zeros((16,16)))
    assert np.array_equal(K2,-1.5*np.eye(16))
    assert np.array_equal(K3,P.T@D@P-1.5*P.T@S@B@P)
    assert np.array_equal(H0@W1+V@P,np.zeros_like(P))
    assert np.array_equal(H0@W2+V@W1-P@K2,np.zeros_like(P))
    assert np.array_equal(H0@W3+V@W2+D@P-W1@K2-P@K3,np.zeros_like(P))
    # Independent word enumeration: two distinct edges cannot return;
    # precisely the six permutations of three different edges survive.
    import itertools
    word_rows=[]
    for word in itertools.product(range(3),repeat=3):
        a=P.T@(-m["T"][word[0]])@R@(-m["T"][word[1]])@R@(-m["T"][word[2]])@P
        nonzero=bool(np.any(a))
        assert nonzero==(len(set(word))==3)
        if nonzero:
            assert np.array_equal(a,-.25*P.T@S@B@P)
        word_rows.append(dict(word=list(word),survives=nonzero))
    Vtree=-(m["T"][0]+m["T"][1])
    assert np.array_equal(P.T@Vtree@R@Vtree@R@Vtree@P,np.zeros((16,16)))
    # Full-basis reconstruction verifies the reduced sector is a physical
    # Gauss=+1 sector, not a reduced matrix picked to fit the answer.
    full_indices=[]
    for k in range(64):
        link=k&7; data=k>>3
        ports=[]
        for i in range(3):
            q=int(m["charges"][i][k,k])
            incoming=(data>>i)&1
            ports.extend((incoming^q,incoming))
        full=sum(bit<<j for j,bit in enumerate(ports)) | (link<<6)
        full_indices.append(full)
    lookup={f:k for k,f in enumerate(full_indices)}
    assert len(lookup)==64
    for k,full in enumerate(full_indices):
        for e,(u,v) in enumerate(m["edges"]):
            dest=full^(1<<(2*u))^(1<<(2*v+1))^(1<<(6+e))
            assert m["T"][e][lookup[dest],k]==1
    # Decisive representation screen: a fixed joint permutation removes
    # the conserved record signs from H, but moves them into the effects.
    F=np.zeros((64,64))
    for k in range(64):
        target=k
        for e,(u,v) in enumerate(m["edges"]):
            if (k>>e)&1:target^=1<<(3+v)
        F[target,k]=1
    assert np.array_equal(F@F,I)
    for e in range(3):
        assert np.array_equal(F@m["T"][e]@F,permutation(1<<e))
    assert np.array_equal(F@H0@F,H0)
    assert np.array_equal(F@D@F,D)
    assert np.array_equal(F@B@F,S@B)
    assert all(np.array_equal(F@a@F,a) for a in m["record"])
    exact_bounds={}
    ub={}
    for key in ("W1","W2","W3","A4","A5","A6"):
        ub[key],exact_bounds[key]=exact_upper(m[key])
    rows=[]
    for epsilon in (.01,.005,.0025):
        e=Fraction(str(epsilon))
        eta=2*sum(e**k*ub["W"+str(k)] for k in (1,2,3))
        eta+=e*ub["A4"]+e**2*ub["A5"]+e**3*ub["A6"]
        H=H0+epsilon*(3*I+V)+epsilon**3*D
        dressed=F@H@F
        assert norm(dressed-np.kron(np.eye(8),dressed[:8,:8]))<1e-14
        K=3*epsilon*np.eye(16)+epsilon**2*K2+epsilon**3*K3
        W=P+epsilon*W1+epsilon**2*W2+epsilon**3*W3
        residual=H@W-W@K
        expansion=epsilon**4*A4+epsilon**5*A5+epsilon**6*A6
        assert norm(residual-expansion)<1e-14
        for ri in m["record"]:
            assert norm(H@ri-ri@H)<1e-14
        Omega=3*math.sqrt(5)/4
        time=math.pi/(2*Omega*epsilon**3)
        assert time<epsilon**-3
        evals,U=np.linalg.eigh(H)
        ek,Uk=np.linalg.eigh(K)
        # Exact spectral propagation numerically; analytic bounds above
        # establish the theorem independently of long-time phase numerics.
        phases=np.exp(-1j*np.remainder(evals*time,2*math.pi))
        Ufull=(U*phases)@U.T
        Ueff=(Uk*np.exp(-1j*np.remainder(ek*time,2*math.pi)))@Uk.T
        operator_error=norm(Ufull@P-P@Ueff)
        assert operator_error<float(eta)
        effect=(I+B)/2
        signal=[]
        energies=[]
        negative_flux=[]
        parts=(H0,epsilon*(3*I+V),epsilon**3*D)
        for sign in (1,-1):
            psi=np.zeros(64)
            for data in range(8):
                psi[8*data]=sign**(data&1)/math.sqrt(8)
            out=Ufull@psi
            probability=float(np.vdot(out,effect@out).real)
            negative_flux.append(float(np.vdot(out,(I-m["Xlink"][0])@out).real/2))
            old_record=float(np.vdot(out,m["record"][0]@out).real)
            leakage=float(np.linalg.norm(Q@out)**2)
            signal.append(dict(record_sign=sign,loop_probability=probability,
                               retained_record=old_record,organization_leakage=leakage))
            before=[float(np.vdot(psi,A@psi).real) for A in parts]
            after=[float(np.vdot(out,A@out).real) for A in parts]
            assert abs(sum(after)-sum(before))<1e-12
            energies.append(dict(record_sign=sign,before=before,after=after,
                                 changes=[b-a for a,b in zip(before,after)]))
        contrast=signal[0]["loop_probability"]-signal[1]["loop_probability"]
        lower=Fraction(4,5)-4*eta
        assert contrast>=float(lower)
        # Equal mixture of both conserved records has <S>=0. Replacing S
        # by that mean freezes initial electric flux; the actual mixture
        # still flips it with leading probability 4/5.
        mixed_flux=sum(negative_flux)/2
        mixed_lower=Fraction(4,5)-2*eta
        assert mixed_flux>=float(mixed_lower)
        source=np.diag([sum(m["charges"][i][k,k] for i in range(3)) for k in range(64)])
        comparison_source=3*I+V+3*epsilon**2*D  # total epsilon derivative
        currents=[1j*(H@A-A@H) for A in parts]
        assert norm(sum(currents))<1e-14
        # Bare electric flux and closed comparison both actually respond.
        flux_current=norm(1j*(H@m["Xlink"][0]-m["Xlink"][0]@H))
        loop_current=norm(1j*(H@B-B@H))
        assert flux_current>0 and loop_current>0
        rows.append(dict(epsilon=epsilon,Delta=1.,mu=epsilon**3/2,
             generated_beta=1.5*epsilon**3,time=time,
             all_input_uniform_error_bound=str(eta),
             all_input_uniform_error_bound_float=float(eta),
             certified_loop_contrast_lower=str(lower),
             certified_loop_contrast_lower_float=float(lower),
             predicted_loop_contrast=.8,actual_loop_contrast=contrast,
             mixed_record_negative_flux_probability=mixed_flux,
             mean_record_replacement_negative_flux_probability=0.,
             certified_mean_replacement_discrepancy_lower=str(mixed_lower),
             direct_uniform_operator_error=operator_error,
             signals=signal,energy_parts=energies,
             full_energy_current_residual=norm(sum(currents)),
             electric_flux_current_norm=flux_current,
             loop_comparison_current_norm=loop_current,
             comparison_parameter_source_norm=norm(comparison_source),
             gap_parameter_source_norm=norm(source)))
    assert rows[1]["certified_loop_contrast_lower_float"]>0
    # Record coupling cannot be replaced by a scalar beta for arbitrary input.
    loop_part=K3-P.T@D@P
    scalar_guess=-1.5*P.T@B@P
    wrong_norm=norm(loop_part-scalar_guess)
    assert abs(wrong_norm-3)<1e-14
    source_files=[ROOT/"research_cognition_physics/archive_429_466/research_note_435.md",
        ROOT/"research_cognition_physics/archive_429_466/research_note_438.md",
        ROOT/"research_cognition_physics/archive_429_466/research_note_449.md",
        ROOT/"research_cognition_physics/archive_429_466/research_note_458.md",
        ROOT/"research_cognition_physics/archive_554_584/research_note_568.md",
        STAGE/"research_note_961.md",STAGE/"research_note_962.md",
        HERE/"drafts/mechanism_priority_entry.md"]
    return dict(round=966,all_scientific_checks_passed=True,
        physical_sector_dimension=64,low_organization_dimension=16,
        record_qubits=3,relation_logical_qubits=1,
        exact_third_order_surviving_words=word_rows,
        coefficients=dict(second_order_constant="-3/2",record_loop="-3/2"),
        dyadic_operator_bounds=exact_bounds,fixed_menu_rows=rows,
        wrong_scalar_loop_generator_norm=wrong_norm,
        fixed_joint_dictionary_factorizes_H_exactly=True,
        same_dictionary_makes_original_loop_effect_record_dependent=True,
        content_dependent_spectrum_generated=False,
        tree_third_order_zero=True,gauss_sector_invariant=True,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_files},
        scope=dict(cognitive_working_candidate=True,
            inherited_exchange_material_435_reused_as_same_H=False,
            independent_pure_loop_term_derived_for_all_records=False,
            local_comparison_generated_record_loop_coupling=True,
            finite_time_all_logical_inputs_and_passive_references=True,
            dimension_group_energy_inputs_removed=False,
            actual_spacetime_propagation_proved=False,
            SM_or_GR_derived=False,full_goal_complete=False,
            visual_checks_performed=False))

def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{i}")
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-6,abs_tol=2e-7),(path,a,b)
    else:assert a==b,(path,a,b)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.write:
        with OUT.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(OUT.read_text("utf-8")))
    print(json.dumps(dict(round=966,checked=True,mode="write" if args.write else "read-only",
        bounds=[dict(e=r["epsilon"],eta=r["all_input_uniform_error_bound_float"],
                     lower=r["certified_loop_contrast_lower_float"],
                     actual=r["actual_loop_contrast"],error=r["direct_uniform_operator_error"])
                for r in result["fixed_menu_rows"]],
        exact_bounds=result["dyadic_operator_bounds"]),ensure_ascii=False))

