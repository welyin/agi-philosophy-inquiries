"""Round 1073 independent finite-matrix audit (NumPy only).

No imports from author code. Default behavior: print JSON to stdout, never write.
Universal statements are proved in independent_review.md; finite matrices audit
the operator contracts and provide an independent Choi-order upper certificate.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np

TOL = 2.0e-11
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = (X, Y, Z)
rng = np.random.default_rng(107310)
checks = []
max_equality_residual = 0.0

def herm(a):
    return (a + a.conj().T) / 2

def op(a):
    return float(np.linalg.norm(a, 2))

def tn(a):
    return float(np.sum(np.abs(np.linalg.eigvalsh(herm(a)))))

def tr(a):
    return float(np.trace(a).real)

def psqrt(a):
    vals, vecs = np.linalg.eigh(herm(a))
    return (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.conj().T

def equal(name, a, b, tol=TOL):
    global max_equality_residual
    res = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    max_equality_residual = max(max_equality_residual, res)
    checks.append((name, res <= tol))
    if res > tol:
        raise AssertionError((name, res, tol))
    return res

def require(name, condition):
    checks.append((name, bool(condition)))
    if not condition:
        raise AssertionError(name)

def leq(name, a, b, tol=TOL):
    v = float(np.max(np.linalg.eigvalsh(herm(a - b))))
    require(name, v <= tol)
    return max(0.0, v)

def effect(kraus):
    return sum((k.conj().T @ k for k in kraus),
               np.zeros((kraus[0].shape[1],) * 2, dtype=complex))

def act(kraus, rho):
    return sum((k @ rho @ k.conj().T for k in kraus),
               np.zeros((kraus[0].shape[0],) * 2, dtype=complex))

def dual(kraus, a):
    return sum((k.conj().T @ a @ k for k in kraus),
               np.zeros((kraus[0].shape[1],) * 2, dtype=complex))

def tagged(branches):
    result = []
    for s, branch in enumerate(branches):
        for k in branch:
            out = np.zeros((4, 2), dtype=complex)
            out[2*s:2*s+2, :] = k
            result.append(out)
    return result

def choi(kraus):
    # Convention: input tensor output; UNNORMALIZED maximally entangled vector.
    dim = kraus[0].size
    out = np.zeros((dim, dim), dtype=complex)
    for k in kraus:
        vec = k.T.reshape(-1)
        out += np.outer(vec, vec.conj())
    return herm(out)

def choi_cb_upper(k1, k2):
    # For Hermiticity-preserving Delta, ||Delta||_diamond <=
    # ||Tr_output |J(Delta)| ||. This follows from the positive/negative
    # Choi decomposition and CP-order domination; no numerical SDP is used.
    j = choi(k1) - choi(k2)
    vals, vecs = np.linalg.eigh(herm(j))
    ja = (vecs * np.abs(vals)) @ vecs.conj().T
    ni = k1[0].shape[1]
    no = k1[0].shape[0]
    red = np.trace(ja.reshape(ni, no, ni, no), axis1=1, axis2=3)
    return 0.5 * op(red)

def projector(n):
    n = np.asarray(n, dtype=float)
    n /= np.linalg.norm(n)
    return herm((I2 + sum(x*s for x,s in zip(n, PAULI))) / 2)

def unitary(n, angle):
    n = np.asarray(n, dtype=float)
    n /= np.linalg.norm(n)
    return math.cos(angle/2)*I2 - 1j*math.sin(angle/2)*sum(x*s for x,s in zip(n,PAULI))

def mp_branch(e, p):
    pv, pu = np.linalg.eigh(herm(p))
    ket = pu[:, int(np.argmax(pv))]
    ev, eu = np.linalg.eigh(herm(e))
    return [math.sqrt(max(float(ev[j]), 0))*np.outer(ket, eu[:,j].conj())
            for j in range(2)]

def mp_instrument(e, p):
    return [mp_branch(e,p), mp_branch(I2-e,I2-p)]

def lueders(p):
    return [[p], [I2-p]]

def tf(e):
    return e - tr(e)*I2/2

def rotated_instrument(eps, n, angle):
    p = projector(n)
    e0 = eps / 4
    e = (1-e0)*p + e0*(I2-p)
    tmax = 0 if eps == 0 else (eps-e0)/(1-2*e0)
    t = 0.9*tmax
    es = [e, I2-e]
    ps = [p, I2-p]
    us = [unitary([1,2,-1], angle), unitary([-2,1,1], -0.7*angle)]
    branches = []
    for s in range(2):
        ks = [math.sqrt(1-t)*k for k in mp_branch(es[s], ps[s])]
        ks.append(math.sqrt(t)*us[s] @ psqrt(es[s]))
        branches.append(ks)
    return branches, e, p, t

# Reference states: entangled states, mixtures, and independently drawn PSD states.
basis = np.eye(4, dtype=complex)
refs = []
for v in [(basis[:,0]+basis[:,3])/math.sqrt(2),
          (basis[:,0]-basis[:,3])/math.sqrt(2),
          (basis[:,1]+1j*basis[:,2])/math.sqrt(2),
          basis[:,0]]:
    refs.append(np.outer(v,v.conj()))
refs.append(np.eye(4, dtype=complex)/4)
for _ in range(5):
    a = rng.normal(size=(4,4)) + 1j*rng.normal(size=(4,4))
    rr = a @ a.conj().T
    refs.append(rr/tr(rr))

families = []
max_reference_gap = 0.0
min_positive_non_mp_witness = 1.0
for eps in [0.0, 0.0001, 0.01, 0.04, 0.10, 0.25, 0.49]:
    for idx, n in enumerate([[0,0,1],[1,2,3],[-2,1,1]]):
        branches,e,p,t = rotated_instrument(eps,n,0.4+idx*0.6)
        name = f"instrument[{eps},{idx}]"
        ps = [p,I2-p]
        es = [e,I2-e]
        equal(name+"/effect",effect(branches[0]),e)
        equal(name+"/total",effect(branches[0])+effect(branches[1]),I2)
        leq(name+"/positive effect",np.zeros((2,2)),e)
        leq(name+"/bounded effect",e,I2)
        for s in range(2):
            # This operator inequality certifies every input and conditioning
            # branch, including reduced states of arbitrary reference inputs.
            wrong_effect = I2-e if s == 0 else e
            leq(name+f"/uniform_recheck_{s}",
                dual(branches[s],wrong_effect),eps*es[s])
        ev = np.linalg.eigvalsh(e)
        require(name+"/spectral_low",ev[0] <= eps+TOL)
        require(name+"/spectral_high",ev[1] >= 1-eps-TOL)
        require(name+"/trace_window",abs(tr(e)-1) <= eps+TOL)
        q = eps/(1-eps)
        b = min(1.0, math.sqrt(q)+q/2+eps)
        leakage = dual(branches[0],I2-p)+dual(branches[1],p)
        leq(name+"/total_leak",leakage,q*I2)
        full = tagged(branches)
        compressed = tagged([[ps[s]@k for k in branches[s]] for s in range(2)])
        mp = tagged(mp_instrument(e,p))
        ideal = tagged(lueders(p))
        jdelta = choi(mp)-choi(compressed)
        leq(name+"/mp_minus_compression_CP",np.zeros_like(jdelta),jdelta)
        dcp = effect(mp)-effect(compressed)
        equal(name+"/CP_difference_effect",dcp,leakage)
        require(name+"/CP_diamond_bound",op(dcp)<=q+TOL)
        require(name+"/effect_to_projection",op(e-p)<=eps+TOL)
        mp_upper = choi_cb_upper(mp,ideal)
        require(name+"/mp_halfdiamond_bound",mp_upper<=eps+TOL)
        actual_upper = choi_cb_upper(full,ideal)
        require(name+"/independent_Choi_upper",actual_upper<=b+TOL)
        sample_max = 0.0
        for ir,rho in enumerate(refs):
            ff = act([np.kron(k,I2) for k in full],rho)
            cc = act([np.kron(k,I2) for k in compressed],rho)
            mm = act([np.kron(k,I2) for k in mp],rho)
            ii = act([np.kron(k,I2) for k in ideal],rho)
            d = 0.5*tn(ff-ii)
            sample_max = max(sample_max,d)
            require(name+f"/reference_trace_{ir}",abs(tr(ff)-1)<=TOL)
            require(name+f"/gentle_reference_{ir}",tn(ff-cc)<=2*math.sqrt(q)+TOL)
            require(name+f"/CP_reference_{ir}",tn(mm-cc)<=q+TOL)
            require(name+f"/full_reference_{ir}",d<=b+TOL)
            require(name+f"/Choi_certifies_reference_{ir}",d<=actual_upper+TOL)
        max_reference_gap = max(max_reference_gap,sample_max)
        # For eps>0 the coherent Kraus term really retains off-diagonal input:
        # unlike a measure-and-prepare branch, it acts nontrivially on |p><p_perp|.
        pv,pu = np.linalg.eigh(p)
        cross = np.outer(pu[:,1],pu[:,0].conj())
        witness = op(act(branches[0],cross))
        if eps>0:
            require(name+"/not_measure_prepare",witness>1e-9)
            min_positive_non_mp_witness=min(min_positive_non_mp_witness,witness)
        else:
            equal(name+"/zero_unique_Lueders",choi(full),choi(ideal))
        families.append({
            "epsilon":eps,"orientation":idx,"coherent_weight":t,
            "universal_halfdiamond_bound":b,
            "independent_choi_halfdiamond_upper":actual_upper,
            "reference_sample_max":sample_max,
            "worst_leakage":op(leakage),
            "offdiagonal_coherent_response":witness
        })

# Sharp sqrt(eps) obstruction: exact effect P, rotate the actual posterior.
lower_cases=[]
p=projector([0,0,1])
for eps in [0,0.0001,0.01,0.09,0.25,0.49]:
    theta=math.asin(math.sqrt(eps))
    v=np.array([[math.cos(theta),-math.sin(theta)],
                [math.sin(theta),math.cos(theta)]],dtype=complex)
    inst=[[v@p],[v@(I2-p)]]
    full=tagged(inst)
    ideal=tagged(lueders(p))
    equal(f"sqrt[{eps}]/same_effect",effect(inst[0]),p)
    equal(f"sqrt[{eps}]/branch_plus",dual(inst[0],I2-p),eps*p)
    equal(f"sqrt[{eps}]/branch_minus",dual(inst[1],p),eps*(I2-p))
    cb=choi_cb_upper(full,ideal)
    equal(f"sqrt[{eps}]/Choi_upper_exact",cb,math.sqrt(eps))
    for ir,rho in enumerate(refs):
        d=0.5*tn(act([np.kron(k,I2) for k in full],rho)
                 -act([np.kron(k,I2) for k in ideal],rho))
        equal(f"sqrt[{eps}]/reference_exact_{ir}",d,math.sqrt(eps))
    lower_cases.append({"epsilon":eps,"halfdiamond_exact":math.sqrt(eps),
                        "independent_choi_upper":cb})

# Strict reverse-bias threshold. The equality example saturates exactly and
# demonstrates why >= epsilon is insufficient.
threshold_cases=[]
for eps in [0,0.001,0.03,0.2,0.49]:
    p=projector([1,-1,2])
    e=(1-eps)*p
    f=p+eps*(I2-p)
    for label,g in [("E",e),("F",f)]:
        inst=mp_instrument(g,p)
        for s in range(2):
            wrong=I2-g if s==0 else g
            es=g if s==0 else I2-g
            leq(f"threshold[{eps}]/{label}/repeat{s}",
                dual(inst[s],wrong),eps*es)
    equal(f"threshold[{eps}]/scalar_gap",f-e,eps*I2)
    equal(f"threshold[{eps}]/tf_collision",tf(e),tf(f))
    equal(f"threshold[{eps}]/gap",op(e-f),eps)
    fstrict=(1-eps)*(I2-p)
    equal(f"threshold[{eps}]/strict_gap",op(e-fstrict),1-eps)
    require(f"threshold[{eps}]/strictly_above",op(e-fstrict)>eps)
    require(f"threshold[{eps}]/strict_tf_distinct",op(tf(e)-tf(fstrict))>0)
    # Independently generated allowed spectral windows, arbitrary orientations.
    for j in range(12):
        n=rng.normal(size=3)
        pp=projector(n)
        aa=eps*rng.random()
        bb=1-eps*rng.random()
        g=aa*(I2-pp)+bb*pp
        scalar=abs(tr(e-g))/2
        require(f"threshold[{eps}]/scalar_window{j}",scalar<=eps+TOL)
        require(f"threshold[{eps}]/quantitative_tf{j}",
                op(tf(e)-tf(g))>=op(e-g)-eps-TOL)
    threshold_cases.append({"epsilon":eps,"bias_only_gap":op(e-f),
                            "bias_only_tf_gap":op(tf(e)-tf(f)),
                            "strict_gap":op(e-fstrict)})

# Same adaptive feedback on both processes, complete classical history retained.
def adaptive_output(rho,eps,steps,ideal=False):
    blocks={():rho}
    for step in range(steps):
        new={}
        for history,state in blocks.items():
            n=[1,2,3] if sum(history)%2 else [0,0,1]
            inst,e,p,t=rotated_instrument(eps,n,0.9)
            if ideal:
                inst=lueders(p)
            for s in range(2):
                branch=act([np.kron(k,I2) for k in inst[s]],state)
                fb=unitary([1,1,0],0.2*(1+s+sum(history)))
                kr=np.kron(fb,I2)
                new[history+(s,)]=kr@branch@kr.conj().T
        blocks=new
    return blocks

adaptive_cases=[]
for eps in [0.0,0.0001,0.01,0.1]:
    q=eps/(1-eps)
    b=min(1,math.sqrt(q)+q/2+eps)
    for steps in [1,2,4]:
        real=adaptive_output(refs[0],eps,steps)
        ideal=adaptive_output(refs[0],eps,steps,True)
        d=0.5*sum(tn(real[h]-ideal[h]) for h in real)
        equal(f"adaptive[{eps},{steps}]/normalization",
              sum(tr(x) for x in real.values()),1)
        require(f"adaptive[{eps},{steps}]/budget",d<=min(1,steps*b)+TOL)
        adaptive_cases.append({"epsilon":eps,"steps":steps,
                               "trace_distance":d,"budget":min(1,steps*b)})

# Hidden writable memory can fake perfect "repeat". The SAME true instrument
# acts on D tensor a three-state memory; its first-use reduction on D is p I.
# It is not a closed qubit instrument on subsequent calls.
coin=0.3
m=np.eye(3,dtype=complex)
ksplus=[math.sqrt(coin)*np.kron(I2,np.outer(m[:,1],m[:,0])),
        np.kron(I2,np.outer(m[:,1],m[:,1]))]
ksminus=[math.sqrt(1-coin)*np.kron(I2,np.outer(m[:,2],m[:,0])),
         np.kron(I2,np.outer(m[:,2],m[:,2]))]
equal("memory/TP",effect(ksplus)+effect(ksminus),np.eye(6))
rho=np.kron(I2/2,np.outer(m[:,0],m[:,0]))
outplus=act(ksplus,rho)
outminus=act(ksminus,rho)
equal("memory/first_plus",tr(outplus),coin)
equal("memory/plus_repeats",tr(act(ksplus,outplus)),tr(outplus))
equal("memory/minus_repeats",tr(act(ksminus,outminus)),tr(outminus))
embedding=np.kron(I2,m[:,[0]])
ereduced=embedding.conj().T@effect(ksplus)@embedding
equal("memory/reduced_effect",ereduced,coin*I2)
require("memory/not_closed_qubit",abs(tr(ereduced)-1)>0.1)
equal("memory/reset_then_plus",tr(outplus)*coin,coin**2)

summary={
    "round":1073,
    "independent":True,
    "author_code_imported":False,
    "checks":len(checks),
    "all_passed":all(ok for _,ok in checks),
    "max_equality_residual":max_equality_residual,
    "tolerance":TOL,
    "instrument_families":families,
    "sqrt_obstruction":lower_cases,
    "reverse_threshold":threshold_cases,
    "adaptive_examples":adaptive_cases,
    "memory_counterexample":{
        "first_qubit_effect_eigenvalues":[coin,coin],
        "actual_instrument_input_dimension":6,
        "retained_memory_conditional_repeat_error":0.0,
        "reset_memory_plus_repeat_error":1-coin,
        "reason":"Retained memory belongs to the recurrent task; deleting it breaks qubit closure."
    },
    "largest_reference_sample_distance":max_reference_gap,
    "smallest_positive_non_measure_prepare_witness":min_positive_non_mp_witness,
    "scope":[
        "Operator inequalities certify all qubit inputs; passive references use the same marginal inequalities.",
        "Choi positive/negative decomposition gives an independent complete-norm upper certificate.",
        "Reference and adaptive samples audit implementation, not all-input or all-history proofs.",
        "No SDP, image checks, author-code imports, or automatic result-file writes.",
        "No generated endpoints, actual reverse operation, or actual spatial stabilizer are claimed."
    ],
    "numpy_version":np.__version__,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
print(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False))

