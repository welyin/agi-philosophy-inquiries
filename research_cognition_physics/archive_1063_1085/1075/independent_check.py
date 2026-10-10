"""Round 1075 independent finite-dimensional checks. stdout JSON; no files written.

Analytic classification and universal bounds require the accompanying proof.
This independent numerical check imports no author code and uses no solver.
"""
import json
import platform
import numpy as np

SEED = 107501
RNG = np.random.default_rng(SEED)
TOL = 3e-11
GROUPS = {}
CURRENT = None

def start(name):
    global CURRENT
    CURRENT = {"assertions": 0, "max_equality_residual": 0.0,
               "minimum_inequality_margin": None}
    GROUPS[name] = CURRENT

def equal(a, b, label):
    residual = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    CURRENT["assertions"] += 1
    CURRENT["max_equality_residual"] = max(
        CURRENT["max_equality_residual"], residual)
    if residual > TOL:
        raise AssertionError((label, residual))

def lower(value, bound, label):
    margin = float(value - bound)
    CURRENT["assertions"] += 1
    old = CURRENT["minimum_inequality_margin"]
    CURRENT["minimum_inequality_margin"] = margin if old is None else min(old, margin)
    if margin < -TOL:
        raise AssertionError((label, value, bound, margin))

def dagger(x):
    return x.conj().T

def ket(d):
    v = RNG.normal(size=d) + 1j * RNG.normal(size=d)
    return v / np.linalg.norm(v)

def proj(v):
    return np.outer(v, v.conj())

def density(d):
    a = RNG.normal(size=(d,d)) + 1j * RNG.normal(size=(d,d))
    a = a @ dagger(a)
    return a / np.trace(a)

def td(a, b):
    return float(np.linalg.svd(a - b, compute_uv=False).sum() / 2)

def minval(a):
    return float(np.linalg.eigvalsh((a + dagger(a))/2)[0])

def apply(kraus, x):
    return sum((k @ x @ dagger(k) for k in kraus),
               np.zeros((kraus[0].shape[0], kraus[0].shape[0]), complex))

def swap(d):
    s = np.zeros((d*d,d*d), complex)
    for i in range(d):
        for j in range(d):
            s[j*d+i,i*d+j] = 1
    return s

def anti_basis(d):
    cols = []
    for i in range(d):
        for j in range(i+1,d):
            v = np.zeros(d*d, complex)
            v[i*d+j], v[j*d+i] = 1/np.sqrt(2), -1/np.sqrt(2)
            cols.append(v)
    return np.column_stack(cols)

def ptr2(x, d):
    return np.einsum("ijkj->ik", x.reshape(d,d,d,d))

def ptr1(x, d):
    return np.einsum("ijil->jl", x.reshape(d,d,d,d))

def random_channel(din, dout, nk):
    a = RNG.normal(size=(dout*nk,din)) + 1j*RNG.normal(size=(dout*nk,din))
    q, _ = np.linalg.qr(a)
    return [q[j*dout:(j+1)*dout,:] for j in range(nk)]

def reset_channel(din, state):
    vals, vecs = np.linalg.eigh(state)
    return [np.sqrt(max(v,0))*np.outer(vecs[:,j], np.eye(din)[i])
            for j,v in enumerate(vals) for i in range(din) if v > 1e-15]

def classified(d, anti_k):
    s = swap(d)
    ps = (np.eye(d*d)+s)/2
    a = anti_basis(d)
    sym_k = [np.kron(np.eye(d), np.eye(d)[j:j+1]) @ ps for j in range(d)]
    ks = sym_k + [k @ dagger(a) for k in anti_k]
    def block(x):
        return ptr2(ps @ x @ ps,d) + apply(anti_k, dagger(a) @ x @ a)
    return ks, block, ps, a

def choi_from_action(fn, din, dout):
    out = np.zeros((din*dout,din*dout), complex)
    for i in range(din):
        for j in range(din):
            e = np.zeros((din,din), complex)
            e[i,j] = 1
            out[i*dout:(i+1)*dout,j*dout:(j+1)*dout] = fn(e)
    return out

def choi_from_kraus(ks):
    return sum((proj(k.T.reshape(-1)) for k in ks))

def orthonormal_pair(d):
    u = ket(d)
    v = ket(d)
    v -= u * np.vdot(u,v)
    v /= np.linalg.norm(v)
    return u,v

def mix_kraus(k1, k2, t):
    return [np.sqrt(1-t)*k for k in k1] + [np.sqrt(t)*k for k in k2]

def closest_pure(state):
    _, vecs = np.linalg.eigh((state+dagger(state))/2)
    return proj(vecs[:,-1])

start("classification_choi_and_full_operator_blocks")
CLASS_CHANNELS = []
for d in (2,3):
    da = d*(d-1)//2
    for case in range(4):
        ak = reset_channel(da, density(d)) if case < 2 else random_channel(da,d,4)
        ks, fn, ps, a = classified(d,ak)
        pa = a @ dagger(a)
        equal(sum(dagger(k)@k for k in ks),np.eye(d*d),"TP on full input")
        cb = choi_from_action(fn,d*d,d)
        ck = choi_from_kraus(ks)
        equal(cb,ck,"Choi via independent matrix-unit action")
        lower(minval(cb),0,"Choi positivity")
        equal(np.einsum("iaja->ij", cb.reshape(d*d,d,d*d,d)),np.eye(d*d),"Choi TP")
        for _ in range(6):
            x = RNG.normal(size=(d*d,d*d))+1j*RNG.normal(size=(d*d,d*d))
            equal(apply(ks,x),fn(x),"arbitrary non-Hermitian block formula")
            equal(fn(swap(d)@x@swap(d)),fn(x),"global role symmetry")
            equal(fn(ps@x@pa),np.zeros((d,d)),"symmetric/antisymmetric cross killed")
            p = proj(ket(d))
            equal(fn(np.kron(p,p)),p,"pure diagonal idempotence")
        CLASS_CHANNELS.append((d,ks,fn,ps,a))

start("forced_spectrum_any_pure_target_and_pairwise_sharpness")
PAIR_CASES = []
for d,ks,fn,ps,a in CLASS_CHANNELS:
    for c in (0.0,0.05,0.25,0.5,0.8,0.98,1.0):
        u,v = orthonormal_pair(d)
        w = c*u+np.sqrt(max(0,1-c*c))*v
        inp = proj(np.kron(u,w))
        sym = ptr2(ps@inp@ps,d)
        am, ap = (1-c)**2/4, (1+c)**2/4
        expected = np.sort(np.array([am,ap]+[0.]*(d-2)))
        equal(np.linalg.eigvalsh(sym),expected,"forced sym-block eigenvalues")
        anti = fn(inp)-sym
        lower(minval(anti),0,"anti contribution positive")
        equal(np.trace(anti).real,(1-c*c)/2,"anti branch weight")
        out = fn(inp)
        qnear = closest_pure(out)
        equal(td(out,qnear),1-np.linalg.eigvalsh(out)[-1],"closest-pure spectral formula")
        lower(td(out,qnear),am,"closest-pure universal lower bound")
        for _ in range(3):
            lower(td(out,proj(ket(d))),am,"finite other pure target")
        plus = (u+w)/np.linalg.norm(u+w)
        _,sharp_fn,_,_ = classified(d,reset_channel(d*(d-1)//2,proj(plus)))
        sharp_out = sharp_fn(inp)
        equal(td(sharp_out,proj(plus)),am,"specified pair reaches sharp bound")
        if c < 1:
            lower(np.linalg.eigvalsh(out)[-2],am,"mixedness retained")
        PAIR_CASES.append({"dimension":d,"c":c,"sharp_bound":am})

start("qubit_covariance_specialization")
kc,fc,ps2,a2 = classified(2,reset_channel(1,np.eye(2)/2))
for _ in range(12):
    x = RNG.normal(size=(4,4))+1j*RNG.normal(size=(4,4))
    equal(fc(x),(ptr1(x,2)+ptr2(x,2))/2,"covariant channel is averaged marginal")
    v = RNG.normal(size=(2,2))+1j*RNG.normal(size=(2,2))
    u,_ = np.linalg.qr(v)
    uu = np.kron(u,u)
    equal(fc(uu@x@dagger(uu)),u@fc(x)@dagger(u),"full simultaneous unitary covariance")
for c in (0.,0.2,0.5,0.9,1.):
    u,v = orthonormal_pair(2)
    w = c*u+np.sqrt(1-c*c)*v
    out = fc(proj(np.kron(u,w)))
    equal(td(out,closest_pure(out)),(1-c)/2,"covariant nearest-pure gap")

z0,z1 = np.eye(2,dtype=complex)
PX = [proj((z0+z1)/np.sqrt(2)),proj((z0-z1)/np.sqrt(2))]
PY = [proj((z0+1j*z1)/np.sqrt(2)),proj((z0-1j*z1)/np.sqrt(2))]
PZ = [proj(z0),proj(z1)]
SIX = PX+PY+PZ
SIGNS = [1,1,1,1,-1,-1]
t0 = proj((np.kron(z0,z1)+np.kron(z1,z0))/np.sqrt(2))
singlet = proj((np.kron(z0,z1)-np.kron(z1,z0))/np.sqrt(2))
in01,in10 = proj(np.kron(z0,z1)),proj(np.kron(z1,z0))
start("six_pauli_tensor_witness")
equal(sum(np.kron(p,p) for p in SIX),2*ps2,"six-state symmetric projector identity")
equal(sum(s*np.kron(p,p) for s,p in zip(SIGNS,SIX))/2,t0,"signed triplet identity")
equal(in01+in10,t0+singlet,"two ordered inputs decomposition")
equal(sum(s*p for s,p in zip(SIGNS,SIX))/2,np.eye(2)/2,"signed ideal outputs")

start("eight_input_robust_bound_random_and_near_contract_channels")
ROBUST_ROWS = []
def robust_case(ks,name):
    outputs = [apply(ks,np.kron(p,p)) for p in SIX]
    delta = max(td(y,p) for y,p in zip(outputs,SIX))
    fwd,rev = apply(ks,in01),apply(ks,in10)
    kappa = td(fwd,rev)
    raw = .25-1.5*delta-.5*kappa
    bound = max(0.,raw)
    trip = apply(ks,t0)
    witness = sum(s*y for s,y in zip(SIGNS,outputs))/2
    equal(witness,trip,"direct map agrees with six-output witness")
    lower(minval(trip),.5-3*delta,"triplet matrix lower bound")
    lower(minval((fwd+rev)/2),.25-1.5*delta,"averaged ordered output matrix bound")
    lower(minval(fwd),raw,"role-defect corrected matrix bound")
    nearest = td(fwd,closest_pure(fwd))
    lower(nearest,bound,"robust closest-pure distance")
    for _ in range(2):
        lower(td(fwd,proj(ket(2))),bound,"robust other pure target")
    ROBUST_ROWS.append({"name":name,"delta":delta,"kappa":kappa,
                        "bound":bound,"closest_pure_distance":nearest,
                        "margin":nearest-bound})

for i in range(32):
    robust_case(random_channel(4,2,4+(i%3)), "generic_CPTP_"+str(i))
partialA = [np.kron(np.eye(2),np.eye(2)[j:j+1]) for j in range(2)]
reset0 = reset_channel(4,proj(z0))
base,_,_,_ = classified(2,reset_channel(1,proj(z0)))
for t in (0.,.001,.01,.03,.08,.15,.3):
    noise = random_channel(4,2,5)
    robust_case(mix_kraus(base,noise,t),"generic_small_contamination_"+str(t))
    robust_case(mix_kraus(base,partialA,t),"role_asymmetric_marginal_"+str(t))
    robust_case(mix_kraus(base,reset0,t),"diagonal_imperfect_reset_"+str(t))
for theta in (.001,.01,.04,.1,.2):
    rot = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    rotated = [rot@k for k in base]
    robust_case(rotated,"rotated_output_"+str(theta))
    robust_case(mix_kraus(rotated,partialA,.07),"rotation_plus_role_asymmetry_"+str(theta))
positive = [r for r in ROBUST_ROWS if r["bound"] > 1e-12]
asymmetric_positive = [r for r in positive if r["kappa"] > 1e-6]
lower(len(positive),20,"nonvacuous robust sample count")
lower(len(asymmetric_positive),10,"nonvacuous non-swap-symmetric sample count")

start("deleted_assumption_witnesses")
equal(apply(partialA,np.kron(proj(z0),proj(z1))),proj(z0),"dropping symmetry permits pure output")
equal(td(apply(partialA,in01),apply(partialA,in10)),1.,"marginal channel violates role symmetry")
equal(apply(reset0,in01),proj(z0),"dropping diagonal preservation permits pure output")
equal(td(apply(reset0,np.kron(proj(z1),proj(z1))),proj(z1)),1.,"reset violates diagonal task")
equal(fc(in01),np.eye(2)/2,"dropping pure-result demand leaves legal average")

summary = {
    "status":"PASS",
    "seed":SEED,
    "runtime":{"python":platform.python_version(),"numpy":np.__version__},
    "tolerance":TOL,
    "assertion_groups":len(GROUPS),
    "assertions":sum(g["assertions"] for g in GROUPS.values()),
    "max_equality_residual":max(g["max_equality_residual"] for g in GROUPS.values()),
    "minimum_inequality_margin":min(g["minimum_inequality_margin"] for g in GROUPS.values()
                                    if g["minimum_inequality_margin"] is not None),
    "groups":GROUPS,
    "coverage":{"classified_channels":len(CLASS_CHANNELS),
                "pure_input_pairs":len(PAIR_CASES),
                "robust_channels":len(ROBUST_ROWS),
                "nonvacuous_robust_channels":len(positive),
                "nonvacuous_role_asymmetric_channels":len(asymmetric_positive)},
    "robust_examples":ROBUST_ROWS,
    "scope":[
        "Matrix and finite-input checks corroborate, but do not prove, universal statements.",
        "Classification tests use d=2 and d=3; the open-patch extension is analytic.",
        "Sharpness is per specified pair, not simultaneous optimization by one fixed antisymmetric channel.",
        "The eight-input bound assumes neither global swap symmetry nor exact diagonal preservation.",
        "No assertion identifies internal pure states with actual spatial endpoints.",
        "No author module imported; script writes no files."
    ]
}
print(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False))