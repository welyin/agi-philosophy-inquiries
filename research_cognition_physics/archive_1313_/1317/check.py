"""Round 1317. Default / --check are read-only; --write-new exclusively initializes JSON.
No scientific module is imported. NumPy diagnostics complement the analytic proof.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent
ATOL = 3e-10
CHECKS = 0

def require(test, label):
    global CHECKS
    CHECKS += 1
    if not bool(test):
        raise AssertionError(label)

def op(a):
    return float(np.linalg.norm(a, 2))

def close(a, b, label, tol=ATOL):
    err = op(np.atleast_2d(a - b))
    require(err <= tol, label + ": " + str(err))
    return err

def iv(x):
    return (F(x), F(x))

def add(x, y):
    return (x[0]+y[0], x[1]+y[1])

def neg(x):
    return (-x[1], -x[0])

def mul(x, y):
    v = [a*b for a in x for b in y]
    return (min(v), max(v))

def div(x, y):
    require(y[0] > 0, "positive interval divisor")
    return mul(x, (1/y[1], 1/y[0]))

def power(x, n):
    out = iv(1)
    for _ in range(n):
        out = mul(out, x)
    return out

def trig_interval(x, sine):
    # Taylor polynomial through degree 49 (sin) or 48 (cos);
    # use Lagrange remainder through total degree 49 for both.
    total = iv(0)
    for k in range(1 if sine else 0, 50, 2):
        sign = (-1)**((k-1)//2 if sine else k//2)
        total = add(total, mul(iv(F(sign, math.factorial(k))), power(x, k)))
    r = max(abs(x[0]), abs(x[1]))**50 / math.factorial(50)
    return (total[0]-r, total[1]+r)

def rational_lower(x, digits=12):
    den = 10**digits
    return F((x.numerator * den)//x.denominator, den)

def rational_upper(x, digits=12):
    return -rational_lower(-x, digits)

def bounds(x):
    return [str(rational_lower(x[0])), str(rational_upper(x[1]))]

def exact_gramian_certificate():
    w = (F(1414213562373095, 10**15), F(1414213562373096, 10**15))
    require(w[0]**2 < 2 < w[1]**2, "sqrt2 rational bracket")
    s2 = trig_interval(mul(iv(2), w), True)
    c2 = trig_interval(mul(iv(2), w), False)
    s4 = trig_interval(mul(iv(4), w), True)
    x = iv(F(4, 3))
    y = add(mul(iv(F(-1, 2)), c2), div(s2, mul(iv(4), w)))
    z = add(iv(F(1, 4)), neg(div(s4, mul(iv(16), w))))
    lam = iv(F(1, 20))
    a, b = add(x, neg(lam)), add(z, neg(lam))
    determinant = add(mul(a, b), neg(mul(y, y)))
    require(a[0] > 0 and determinant[0] > 0, "W-I/20 Sylvester certificate")
    return {
        "method": "rational sqrt bracket, degree49 Taylor with exact Lagrange remainder",
        "W_pm_entries": {"++": bounds(x), "+-": bounds(y), "--": bounds(z)},
        "first_minor_lower": str(rational_lower(a[0])),
        "determinant_lower": str(rational_lower(determinant[0])),
        "lambda_certified": "1/20",
    }

def chain(N, n, m=1., q=1.):
    Q = np.zeros((N, N+1))
    for j in range(N):
        Q[j, j], Q[j, j+1] = 1., -q
    u = q**(-np.arange(N+1, dtype=float))
    u /= np.linalg.norm(u)
    return m*np.kron(Q, np.eye(n)), u

def block_h(Z, T, rho=1., cut=None):
    Tc = T.copy()
    n = Z.shape[0]
    if cut is not None:
        Tc[cut*n:(cut+1)*n] = 0.
    M = np.vstack((rho*Z, Tc))
    h = np.block([[np.zeros((M.shape[0], M.shape[0])), M],
                  [M.conj().T, np.zeros((M.shape[1], M.shape[1]))]])
    return h, M, Tc

def embed(n, N):
    size = 2*(N+1)*n
    E0 = np.zeros((size, n), complex)
    E0[(N+1)*n:(N+2)*n] = np.eye(n)
    return E0

def evolve(h, t):
    w, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*w)) @ v.conj().T

def sinfunc(B, t):
    w, v = np.linalg.eigh(B)
    val = np.full(len(w), t)
    mask = w > 1e-12
    val[mask] = np.sin(t*np.sqrt(w[mask]))/np.sqrt(w[mask])
    return (v*val)@v.conj().T

def heavy(T):
    G = T@T.conj().T
    return T.conj().T@np.linalg.solve(G, np.linalg.solve(G, T))

def fixed_examples():
    alpha, zmax, eps, tau = 1e-6, .6, 1e-4, 2.
    a_cert = eps*(1+zmax*zmax*tau*tau/2)*math.sqrt(tau/.05)
    d_cert = math.sqrt(3)*a_cert
    c_cert = 2*zmax*d_cert+d_cert*d_cert
    require(c_cert < .00227 and c_cert/.18**2 < .071, "nonvacuous declared bound")
    rows = []
    for n in (3, 4, 5):
        N=2
        j=np.arange(n)
        U=np.exp(2j*np.pi*np.outer(j,j)/n)/math.sqrt(n)
        Y=.18*U@np.diag(1+j/20)
        Z0=alpha*np.diag((j+1)/n)
        Z1=1j*alpha*np.roll(np.eye(n), 1, axis=0)
        Z2=math.sqrt(3)*Y-Z0-Z1
        Z=np.hstack((Z0,Z1,Z2))
        T,u=chain(N,n)
        f=u[-1]
        Zpre=Z[:,:N*n]
        Zref=np.hstack((np.zeros((n,N*n)),Y/f))
        require(op(Z)<zmax, "registered portal bound")
        require(op(Z0@Z1-Z1@Z0)>alpha**2/10, "noncommuting blocks")
        close(Z@np.kron(u[:,None],np.eye(n)),Y,"same Y")
        require(np.min(np.linalg.svd(Y,compute_uv=False))>.179999999, "Y lower singular value")
        # Entire time window uses Duhamel, not this diagnostic grid.
        p_analytic=tau*tau*2*alpha*alpha
        require(p_analytic<eps*eps, "analytic whole window")
        hc,Mc,Tc=block_h(Z,T,cut=N-1)
        E0=embed(n,N)
        grid=np.linspace(0,tau,81)
        maxima=max(op((evolve(hc,t)@E0)[:n]) for t in grid)
        require(maxima<=math.sqrt(p_analytic)+ATOL,"grid diagnostic inside analytic bound")
        R=heavy(T)
        g2=float(np.min(np.linalg.eigvalsh(T@T.conj().T)))
        C=Z@R@Z.conj().T
        Cr=Zref@R@Zref.conj().T
        close(Cr,5/3*Y@Y.conj().T,"endpoint identity")
        delta=op(Z-Zref)
        require(delta<=op(Zpre)/f+ATOL,"same-Y portal factor")
        require(op(Zpre)<=a_cert and delta<d_cert,"finite certificate bound")
        actual=op(C-Cr)
        require(actual<=c_cert,"heavy response bound")
        D=Y@Y.conj().T
        w,v=np.linalg.eigh(D)
        inv=(v*(w**-.5))@v.conj().T
        spectrum=np.linalg.eigvalsh(inv@C@inv)
        require(max(abs(spectrum-5/3))<c_cert/.18**2,"normalized response")
        # Arbitrary passive partner: both outcomes, no normalized postselection.
        rng=np.random.default_rng(131700+n)
        psi=rng.normal(size=(n,n+1))+1j*rng.normal(size=(n,n+1))
        psi/=np.linalg.norm(psi)
        t=1.13
        V=evolve(hc,t)@E0
        Ks=V[:n]
        Kf=V[n:]
        close(Ks.conj().T@Ks+Kf.conj().T@Kf,np.eye(n),"complete instrument")
        success=Ks@psi
        failure=Kf@psi
        rerr=close(success.conj().T@success+failure.conj().T@failure,
                   psi.conj().T@psi,"passive partner retained")
        require(abs(np.linalg.norm(success)**2+np.linalg.norm(failure)**2-1)<ATOL,"full probability")
        # Finite selected parameter derivatives; source identity precedes constraints.
        dZ=np.zeros_like(Z,dtype=complex); dZ[0,0]=1j
        dh,_ ,_=block_h(dZ,np.zeros_like(T),cut=N-1)
        step=1e-6
        hp=block_h(Z+step*dZ,T,cut=N-1)[0]
        hm=block_h(Z-step*dZ,T,cut=N-1)[0]
        close((hp-hm)/(2*step),dh,"imag portal source")
        dM=np.vstack((dZ,np.zeros_like(T)))
        close((Mc+step*dM).conj().T@(Mc+step*dM)-
              (Mc-step*dM).conj().T@(Mc-step*dM),
              2*step*(dM.conj().T@Mc+Mc.conj().T@dM),"complete quadratic source",1e-12)
        rows.append({"n":n,"portal_norm":op(Z),"prefix_norm":op(Zpre),
                     "commutator_norm":op(Z0@Z1-Z1@Z0),"same_Y_delta":delta,
                     "probability_entire_window_bound":p_analytic,
                     "sample_max_amplitude_diagnostic":maxima,
                     "response_difference":actual,"generalized_spectrum":spectrum.tolist(),
                     "partner_error":rerr,"gap_squared":g2})
    return {"fixed_contract":{"N":2,"m":1,"q":1,"rho":1,"tau":2,"zmax":.6,
                              "probability_tolerance":1e-8,"alpha":alpha},
            "certified_bounds":{"prefix":a_cert,"same_Y_delta":d_cert,
                                "response":c_cert,"generalized":c_cert/.18**2,
                                "relative_generalized":c_cert/.18**2/(5/3)},
            "examples":rows}

def volterra_diagnostics():
    rng=np.random.default_rng(131703)
    N,n,m,q,rho=3,2,.7,1.3,.9
    T,u=chain(N,n,m,q)
    Z=.18*(rng.normal(size=(n,(N+1)*n))+1j*rng.normal(size=(n,(N+1)*n)))
    h,M,Tc=block_h(Z,T,rho,cut=N-1)
    B=Tc.conj().T@Tc
    E0=embed(n,N)
    Er=np.zeros(((N+1)*n,n)); Er[:n]=np.eye(n)
    nodes,weights=np.polynomial.legendre.leggauss(64)
    errors=[]
    for t in (.17,.71,1.6):
        A=(evolve(h,t)@E0)[:n]
        integ=np.zeros((n,n),complex)
        for z,w in zip(nodes,weights):
            s=t*(z+1)/2
            As=(evolve(h,s)@E0)[:n]
            integ+=w*t/2*Z@sinfunc(B,t-s)@Z.conj().T@As
        rhs=-1j*rho*Z@sinfunc(B,t)@Er-rho*rho*integ
        errors.append(close(A,rhs,"Volterra exact identity diagnostic"))
    # Fixed grid plus derivative is a separate, coarse but valid enclosure.
    tau=1.6
    grid=np.linspace(0,tau,81)
    gridmax=max(op((evolve(h,t)@E0)[:n]) for t in grid)
    enclosure=gridmax+rho*op(Z)*(grid[1]-grid[0])/2
    densemax=max(op((evolve(h,t)@E0)[:n]) for t in np.linspace(0,tau,801))
    require(densemax<=enclosure+ATOL,"grid derivative upper bound diagnostic")
    return {"parameters":[N,n,m,q,rho],"identity_errors":errors,
            "grid_derivative_enclosure":enclosure,"dense_sample_max":densemax,
            "scope":"Quadrature and dense values are diagnostics; proof supplies all-time bounds."}

def deletion_examples():
    simple=[]
    for K in (1.,10.,100.):
        T,u=chain(1,1)
        Z=np.array([[1.,K]])
        h,M,Tc=block_h(Z,T,cut=0)
        t=math.pi/(2*math.sqrt(1+K*K))
        amp=(evolve(h,t)@embed(1,1))[0,0]
        require(abs(abs(amp)-1/math.sqrt(1+K*K))<ATOL,"N1 shielding")
        simple.append({"K":K,"prefix":1.,"global_probability_max":1/(1+K*K)})
    fixed=[]
    for K in (4.,16.,64.,256.):
        T,u=chain(2,1)
        Z=np.array([[0.,K,1-K]])
        close(Z@u[:,None],np.array([[1/math.sqrt(3)]]),"fixed Y shielding")
        h,M,Tc=block_h(Z,T,cut=1)
        G=np.array([[K*K+(1-K)**2,-K],[-K,2.]])
        vals=np.linalg.eigvalsh(G)
        lm,lp=vals
        bound=K/(lp-lm)*(1/math.sqrt(lp)+1/math.sqrt(lm))
        maxima=max(abs((evolve(h,t)@embed(1,2))[0,0]) for t in np.linspace(0,6,61))
        require(maxima<=bound+ATOL,"fixed Y full time upper bound")
        require(abs(np.linalg.det(G)-(3*K*K-4*K+2))<1e-4,"G determinant")
        fixed.append({"K":K,"prefix":K,"Y":1/math.sqrt(3),
                      "global_amplitude_bound":bound,"lambda_minus":lm,
                      "sample_max_diagnostic":maxima})
    require(fixed[-1]["global_amplitude_bound"]<fixed[0]["global_amplitude_bound"]/50,
            "shielding trend")
    # AND is a changed permission class; same-Y all-yes family, all-no portal zero.
    T,u=chain(2,3)
    Y=np.diag([.12,.15,.18]).astype(complex)
    R=heavy(T); f=u[-1]
    androws=[]
    for t in (0.,.5,1.):
        p=t*u.copy(); p[-1]+=(1-t)/f
        Z=np.hstack([p[j]*Y for j in range(3)])
        close(Z@np.kron(u[:,None],np.eye(3)),Y,"AND same Y")
        C=Z@R@Z.conj().T
        close(C,(1-t)**2*5/3*Y@Y.conj().T,"AND freely variable response")
        h=block_h(np.zeros_like(Z),T,cut=1)[0]
        leak=op((evolve(h,.83)@embed(3,2))[:3])
        require(leak<ATOL,"AND veto")
        androws.append({"t":t,"response_norm":op(C),"refusal_amplitude":leak})
    return {"N1_without_fixed_Y":simple,"N2_fixed_Y":fixed,
            "AND_changed_control_class":androws}

def car_instrument():
    modes=4; dim=1<<modes
    creation=[]
    for k in range(modes):
        a=np.zeros((dim,dim),complex)
        for state in range(dim):
            if not ((state>>k)&1):
                a[state|(1<<k),state]=(-1)**((state&((1<<k)-1)).bit_count())
        creation.append(a)
    for i in range(modes):
        for j in range(modes):
            close(creation[i].conj().T@creation[j]+creation[j]@creation[i].conj().T,
                  np.eye(dim)*(i==j),"CAR")
    def second(h):
        out=np.zeros((dim,dim),complex)
        for i in range(modes):
            for j in range(modes):
                out+=h[i,j]*creation[i]@creation[j].conj().T
        return out
    T,u=chain(1,1,.7,1.2)
    Z=np.array([[.1+.05j,.2-.03j]])
    HC=np.diag([0.,.3])
    blocks=[]
    single=np.eye(dim)[:,[1<<k for k in range(modes)]]
    number=np.diag([state.bit_count() for state in range(dim)])
    for c in (0,1):
        h=block_h(Z,T,.8,cut=0 if c==0 else None)[0]
        dg=second(h)
        close(single.conj().T@dg@single,h,"single sector exact")
        close(dg@number-number@dg,np.zeros_like(dg),"number conservation")
        # Offset is full Fock identity; NOT chemical potential.
        blocks.append(dg+(10.+HC[c,c])*np.eye(dim))
    full=np.block([[blocks[0],np.zeros((dim,dim))],
                   [np.zeros((dim,dim)),blocks[1]]])
    perm=np.kron(np.diag([0.,1.]),np.eye(dim))
    close(full@perm-perm@full,np.zeros_like(full),"orthogonal permission record retained")
    require(np.min(np.linalg.eigvalsh(full))>0,"fixed global offset positive")
    t=.61
    U=evolve(full,t)
    close(U.conj().T@full@U,full,"static complete energy conservation")
    return {"CAR_modes":modes,"Fock_dimension":dim,"with_permission_dimension":2*dim,
            "fixed_full_Fock_offset":10.,"minimum_energy":float(np.min(np.linalg.eigvalsh(full))),
            "scope":"Finite exact CAR parent; no spatial field-device or autonomous withdrawal claim."}

def provenance():
    paths=[
        HERE/"proof.md", HERE/"NEXT.md", STAGE/"research_note_1317.md",
        STAGE/"_admission/round1317_admission.md",
        STAGE/"1316/proof.md", STAGE/"1316/acceptance.json", STAGE/"1313/proof.md",
        STAGE/"_admission/refusal_chain_screen_after1316.md",
        STAGE/"_admission/joint_consent_scope_screen_after1316.md",
        STAGE/"_admission/joint_mechanism_screen_after1316.md",
        ROOT/"archive_956_989/research_note_988.md",
        ROOT/"archive_923_934/research_note_925.md",
        ROOT/"archive_923_934/research_note_933.md",
        ROOT/"archive_935_955/research_note_954.md",
        ROOT/"archive_956_989/research_note_968.md",
    ]
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def run():
    global CHECKS
    CHECKS=0
    out={"round":1317,"status":"PASS","new_axioms":0,
         "exact_gramian":exact_gramian_certificate(),
         "finite_positive_examples":fixed_examples(),
         "volterra":volterra_diagnostics(),
         "deletions":deletion_examples(),
         "finite_CAR":car_instrument()}
    out["sources_sha256"]=provenance()
    out["code_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out["assertions"]=CHECKS
    return out

def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    path=HERE/"results.json"
    serialized=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n"
    if args.write_new:
        with path.open("x",encoding="utf-8",newline="\n") as f:
            f.write(serialized)
    else:
        stored=json.loads(path.read_text(encoding="utf-8"))
        if result!=stored:
            raise AssertionError("Saved result differs; read-only mode will not rewrite it.")
    print(json.dumps({"status":"PASS","round":1317,"assertions":result["assertions"],
                      "mode":"exclusive-init" if args.write_new else "read-only",
                      "response_bound":result["finite_positive_examples"]["certified_bounds"]["response"],
                      "source_hashes":len(result["sources_sha256"])},ensure_ascii=False))
if __name__=="__main__":
    main()
