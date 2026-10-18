"""1295 independent checks; no author/historical computational imports.
Exact rational structural certificates plus small finite CAR/Fourier diagnostics.
Default and --check are read-only; --write creates results exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/"research_cognition_physics").is_dir())
OUT=HERE/"independent_results.json"
EPS=F(31,2**34)
ALPHA=F(31,2**27)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def add(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(k,a): return [[k*x for x in r] for r in a]
def mul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def transpose(a): return list(map(list,zip(*a)))
def rank(a):
    a=[list(map(F,r)) for r in a]; n=len(a[0]); k=0
    for j in range(n):
        pivot=next((i for i in range(k,len(a)) if a[i][j]),None)
        if pivot is None: continue
        a[k],a[pivot]=a[pivot],a[k]
        q=a[k][j]; a[k]=[x/q for x in a[k]]
        for i in range(len(a)):
            if i!=k and a[i][j]:
                q=a[i][j]; a[i]=[x-q*y for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a): break
    return k

def commutant_rank(mats):
    n=len(mats[0]); columns=[]
    for i in range(n):
        for j in range(n):
            e=[[F(0) for _ in range(n)] for _ in range(n)]; e[i][j]=1
            col=[]
            for a in mats:
                c=add(mul(a,e),scale(-1,mul(e,a)))
                col.extend(x for r in c for x in r)
            columns.append(col)
    return rank(transpose(columns))

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0: p.pop()
    return p
def padd(a,b):
    z=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): z[i]+=x
    for i,x in enumerate(b): z[i]+=x
    return trim(z)
def pdivide(a,b):
    a=trim(a); b=trim(b)
    while len(a)>=len(b) and any(a):
        q=a[-1]/b[-1]; k=len(a)-len(b)
        for j,x in enumerate(b): a[j+k]-=q*x
        a=trim(a)
    return a
def pgcd(a,b):
    while any(b): a,b=b,pdivide(a,b)
    return [x/a[-1] for x in a]

def structures():
    out=[]
    for n in (3,4):
        I=eye(n)
        D=[[F(i+1 if i==j else 0) for j in range(n)] for i in range(n)]
        P=[[F(abs(i-j)==1) for j in range(n)] for i in range(n)]
        A=add(D,P)
        yn=add(I,scale(EPS,A)); ye=add(scale(2,I),scale(EPS,D))
        hn=mul(yn,transpose(yn)); he=mul(ye,transpose(ye))
        # Remove scalar parts and divide by epsilon exactly, before any rank test.
        hn_scaled=scale(1/EPS,add(hn,scale(-1,I)))
        he_scaled=scale(1/EPS,add(he,scale(-4,I)))
        assert hn_scaled==add(scale(2,A),scale(EPS,mul(A,A)))
        assert he_scaled==add(scale(4,D),scale(EPS,mul(D,D)))
        neighbors=[hn[i][i+1] for i in range(n-1)]
        assert all(x==2*EPS+EPS**2*(2*i+3)>0 for i,x in enumerate(neighbors))
        cr=commutant_rank([hn_scaled,he_scaled])
        assert cr==n*n-1
        dcr=commutant_rank([D,mul(D,D)])
        assert dcr==n*n-n
        # Integer Jacobi recurrence for det(x I-A); gcd with derivative proves simple roots.
        pm2=[F(1)]; pm1=[F(-1),F(1)]
        for j in range(2,n+1):
            q=padd([0]+pm1,[-j*x for x in pm1])
            q=padd(q,[-x for x in pm2])
            pm2,pm1=pm1,q
        deriv=[i*pm1[i] for i in range(1,len(pm1))]
        assert len(pgcd(pm1,deriv))==1
        mass_lo=1-EPS; mass_hi=1+(n+2)*EPS
        assert mass_lo>0
        out.append(dict(n=n,epsilon=EPS,A_characteristic_coefficients_ascending=pm1,
          characteristic_gcd_degree=0,nu_mass_interval=[mass_lo,mass_hi],
          electron_masses=[2+EPS*j for j in range(1,n+1)],
          adjacent_nu_squared_entries=neighbors,commutant_complex_dimension=n*n-cr,
          path_deleted_commutant_dimension=n*n-dcr,
          epsilon_deleted_commutant_dimension=n*n,
          normalized_nu_squared=hn_scaled,normalized_e_squared=he_scaled))
    d=3*EPS
    assert d/ALPHA==F(3,128)
    gap=ALPHA**2-2*ALPHA*d
    assert gap==F(61,64)*ALPHA**2
    singular=ALPHA-d
    assert singular==F(125,128)*ALPHA
    energy_squared=4*(F(524288)+(1+6*EPS)**2)
    assert energy_squared<1450**2
    # An explicit additional mass norm margin; irreducibility/spectral openness is separate.
    open_d=4*EPS
    return dict(families=out,uniform_operator_difference_upper=d,
      whole_window_probability_gap_lower=gap,receiver_singular_value_lower=singular,
      energy_squared_strict_upper=energy_squared,energy_strict_upper=1450,
      additional_nu_perturbation_norm_less_than=EPS,
      enlarged_record_difference_upper=open_d,
      enlarged_probability_lower=ALPHA**2-2*ALPHA*open_d,
      enlarged_receiver_singular_lower=ALPHA-open_d,
      open_irreducibility_requires_intersection_with_positive_rank_and_simple_spectrum_neighborhood=True)

def nparr(a): return np.array([[float(x) for x in row] for row in a],complex)
def spectral(H,fn):
    v,U=np.linalg.eigh(H); return (U*fn(v))@U.conj().T
def unitary(H,t): return spectral(H,lambda v:np.exp(-1j*t*v))
def masses(n):
    D=np.diag(np.arange(1,n+1,dtype=float)); P=np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
    return np.eye(n)+float(EPS)*(D+P),2*np.eye(n)+float(EPS)*D
def blockdiag(a,b):
    z=np.zeros((a.shape[0],b.shape[0]),complex)
    return np.block([[a,z],[z.conj().T,b]])

def fourier_diagnostics():
    sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.])
    p=(.37,-.29,.11); spin=p[0]*sx+p[1]*sy+p[2]*sz
    zeta=np.array([1.,1j])/math.sqrt(2)
    out=[]; maximum=0.
    for n in (3,4):
        M,_=masses(n); S=np.kron(spin,np.eye(n)); B=np.kron(np.eye(2),M)
        h=np.block([[-S,B],[B,S]])
        h0=np.block([[-S,np.eye(2*n)],[np.eye(2*n),S]])
        hm=np.block([[S,B],[B,-S]])
        hp=blockdiag(h,hm)
        V=np.zeros((8*n,2),complex)
        for j in range(2):
            vec=np.kron(zeta,np.eye(n)[:,j])
            V[2*n:4*n,j]=vec/math.sqrt(2)
            V[6*n:8*n,j]=vec/math.sqrt(2)
        FF=V@V.conj().T; QQ=np.eye(8*n)-FF
        Ph=spectral(hp,lambda x:(x<0).astype(float))
        absH=spectral(hp,np.abs)
        e_expected=float(np.trace(FF@absH).real)
        fhf=float(np.linalg.norm(V.conj().T@hp@V,2))
        identity_error=0.; energy_values=[]
        for rho in [np.diag([1.,0.]),np.eye(2)/2,np.array([[.4,.2+.1j],[.2-.1j,.6]])]:
            C=QQ@Ph@QQ+V@rho@V.conj().T
            e=float(np.trace(hp@(C-Ph)).real)
            energy_values.append(e); identity_error=max(identity_error,abs(e-e_expected))
        times=[]
        for t in (7/16,1/2):
            U=unitary(h,t); U0=unitary(h0,t)
            diff=float(np.linalg.norm(U-U0,2))
            cap=t*float(EPS)*float(np.linalg.norm(np.diag(np.arange(1,n+1))+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1),2))
            assert diff<=cap+3e-14
            times.append(dict(time=t,full_operator_difference=diff,duhamel_upper=cap,
                              floating_diagnostic_allowance=3e-14))
        # A genuine Higgs-phase tangent in this neutrino block; full other-source completion is analytic.
        phase=np.block([[np.zeros_like(B),1j*B],[-1j*B,np.zeros_like(B)]])
        J=blockdiag(phase,phase)
        Up=unitary(hp,1/2)
        logical=V.conj().T@Up.conj().T@J@Up@V
        traceless=logical-np.trace(logical)*np.eye(2)/2
        nonscalar=float(np.linalg.norm(traceless,2))
        assert nonscalar>float(EPS)*.01
        maximum=max(maximum,fhf,identity_error)
        out.append(dict(n=n,times=times,FhF_norm=fhf,energy_values=energy_values,
             energy_trace_F_abs_h=e_expected,energy_identity_residual=identity_error,
             phase_source_code_traceless_norm=nonscalar,
             phase_source_traceless_divided_by_epsilon=nonscalar/float(EPS)))
    assert maximum<2e-13
    return dict(cases=out,maximum_algebra_residual=maximum,
      finite_momentum_samples_are_not_continuum_record_or_Hadamard_proofs=True)

def annihilator(m,j):
    a=np.zeros((2**m,2**m),complex)
    for b in range(2**m):
        if (b>>j)&1: a[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
    return a
def dq(K,ops):
    out=np.zeros_like(ops[0])
    for i in range(len(ops)):
        for j in range(len(ops)):
            if abs(K[i,j])>1e-16: out+=K[i,j]*(ops[i].conj().T@ops[j])
    return out
def fock_density(C):
    # Exterior powers of the diagonalizing unitary, distinct from Gibbs-log methods.
    v,U=np.linalg.eigh(C); m=len(v); v=np.clip(v,0,1)
    sets=[[i for i in range(m) if (b>>i)&1] for b in range(2**m)]
    out=np.zeros((2**m,2**m),complex)
    for occ in sets:
        weight=float(np.prod([v[i] if i in occ else 1-v[i] for i in range(m)]))
        if weight<1e-18: continue
        col=np.zeros(2**m,complex)
        for b,rows in enumerate(sets):
            if len(rows)==len(occ):
                col[b]=np.linalg.det(U[np.ix_(rows,occ)]) if occ else 1.
        out+=weight*np.outer(col,col.conj())
    return out

def finite_car_diagnostic():
    n=3; M,_=masses(n); k=.7; m=2*n
    h0=np.block([[-k*np.eye(n),M],[M,k*np.eye(n)]])
    order=[n,n+1]+list(range(n))+[n+2]
    h=h0[np.ix_(order,order)]
    Ph=spectral(h,lambda x:(x<0).astype(float))
    env=fock_density(Ph[2:,2:])
    ops=[annihilator(m,j) for j in range(m)]; size=2**m
    T=[]
    for occ in (1,2):
        E=np.zeros((size,2**(m-2)),complex)
        for b in range(2**(m-2)): E[4*b+occ,b]=1
        T.append(E)
    rho=np.array([[.4,.2+.1j],[.2-.1j,.6]],complex)
    state=sum((rho[i,j]*T[i]@env@T[j].conj().T for i in range(2) for j in range(2)),
              np.zeros((size,size),complex))
    C=np.zeros((m,m),complex)
    for i in range(m):
        for j in range(m): C[i,j]=np.trace(ops[j].conj().T@ops[i]@state)
    target=np.zeros((m,m),complex); target[2:,2:]=Ph[2:,2:]; target[:2,:2]=rho
    covariance_error=float(np.linalg.norm(C-target))
    N0=ops[0].conj().T@ops[0]; N1=ops[1].conj().T@ops[1]
    occupation_product=float(np.trace(N0@N1@state).real)
    assert abs(occupation_product)<1e-13 and np.linalg.det(rho).real>.1
    # Full operator identities on input blocks imply the arbitrary passive-partner statement.
    receiver_indices=[order.index(0),order.index(1)]
    words=[np.eye(size)]
    for j in receiver_indices:
        words=words+[x@ops[j] for x in words]+[x@ops[j].conj().T for x in words]
    blind=0.
    for O in words:
        common=np.trace(T[0].conj().T@O@T[0]@env)
        for i in range(2):
            for j in range(2):
                val=np.trace(T[j].conj().T@O@T[i]@env)
                blind=max(blind,abs(val-(common if i==j else 0)))
    # After propagation, all source expectations retain actual logical matrices.
    U=unitary(h,.49)
    tangent0=np.block([[np.zeros((n,n)),1j*M],[-1j*M,np.zeros((n,n))]])
    K=tangent0[np.ix_(order,order)]
    Kt=U.conj().T@K@U; O=dq(Kt,ops)
    baseline=np.trace(Kt[2:,2:]@Ph[2:,2:])
    source_error=0.
    for i in range(2):
        for j in range(2):
            actual=np.trace(T[j].conj().T@O@T[i]@env)
            expected=(baseline if i==j else 0)+Kt[j,i]
            source_error=max(source_error,abs(actual-expected))
    # Input/partner expansion coefficients tested directly, not inferred from a Gaussian replacement.
    bell=np.array([[1,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,1]],complex)/2
    joint=sum((np.kron(T[i]@env@T[j].conj().T,bell[2*i:2*i+2,2*j:2*j+2])
               for i in range(2) for j in range(2)),np.zeros((2*size,2*size),complex))
    observed=np.array([[np.trace(O@joint[a::2,b::2]) for b in range(2)] for a in range(2)])
    expected=baseline*np.eye(2)/2+Kt[:2,:2].T/2
    partner_error=float(np.linalg.norm(observed-expected))
    assert max(covariance_error,blind,source_error,partner_error)<5e-13
    return dict(flavor_number=n,CAR_mode_count=m,
      covariance_residual=covariance_error,initial_receiver_full_algebra_blind_residual=blind,
      true_double_occupation=occupation_product,Gaussian_surrogate_double_occupation=float(np.linalg.det(rho).real),
      source_affine_matrix_residual=source_error,source_partner_matrix_residual=partner_error,
      source_logical_traceless_norm=float(np.linalg.norm(Kt[:2,:2]-np.trace(Kt[:2,:2])*np.eye(2)/2,2)),
      actual_reference_is_own_negative_spectrum=True,
      source_matrix_may_correlate_with_partner=True,
      finite_CAR_diagnostic_does_not_replace_continuum_source_domain_proof=True)

def pack(x):
    if isinstance(x,F): return {"numerator":x.numerator,"denominator":x.denominator,"decimal":float(x)}
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [pack(v) for v in x]
    if isinstance(x,np.generic): return x.item()
    return x
def calculate():
    return pack(dict(status="PASS_WITH_SCOPE",independent_scientific_groups=1,
        exact_rational=structures(),fourier_diagnostics=fourier_diagnostics(),
        finite_CAR=finite_car_diagnostic(),
        scope={"quadratic_original_field_models_only":True,"unknown_input_and_passive_partner":True,
               "no_actual_readout_instrument_or_noiseless_recovery_claim":True,
               "no_interacting_SM_or_shared_self_gravitating_solution_claim":True}))
def bindings():
    paths=[HERE/"proof.md",HERE/"check.py",HERE/"results.json",
           HERE.parent/"_shared"/"round1295_admission.md",Path(__file__).resolve()]
    return {p.relative_to(ROOT).as_posix():sha(p) for p in paths}
def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true"); group.add_argument("--check",action="store_true")
    group.add_argument("--diagnose",action="store_true")
    args=parser.parse_args()
    scientific=calculate()
    if args.diagnose:
        print(json.dumps(scientific,ensure_ascii=False,indent=2)); return
    result={"scientific":scientific,"source_bindings":bindings()}
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
        print("CREATED PASS_WITH_SCOPE "+sha(OUT)); return
    saved=json.loads(OUT.read_text(encoding="utf-8"))
    assert saved==result,"READ_ONLY_MISMATCH"
    print("READ_ONLY_MATCH PASS_WITH_SCOPE "+sha(OUT))
if __name__=="__main__": main()
