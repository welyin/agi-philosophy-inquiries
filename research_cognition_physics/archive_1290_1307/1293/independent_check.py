"""1293 independent checks: analytic rational guardrails and finite CAR diagnostics.
No author/historical calculation imports. Default/--check are read-only; --write exclusive.
Finite matrices do not prove continuum Hadamard or interacting-SM statements.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir() and (p/'scripts').is_dir())
OUT=HERE/'independent_results.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rational_guardrails():
    a,L,t0,t1=F(1,32),F(1,4),F(7,16),F(1,2)
    rlo,rhi=L-2*a,L+2*a
    zlo,zhi=t0*t0-rhi*rhi,t1*t1-rlo*rlo
    assert rhi<t0 and zlo==F(3,32) and zhi==F(55,256)<F(1,4)
    # J1(s)/s = sum (-1)^j z^j/[2^(2j+1)j!(j+1)!].
    # Successive absolute-term ratio <= zhi/8 <1 at every j>=0.
    bessel_lo=F(1,2)-zhi/16
    assert bessel_lo==F(1993,4096)>F(31,64)
    ratio0=zhi/8
    assert ratio0<1
    # b>1/4 on radius-a/2 ball, sup b<1/2. pi>3 => L1-normalized norm^2<4/a^3.
    density_norm2_upper=4/a**3
    assert density_norm2_upper==2**17
    product_L1_lower=1/density_norm2_upper
    kernel_lower_without_pi=F(31,256)
    # pi<4, strict inequalities inherited from bump and Bessel steps.
    amplitude_lower=kernel_lower_without_pi/F(4)*product_L1_lower
    probability_lower=amplitude_lower**2
    assert amplitude_lower==F(31,2**27)
    assert probability_lower==F(961,2**54)
    # exp(4/3)<4: Taylor through n=12 plus a geometric upper tail.
    x=F(4,3)
    terms=[F(1)]
    for j in range(1,14): terms.append(terms[-1]*x/j)
    exp_upper=sum(terms[:13])+terms[13]/(1-x/14)
    assert exp_upper<4
    # |grad b|<2/a follows max_{z>=1} z^2 exp(-z)=4/e^2<1.
    # Inner half-ball supplies int b^2>pi*a^3/96.
    gradient_norm_squared_upper=F(512)/a**2
    energy_squared_upper=4*(1+gradient_norm_squared_upper)
    assert gradient_norm_squared_upper==524288 and energy_squared_upper<1450**2
    # exp(1)>2 gives sup(b)<1/2, without decimal transcendental assumptions.
    assert sum(F(1,math.factorial(j)) for j in range(3))>2
    return dict(a=a,L=L,time_window=[t0,t1],distance_window=[rlo,rhi],
        timelike_s_squared=[zlo,zhi],bessel_lower=bessel_lo,
        alternating_ratio_upper=ratio0,bump_exp_upper=exp_upper,
        L1_density_L2_squared_strict_upper=density_norm2_upper,
        normalized_bump_L1_product_strict_lower=product_L1_lower,
        amplitude_strict_lower=amplitude_lower,probability_strict_lower=probability_lower,
        normalized_gradient_squared_strict_upper=gradient_norm_squared_upper,
        relative_energy_squared_strict_upper=energy_squared_upper,relative_energy_strict_upper=1450,
        continuous_proof_required=True)

def continuum_diagnostics():
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    eye=np.eye(2)
    maximum=0.
    samples=[]
    for mass,p,t in [(1.,(.17,-.21,.31),.44),(1.,(-.4,.3,.2),.49),(0.,(.3,.4,.1),.47)]:
        hp=p[0]*sx+p[1]*sy+p[2]*sz
        H=np.block([[-hp,mass*eye],[mass*eye,hp]])
        E=math.sqrt(mass*mass+sum(x*x for x in p))
        ev,vec=np.linalg.eigh(H)
        U=(vec*np.exp(-1j*t*ev))@vec.conj().T
        formula=math.cos(t*E)*np.eye(4)-1j*math.sin(t*E)/E*H
        err=max(np.linalg.norm(H@H-E*E*np.eye(4)),np.linalg.norm(U-formula),
                np.linalg.norm(U[:2,2:]+1j*mass*math.sin(t*E)/E*eye))
        maximum=max(maximum,float(err))
        samples.append(dict(mass=mass,momentum=p,time=t,residual=float(err)))
    # Radial bump integrals are diagnostics only; strict bound is in rational_guardrails.
    bump=[]
    a=1/32
    for n in (48,96):
        z,w=np.polynomial.legendre.leggauss(n)
        r=(z+1)/2; w=w/2
        b=np.exp(-1/(1-r*r))
        i1=4*math.pi*a**3*float(np.dot(w,r*r*b))
        i2=4*math.pi*a**3*float(np.dot(w,r*r*b*b))
        bump.append(dict(order=n,L1=i1,L2_squared=i2,
                         L1_normalized_L2_squared=i2/i1**2,
                         L2_normalized_L1_squared=i1*i1/i2))
    assert maximum<1e-12
    assert all(d['L1_normalized_L2_squared']<2**17 for d in bump)
    return dict(fourier_dirac=samples,maximum_residual=maximum,bump_radial_diagnostics=bump,
                quadrature_is_not_certified_tail_or_PDE_error=True)

def car(n,j):
    out=np.zeros((2**n,2**n),complex)
    for word in range(2**n):
        if (word>>j)&1:
            out[word^(1<<j),word]=(-1)**((word&((1<<j)-1)).bit_count())
    return out

def second_quantized(K,c):
    return sum((K[i,j]*(c[i].conj().T@c[j]) for i in range(len(c)) for j in range(len(c))),
               np.zeros_like(c[0]))

def unitary(H,t):
    e,v=np.linalg.eigh(H)
    return (v*np.exp(-1j*t*e))@v.conj().T

def finite_car():
    # Four actual CAR modes: two flavor copies at input f and at its orthogonal complement.
    # Index order (f0,f1,e0,e1), Fock bits little-endian; no author ordering/code reused.
    c=[car(4,j) for j in range(4)]; dim=16
    ident=np.eye(dim); parity=np.diag([(-1)**i.bit_count() for i in range(dim)])
    residuals={}
    residuals['CAR']=max(float(np.linalg.norm(c[i]@c[j].conj().T+c[j].conj().T@c[i]-(ident if i==j else 0)))
                         for i in range(4) for j in range(4))
    # A correlated faithful, gauge-invariant Gaussian reference. Reset need not assume its factorization.
    hspace=np.array([[.7,.31-.19j],[.31+.19j,1.3]],complex)
    h=np.kron(hspace,np.eye(2))
    H=second_quantized(h,c)
    e,v=np.linalg.eigh(H); weights=np.exp(-e); weights/=weights.sum()
    old=(v*weights)@v.conj().T
    sigma=np.zeros((4,4),complex)
    for i in range(4):
        for j in range(4): sigma[i,j]=sum(old[4*i+a,4*j+a] for a in range(4))
    embeddings=[]
    for occ in (1,2):
        E=np.zeros((16,4),complex)
        for env in range(4): E[4*env+occ,env]=1
        embeddings.append(E)
    E0,E1=embeddings
    # The maximally mixed single-occupation input; reset is exact finite CP and parity covariant.
    reset=np.zeros_like(old); completeness=np.zeros_like(old)
    for E in embeddings:
        for a in range(4):
            extraction=np.zeros((4,16),complex)
            for env in range(4): extraction[env,4*env+a]=1
            K=E@extraction/math.sqrt(2)
            reset+=K@old@K.conj().T; completeness+=K.conj().T@K
    expected=(E0@sigma@E0.conj().T+E1@sigma@E1.conj().T)/2
    residuals['reset_completeness']=float(np.linalg.norm(completeness-ident))
    residuals['reset_state']=float(np.linalg.norm(reset-expected))
    residuals['reset_parity']=float(np.linalg.norm(reset@parity-parity@reset))
    correlation_removed=float(np.linalg.norm(old-reset))
    # Actual single occupancy is non-Gaussian even when its one-body matrix is I/2.
    occupation_product=float(np.trace((c[0].conj().T@c[0])@(c[1].conj().T@c[1])@reset).real)
    assert abs(occupation_product)<1e-14
    # Receiver's full two-mode algebra, implemented with its actual CAR generators.
    receiver_words=[ident]
    for op in (c[2],c[2].conj().T,c[3],c[3].conj().T):
        receiver_words += [x@op for x in list(receiver_words)]
    blind=0.
    for O in receiver_words:
        common=E0.conj().T@O@E0
        blind=max(blind,float(np.linalg.norm(E1.conj().T@O@E1-common)),
                  float(np.linalg.norm(E0.conj().T@O@E1)),float(np.linalg.norm(E1.conj().T@O@E0)))
    residuals['receiver_code_matrix_identity']=blind
    # Unknown reference diagnostic. Operator block identity above, not this sample, gives arbitrary R.
    rhoQR=np.array([[.37,0,0,.12+.08j],[0,.13,.04j,0],[0,-.04j,.19,0],[.12-.08j,0,0,.31]],complex)
    assert np.linalg.eigvalsh(rhoQR).min()>0
    joint=np.zeros((32,32),complex)
    for i,Ei in enumerate(embeddings):
        for j,Ej in enumerate(embeddings): joint+=np.kron(Ei@sigma@Ej.conj().T,rhoQR[2*i:2*i+2,2*j:2*j+2])
    rhoR=rhoQR[:2,:2]+rhoQR[2:,2:]
    partner=0.
    for O in receiver_words:
        reduced=np.zeros((2,2),complex)
        for a in range(2):
            for b in range(2): reduced[a,b]=np.trace(O@joint[a::2,b::2])
        target=np.trace(E0.conj().T@O@E0@sigma)*rhoR
        partner=max(partner,float(np.linalg.norm(reduced-target)))
    residuals['initial_receiver_R_factorization']=partner
    # Common flavor-blind propagation. Probability difference equals actual one-particle overlap squared.
    Gspace=np.array([[0,-1j],[1j,0]],complex)
    G=np.kron(Gspace,np.eye(2)); angle=math.pi/6
    U=unitary(second_quantized(G,c),angle); u=unitary(G,angle)
    n2=c[2].conj().T@c[2]
    states=[E@sigma@E.conj().T for E in embeddings]
    probs=[float(np.trace(n2@U@rho@U.conj().T).real) for rho in states]
    alpha2=float(abs(u[2,0])**2)
    residuals['probability_transfer']=abs((probs[0]-probs[1])-alpha2)
    residuals['full_probability_normalization']=max(abs(float(np.trace(U@rho@U.conj().T).real)-1) for rho in states)
    # All flavor-singlet quadratic jets reduce to one arbitrary spatial-spin matrix K.
    source_errors=[]
    Klist=[np.eye(2),np.array([[2,1+.3j],[1-.3j,-.4]]),np.array([[0,1j],[-1j,0]]),hspace]
    for Kspace in Klist:
        O=U.conj().T@second_quantized(np.kron(Kspace,np.eye(2)),c)@U
        x=E0.conj().T@O@E0
        source_errors.append(max(float(np.linalg.norm(E1.conj().T@O@E1-x)),
            float(np.linalg.norm(E0.conj().T@O@E1)),float(np.linalg.norm(E1.conj().T@O@E0))))
    residuals['all_tested_source_code_scalar']=max(source_errors)
    # A flavor spurion is deliberately not a permitted background singlet source.
    spurion=c[0].conj().T@c[0]
    spurion_difference=float(np.trace(spurion@(states[0]-states[1])).real)
    assert abs(spurion_difference-1)<1e-12
    assert max(residuals.values())<1e-11
    assert correlation_removed>0.1 and alpha2>0
    return dict(residuals=residuals,reference_reset_distance_diagnostic=correlation_removed,
        receiver_probabilities=probs,overlap_squared=alpha2,
        flavor_spurion_expectation_difference=spurion_difference,
        actual_input_occupation_product=occupation_product,Gaussian_substitute_occupation_product=0.25,
        finite_modes=4,Fock_dimension=16,source_matrices=len(Klist),receiver_words=len(receiver_words),
        arbitrary_R_from_operator_identity_not_samples=True,
        reset_is_initial_state_contract_not_autonomous_preparation=True,
        source_scalar_does_not_mean_higher_moments_or_conditional_sources_equal=True)

def energy_identity_diagnostics():
    # Independent two opposite Fourier cells: real-even pure-R input, both flavors.
    # These finite matrices check the algebra of the reset energy, not a continuum cutoff.
    pauli=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.])]
    rows=[]
    for momentum in ((.23,-.17,.41),(-.31,.27,.09)):
        sp=sum(p*s for p,s in zip(momentum,pauli))
        hp=np.block([[-sp,np.eye(2)],[np.eye(2),sp]])
        hm=np.block([[sp,np.eye(2)],[np.eye(2),-sp]])
        hs=np.block([[hp,np.zeros((4,4))],[np.zeros((4,4)),hm]])
        H=np.kron(hs,np.eye(2)); ev,vec=np.linalg.eigh(H)
        P=(vec*(ev<0))@vec.conj().T
        absolute=(vec*np.abs(ev))@vec.conj().T
        wave=np.zeros(8,complex);wave[2]=wave[6]=1/math.sqrt(2)
        V=np.stack([np.kron(wave,np.eye(2)[:,j]) for j in range(2)],axis=1)
        Fmat=V@V.conj().T;Q=np.eye(16)-Fmat
        target=float(np.trace(Fmat@absolute).real)
        errs=[float(np.linalg.norm(V.conj().T@H@V))]
        energies=[]
        for rho in (np.diag([1.,0.]),np.diag([0.,1.]),np.array([[.6,.2+.1j],[.2-.1j,.4]])):
            delta=Q@P@Q-P+V@rho@V.conj().T
            energy=float(np.trace(H@delta).real);energies.append(energy)
            errs.append(abs(energy-target))
        analytic=2*math.sqrt(1+sum(p*p for p in momentum))
        errs.append(abs(target-analytic))
        assert max(errs)<1e-12
        rows.append(dict(momentum=momentum,energies=energies,target=target,max_residual=max(errs)))
    return dict(samples=rows,finite_algebra_only=True)

def yukawa_spurion_diagnostic():
    # A genuine complex Yukawa variation, rather than a background Higgs/gauge variation.
    # For H_mass=[[0,Y],[Y^dagger,0]], delta Y=i diag(1,-1) gives -sigma_y tensor Z.
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    z=np.diag([1.,-1.]); t=.37
    H=np.kron(sx,np.eye(2)); J=np.vstack([np.zeros((2,2)),np.eye(2)])
    V=unitary(H,t)@J
    spurion=-np.kron(sy,z)
    compressed=V.conj().T@spurion@V
    predicted=-math.sin(2*t)*z
    error=float(np.linalg.norm(compressed-predicted))
    scalar_defect=float(np.linalg.norm(compressed-np.trace(compressed)/2*np.eye(2)))
    assert error<1e-12 and scalar_defect>.5
    # The corresponding common Higgs phase vertex has identity flavor and stays scalar.
    common=V.conj().T@(-np.kron(sy,np.eye(2)))@V
    common_error=float(np.linalg.norm(common+math.sin(2*t)*np.eye(2)))
    assert common_error<1e-12
    return dict(time=t,complex_Yukawa_scalar_defect=scalar_defect,residual=error,
        common_Higgs_phase_residual=common_error,
        arbitrary_coupling_variations_are_outside_source_menu=True)

def norm(v):
    if isinstance(v,F): return {'numerator':v.numerator,'denominator':v.denominator}
    if isinstance(v,dict): return {k:norm(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [norm(x) for x in v]
    if isinstance(v,(np.floating,np.integer)): return v.item()
    return v

def calculate():
    science=norm(dict(rational=rational_guardrails(),continuum_diagnostics=continuum_diagnostics(),finite_CAR=finite_car(),energy_identity=energy_identity_diagnostics(),Yukawa_spurion=yukawa_spurion_diagnostic()))
    bindings={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [
        Path(__file__).resolve(),HERE/'proof.md',HERE/'check.py',HERE/'results.json',
        HERE.parent/'_shared'/'round1293_admission.md']}
    return dict(round=1293,status='PASS_WITH_SCOPE',independent=True,
        scientific_payload_sha256=hashlib.sha256(json.dumps(science,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode()).hexdigest(),
        source_bindings=bindings,science=science,
        scope=['exact quadratic Dirac/CAR object only','continuum support and unknown-reference quantifiers are analytic',
               'finite diagnostics do not certify full interacting SM, devices, self-gravity or all source moments'])

def main():
    parser=argparse.ArgumentParser(); group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=calculate()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f: json.dump(result,f,ensure_ascii=False,sort_keys=True,indent=2);f.write('\n')
        print('CREATED independent_results.json')
    else:
        previous=json.loads(OUT.read_text(encoding='utf-8'))
        assert result==previous,'independent results or source binding changed'
        print('READ_ONLY_MATCH')
    print(json.dumps({'status':result['status'],'scientific_payload_sha256':result['scientific_payload_sha256'],
        'largest_CAR_residual':max(result['science']['finite_CAR']['residuals'].values())},sort_keys=True))

if __name__=='__main__': main()
