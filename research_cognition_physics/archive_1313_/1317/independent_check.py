"""Independent checks of round 1317. No author module is loaded.

The all-window diagnostic uses an analytic derivative bound between grid points;
the analytic theorem, not the floating-point grid, supplies the universal claim.
Default/--check read only; --write-new exclusively creates this agent's result.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
ADMISSION = HERE.parent / "_admission" / "round1317_admission.md"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def herm(a):
    return (a+a.conj().T)/2


def norm(a):
    return float(np.linalg.norm(a, 2))


class Audit:
    def __init__(self):
        self.errors, self.count = {}, 0

    def close(self, label, a, b, tol=2e-9):
        a, b = np.atleast_2d(a), np.atleast_2d(b)
        err = norm(a-b)/max(1., norm(b))
        self.errors[label] = max(self.errors.get(label, 0.), err)
        self.count += 1
        if err > tol:
            raise AssertionError((label, err, tol))

    def truth(self, label, assertion):
        self.count += 1
        self.errors.setdefault(label, 0.)
        if not bool(assertion):
            raise AssertionError(label)


def qmatrix(N, q):
    mat = np.zeros((N, N+1))
    for j in range(N):
        mat[j, j] = 1.
        mat[j, j+1] = -q
    return mat


def basis(N, q):
    v = q**(-np.arange(N+1, dtype=float))
    return v/np.linalg.norm(v)


def mass_and_hamiltonian(N, n, m, q, rho, Z, cut=False):
    chain = m*qmatrix(N, q)
    if cut:
        chain[-1] = 0.
    t = np.kron(chain, np.eye(n))
    M = np.vstack((rho*Z, t))
    zero = np.zeros_like(M)
    H = np.block([[zero, M.conj().T], [M, zero]])
    e0 = np.eye(len(H), dtype=complex)[:, :n]
    el = np.eye(len(H), dtype=complex)[:, M.shape[1]:M.shape[1]+n]
    return M, H, e0, el


def evolve(H, t):
    w, v = np.linalg.eigh(H)
    return (v*np.exp(-1j*t*w))@v.conj().T


def beta_spectrum(N, m, q):
    chain = m*qmatrix(N, q)
    chain[-1] = 0.
    B = chain.T@chain
    w, v = np.linalg.eigh(B[:N, :N])
    w[np.abs(w) < 1e-12] = 0.
    return np.sqrt(np.maximum(0., w)), v, B


def gramian(N, m, q, T):
    omega, v, B = beta_spectrum(N,m,q)
    integ = np.zeros((N,N))
    for i, a in enumerate(omega):
        for j, b in enumerate(omega):
            if a == 0 and b == 0:
                integ[i,j] = T**3/3
            elif a == 0 or b == 0:
                s = max(a,b)
                integ[i,j] = (np.sin(s*T)-s*T*np.cos(s*T))/s**3
            else:
                integ[i,j] = T*(np.sinc((a-b)*T/np.pi)-np.sinc((a+b)*T/np.pi))/(2*a*b)
    coeff = v[0,:]
    W = (v*coeff)@integ@(v*coeff).T
    return herm(W), omega, v, B


def beta_at(t, omega, v):
    factors = t*np.sinc(omega*t/np.pi)
    return (v*factors)@v[0,:]


def response_kernel(N, m, q):
    t = m*qmatrix(N,q)
    G = t@t.T
    inv = np.linalg.inv(G)
    return t.T@inv@inv@t, float(np.linalg.eigvalsh(G)[0])


def random_unitary(rng,n):
    a = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    return np.linalg.qr(a)[0]


def finite_case(n, rng, audit):
    # Fixed audit parameters, registered before computing the response.
    N,m,q,rho,T,Zmax = 2,1.2,1.25,0.4,3.,0.5
    blocks = [0.0002*random_unitary(rng,n),
              0.0003*random_unitary(rng,n),
              random_unitary(rng,n)@np.diag(np.linspace(.25,.35,n))@random_unitary(rng,n)]
    Z = np.hstack(blocks)
    _,H,e0,el = mass_and_hamiltonian(N,n,m,q,rho,Z,cut=True)
    W,omega,v,B = gramian(N,m,q,T)
    lmin = float(np.linalg.eigvalsh(W)[0])
    audit.truth("Gramian_positive", lmin>0)
    audit.truth("registered_portal_bound", norm(Z)<Zmax)
    nodes,weights = np.polynomial.legendre.leggauss(96)
    Wnum = np.zeros_like(W)
    for node,weight in zip(nodes,weights):
        b=beta_at(T*(node+1)/2,omega,v)
        Wnum += T/2*weight*np.outer(b,b)
    audit.close("analytic_Gramian_vs_quadrature",W,Wnum)

    # Eigen-expansion of the actual Hermitian H, with a rigorous derivative
    # envelope ||A'(t)|| <= rho*Zmax between adjacent sampled times.
    w,uu=np.linalg.eigh(H)
    left,right=el.conj().T@uu,uu.conj().T@e0
    times=np.linspace(0,T,5001)
    amps=np.einsum('ak,tk,kb->tab',left,np.exp(-1j*np.outer(times,w)),right,optimize=True)
    sampled_max=float(np.max(np.linalg.svd(amps,compute_uv=False)[:,0]))
    step=T/(len(times)-1)
    epsilon=sampled_max+rho*Zmax*step/2+2e-12
    R=1+rho*rho*Zmax*Zmax*T*T/2
    entry_bound=epsilon/rho*np.sqrt(T/lmin)*R
    pre=Z[:,:N*n]
    audit.truth("finite_refusal_to_entry_bound",norm(pre)<=entry_bound+1e-10)
    audit.truth("bound_is_not_vacuous_for_audit_case",entry_bound<Zmax/10)

    u=basis(N,q)
    V0=np.kron(u.reshape(-1,1),np.eye(n))
    Y=Z@V0
    Zref=np.hstack([np.zeros((n,n)) for _ in range(N)]+[Y/u[-1]])
    delta=norm(Z-Zref)
    dmax=entry_bound/u[-1]
    audit.close("same_light_Y",Zref@V0,Y)
    audit.truth("delta_sharp_norm_bound",delta<=norm(pre)/u[-1]+1e-10)
    kernel,gap=response_kernel(N,m,q)
    K=np.kron(kernel,np.eye(n))
    C=herm(Z@K@Z.conj().T)
    Cref=herm(Zref@K@Zref.conj().T)
    cbound=(2*Zmax*dmax+dmax*dmax)/gap
    audit.truth("same_Y_heavy_response_bound",norm(C-Cref)<=cbound+1e-10)
    audit.truth("response_bound_smaller_than_reference",cbound<norm(Cref))

    # Direct Volterra identity in the full matrix order.
    tfinal=.83*T
    actual=el.conj().T@evolve(H,tfinal)@e0
    integral=np.zeros((n,n),complex)
    wb,vb=np.linalg.eigh(B)
    wb[np.abs(wb)<1e-12]=0.
    sb=np.sqrt(np.maximum(0.,wb))
    for node,weight in zip(nodes,weights):
        s=tfinal*(node+1)/2
        st=(vb*((tfinal-s)*np.sinc(sb*(tfinal-s)/np.pi)))@vb.T
        aa=el.conj().T@evolve(H,s)@e0
        integral += tfinal/2*weight*(Z@np.kron(st,np.eye(n))@Z.conj().T@aa)
    bfull=np.r_[beta_at(tfinal,omega,v),0.]
    audit.close("noncommuting_Volterra_order",rho*Z@np.kron(bfull[:,None],np.eye(n)),1j*(actual+rho*rho*integral))

    # Units: energy rescaling by k and time by 1/k leaves the entry bound.
    scale=2.7
    W2,_,_,_=gramian(N,scale*m,q,T/scale)
    audit.close("Gramian_time_units",W2,W/scale**3)
    eb2=epsilon/(rho*scale)*np.sqrt((T/scale)/np.linalg.eigvalsh(W2)[0])*R
    audit.close("entry_bound_dimensionless",np.array([eb2]),np.array([entry_bound]))
    K2,_=response_kernel(N,scale*m,q)
    audit.close("response_inverse_energy_squared",K2,kernel/scale**2)

    # Complete success/failure instrument with an arbitrary maximally entangled
    # partner; keeping all output modes preserves the whole input isometry.
    U=evolve(H,T)
    isom=U@e0
    P=el@el.conj().T
    success=P@isom
    failure=(np.eye(len(H))-P)@isom
    audit.close("instrument_complete",success.conj().T@success+failure.conj().T@failure,np.eye(n))
    bell=isom/np.sqrt(n)
    audit.close("unvisited_partner_preserved",bell.conj().T@bell,np.eye(n)/n)
    audit.truth("worst_input_probability_certificate",norm(el.conj().T@isom)**2<=epsilon**2+1e-10)
    _,Href,_,_=mass_and_hamiltonian(N,n,m,q,rho,Zref,cut=True)
    for t in [0.,.2,T]:
        audit.close("exact_endpoint_refusal_tool",el.conj().T@evolve(Href,t)@e0,np.zeros((n,n)))
    return dict(n=n,N=N,m=m,q=q,rho=rho,T=T,Zmax=Zmax,
                minimum_Gramian_eigenvalue=lmin,sampled_amplitude_max=sampled_max,
                whole_window_amplitude_upper=epsilon,probability_upper=epsilon**2,
                actual_nonterminal_norm=norm(pre),certified_nonterminal_norm=entry_bound,
                actual_portal_delta=delta,certified_portal_delta=dmax,
                heavy_response_difference=norm(C-Cref),heavy_response_bound=cbound,
                reference_response_norm=norm(Cref))


def strong_family(audit):
    N,n,m,q,rho,Y=2,1,1.,1.4,.7,.3
    u=basis(N,q)
    rows=[]
    kernel,_=response_kernel(N,m,q)
    for strength in [10.,30.,100.,300.]:
        a=strength
        b=Y/u[-1]-q*strength
        Z=np.array([[0.,a,b]])
        M,H,e0,el=mass_and_hamiltonian(N,n,m,q,rho,Z,cut=True)
        G=M[:2]@M[:2].T
        lam=np.linalg.eigvalsh(G)
        bound=rho*abs(a)*m*m*q/(lam[1]-lam[0])*(1/np.sqrt(lam[1])+1/np.sqrt(lam[0]))
        audit.close("strong_family_fixed_Y",Z@u[:,None],np.array([[Y]]))
        maxamp=0.
        for t in np.linspace(0,5,103):
            amp=(el.conj().T@evolve(H,t)@e0)[0,0]
            f=(np.sin(t*np.sqrt(lam[1]))/np.sqrt(lam[1])-np.sin(t*np.sqrt(lam[0]))/np.sqrt(lam[0]))/(lam[1]-lam[0])
            exact=-1j*m*G[0,1]*f
            audit.close("strong_exact_two_mode_formula",np.array([amp]),np.array([exact]),2e-8)
            maxamp=max(maxamp,abs(amp))
        audit.truth("strong_all_time_envelope",maxamp<=bound+1e-9)
        rows.append(dict(strength=strength,light_Y=Y,nonterminal_norm=norm(Z[:,:2]),
                         total_portal_norm=norm(Z),all_time_amplitude_bound=float(bound),
                         heavy_response=float((Z@kernel@Z.T)[0,0])))
    audit.truth("strong_hiding_envelope_decreases",all(rows[j+1]['all_time_amplitude_bound']<rows[j]['all_time_amplitude_bound'] for j in range(3)))
    audit.truth("strong_hiding_response_grows",all(rows[j+1]['heavy_response']>rows[j]['heavy_response'] for j in range(3)))
    # General AND control is a changed control-scope contract, not a counterexample
    # to the fixed-portal theorem: each refusal sets every portal to zero.
    for flags in [(0,0),(0,1),(1,0)]:
        chain=m*qmatrix(N,q)
        for j,c in enumerate(flags):
            chain[j]*=c
        M=np.vstack((np.prod(flags)*rho*Z,chain))
        H=np.block([[np.zeros_like(M),M.T],[M,np.zeros_like(M)]])
        audit.close("AND_exactly_decouples_receiver",el.conj().T@H,np.zeros((1,len(H))))
        audit.close("AND_nonterminal_still_registered",Z[:,:2],np.array([[0.,300.]]))
    return rows


def car_lift(h):
    modes=len(h)
    size=1<<modes
    out=np.zeros((size,size),complex)
    for state in range(size):
        for j in range(modes):
            if not (state>>j)&1:
                continue
            intermediate=state^(1<<j)
            signj=(-1)**((state&((1<<j)-1)).bit_count())
            for i in range(modes):
                if (intermediate>>i)&1 or h[i,j]==0:
                    continue
                final=intermediate|(1<<i)
                signi=(-1)**((intermediate&((1<<i)-1)).bit_count())
                out[final,state]+=signi*signj*h[i,j]
    return out


def car_and_record(audit):
    N,n,m,q,rho,Zmax=1,2,1.,1.2,.4,.6
    Z=np.hstack((np.array([[.1,.03j],[.02,.12]]),np.array([[.3,.02],[-.01j,.35]])))
    _,h_on,_,_=mass_and_hamiltonian(N,n,m,q,rho,Z,False)
    _,h_off,_,_=mass_and_hamiltonian(N,n,m,q,rho,Z,True)
    audit.truth("CAR_registered_portal_bound",norm(Z)<Zmax)
    lifted=[]
    modes=len(h_on)
    indices=[1<<j for j in range(modes)]
    number=np.diag([s.bit_count() for s in range(1<<modes)])
    # Same constant for both permission sectors, chosen from declared bounds;
    # it is a full-Fock scalar, not a number-dependent chemical potential.
    shift=modes*(m*np.sqrt(1+q*q)+rho*Zmax)
    for label,h in [('yes',h_on),('no',h_off)]:
        fock=car_lift(h)
        audit.close("CAR_Hermitian_"+label,fock,fock.conj().T)
        audit.close("CAR_one_particle_dictionary_"+label,fock[np.ix_(indices,indices)],h)
        audit.close("CAR_number_conservation_"+label,fock@number,number@fock)
        audit.truth("CAR_common_scalar_offset_bounded_below_"+label,np.linalg.eigvalsh(fock)[0]+shift>=0)
        lifted.append(fock)
    # An explicit source identity for the changed chain permissions.
    audit.close("CAR_same_source_lift",lifted[0]-lifted[1],car_lift(h_on-h_off))
    # The block control is diagonal in the prepared permission record, so its
    # commutator with that record projector vanishes for every Fock input.
    z=np.zeros_like(lifted[0])
    total=np.block([[lifted[1],z],[z,lifted[0]]])
    permit=np.diag(np.r_[np.zeros(len(z)),np.ones(len(z))])
    audit.close("permission_record_preserved",total@permit,permit@total)
    return dict(modes=modes,Fock_dimension=1<<modes,permission_Fock_dimension=len(total),
                common_full_Fock_scalar_offset=shift,
                scope="Finite CAR bilinear; not a Lorentz mass or full-field implementation certificate")


def rational_gramian_certificate(audit):
    # Independent exact integral route for N=2, q=m=1, T=2. In the eigenbasis
    # of B the entries are 4/3, (sin(2sqrt2)/sqrt2-2cos(2sqrt2))/4,
    # and 1/4-sin(4sqrt2)/(16sqrt2). Bound each trig value rationally.
    order=18
    def series(power, sine):
        center=sum((Fraction((-power)**k,math.factorial(2*k+1 if sine else 2*k))
                    for k in range(order+1)),Fraction(0))
        denom=math.factorial(2*order+3 if sine else 2*order+2)
        radius=Fraction(power**(order+1),denom)
        audit.truth("alternating_tail_decreases",power<(2*order+2)*(2*order+3))
        return center,radius
    s8,es8=series(8,True)
    c8,ec8=series(8,False)
    s32,es32=series(32,True)
    off=(s8-c8)/2
    eoff=(es8+ec8)/2
    diagonal=Fraction(1,4)-s32/4
    ediag=es32/4
    a=Fraction(4,3)-Fraction(1,20)
    bmin=diagonal-ediag-Fraction(1,20)
    cabs=abs(off)+eoff
    determinant_lower=a*bmin-cabs*cabs
    audit.truth("exact_rational_Gramian_diagonal_margin",a>0 and bmin>0)
    audit.truth("exact_rational_Gramian_determinant_margin",determinant_lower>0)
    W,*_=gramian(2,1.,1.,2.)
    audit.truth("numeric_Gramian_above_certified_lower",np.linalg.eigvalsh(W)[0]>.05)
    return dict(N=2,m=1,q=1,T=2,certified_lambda_min_lower="1/20",
                order=order,exact_determinant_lower=str(determinant_lower),
                determinant_lower_decimal=float(determinant_lower),
                route="Closed spectral integrals and independent rational alternating-series enclosures")


def run():
    audit=Audit()
    rng=np.random.default_rng(131739)
    cases=[finite_case(n,rng,audit) for n in [2,3,4,5]]
    strong=strong_family(audit)
    car=car_and_record(audit)
    rational=rational_gramian_certificate(audit)
    return dict(status="PASS_WITH_SCOPE",round=1317,seed=131739,
                independent_of_author_code=True,assertion_count=audit.count,
                check_categories=len(audit.errors),max_normalized_error=max(audit.errors.values()),
                errors=audit.errors,finite_cases=cases,strong_fixed_Y_family=strong,CAR=car,
                rational_Gramian_certificate=rational,
                bindings={"independent_check.py":sha(__file__),"admission":sha(ADMISSION)},
                scope="Analytic finite-window inference with registered portal bound; no general consensus axiom, autonomous withdrawal or full QFT certification")


def main():
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--write-new',action='store_true')
    args=ap.parse_args()
    result=run()
    if args.write_new:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
            f.write('\n')
    elif OUT.exists():
        old=json.loads(OUT.read_text(encoding='utf-8'))
        if old!=result:
            raise AssertionError('Independent saved result differs from recomputation')
    print(json.dumps({key:result[key] for key in ['status','assertion_count','check_categories','max_normalized_error']},ensure_ascii=False))


if __name__=='__main__':
    main()
