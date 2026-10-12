"""Round 1141: exact budgets and finite, noncommuting column diagnostics.
This is not a continuous-field or actual BT simulation.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESULT=HERE/"results.json"

def norm(x): return float(np.linalg.norm(x))
def compare(a,b,p="root"):
    if isinstance(a,dict):
        assert set(a)==set(b),(p,"keys")
        for k in a: compare(a[k],b[k],p+"."+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(p,"length")
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,p+f"[{i}]")
    elif isinstance(a,float):
        assert np.isclose(a,b,rtol=1e-8,atol=3e-9),(p,a,b)
    else: assert a==b,(p,a,b)
def frac(x):
    return dict(numerator=str(x.numerator),denominator=str(x.denominator),value=float(x))
def fingerprint_group():
    names=["1129/proof.md","1133/proof.md","1137/proof.md",
           "_shared/after1140_direction_decision.md"]
    return dict(status="PASS",sha256={n:hashlib.sha256((STAGE/n).read_bytes()).hexdigest() for n in names},
       preserved=["original epsilon/A/Y/O/Z/source","full normal CAR and Q/R","all quantum receiver/hold tails"],
       new_cognitive_axioms=0,scientific_calibrations=0,
       actual_BT_migration=False,complete_local_energy_balance=False,
       manufactured_energy_instrument=False)
def budget_group():
    eps,dc,A,t=F(1,10**5),F(1,10**6),F(10**34),F(5,4)
    lam,kap=F(11,7),F(11,14)
    aa=lam*512*64/eps**3
    bb=42*8/dc+(lam+kap)*512
    D=23*aa+3*bb
    assert F(524)<23**2 and F(8)<3**2
    # nu <= ||rho||_2 ||grad rho||_2; squared constant < 100^2.
    nu_squared_bound=F(8)*16**2*F(88,21)
    assert nu_squared_bound<100**2
    assert 524*100<230**2
    assert F(2146)<47**2
    Cw=47*A+230/eps**2
    d=t*D/A
    energy=Cw*d+d*d/2
    qforce=D/A
    work=energy
    flow=16*(Cw*qforce*t*t/2+qforce*qforce*t*t*t/6)
    assert work<F(71,10)*10**22<F(10**23)
    assert flow<F(7)*10**23<F(10**24)
    assert D<F(12,10)*10**21
    assert energy<F(71,10)*10**22<F(10**23)
    return dict(status="PASS",epsilon=frac(eps),A=frac(A),t_max=frac(t),
        handoff_density_bound=frac(aa),receiver_derivative_bound=frac(bb),
        force_column_D=frac(D),relative_wave_column=frac(d),
        free_column_Cw=frac(Cw),relative_energy_bound=frac(energy),
        final_certificate=10**23,nu_squared_constant_bound=frac(nu_squared_bound),
        region_gradient_bound=16,absolute_work_integral_bound=frac(work),
        absolute_relative_flux_integral_bound=frac(flow),
        work_certificate=10**23,relative_flux_certificate=10**24,
        fixed_parameters_changed=False)
def column_and_form_groups():
    # Finite diagnostic only: 17 periodic locations, 2-position support, 2-state material.
    n,dim=17,4
    x=np.arange(n)
    circ=lambda center: np.minimum((x-center)%n,(center-x)%n)
    eta=np.maximum(1-circ(3)/5,0)
    chi=np.maximum(1-circ(3)/3,0)
    kernel=np.zeros(n)
    kernel[0],kernel[1],kernel[-1]=.5,.25,.25
    mult=np.fft.fft(kernel)
    assert np.max(np.abs(mult))<=1+1e-14
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]).astype(complex)
    eye=np.eye(dim,dtype=complex)
    alpha,beta=1.3,.7
    raw=np.zeros((n,dim,dim),complex)
    for y,shift in enumerate((0,6)):
        ey=np.roll(eta,shift)
        cy=np.roll(chi,shift)
        for j in range(n):
            angle=.31*j
            a=np.cos(angle)*sx+np.sin(angle)*sz
            b=np.cos(.27*j)*sy+np.sin(.27*j)*sz
            raw[j,2*y:2*y+2,2*y:2*y+2]=alpha*ey[j]*a+beta*cy[j]*b
    def spatial(a,m):
        return np.fft.ifft(np.fft.fft(a,axis=0)*m[:,None,None],axis=0)
    force=spatial(raw,mult)
    def op_column(a):
        gram=sum((v.conj().T@v for v in a),np.zeros((dim,dim),complex))
        return float(np.sqrt(max(np.linalg.eigvalsh(gram).max(),0)))
    D=alpha*norm(eta)+beta*norm(chi)
    raw_norm,force_norm=op_column(raw),op_column(force)
    assert force_norm<=raw_norm+1e-12 and raw_norm<=D+1e-12
    H=.31*np.kron(sx,np.eye(2))+.47*np.kron(sz,sx)+.29*np.kron(np.eye(2),sz)
    vals,vecs=np.linalg.eigh(H)
    def unitary(t): return (vecs*np.exp(-1j*t*vals))@vecs.conj().T
    k=2*np.pi*np.fft.fftfreq(n)
    omega=np.abs(k)
    t,steps,Atoy=.25,12,1000.
    dt=t/steps
    dpi=np.zeros_like(force)
    dgrad=np.zeros_like(force)
    max_force_norm=0.
    for j in range(steps):
        s=(j+.5)*dt
        u=unitary(s)
        fs=np.einsum("ab,xbc,cd->xad",u.conj().T,force,u)
        max_force_norm=max(max_force_norm,op_column(fs))
        tau=t-s
        dpi-=dt/Atoy*spatial(fs,mult*np.cos(omega*tau))
        dgrad-=dt/Atoy*spatial(fs,mult*1j*np.sign(k)*np.sin(omega*tau))
    delta=np.stack([dpi,dgrad],axis=0)
    delta_column=op_column(delta.reshape(-1,dim,dim))
    column_bound=t*D/Atoy
    assert delta_column<=column_bound+1e-12
    herm=max(norm(a-a.conj().T) for a in delta.reshape(-1,dim,dim))
    assert herm<1e-12
    # A common free history of already energy-smoothed fields; matrices do not commute.
    beta0=np.exp(-(circ(8)/4)**2)
    pi_init=np.array([Atoy*beta0[j]*eye+.3*np.cos(.4*j)*np.kron(sz,sx) for j in range(n)])
    phi_init=np.array([.2*np.sin(.3*j)*np.kron(sx,sy) for j in range(n)])
    pi_free=spatial(pi_init,np.cos(omega*t))-spatial(phi_init,omega*np.sin(omega*t))
    sinc=np.full(n,t)
    nonzero=omega>0
    sinc[nonzero]=np.sin(omega[nonzero]*t)/omega[nonzero]
    phi_free=spatial(phi_init,np.cos(omega*t))+spatial(pi_init,sinc)
    grad_free=spatial(phi_free,1j*k)
    z0=np.stack([pi_free,grad_free],axis=0)
    w=np.maximum(1-circ(5)/6,0)
    weights=np.tile(w,2)
    zflat=z0.reshape(-1,dim,dim)
    dflat=delta.reshape(-1,dim,dim)
    def energy(col):
        return .5*sum((wj*(a.conj().T@a) for wj,a in zip(weights,col)),
                      np.zeros((dim,dim),complex))
    e0=energy(zflat)
    e=energy(zflat+dflat)
    direct=e-e0
    expanded=.5*sum((wj*(a.conj().T@d+d.conj().T@a+d.conj().T@d)
                      for wj,a,d in zip(weights,zflat,dflat)),np.zeros((dim,dim),complex))
    residual=norm(direct-expanded)
    assert residual<3e-9
    psi=np.array([1.,.3j,-.4,.2+.1j],complex)
    psi/=norm(psi)
    free_vec=np.einsum("xab,b->xa",zflat,psi)*np.sqrt(weights)[:,None]
    delta_vec=np.einsum("xab,b->xa",dflat,psi)*np.sqrt(weights)[:,None]
    cross=float(np.vdot(free_vec,delta_vec).real)
    square=norm(delta_vec)**2/2
    expval=float(np.vdot(psi,expanded@psi).real)
    cauchy=norm(free_vec)*column_bound+column_bound**2/2
    assert abs(expval-(cross+square))<1e-10
    assert abs(expval)<=cauchy+1e-12
    noncomm=max(norm(a@d-d@a) for a,d in zip(zflat,dflat))
    assert noncomm>1e-8
    column=dict(status="PASS",scope="finite matrix-valued periodic wave quadrature; not actual fields",
        locations=n,hilbert_dimension=dim,clock_or_actual_BT_simulated=False,
        raw_force_column=raw_norm,smoothed_force_column=force_norm,
        fiber_Young_bound=D,max_Heisenberg_force_column=max_force_norm,
        joint_retarded_column=delta_column,triangle_certificate=column_bound,
        Hermitian_residual=herm,all_bounds_pass=True)
    form=dict(status="PASS",scope="noncommuting finite quadratic forms, same initial vector and free history",
        operator_identity_residual=residual,noncommuting_column_witness=noncomm,
        relative_energy_expectation=expval,cross_term=cross,
        response_square_term=square,Cauchy_certificate=cauchy,
        free_energy_expectation=float(np.vdot(psi,e0@psi).real),
        common_free_term_removed=True,not_a_positive_work_claim=True)
    return column,form
def run():
    col,form=column_and_form_groups()
    return dict(round=1141,status="PASS",groups={
       "frozen_source_and_scope":fingerprint_group(),
       "exact_rational_budget":budget_group(),
       "force_and_retarded_columns":col,
       "noncommuting_relative_energy_form":form})
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    data=run()
    if args.write_results:
        RESULT.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("PASS: wrote round 1141 results, four groups.")
    else:
        compare(data,json.loads(RESULT.read_text(encoding="utf-8")))
        print("PASS: read-only round 1141 replay, four groups.")
    print(json.dumps({name:g["status"] for name,g in data["groups"].items()},ensure_ascii=False))
if __name__=="__main__": main()
