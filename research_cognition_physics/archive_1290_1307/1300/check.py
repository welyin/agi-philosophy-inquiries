"""1300 finite diagnostics; continuum PDE claims are not numerical outputs.
Default / --check recompute read-only; --write creates results.json exclusively.
Uses only standard Python and the existing NumPy runtime.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'


def unitary(h,t):
    v,w=np.linalg.eigh(h)
    return (w*np.exp(-1j*t*v))@w.conj().T


def annihilator(j,n):
    a=np.zeros((1<<n,1<<n),complex)
    for s in range(1<<n):
        if s&(1<<j):
            a[s^(1<<j),s]=(-1)**((s&((1<<j)-1)).bit_count())
    return a


def finite_cp_source():
    # Three-mode example, independently specified; not a discretization of 1298.
    aa=[annihilator(j,3) for j in range(3)]
    number=[a.conj().T@a for a in aa]
    qmap=lambda a:sum(a[i,j]*aa[i].conj().T@aa[j] for i in range(3) for j in range(3))
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    eye2=np.eye(2)
    time=.8
    gen=np.zeros((3,3),complex);gen[:2,:2]=sy;gen[0,2]=gen[2,0]=.4
    pref=np.diag([.2,.5,.8])
    code_basis=np.eye(4)[:,[1,2]]

    def objects(k):
        w=unitary(gen,k)
        p=w@pref@w.conj().T
        q=float(p[2,2].real)
        comp=np.diag([1-q,q])
        h=np.array([[k,.2,1],[.2,-.3*k,.1],[1,.1,.4+k]],complex)
        hf=qmap(h)
        u=unitary(hf,time)
        return p,comp,h,u

    def prepared(k,rho):
        _,comp,_,_=objects(k)
        return np.kron(comp,code_basis@rho@code_basis.conj().T)

    k=.17
    p,comp,h,u=objects(k)
    bell=np.zeros((16,16),complex)
    for i in range(2):
        for j in range(2):
            bell+=np.kron(np.kron(comp,np.outer(code_basis[:,i],code_basis[:,j])),np.outer(eye2[:,i],eye2[:,j]))/2
    initial_receiver_partner=np.trace(bell.reshape(2,4,2,2,4,2),axis1=1,axis2=4).reshape(4,4)
    residual=max(np.linalg.norm(initial_receiver_partner-np.kron(comp,eye2/2)),abs(np.trace(bell)-1))
    full_u=np.kron(u,eye2)
    evolved=full_u@bell@full_u.conj().T
    effect=np.kron(number[2],eye2)
    post=effect@evolved@effect+(np.eye(16)-effect)@evolved@(np.eye(16)-effect)
    partner=np.trace(post.reshape(8,2,8,2),axis1=0,axis2=2)
    residual=max(residual,np.linalg.norm(partner-eye2/2),abs(np.trace(post)-1))
    probs=[]
    for i in range(2):
        rho=np.outer(eye2[:,i],eye2[:,i])
        state=prepared(k,rho)
        probs.append(float(np.trace(number[2]@u@state@u.conj().T).real))
    assert probs[0]-probs[1]>.4
    # Explicit operator-source dependence, including a fixed nonzero baseline.
    rho=(eye2+.4*sx)/2
    def source(k):
        _,_,_,uk=objects(k)
        state=uk@prepared(k,rho)@uk.conj().T
        insertion=qmap(np.diag([.3+k,-.4,1.2]))
        return state,insertion,float(np.trace(insertion@state).real)+.7
    eps=1e-5
    st,ins,val=source(k)
    sm,_,vm=source(k-eps);sp,_,vp=source(k+eps)
    state_part=float(np.trace(ins@(sp-sm)/(2*eps)).real)
    contact=float(np.trace(number[0]@st).real)
    full=(vp-vm)/(2*eps)
    assert abs(full-state_part-contact)<1e-9 and contact>.1
    # A nonzero baseline .7 is not deleted by transporting the difference.
    assert abs(val-float(np.trace(ins@st).real)-.7)<1e-12
    assert residual<1e-12
    return {'scope':'three-mode CAR diagnostic, not the continuous Dirac model',
            'record_gap':round(probs[0]-probs[1],12),
            'source_full_derivative':round(full,10),
            'missed_insertion_contact':round(contact,10),
            'preserved_source_baseline':.7,
            'receiver_partner_initial_factorization':True,
            'partner_marginal_after_nonselective_record':True,
            'global_correlations_not_claimed_unchanged':True,
            'full_state_not_replaced_by_covariance':True,
            'finite_CP_residual_below_1e_minus_12':True}


def bounded_thermal():
    # An exactly specified finite thermal bath checks ordering and error logic;
    # it supplies no spectral mass data for the continuous Dirac task.
    sx=np.array([[0,1],[1,0]],complex)
    sz=np.diag([1.,-1.])
    nu=1.1;E=nu;beta=.7;T=6.;center=1.2;lam=.015
    h=.5*nu*sz
    vals=np.exp(-beta*np.diag(h));rho=np.diag(vals/sum(vals))
    def chit(omega):return T*np.sinc(omega*T/(2*np.pi))**2*np.exp(1j*omega*center)
    def amplitude(sign):
        a=np.zeros((2,2),complex)
        a[0,1]=.5*chit(sign*E+nu)
        a[1,0]=.5*chit(sign*E-nu)
        return a
    response=[];wrong=[]
    for sign in [1,-1]:
        a=amplitude(sign)
        response.append(float(np.trace(rho@a.conj().T@a).real))
        wrong.append(float(np.trace(rho@a@a.conj().T).real))
    assert min(response)>0
    assert np.allclose(response,wrong[::-1],atol=1e-13)
    logratio=np.log(response[1]/response[0])
    assert logratio>0 and abs(logratio-beta*E)<.005
    # Full finite noncommuting evolution; midpoint resolution is a diagnostic,
    # while the analytic norm inequalities in proof.md supply the guarantee.
    def evolve(k,steps):
        dt=2*T/steps;u=np.eye(4,dtype=complex)
        for t in center-T+(np.arange(steps)+.5)*dt:
            b=.5*np.array([[0,np.exp(1j*(nu+.3*k)*t)],[np.exp(-1j*(nu+.3*k)*t),0]],complex)
            flip=np.array([[0,np.exp(1j*E*t)],[np.exp(-1j*E*t),0]],complex)
            v=lam*(1-abs(t-center)/T)*np.kron(flip,b)
            u=unitary(v,dt)@u
        return u
    def probabilities(u):
        out=[]
        for i in [1,0]: # ground -> excited is F+, reverse is F-
            e=np.eye(2)[:,i];r=np.kron(np.outer(e,e),rho)
            target=np.eye(2)[:,1-i];project=np.kron(np.outer(target,target),np.eye(2))
            out.append(float(np.trace(project@u@r@u.conj().T).real))
        return np.array(out)
    u=evolve(0,2048);u2=evolve(0,4096);p=probabilities(u2)
    p_error=float(np.max(np.abs(p-probabilities(u))))
    x=lam*T;dyson_bound=x**4/2
    leading=lam**2*np.array(response)
    observed=float(np.max(np.abs(p-leading)))
    assert observed<dyson_bound and p_error<2e-9
    delta=1e-3;uk=evolve(delta,4096)
    LB=.15*(T+abs(center))
    du=float(np.linalg.norm(uk-u2,2))
    assert du<delta*abs(lam)*T*LB
    return {'scope':'finite two-level bath benchmark; no Dirac spectral claim',
            'response_log_ratio':round(float(logratio),12),
            'reversed_operator_order_log_ratio':round(float(-logratio),12),
            'full_flip_probabilities':[round(float(v),12) for v in p],
            'observed_leading_error':round(observed,12),
            'analytic_Dyson_bound':dyson_bound,
            'midpoint_resolution_difference':round(p_error,16),
            'finite_Duhamel_difference':round(du,12),
            'finite_Duhamel_bound':delta*abs(lam)*T*LB,
            'positive_responses_and_full_probabilities':True}


def rational_budget():
    # Normalize m=1 only for this algebraic budget certificate.
    V=[[1,2,-1,0,1],[2,-1,1,1,0],[0,1,2,-1,1]]
    a=[F(n,8192) for n in [1,-1,2,-1,1]]
    c=[F(n,4096) for n in [1,-1,2]]
    va=[sum(F(v)*z for v,z in zip(row,a)) for row in V]
    assert va==[F(-2,8192),F(4,8192),F(5,8192)]
    assert len(set(va))==3 and all(z!=0 for z in va)
    assert max(abs(z) for z in a)<=F(1,2048)
    assert max(abs(z) for z in c)<=F(1,1024)
    bound=F(1,64)+F(5,2048)+F(3,1024)
    assert bound==F(43,2048)<F(1,32)
    assert 2*bound==F(86,2048)<F(1,8)
    assert 2*bound+F(1,32)==F(150,2048)<F(1,8)
    # Exact probabilities of rounded uniform-noise records, independent inputs.
    means=[F(9,8)+F(1,16384),F(9,8)+F(1,256),F(9,8)-F(1,128)]
    radius=F(1,1024);step=F(1,512)
    channels=[]
    for mu,bias,cal in zip(means,va,c):
        base=mu+bias+cal;lo=base-radius;hi=base+radius
        lower=int(lo//step)-2;upper=int(hi//step)+3
        probs=[]
        for j in range(lower,upper+1):
            left=max(lo,(F(j)-F(1,2))*step)
            right=min(hi,(F(j)+F(1,2))*step)
            p=max(F(0),right-left)/(2*radius)
            if p:probs.append((j,p))
        assert sum(p for _,p in probs)==1 and all(p>0 for _,p in probs)
        channels.append([[j,str(p)] for j,p in probs])
    assert F(2)+F(2)+F(4)+F(1)+F(1)==10
    assert F(2)+F(19,10)+F(37,10)+F(11,10)+F(6,5)+F(1,10)==10
    # Fixed apparatus source is not quadratic in kappa nor removed by s~1/k.
    src=F(1,2048);k=F(1,65536);field=F(3)
    assert (field*k*k+src*k)/k==field*k+src
    return {'scope':'exact normalized joint effective input model, not hardware',
            'source_projections':[str(v) for v in va],
            'single_output_bound':str(bound),
            'geometry_source_bound':str(2*bound),
            'thermal_source_bound_including_intrinsic_bias':str(2*bound+F(1,32)),
            'independent_channel_probabilities':channels,
            'noise_variance':str(radius**2/3),
            'record_grid_spacing':str(step),
            'positive_source_margin':str(1-bound),
            'ordinary_total_energy':10,
            'source_remaining_after_dividing_by_kappa':str(src),
            'energy_table_does_not_determine_V':True}


def payload():
    paths=['proof.md','../_shared/round1300_admission.md','../1298/proof.md',
           '../../archive_1219_1233/1232/proof.md']
    return {'status':'PASS_AUTHOR_FINITE_DIAGNOSTICS',
            'finite_CP_source':finite_cp_source(),
            'finite_thermal_order_and_Duhamel':bounded_thermal(),
            'joint_effective_budget':rational_budget(),
            'continuum_PDE_proved_by_code':False,
            'Dirac_spectral_density_computed':False,
            'autonomous_instrument_factory_proved':False,
            'full_SM_or_M4_completed':False,
            'new_axioms':0,
            'source_bindings':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
            'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args();out=payload()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
        print('CREATED',out['status'])
    else:
        old=json.loads(OUT.read_text(encoding='utf-8'))
        if old!=out:raise SystemExit('MISMATCH: author result or input binding changed')
        print('READ_ONLY_MATCH',out['status'])

if __name__=='__main__':main()
