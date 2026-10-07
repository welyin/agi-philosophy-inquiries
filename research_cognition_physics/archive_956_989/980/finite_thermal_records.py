"""980: finite native thermal copies, exact collision channel, bounded averaging.
No large reservoir state or infinite bath is numerically constructed.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
from decimal import Decimal, localcontext
from fractions import Fraction as Fr
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'finite_thermal_records_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def op(a):return float(np.linalg.norm(a,2))
def entropy(a):
    ev=np.linalg.eigvalsh((a+a.conj().T)/2);ev=ev[ev>1e-15]
    return float(-np.sum(ev*np.log(ev)))
def ptrace(a,part):
    t=a.reshape(4,4,4,4)
    return np.einsum('abcb->ac',t) if part==0 else np.einsum('abad->bd',t)
def unitary(h,t):
    e,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*e))@v.conj().T
def coarse_certificate():
    # Three rational sign brackets exhaust the roots of the exact cubic.
    poly=lambda x:x**3-Fr('0.975')*x*x-Fr('0.06546875')*x-Fr('0.00053125')
    brackets=[(Fr('-.054072'),Fr('-.054071')),(Fr('-.009461'),Fr('-.009460')),
              (Fr('1.038532'),Fr('1.038533'))]
    assert all(poly(a)*poly(b)<0 for a,b in brackets)
    def add(a,b):return a[0]+b[0],a[1]+b[1]
    def mul(a,b):
        v=[x*y for x in a for y in b];return min(v),max(v)
    def const(x):v=Fr(x);return v,v
    def div(a,b):
        assert b[0]*b[1]>0
        return mul(a,(1/b[1],1/b[0]))
    weights=[]
    for x in brackets:
        num=add(mul(x,add(x,const('.025'))),const('-.00046875'))
        den=add(add(mul(const(3),mul(x,x)),mul(const('-1.95'),x)),const('-.06546875'))
        weights.append(div(num,den))
    assert weights[0][0]>Fr('.022') and weights[0][1]<Fr('.023')
    assert weights[1][0]>Fr('.013') and weights[1][1]<Fr('.014')
    assert weights[2][0]>Fr('.96') and weights[2][1]<Fr('.97')
    # Positive Gibbs weights are bounded at interval endpoints. 65-digit
    # transcendental arithmetic is inflated by 1e-40 before each use.
    with localcontext() as ctx:
        ctx.prec=65;D=Decimal;pad=D('1e-40');b=2*D(2).ln()
        def dec(x):return D(x.numerator)/D(x.denominator)
        boxes=[brackets[0],brackets[1],(Fr(1),Fr(1)),brackets[2]]
        raw=[]
        for x in boxes:
            vals=[-(b+sign*pad)*dec(y) for sign in (-1,1) for y in x]
            raw.append((min(vals).exp()-pad,max(vals).exp()+pad))
        denlo=sum(x[0] for x in raw);denhi=sum(x[1] for x in raw)
        probs=[(a/denhi-pad,z/denlo+pad) for a,z in raw]
        assert all(a>D('.09') and z<D('.42') for a,z in probs)
        assert probs[2][0]>D('.096') and probs[2][1]<D('.1')
        assert probs[3][0]>D('.09') and probs[3][1]<D('.1')
        assert probs[0][1]+probs[1][1]<D('.83')
    energy_boxes=[brackets[0],brackets[1],(Fr(1),Fr(1)),brackets[2]]
    q_edges=[(i,2) for i in (0,1,3)]+[(2,i) for i in (0,1,3)]
    off_gap_bounds=[]
    for i,j in q_edges:
        for k,l in q_edges:
            if sorted((i,k))==sorted((j,l)):continue
            gaplo=energy_boxes[i][0]+energy_boxes[k][0]-energy_boxes[j][1]-energy_boxes[l][1]
            gaphi=energy_boxes[i][1]+energy_boxes[k][1]-energy_boxes[j][0]-energy_boxes[l][0]
            assert gaplo*gaphi>0
            off_gap_bounds.append(min(abs(gaplo),abs(gaphi)))
    assert min(off_gap_bounds)>Fr('.03')
    x=Fr('.013');gap=Fr('.096')*(x-x**3/6)**2
    y=Fr('.96');central_gap=Fr('.09')*(y*y/2-y**4/24)
    assert gap>Fr(1,62500) and central_gap>Fr(1,62500)
    assert Fr('.1')+Fr('.83')*Fr('.023')**2<Fr('.101')
    # r^600000 <= exp(-9.6); certify exp(9.6)>14000 by a positive Taylor sum.
    expsum=sum((Fr('9.6')**j)/math.factorial(j) for j in range(33))
    assert expsum>14000
    return dict(root_brackets=[[str(a),str(b)] for a,b in brackets],
        weight_intervals=[[float(a),float(b)] for a,b in weights],
        thermal_intervals=[[float(a),float(b)] for a,b in probs],
        nonresonant_gap_rational_lower=float(min(off_gap_bounds)),
        analytic_spectral_gap_lower=float(gap),uniform_rate='1 - 1/62500',
        minimum_population_diagonal_lower=.899,exp_9_6_rational_lower=float(expsum))
def run():
    _,h,q,F,_,_=load('native968',STAGE/'968/internal_relay.py').material()
    eta=.05;h=h+eta*F;e,R=np.linalg.eigh(h);q=R.T@q@R
    center=int(np.argmin(abs(e-1.)));outer=[i for i in range(4) if i!=center]
    beta=2*math.log(2);tau=np.exp(-beta*(e-e[0]));tau/=sum(tau)
    theta=1.;weights=np.array([q[i,center]**2 for i in range(4)])
    assert center==2 and min(weights[outer])>.01 and abs(sum(weights)-1)<1e-13
    # Independently evaluate the small eigensystem residual using 65-digit input.
    with localcontext() as ctx:
        ctx.prec=65;D=Decimal
        hh=[[D(0),D('-.2'),D(0),-D('.05')*D(3).sqrt()/4],
            [D('-.2'),D(1),D(0),D(0)],[D(0),D(0),D(1),D(0)],
            [-D('.05')*D(3).sqrt()/4,D(0),D(0),D('-.025')]]
        rr=[[D.from_float(float(R[i,j])) for j in range(4)] for i in range(4)]
        ee=[D.from_float(float(x)) for x in e]
        residual=sum((sum(hh[i][k]*rr[k][j] for k in range(4))-rr[i][j]*ee[j])**2
                     for i in range(4) for j in range(4)).sqrt()
    gram=op(R.T@R-np.eye(4));res=float(residual)
    assert res<1e-13 and gram<1e-13 and min(np.diff(e))>.038
    H0=np.diag((e[:,None]+e[None,:]).ravel());V=np.kron(q,q);Vbar=np.zeros((16,16))
    for i in outer:
        a=4*i+center;b=4*center+i;Vbar[a,b]=Vbar[b,a]=weights[i]
    diff=np.diag(H0)[:,None]-np.diag(H0)[None,:]
    off=V-Vbar;mask=abs(off)>1e-12
    min_off=float(np.min(abs(diff[mask])))
    assert min_off>.03
    Cactual=float(2*np.sum(abs(off[mask])/abs(diff[mask])))
    # ||V||_F=2; sum of absolute entries <= 16*2. Thus C<=64/.03<4000.
    C=4000.;assert Cactual<C and op(V)<1+1e-13 and op(Vbar)<1+1e-13
    cosines=np.ones((4,4));P=np.eye(4)
    for i in outer:
        s=math.sin(theta*weights[i])**2;c=math.cos(theta*weights[i])
        P[center,i]=tau[center]*s;P[i,center]=tau[i]*s
        P[i,i]-=tau[center]*s;P[center,center]-=tau[i]*s
        cosines[i,center]=cosines[center,i]=c
    alpha=(cosines*tau)@cosines.T
    def channel(rho):
        out=alpha*rho;np.fill_diagonal(out,P@np.diag(rho));return out
    sqrt=np.sqrt(tau);B=P*sqrt[None,:]/sqrt[:,None]
    assert op(B-B.T)<1e-14 and op(P@tau-tau)<1e-14
    ev=np.linalg.eigvalsh(B);raw_pop=max(abs(ev[:-1]));raw_coh=max(abs(alpha[i,j]) for i in range(4) for j in range(4) if i!=j)
    raw_rate=max(raw_pop,raw_coh)
    certificate=coarse_certificate();rate=1-1/62500
    assert raw_rate<rate<1
    # Norm conversion: half diamond <= d/2 * Frobenius(Choi) for d=4.
    # Pop part <= kappa*sqrt(3)*r^N; offdiagonals <=sqrt(12)*r^N.
    kappa=math.sqrt(max(tau)/min(tau));kappa_upper=3.
    assert kappa<kappa_upper
    prefactor=13.;assert 2*math.sqrt(3*kappa_upper**2+12)<prefactor
    effective_budget=.001;parent_budget=.0001
    N=100000*math.ceil(62500*math.log(prefactor/effective_budget)/100000)
    assert N==600000
    effective_bound=13/14000
    g=parent_budget/(N*C*(2+2*theta))
    total_bound=effective_bound+parent_budget
    assert total_bound<=.0011 and g>0
    # Independent exact 16x16 dilation of ONE comparison contact.
    U=unitary(Vbar,theta);dilation_error=0.
    for i in range(4):
        for j in range(4):
            rho=np.zeros((4,4),complex);rho[i,j]=1
            actual=ptrace(U@np.kron(rho,np.diag(tau))@U.conj().T,0)
            dilation_error=max(dilation_error,op(actual-channel(rho)))
    assert dilation_error<1e-12
    # Channel power and Choi check at the final N, using the reversible eigensystem.
    vals,z=np.linalg.eigh(B)
    correction=(z[:,:-1]*(vals[:-1]**N))@z[:,:-1].T
    PN=np.outer(tau,np.ones(4))+sqrt[:,None]*correction/sqrt[None,:]
    choi=np.zeros((16,16),complex)
    for i in range(4):
        for j in range(4):
            block=np.zeros((4,4),complex)
            if i==j:np.fill_diagonal(block,PN[:,i]-tau)
            else:block[i,j]=alpha[i,j]**N
            choi[4*i:4*i+4,4*j:4*j+4]=block
    choi_bound=.5*float(np.sum(abs(np.linalg.eigvalsh(choi))))
    assert choi_bound<=effective_bound
    # Thermodynamic identity on a genuinely coherent initial material state.
    psi=np.ones(4,dtype=complex)/2;rho=np.outer(psi,psi.conj())
    joint=U@np.kron(rho,np.diag(tau))@U.conj().T
    final=ptrace(joint,0);bath=ptrace(joint,1)
    heat=float(e@np.real(np.diag(bath)-tau))
    material_change=float(e@np.real(np.diag(final-rho)))
    mutual=entropy(final)+entropy(bath)-entropy(joint)
    rel_bath=-entropy(bath)-float(np.real(np.diag(bath))@np.log(tau))
    sigma=entropy(final)-entropy(rho)+beta*heat
    landauer_error=abs(sigma-mutual-rel_bath)
    assert abs(heat+material_change)<1e-12 and landauer_error<1e-12 and sigma>=0
    # Actual Q*Q contact for independent moderate-g diagnostic (not large-time simulation).
    diagnostics=[]
    for test_g in (1e-3,5e-4):
        t=theta/test_g
        actual=unitary(H0+test_g*V,t);reference=unitary(H0+test_g*Vbar,t)
        err=op(actual-reference);bound=test_g*Cactual*(2+theta*(op(V)+op(Vbar)))
        assert err<bound
        diagnostics.append(dict(g=test_g,time=t,full_unitary_error=err,averaging_bound=bound))
    S_tau=-float(tau@np.log(tau));mean_above_ground=float((e-e[0])@tau)
    info_upper=2*(-total_bound*math.log(total_bound)-(1-total_bound)*math.log(1-total_bound)+total_bound*math.log(3))
    # This thermal reset state is not a pure blank.
    pure_blank_distance=1-max(tau)
    assert pure_blank_distance>.58
    files=[Path(__file__),HERE/'drafts/finite_thermal_decision.md',STAGE/'968/internal_relay.py',
        STAGE/'research_note_972.md',STAGE/'975/drafts/record_arrow_adoption.md',
        STAGE/'979/spectral_mass_bridge_results.json']
    return dict(round=980,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(eta=eta,beta=beta,theta=theta,energies=e.tolist(),thermal_populations=tau.tolist(),
            Q_star_center=center,resonant_weights=weights.tolist()),
        eigensystem_certificate=dict(decimal_residual=res,orthogonality_error=gram,
            minimum_single_gap=float(min(np.diff(e))),coarse_analytic_certificate=certificate,
            minimum_nonresonant_V_gap=min_off,primitive_integral_bound=C,
            diagnostic_integral_bound=Cactual),
        channel=dict(population_matrix=P.tolist(),coherence_multipliers=alpha.tolist(),
            reversible_population_spectrum=ev.tolist(),population_rate=float(raw_pop),
            coherence_rate=float(raw_coh),certified_rate=float(rate),norm_prefactor=prefactor,
            stationary_error=op(P@tau-tau),one_contact_dilation_error=dilation_error,
            contacts=N,comparison_thermal_reset_bound=effective_bound,final_choi_reset_bound=choi_bound,
            full_charge_reset_bound=total_bound,arbitrary_unknown_input_and_passive_reference=True,
            final_label_information_upper_nats=info_upper),
        finite_resources=dict(coupling_g=g,each_contact_time=theta/g,total_contact_time=N*theta/g,
            prepared_thermal_copies=N,initial_bath_energy_above_ground=N*mean_above_ground,
            initial_bath_entropy=N*S_tau,thermal_state_to_best_pure_blank_distance=pure_blank_distance,
            switching_work_absolute_bound=N*g,complete_bath_bare_energy_change_bound=float(e[-1]-e[0])+N*g,
            no_bath_reused_or_reset=True,controller_preparation_not_constructed=True),
        thermodynamic_check=dict(material_energy_change=material_change,bath_heat=heat,
            energy_balance_error=abs(heat+material_change),endpoint_entropy_production=sigma,
            mutual_information=mutual,bath_relative_entropy=rel_bath,landauer_identity_error=landauer_error),
        averaging_diagnostics=diagnostics,
        scope=dict(finite_positive_full_charge_collision_sequence=True,infinite_bath_required=False,
            new_primitive_reset_operator_inserted=False,constant_autonomous_controller_constructed=False,
            pure_blank_reset_proved=False,universe_permanent_arrow_proved=False,
            actual_geometric_feedback_computed=False,whole_lifecycle_with_976_verified=False,
            full_SM_matching_completed=False,full_goal_completed=False),
        references=['https://arxiv.org/abs/2106.11974','https://arxiv.org/abs/1306.4352'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-6,abs_tol=1e-11),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not OUT.exists()
    out=run()
    if a.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(OUT.read_text('utf-8-sig')))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
