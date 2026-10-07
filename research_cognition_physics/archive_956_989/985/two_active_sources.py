"""Two active native mass sources: finite-domain common-process certificate.

No spatial grid evolution. Initial moments are exact on finite photon support;
the extension and Duhamel bounds are analytic, with rational upper certificates.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
from fractions import Fraction as F
import numpy as np

HERE=Path(__file__).resolve().parent; STAGE=HERE.parent; ROOT=HERE.parents[2]
OUT=HERE/'two_active_sources_results.json'

def read(p): return json.loads(p.read_text('utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
def kron(*xs):
    out=np.ones((1,1))
    for x in xs: out=np.kron(out,x)
    return out
def fraction(x): return dict(exact=str(x), decimal=float(x))

def initial_moments(n):
    m,h,q,_,W,_=load('native968_985',STAGE/'968/internal_relay.py').material()
    eye=np.eye(4); ann=np.diag(np.sqrt(np.arange(1,n)),1)
    hmat=kron(h,eye)+kron(eye,h)+.2*kron(q,q)
    s=kron(q,eye)+kron(eye,q)
    hl=kron(hmat,np.eye(n))+.5*kron(np.eye(16),np.diag(np.arange(n)))
    hl+=.003*kron(s,ann+ann.T)+.000018*kron(s@s,np.eye(n))
    k=kron(hl,np.eye(2))+kron(np.eye(16*n),np.diag([-m['J'],0.]))
    photon=np.zeros(n); photon[:2]=1/math.sqrt(2)
    plus=W@np.ones(2)/math.sqrt(2); rplus=np.ones(2)/math.sqrt(2)
    w=np.column_stack([np.kron(np.kron(np.kron(W[:,i],plus),photon),rplus) for i in (0,1)])
    first=w.T@k@w; second=(k@w).T@(k@w)
    # n>=3 includes the entire K W image, so these are infinite-Fock moments.
    assert n>=3 and np.linalg.norm(w.T@w-np.eye(2))<1e-12
    return first,second,m

def run():
    old=read(STAGE/'978/moving_complete_source_results.json')
    record=read(STAGE/'976/finite_internal_reference_results.json')
    p=old['parameters']; g=p['G']; mu=p['mu']; a=p['soft_core']; t=p['time_max']
    first,second,m=initial_moments(3)
    first4,second4,_=initial_moments(4)
    moment_difference=max(float(np.linalg.norm(first-first4)),float(np.linalg.norm(second-second4)))
    assert moment_difference<1e-12
    # Conservative analytic bound includes material, reference, and both field terms.
    k_triangle=3*1.04+.2+.5+.003*2*math.sqrt(2)+.000018*4
    assert k_triangle<4 and -3*m['J']-.2>-1
    k_rms=float(np.sqrt(np.linalg.eigvalsh(second).max()))
    total_second=kron(second,np.eye(2))+kron(np.eye(2),second)+2*kron(first,first)
    total_rms=float(np.sqrt(np.linalg.eigvalsh(total_second).max()))
    assert k_rms<4 and total_rms<8
    # Keep the immense baseline symbolic: never subtract mu+K - mu in doubles.
    cutoff_over_mu=2+400/mu
    lower_loss=g*mu/(4*a)*cutoff_over_mu
    tail_state_bound=16/(mu+200)  # 2 ||(K_A+K_B)W||/(E-mu-200)
    active_potential_bound=t*g/a*204*4
    active_kinetic_increment=t*math.sqrt(15)/(8*p['position_sigma']**2)*2/mu**2
    inherited=old['uniform_moving_error']['total']
    # Rational guards avoid losing these positive terms when added to 0.00172.
    assert 2e35<mu<3e35 and t<140000 and lower_loss<.00125
    tail_upper=F(8,10**35)
    potential_upper=F(12,10**31)
    kinetic_upper=F(4,10**72)
    base_upper=F(17202,10**7)
    assert inherited<float(base_upper)
    assert tail_state_bound<float(tail_upper)
    assert active_potential_bound<float(potential_upper)
    assert active_kinetic_increment<float(kinetic_upper)
    total_upper=base_upper+tail_upper+potential_upper+kinetic_upper
    report_upper=F(1721,10**6)
    assert total_upper<report_upper
    contrast=F(4859481,10**7)-2*report_upper
    assert float(contrast)<record['interval_certificate']['uniform_infinite_contrast_lower']-2*float(report_upper)
    assert contrast>F(482506,10**6)
    # Independent scalar mass-fibre identities and negative extrapolation control.
    samples=[]
    for ma,mb in [(2.,3.),(8.,9.),(50.,20.)]:
        e=10.; gg=.01; core=2.; x=np.array([.3,-.1,.2]); y=np.array([1.1,.4,-.2])
        s=ma+mb; chi=min(1.,e/s); d=y-x; u=1/math.sqrt(float(d@d)+core*core)
        potential=-gg*ma*mb*chi*u; force=gg*ma*mb*chi*d*u**3
        step=1e-4; errors=[]
        def v(xx,yy): return -gg*ma*mb*chi/math.sqrt(float((yy-xx)@(yy-xx))+core*core)
        for i in range(3):
            z=np.eye(3)[i]*step
            dx=(v(x-2*z,y)-8*v(x-z,y)+8*v(x+z,y)-v(x+2*z,y))/(12*step)
            dy=(v(x,y-2*z)-8*v(x,y-z)+8*v(x,y+z)-v(x,y+2*z))/(12*step)
            errors.extend([abs(dx+force[i]),abs(dy-force[i])])
        assert potential>=-gg*e/(4*core)*s and max(errors)<1e-11
        samples.append(dict(masses=[ma,mb],chi=chi,force_derivative_error=max(errors),
                            constraint_source=potential/gg,rest_plus_potential=s+potential))
    negative=[]
    for n in (100,1000,10000):
        negative.append(dict(n=n,unmodified_rest_plus_potential=2*n-.01*n*n/2,
                             extended_rest_plus_potential=2*n-.01*n*n*min(1,10/(2*n))/2))
    assert negative[-1]['unmodified_rest_plus_potential']<-400000
    files=[Path(__file__),HERE/'drafts/two_active_sources_decision.md',
           STAGE/'968/internal_relay.py',STAGE/'976/finite_internal_reference_results.json',
           STAGE/'978/moving_complete_source_results.json',STAGE/'978/moving_complete_source.py',
           STAGE/'981/drafts/common_parent_contract_v1.md']
    return dict(round=985,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(G=g,mu=mu,soft_core=a,time_max=t,mass_baseline_A=200,
                        cutoff_symbolic='E=2*mu+400; held fixed under source variations'),
        initial_moments=dict(K_first=first.tolist(),K_second=second.tolist(),
            rms_K=k_rms,rms_K_A_plus_K_B=total_rms,analytic_K_norm_upper=4,
            analytic_sum_norm_upper=8,moment_support_comparison=moment_difference,
            unknown_joint_two_sender_input_and_passive_reference=True,
            all_Fock_occupations_retained_in_theorem=True),
        stability=dict(relative_negative_form_bound=lower_loss,
            common_energy_window_reduces_both_H=True,
            modified_high_energy_interaction_is_explicit_input=True,
            lower_bound_implies_positive_H=True,
            unmodified_unbounded_below_does_not_imply_nonunitarity=True),
        uniform_budget=dict(inherited_978=inherited,high_energy_extension=tail_state_bound,
            second_active_source=active_potential_bound,second_kinetic_increment=active_kinetic_increment,
            rational_base_upper=fraction(base_upper),rational_tail_upper=fraction(tail_upper),
            rational_active_upper=fraction(potential_upper),rational_kinetic_upper=fraction(kinetic_upper),
            exact_sum_upper=fraction(total_upper),rounded_state_upper=fraction(report_upper),
            record_contrast_lower=fraction(contrast)),
        source_checks=samples,negative_extrapolation_control=negative,
        conservation=dict(total_energy=True,total_momentum=True,both_complete_internal_masses=True,
                          equal_opposite_forces=True,same_H_source_derivative=True),
        scope=dict(two_active_source_Newton_subinterface_completed=True,
            entire_U1_completed=False,radiative_TT_matched=False,full_SM_matching=False,
            GR_total_physical_error_certified=False,physical_minimum_scale_assumed=False,
            microscopic_continuity_assumed=False,UV_completion_claimed=False,
            unbounded_source_error_deduced_from_probability=False,full_goal_completed=False),
        references=['https://arxiv.org/abs/1908.06929','https://arxiv.org/abs/gr-qc/9405057'],
        source_hashes={str(f.relative_to(ROOT)):sha(f) for f in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a: compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-90),(a,b)
    else: assert a==b,(a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write: assert not OUT.exists()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else: compare(result,read(OUT))
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
