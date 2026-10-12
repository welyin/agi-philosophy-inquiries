"""1155: normal-source analytic-tool calibrations; no full field or large-T simulation."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math
import numpy as np
HERE=Path(__file__).resolve().parent

def norm(x):return float(np.linalg.norm(x,2))
def exp_h(h,t=1):
    e,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*e))@v.conj().T
def outer(v):return np.outer(v,v.conj())
def rat(x):return dict(numerator=x.numerator,denominator=x.denominator,float=float(x))
def car(j):
    a=np.zeros((16,16),complex)
    for m in range(16):
        if m&(1<<j):
            a[m^(1<<j),m]=(-1)**((m&((1<<j)-1)).bit_count())
    return a
def car_data():
    aa=[car(j) for j in range(4)]
    nn=[a.conj().T@a for a in aa]
    na,nb=nn[1],nn[3]
    xa=aa[0].conj().T@aa[1]+aa[1].conj().T@aa[0]
    return aa,nn,na,nb,xa

def budgets_and_CAR():
    T=F(10**10);alpha2=F(24000);q2=F(3);fpp=F(24);F2=F(10,7)
    delta=alpha2*q2*fpp/T
    eself=alpha2*q2*fpp*fpp/(2*T**3)
    h=math.sqrt(float(eself))
    record=float(delta)+2*h
    gamma=F(5,4);kappa=F(3,4)
    energy=gamma*2*eself+kappa*delta
    momentum=kappa*2*eself+gamma*delta
    assert delta==F(1728,10**7)
    assert record<.000172801
    assert gamma*gamma-kappa*kappa==1
    assert float(energy)<.000129601 and float(momentum)<.000216001
    aa,nn,na,nb,xa=car_data()
    I=np.eye(16)
    anti=max(norm(aa[i]@aa[j].conj().T+aa[j].conj().T@aa[i]-(I if i==j else 0))
             for i in range(4) for j in range(4))
    assert anti==0 and norm(na@nb-nb@na)==0
    thetaA=.271;thetaB=-.419
    # Exact self phases and finite material clock phases are retained, not deleted.
    unitary=exp_h(-(thetaA*na@na+thetaB*nb@nb+math.pi*na@nb))
    W=np.zeros((16,2),complex)
    for a in range(2):
        for b in range(2):W[(1<<a)+(1<<(2+b)),a]=1/math.sqrt(2)
    EB=np.zeros((16,16),complex)
    for m in range(4):
        beta=np.zeros(16,complex)
        beta[m+4]=1/math.sqrt(2);beta[m+8]=np.exp(1j*thetaB)/math.sqrt(2)
        EB+=outer(beta)
    pulled=W.conj().T@unitary.conj().T@EB@unitary@W
    error=norm(pulled-np.diag([1.,0.]))
    assert error<1e-14 and norm(unitary.conj().T@unitary-I)<1e-14
    assert norm((na@na)-na)==0 and norm((nb@nb)-nb)==0
    insertion=exp_h(xa,math.pi/4)
    assert norm(insertion@na-na@insertion)>.7
    # The column-operator identity, not sampled unknown states, proves the reference extension.
    referror=norm(np.kron(pulled-np.diag([1.,0.]),np.eye(3)))
    assert abs(referror-error)<1e-14
    return dict(name="prior_budgets_and_complete_CAR_record",passed=True,
      duration=int(T),alpha_squared_upper=int(alpha2),q_norm_squared_upper=int(q2),
      f_second_derivative_L1_upper=int(fpp),F2=rat(F2),
      phase_and_spatial_impulse_error=rat(delta),single_source_radiation_energy_bound=rat(eself),
      single_source_displacement_bound=h,complete_record_error_bound=record,
      boosted_station_energy_error=rat(energy),boosted_station_momentum_error=rat(momentum),
      complete_CAR_dimension=16,CAR_residual=anti,bare_pointer_effect_error=error,
      self_phases_retained=True,arbitrary_reference_supported_by_operator_identity=True,
      large_T_field_numerically_simulated=False)

def frequency_diagnostic():
    # Gaussian spatial profile is only an independently integrable diagnostic.
    # It is not the compact support kernel and cannot certify the physical causality menu.
    T=8.;alpha=.8;sigma=.4;d=3.
    def f(u):return 30*u*u*(1-u)*(1-u)
    kn,kw=np.polynomial.legendre.leggauss(120)
    k=(kn+1)*8.;wk=kw*8.;omega=np.sqrt(1+k*k)
    measure=wk*k*k*np.exp(-sigma*sigma*k*k)/(2*math.pi**2)
    nt=96
    xn,xw=np.polynomial.legendre.leggauss(nt)
    t=(xn+1)*T/2;tw=xw*T/2;s=alpha/math.sqrt(T)*f(t/T)
    # Independent triangular time quadrature of the actual retarded Green function.
    u=t[:,None]*(xn[None,:]+1)/2
    uw=t[:,None]*xw[None,:]/2
    su=alpha/math.sqrt(T)*f(u/T)
    delay=t[:,None]-u
    argument=omega[:,None,None]*delay[None,:,:]
    amplitude=tw[:,None]*s[:,None]*uw*su
    R=np.einsum("ktu,tu->k",np.sin(argument)/omega[:,None,None],amplitude)
    time_work=np.einsum("ktu,tu->k",np.cos(argument),amplitude)
    temporal=np.exp(1j*omega[:,None]*t[None,:])@(tw*s)
    energy_frequency=.5*np.abs(temporal)**2
    assert np.max(np.abs(time_work-energy_frequency))<1e-12
    # Independent exact forced-oscillator update for midpoint piecewise-constant sources.
    steps=4096;dt=T/steps
    y=np.zeros_like(k);v=np.zeros_like(k);Rstep=np.zeros_like(k);Qstep=np.zeros_like(k)
    c=np.cos(omega*dt);si=np.sin(omega*dt)
    for n in range(steps):
        sn=alpha/math.sqrt(T)*f((n+.5)/steps)
        next_y=c*y+si*v/omega+sn*(1-c)/omega**2
        next_v=-omega*si*y+c*v+sn*si/omega
        Rstep+=sn*(y*si/omega+v*(1-c)/omega**2+
                    sn*(dt/omega**2-si/omega**3))
        Qstep+=sn*(next_y-y)
        y,v=next_y,next_v
    energy_step=.5*(v*v+omega*omega*y*y)
    assert np.max(np.abs(Qstep-energy_step))<1e-10
    assert np.max(np.abs(Rstep-R))<2e-6
    assert np.max(np.abs(energy_step-energy_frequency))<2e-6
    sinc=np.sinc(k*d/math.pi)
    phase=float(np.dot(measure,R*sinc))
    hstep=1e-4
    def phase_at(distance):
        return float(np.dot(measure,R*np.sinc(k*distance/math.pi)))
    derivative=(phase_at(d+hstep)-phase_at(d-hstep))/(2*hstep)
    # Actual spatial derivative under the original angular source integral.
    mu,mw=np.polynomial.legendre.leggauss(80)
    angular=.5*np.sum(mw[None,:]*(-k[:,None]*mu[None,:])*
                       np.sin(k[:,None]*d*mu[None,:]),axis=1)
    impulse_B=float(np.dot(measure,R*angular))
    impulse_A=-impulse_B
    Es=float(np.dot(measure,time_work))
    Ec=float(np.dot(measure,time_work*sinc))
    source_B_energy=-(Es+Ec);source_A_energy=source_B_energy
    # Independent total final classical/normal coherent field energy, with both source profiles.
    interference=.5*np.sum(mw[None,:]*
                  np.abs(1+np.exp(1j*k[:,None]*d*mu[None,:]))**2,axis=1)
    final_energy=float(np.dot(measure,energy_frequency*interference))
    assert Es>1e-8 and abs(source_B_energy)>1e-8
    assert abs(Ec)<=Es+1e-14
    assert abs(impulse_B-derivative)<1e-9
    assert impulse_B<0 and abs(impulse_A+impulse_B)<1e-15
    assert abs(source_A_energy+source_B_energy+final_energy)<1e-12
    return dict(name="finite_pulse_retarded_phase_and_station_exchange",passed=True,
      diagnostic_duration=T,diagnostic_alpha=alpha,Gaussian_sigma=sigma,separation=d,
      spatial_k_cutoff=16.,spatial_nodes=120,triangular_time_nodes=nt,
      forced_oscillator_steps=steps,phase=phase,phase_derivative=derivative,
      station_B_spatial_impulse=impulse_B,station_A_spatial_impulse=impulse_A,
      self_energy_from_retarded_work=Es,cross_energy_from_retarded_work=Ec,
      station_B_time_exchange=source_B_energy,station_A_time_exchange=source_A_energy,
      final_field_energy_independent=final_energy,
      energy_balance_residual=abs(source_A_energy+source_B_energy+final_energy),
      retarded_R_crosscheck=float(np.max(np.abs(Rstep-R))),
      radiation_energy_crosscheck=float(np.max(np.abs(energy_step-energy_frequency))),
      per_frequency_work_identity=float(np.max(np.abs(time_work-energy_frequency))),
      Gaussian_is_not_actual_compact_profile=True,
      large_T_not_numerically_simulated=True)

def coherent_overlap(a,b):
    return np.exp(-.5*np.vdot(a,a).real-.5*np.vdot(b,b).real+np.vdot(a,b))

def coherent_history():
    aa,nn,na,nb,xa=car_data()
    av=np.diag(na).real;bv=np.diag(nb).real
    insertion=exp_h(xa,math.pi/4)
    # Two degenerate unit-frequency normal oscillator modes calibrate the exact
    # Weyl grammar; no Fock truncation. Field number equals energy in these units.
    a1=np.array([.07+.02j,.04j]);b1=np.array([-.01+.03j,.05])
    a2=np.array([-.02+.01j,.06-.02j]);b2=np.array([.03+.02j,-.03j])
    h1=np.array([av[m]*a1+bv[m]*b1 for m in range(16)])
    h2=np.array([av[m]*a2+bv[m]*b2 for m in range(16)])
    ph1=.23*av**2-.17*bv**2+math.pi*av*bv
    ph2=-.31*av**2+.29*bv**2+.71*av*bv
    coefficient=np.zeros((16,16),complex)
    displacement=np.zeros((16,16,2),complex)
    target=np.diag(np.exp(1j*ph2))@insertion@np.diag(np.exp(1j*ph1))
    for n in range(16):
        for m in range(16):
            coefficient[n,m]=target[n,m]*np.exp(1j*np.imag(np.vdot(h1[m],h2[n])))
            displacement[n,m]=h1[m]+h2[n]
    actual_gram=np.zeros((16,16),complex)
    overlap_target=np.zeros((16,16),complex)
    number=np.zeros((16,16),complex)
    for m in range(16):
        for l in range(16):
            for n in range(16):
                g=coherent_overlap(displacement[n,m],displacement[n,l])
                product=coefficient[n,m].conj()*coefficient[n,l]
                actual_gram[m,l]+=product*g
                overlap_target[m,l]+=coefficient[n,m].conj()*target[n,l]*coherent_overlap(displacement[n,m],np.zeros(2))
                number[m,l]+=product*np.vdot(displacement[n,m],displacement[n,l])*g
    isometry_error=norm(actual_gram-np.eye(16))
    assert isometry_error<1e-13
    diff=actual_gram+target.conj().T@target-overlap_target-overlap_target.conj().T
    diff=(diff+diff.conj().T)/2
    error=math.sqrt(max(0.,float(np.linalg.eigvalsh(diff)[-1])))
    first=max(float(np.linalg.norm(x)) for x in h1)
    second=max(float(np.linalg.norm(x)) for x in h2)
    assert error<=first+second+1e-13
    assert error>.01
    # Comparison only: a genuine CP trace-and-reprepare field reset retains
    # first-window material dephasing. The actual history above never resets.
    reset_number=target.conj().T@np.diag([np.vdot(x,x).real for x in h2])@target
    first_gram=np.array([[coherent_overlap(h1[m],h1[l]) for l in range(16)] for m in range(16)])
    reset_number=reset_number*first_gram
    assert np.linalg.eigvalsh((reset_number+reset_number.conj().T)/2)[0]>-1e-13
    psi=np.zeros(16,complex)
    for a in range(2):
        for b in range(2):psi[(1<<a)+(1<<(2+b))]=.5
    actual_number=float(np.vdot(psi,number@psi).real)
    erased_number=float(np.vdot(psi,reset_number@psi).real)
    assert abs(actual_number-erased_number)>1e-4
    # Retain every matrix entry rather than testing a few unknown states.
    assert np.linalg.eigvalsh((number+number.conj().T)/2)[0]>-1e-13
    referror=norm(np.kron(diff,np.eye(3)))
    assert abs(referror-norm(diff))<1e-12
    return dict(name="two_exact_Weyl_windows_with_noncommuting_CAR_and_no_reset",passed=True,
      complete_material_dimension=16,normal_oscillator_modes=2,
      full_history_isometry_error=isometry_error,
      complete_output_error_to_retained_phase_vacuum_target=error,
      first_window_displacement_bound=first,second_window_displacement_bound=second,
      sequential_full_output_bound=first+second,
      actual_final_field_number=actual_number,
      CP_reset_control_final_field_number=erased_number,
      CP_reset_control_keeps_first_window_material_dephasing=True,
      field_memory_changes_later_energy=abs(actual_number-erased_number),
      passive_reference_squared_error_norm_difference=abs(referror-norm(diff)),
      Fock_cutoff_used=False,first_window_field_reset=False,
      second_window_vacuum_exchange_not_assumed=True,
      source_causality_and_boost_not_proved_by_this_two_mode_diagnostic=True)

def compare(a,b,where="root"):
    if isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=2e-8,abs_tol=1e-11),(where,a,b)
    elif isinstance(a,dict):
        assert a.keys()==b.keys(),where
        for k in a:compare(a[k],b[k],where+"."+k)
    elif isinstance(a,list):
        assert len(a)==len(b),where
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,where+"."+str(i))
    else:assert a==b,(where,a,b)
def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    src=HERE/"source_hashes.json"
    if src.exists():
        for name,digest in json.loads(src.read_text(encoding="utf-8")).items():
            assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    result=dict(round=1155,all_passed=True,numerical_check_groups=3,
      groups=[budgets_and_CAR(),frequency_diagnostic(),coherent_history()],
      interpretation="Analytic-tool calibrations; the actual compact 3D parent and Lorentz transport require the proof.",
      physical_new_experiments=0,new_adopted_axioms=0,full_goal_complete=False)
    dest=HERE/"results.json"
    if args.write_results and not dest.exists():
        dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        assert dest.exists(),"Run --write-results once."
        compare(result,json.loads(dest.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="PASS",round=1155,groups=3,
        source_hashes_checked=src.exists(),
        record_bound=result["groups"][0]["complete_record_error_bound"],
        diagnostic_energy_balance=result["groups"][1]["energy_balance_residual"],
        Weyl_history_error=result["groups"][2]["full_history_isometry_error"])))
if __name__=="__main__":main()
