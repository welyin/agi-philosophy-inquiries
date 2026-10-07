"""Native CAR material: the same proper-time kernel and lapse insertion.
The numerical witness is a response test on a prescribed background, not GR.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'native_response_source_results.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def unitary(H,t):
    e,v=np.linalg.eigh(H)
    return (v*np.exp(-1j*t*e))@v.conj().T
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-8,abs_tol=3e-12),(a,b)
    else:assert a==b,(a,b)

def run():
    spec=importlib.util.spec_from_file_location('native956',STAGE/'956/native_material_interface.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    m=old.material(1.,.1);h=m['H'];J=m['J'];d=m['occupancy'];Delta=1+J
    cs=[old.annihilate(i) for i in range(4)];ns=[x.T@x for x in cs]
    sec=np.eye(16)[:,[i for i in range(16) if i.bit_count()==2]]
    Q=sec.T@(ns[0]+ns[1]-ns[2]-ns[3])@sec/2
    g=math.sqrt(1-d)*m['s']+math.sqrt(d)*m['d'];exc=Q@m['d']
    assert np.linalg.norm(h@g+J*g)<1e-13
    assert np.linalg.norm(Q@g-math.sqrt(d)*exc)<1e-13
    assert np.linalg.norm(h@exc-exc)<1e-13
    assert abs(np.linalg.norm(Q,2)-1)<1e-14
    parity=np.eye(6)-2*np.outer(exc,exc)
    assert np.linalg.norm(parity@h-h@parity)<1e-14
    assert np.linalg.norm(parity@Q+Q@parity)<1e-14
    evals,ev=np.linalg.eigh(h);gaps=evals+J
    weights=abs(ev.conj().T@Q@g)**2
    def coeff(N,w):return 2*d*N*Delta/((N*Delta)**2-w*w)
    def alpha(w):return coeff(1.,w)
    def beta(w):return 2*d*Delta*(Delta*Delta+w*w)/(Delta*Delta-w*w)**2
    rows=[];spectral_error=0.;derivative_error=0.
    for N in (.8,1.,1.2):
        for ratio in (0.,.1,.3,.6):
            w=ratio*N*Delta
            good=gaps>1e-9
            spectral=float(np.sum(2*N*gaps[good]*weights[good]/((N*gaps[good])**2-w*w)))
            exact=coeff(N,w);spectral_error=max(spectral_error,abs(spectral-exact))
            step=2e-4
            derivative=(coeff(N-2*step,w)-8*coeff(N-step,w)+8*coeff(N+step,w)-coeff(N+2*step,w))/(12*step)
            target=-beta(w/N)/(N*N)
            derivative_error=max(derivative_error,abs(derivative-target))
            rows.append(dict(lapse=N,proper_frequency_ratio=ratio,kernel=exact,
                             mixed_lapse_derivative=target,finite_difference=derivative))
    assert spectral_error<1e-13 and derivative_error<1e-10
    a0=2*d/Delta
    expansion=[]
    for p in (0,1,2):
        z=.01;r=.1
        ap=a0*sum(z**j for j in range(p+1))
        bp=a0*sum((2*j+1)*z**j for j in range(p+1))
        ea=z**(p+1)
        eb=z**(p+1)*((2*p+3)-(2*p+1)*z)/(1+z)
        assert abs((alpha(r*Delta)-ap)/alpha(r*Delta)-ea)<1e-14
        assert abs((beta(r*Delta)-bp)/beta(r*Delta)-eb)<1e-14
        expansion.append(dict(derivative_order=2*p,alpha_relative_remainder=ea,
                              lapse_relative_remainder=eb))
    wrong_ratio=.6
    missing=1-alpha(wrong_ratio*Delta)/beta(wrong_ratio*Delta)
    assert abs(missing-2*wrong_ratio**2/(1+wrong_ratio**2))<1e-14

    # Finite square pulses are an exact bounded Hamiltonian test, not a new device.
    area=.05;width=.1/Delta;gap=(math.pi/2-.1)/Delta
    f=area/width;U_p=unitary(h-f*Q,width)
    total_area=2*area;R=math.sinh(total_area)-total_area
    amp_bound=math.sqrt(d)*total_area
    probability_error=2*amp_bound*R+R*R
    pulse=[]
    for lam in (0.,.25):
        U=U_p@unitary((1+lam)*h,gap)@U_p
        state=U@g;prob=float(abs(exc@state)**2)
        phase=Delta*(width+(1+lam)*gap)
        born=2*d*area*area*np.sinc(Delta*width/(2*math.pi))**2*(1+math.cos(phase))
        assert abs(prob-born)<probability_error
        assert np.linalg.norm(U.conj().T@U-np.eye(6))<1e-13
        pulse.append(dict(lapse_perturbation=lam,exact_probability=prob,
                          second_order_probability=born,absolute_error=abs(prob-born),
                          energy_after=float(np.vdot(state,h@state).real)))
    contrast=abs(pulse[0]['exact_probability']-pulse[1]['exact_probability'])
    born_contrast=abs(pulse[0]['second_order_probability']-pulse[1]['second_order_probability'])
    lower=born_contrast-2*probability_error
    assert lower>5e-5 and contrast>=lower
    # Tangent of one lapse interval, from the same actual U, compared with recomputation.
    U_gap=unitary(h,gap);Ug_prime=-1j*gap*h@U_gap
    v=U_p@U_gap@U_p@g;dv=U_p@Ug_prime@U_p@g
    exact_slope=float(2*np.real(np.vdot(exc,v).conjugate()*np.vdot(exc,dv)))
    def p(lam):return float(abs(exc@U_p@unitary((1+lam)*h,gap)@U_p@g)**2)
    step=1e-4
    fd=(p(-2*step)-8*p(-step)+8*p(step)-p(2*step))/(12*step)
    assert abs(exact_slope-fd)<1e-12
    # Static Hellmann-Feynman check at fixed coordinate field f.
    N=1.1;fs=.002
    def ground(Nx):
        e,V=np.linalg.eigh(Nx*h-fs*Q)
        return float(e[0]),V[:,0]
    energy,gs=ground(N);actual_source=float(gs@h@gs)
    step=2e-4
    derivative=(ground(N-2*step)[0]-8*ground(N-step)[0]+8*ground(N+step)[0]-ground(N+2*step)[0])/(12*step)
    polarization_source=actual_source+J
    source_leading=d*fs*fs/(N*N*Delta)
    assert abs(actual_source-derivative)<1e-11
    # No inference of the full Maxwell or Einstein source from this partial derivative.
    hashes={str((STAGE/r).relative_to(ROOT)):sha(STAGE/r) for r in (
        '956/native_material_interface.py','research_note_952.md','research_note_960.md',
        'research_note_965.md','research_note_971.md','981/drafts/parent_matching_decision.md',
        '981/drafts/common_parent_contract_v1.md')}
    return dict(round=981,all_scientific_checks_passed=True,
        parameters=dict(U=1.,v=.1,eta=0.,native_dimension=6,J=J,d=d,gap=Delta,alpha_zero=a0),
        source_hashes=hashes,response=dict(spectral_error=spectral_error,
            lapse_five_point_error=derivative_error,calibrations=rows,
            low_frequency_ratio_bound=.1,derivative_expansion=expansion,
            frozen_proper_frequency_missing_fraction_at_point6=missing),
        finite_pulse=dict(area_each=area,width=width,gap_duration=gap,
            total_time=2*width+gap,total_absolute_field_area=total_area,
            odd_dyson_amplitude_remainder=R,probability_remainder_bound=probability_error,
            rows=pulse,exact_contrast=contrast,certified_contrast_lower=lower,
            lapse_slope=exact_slope,lapse_slope_difference=abs(exact_slope-fd)),
        static_source=dict(lapse=N,coordinate_field=fs,total_material_source=actual_source,
            hellmann_feynman_error=abs(actual_source-derivative),
            polarization_piece=polarization_source,leading_piece=source_leading),
        scope=dict(common_parent_fields_specified=True,one_material_mixed_response_verified=True,
            full_SM_matching_completed=False,Einstein_background_solved=False,
            dynamic_quantum_geometry_constructed=False,all_unknown_input_pulse_contrast=False,
            finite_pulse_is_autonomous_instrument=False,full_goal_completed=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    value=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as dest:json.dump(value,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(value,json.loads(OUT.read_text('utf-8')))
    print(json.dumps({k:v for k,v in value.items() if k not in ('source_hashes','response')},ensure_ascii=False,indent=2))
