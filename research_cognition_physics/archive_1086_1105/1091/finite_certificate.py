"""Round1091 finite candidate-c certificate; no empirical data or all-c exclusion."""
from fractions import Fraction as Q
from pathlib import Path
import json,math,sys
HERE=Path(__file__).resolve().parent

def run():
    umin,umax,Lmin,Lmax=Q(499,1000),Q(501,1000),Q(999,1000),Q(1001,1000)
    gamma_low,gamma_high=Q(23,20),Q(6,5)
    assert gamma_low**2*(1-umin**2)<1
    assert gamma_high**2*(1-umax**2)>1
    omega,eps=Q(1,5),Q(1,1000)
    pulse=2*omega*eps/umin
    pulse_upper=Q(802,1000000)
    assert pulse<pulse_upper
    # pi < 22/7, and cos(.5) > 1-(.5)^2/2 = 7/8.
    midpoint_upper=3*Q(22,7)/(2*1024)
    stat=Q(1,200)
    ep=pulse_upper+midpoint_upper+stat
    et=ep/Q(7,80)
    lorentz_lower=gamma_low*umin*Lmin
    all_four_time_upper=(1+gamma_high)*2*et
    margin=lorentz_lower-all_four_time_upper
    assert margin>Q(499,10000)
    samples=140000
    exponent=2*samples*stat**2
    assert exponent==7
    # e**7 > sum_{k=0}^{10} 7**k/k! > 800, hence 8*exp(-7)<.01.
    partial=sum(Q(7**k,math.factorial(k)) for k in range(11))
    assert partial>800
    return {'schema':'round1091_finite_certificate_v1','status':'PASS','candidate_c':1,
      'declared_u_interval':[float(umin),float(umax)],'declared_L_interval':[float(Lmin),float(Lmax)],
      'clock_window':[1.5,2.5],'omega':float(omega),'contact_width':float(eps),'integration_steps':1024,
      'clock_probability_error_upper':float(ep),'per_clock_time_error_upper':float(et),
      'lorentz_term_lower':float(lorentz_lower),'all_four_time_error_upper':float(all_four_time_upper),
      'positive_margin_lower':float(margin),'exact_margin_fraction':str(margin),
      'samples_per_record':samples,'record_means':4,'total_record_samples':4*samples,
      'joint_sampling_failure_upper':8*math.exp(-7),'failure_bound_exact_series_verified':True,
      'length_velocity_calibration_bounds_assumed_not_empirically_measured':True,
      'independent_same_preparation_repetitions_assumed':True,'empirical_data':False,
      'uniform_exclusion_of_all_finite_c_at_fixed_error':False}

def main():
    result=run();p=HERE/'finite_certificate_results.json'
    if '--write' in sys.argv:
        with p.open('x',encoding='utf8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        saved=json.loads(p.read_text(encoding='utf8'))
        assert saved==result,'stored finite-certificate result changed'
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
