"""Round 1261: finite joint 1PM propagation and all-menu retiming deletion.
Default is read-only verification; --write creates results.json exclusively.
"""
from fractions import Fraction as F
from math import isqrt, pi, sin
import argparse
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def su(q, digits=100):
    q=F(q)
    if q < 0:
        raise ValueError('Negative radicand')
    scale=10**digits
    n=isqrt(q.numerator*scale**2//q.denominator)
    if F(n*n,scale**2)<q:
        n+=1
    result=F(n,scale)
    assert result*result>=q
    return result

def rat(q):
    q=F(q)
    return str(q.numerator)+'/'+str(q.denominator)

def cap_mass(q,d):
    q=F(q)
    return F(1) if q<=d else max(F(0),2-q/d)

def cap_photon(e,k):
    e=F(e)
    if e<=k/2 or e>=2*k: return F(0)
    if e<2*k/3: return 6*e/k-3
    if e<=4*k/3: return F(1)
    return 3-3*e/(2*k)

def calculate():
    R,mu,k,sig=F(10**40),F(10**34),F(10**12),F(10**12)
    ph,d,om,tc,n=F(1,10**6),F(1,10**4),F(1000),F(300000),40000
    lam=1-ph
    pil,piu=F(333,106),F(355,113)
    tl=(2*n+F(1,2))*pil/lam
    tu=(2*n+F(1,2))*piu/lam
    assert tu<tc and ph*tl>F(1,4) and ph*tu<F(13,50)<pil/12
    # The entire good region lies in the photon plateau, including transverse momentum.
    assert 3*k/4>=2*k/3
    assert (5*k/4)**2+1<(4*k/3)**2
    for e,expected in ((k/2,0),(2*k/3,1),(4*k/3,1),(2*k,0)):
        assert cap_photon(e,k)==expected
    assert cap_mass(d,d)==1 and cap_mass(2*d,d)==0
    BV=1000*ph*k
    assert mu+26-BV>0
    D2=12*sig**2+tc**2+2
    mass8=F(945,16)/sig**8
    trans8=F(24)/sig**8
    p_other=768*D2/R**2+2*mass8/d**8+trans8
    long8=128*((om/2)**8+F(105,256))
    a=k/4
    p_long=long8/a**8
    tail_other=(BV+ph*(32+5*k/4))*su(p_other)
    tail_long=(BV+32*ph)*su(p_long)+ph*(k*su(p_long)+su(long8/a**6))
    good_const=(8*ph*d*d/26+4*ph/k+1000*k/R+
                32*d*k/R+1000*d*d*k/(mu*R)+1000*d*d/(mu*R))
    spatial=(64*ph+2*ph*k)/R*su(D2)
    V_error=tc*(good_const+spatial+tail_other+tail_long)
    free_error=tc*((F(1,26)+1/mu)/(2*sig**2)+1/((k-om/2)*sig**2))
    eps=V_error+free_error
    assert eps<F(1,10**8)
    tail_signal=640/om**2
    joint_lower=F(3,5)*(1-tail_signal)-tail_signal-2*eps
    clock_lower=F(3,32)-eps-free_error
    signal_lower=1-2*tail_signal-eps-free_error
    assert joint_lower>F(59,100)
    assert clock_lower>F(9,100)
    assert signal_lower>F(998,1000)
    # Exact generator coefficient identity proves a whole-operator statement.
    # No large global phase is evaluated numerically.
    hstar={'I':mu,'M':lam,'Pz':lam}
    free_retimed={'I':lam*mu,'M':lam,'Pz':lam}
    phase_generator=mu*(1-lam)
    assert hstar['I']==free_retimed['I']+phase_generator
    assert hstar['M']==free_retimed['M']
    assert hstar['Pz']==free_retimed['Pz']
    assert lam*tu<tc
    retime_error=eps+free_error
    assert retime_error<F(1,10**9)
    T=(2*pi*n+pi/2)/float(lam)
    out={
      'round':1261,'status':'PASS','author_scientific_groups':2,
      'groups':[
        {'id':'joint_finite_budget','status':'PASS',
         'meaning':'Full three-body completed 1PM quantum propagation, bounded energy exchange, fixed-coordinate joint effects.'},
        {'id':'all_menu_common_retiming_deletion','status':'PASS',
         'meaning':'Exact common-comparator identity plus a uniform finite bound against true free propagation at lambda*T.'}],
      'parameters':{'R':rat(R),'mu':rat(mu),'k':rat(k),'sigma_massive_transverse':rat(sig),
        'massive_momentum_good_radius':rat(d),'longitudinal_bandwidth':rat(om),
        'T_definition':'(2*pi*40000+pi/2)/(1-1e-6)',
        'fixed_constants':'hbar=c=G=1 in one fixed unit system',
        'completion':'u=1/max(r,R/4); massive linear hats d to 2d; photon trapezoid k/2,2k/3,4k/3,2k'},
      'strict_rational':{'potential_norm_bound':rat(BV),
        'Husimi_position_second_bound':rat(D2),'other_bad_probability_bound':rat(p_other),
        'long_bad_probability_bound':rat(p_long),'good_velocity_and_pair_bound':rat(good_const),
        'spatial_L2_bound':rat(spatial),'potential_evolution_bound':rat(V_error),
        'free_evolution_bound':rat(free_error),'total_state_bound':rat(eps),
        'midplane_signal_tail_bound':rat(tail_signal),
        'complete_free_energy_exchange_bound':rat(2*BV),
        'energy_exchange_relative_source_baseline':rat(2*BV/mu),
        'joint_probability_gap_lower':rat(joint_lower),'clock_gap_lower':rat(clock_lower),
        'signal_gap_lower':rat(signal_lower),'all_menu_retiming_probability_upper':rat(retime_error)},
      'retiming':{'lambda':rat(lam),'central_phase_generator':rat(phase_generator),
        'hstar_coefficients':{z:rat(v) for z,v in hstar.items()},
        'free_retimed_coefficients':{z:rat(v) for z,v in free_retimed.items()},
        'identity':'U_star(T)=exp(-i*mu*(1-lambda)*T)*U_free_comparator(lambda*T)',
        'uniform_over_all_bounded_joint_effects_and_reference':True,
        'actual_retimed_trigger_implemented':False},
      'display_only':{'T':T,'delta':float(ph)*T,'total_bound':float(eps),
        'joint_lower':float(joint_lower),'clock_lower':float(clock_lower),
        'signal_lower':float(signal_lower),'retiming_upper':float(retime_error),
        'clock_reference_probability':sin(3*pi/8)**6},
      'scope':{'anti_wick_quantization_input':True,'final_canonical_material_dictionary_input':True,
        'original_material_identity_inherited_from_1253_not_new_group':True,
        'all_three_positions_and_recoil_retained':True,'all_internal_reference_states_uniform':True,
        'actual_arrival_instrument':False,'pre_reduction_dictionary_transport':False,
        'independent_geometry_or_curvature_witness':False,
        'full_EH_source_insertion_error':False,'beyond_1PM_error':False,'GR_complete':False},
      'source_hashes':{'proof.md':sha(HERE/'proof.md'),'check.py':sha(Path(__file__)),
        '../1253/proof.md':sha(HERE.parent/'1253'/'proof.md'),
        '../1241/proof.md':sha(HERE.parent/'1241'/'proof.md'),
        '../_admission/same_source_signal_bridge_after1260/README.md':sha(HERE.parent/'_admission'/'same_source_signal_bridge_after1260'/'README.md')}
    }
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true',help='Create results.json exclusively; refuse overwrite.')
    args=parser.parse_args()
    actual=calculate()
    result_path=HERE/'results.json'
    if args.write:
        with result_path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(actual,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('WRITE_ONCE_PASS_1261: 2 author groups; exclusive results.json created')
    else:
        expected=json.loads(result_path.read_text(encoding='utf-8'))
        if actual!=expected:
            raise AssertionError('Saved 1261 result or source hashes differ; no files changed')
        print('READ_ONLY_MATCH_1261: 2 author groups; no files changed')
    print(json.dumps(actual['display_only'],ensure_ascii=False,sort_keys=True))

if __name__=='__main__':
    main()