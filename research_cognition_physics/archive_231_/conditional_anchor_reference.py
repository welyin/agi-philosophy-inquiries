"""Round 531: genuine finite anchor records in the nonlinear quantum model.

The 27-site location guarantee is analytic, using the real relative Gibbs state;
the numerical quantum posterior diagnostic is explicitly one cell only. Geometry,
thermal preparation, marked anchors and instruments are supplied. The successful
source is extremely rare, and no macroscopic geometry or dimension is generated.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from nonlinear_thermal_reference import quantum_spectrum
from thermal_reference_dimension_model import thermal_q

HERE = Path(__file__).resolve().parent
TARGET = HERE/'conditional_anchor_reference_results.json'
THETA, GAIN = 2**20, 2**22
TAU2 = F(3, 8)


def fraction(x):
    x = F(x)
    return dict(numerator=x.numerator, denominator=x.denominator, decimal=float(x))


def graph():
    sites = list(itertools.product((-1, 0, 1), repeat=3))
    index = {x:i for i,x in enumerate(sites)}
    lap = 6*np.eye(27, dtype=int)
    for x,i in index.items():
        for axis in range(3):
            for direction in (-1, 1):
                y = list(x)
                y[axis] = (y[axis]+direction+1) % 3-1
                lap[i,index[tuple(y)]] -= 1
    signatures = np.zeros((27,3), dtype=int)
    for i,x in enumerate(sites):
        for axis in range(3):
            other_nonzero = sum(x[j]!=0 for j in range(3) if j!=axis)
            signatures[i,axis] = x[axis]*(26,5,2)[other_nonzero]
    anchors = []
    for axis in range(3):
        plus,minus = [0,0,0],[0,0,0]
        plus[axis],minus[axis] = 1,-1
        anchors.append((index[tuple(plus)],index[tuple(minus)]))
        rhs = np.zeros(27,dtype=int)
        rhs[anchors[-1][0]],rhs[anchors[-1][1]] = 162,-162
        assert np.array_equal(lap@signatures[:,axis],rhs)
    assert len(set(map(tuple,signatures))) == 27
    separation = min(int(np.max(abs(signatures[i]-signatures[j])))
                     for i in range(27) for j in range(i))
    assert separation == 2
    return sites,index,lap,signatures,anchors


def exact_geometry():
    sites,index,lap,signatures,anchors = graph()
    reflection = np.array([index[(x[0],-x[1],x[2])] for x in sites])
    assert np.array_equal(lap[np.ix_(reflection,reflection)],lap)
    assert all(reflection[a]==a for a in anchors[0])
    assert np.array_equal(signatures[reflection,0],signatures[:,0])
    return dict(sites=sites,signature_numerators=signatures.tolist(),denominator=162,
        anchors=anchors,minimum_linf_separation=fraction(F(1,81)),
        unique_three_component_signatures=27,
        unique_one_component_signatures=len(set(signatures[:,0])),
        single_pair_reflection_preserves_graph_and_anchors=True,
        dipole_equation_checked_in_exact_integers=True)


def common_covariance_check():
    _,_,lap,signatures,anchors = graph()
    eigen,rotation = np.linalg.eigh(lap.astype(float))
    eigen = np.rint(eigen).astype(int)
    assert set(eigen) == {0,3,6,9}
    active = eigen > 0
    classical = np.ones(27)
    classical[active] = 1/eigen[active]
    spectral = np.full(27,float(THETA))
    spectral[active] = thermal_q(eigen[active],THETA)
    covariance = (rotation*spectral)@rotation.T
    residual = spectral-THETA*classical
    assert np.min(residual)>-1e-9 and np.max(residual)<.5
    y0 = float(F(13*GAIN,81))*np.array([1.,-1.])
    rows=[]
    for axis,selected in enumerate(anchors):
        s = covariance[np.ix_(selected,selected)] + float(TAU2)*np.eye(2)
        regression = covariance[:,selected]@np.linalg.inv(s)
        drift = float(np.max(abs(regression@y0/GAIN-signatures[:,axis]/162)))
        assert np.linalg.eigvalsh(s)[0] >= THETA/9
        assert np.linalg.norm(regression,2)<10
        assert drift < 2/THETA
        rows.append(dict(axis=axis,S_min_eigenvalue=float(np.linalg.eigvalsh(s)[0]),
            regression_operator_norm=float(np.linalg.norm(regression,2)),
            ideal_record_common_mean_error=drift))
    return dict(nonzero_eigenvalues=[3,6,9],maximum_quantum_remainder=float(np.max(residual)),
        comparisons=rows,quantum_common_covariance=covariance.tolist())


def one_cell_quantum_posterior():
    grid,w,_,energies,vectors = quantum_spectrum(1.,size=32,halfwidth=8.)
    weights=np.exp(-(energies-energies[0])); weights/=sum(weights)
    diagonal=(vectors*vectors)@weights
    # This is the diagonal of exp(-H_rel), not a classical exp(-V) density.
    variable=np.sqrt(5)*w
    common_var,tau2,gamma=4.,float(TAU2),25.
    s=common_var+tau2
    k=common_var/s
    conditional_gaussian_var=common_var*tau2/s
    precision=np.diag([1/common_var,1/gamma])+np.ones((2,2))/tau2
    comparison=np.linalg.inv(precision)

    def posterior(y):
        likelihood=np.exp(-.5*(y-variable)**2/s)/math.sqrt(2*math.pi*s)
        probability=float(diagonal@likelihood)
        mass=diagonal*likelihood/probability
        mv=float(mass@variable)
        vv=float(mass@((variable-mv)**2))
        mu=k*(y-mv)
        cov=np.array([[conditional_gaussian_var+k*k*vv,-k*vv],[-k*vv,vv]])
        mx=mu+mv
        actual=conditional_gaussian_var+(1-k)**2*vv
        mgfs=[]
        for tilt in (-.4,.4):
            logmgf=(.5*tilt*tilt*conditional_gaussian_var+
                math.log(float(mass@np.exp(tilt*((1-k)*variable+k*y-mx)))))
            proxy=float(np.ones(2)@comparison@np.ones(2))
            assert logmgf<=.5*tilt*tilt*proxy+1e-5
            mgfs.append(dict(tilt=tilt,log_mgf=logmgf,upper=.5*tilt*tilt*proxy))
        assert np.linalg.eigvalsh(comparison-cov)[0]>-1e-5
        assert abs(mu-common_var/s*(y-mv))<1e-14
        return dict(y=float(y),density=probability,mean_U=mu,mean_V=mv,
            covariance=cov.tolist(),mean_X=mx,variance_X=actual,
            wrongly_dropped_cross_covariance_variance=float(np.trace(cov)),
            conditional_V_fourth_cumulant=float(mass@((variable-mv)**4)-3*vv*vv),
            centered_mgf=mgfs)

    rows=[posterior(y) for y in (-4.,-2.,0.,2.,4.)]
    assert all(row['covariance'][0][1]<-.1 for row in rows)
    nodes,quadrature=np.polynomial.legendre.leggauss(40)
    low,high=1.75,2.25
    samples=[posterior((low+high)/2+(high-low)*node/2) for node in nodes]
    masses=np.array([sample['density'] for sample in samples])*quadrature*(high-low)/2
    p=float(sum(masses)); conditional_weights=masses/p
    mean=float(conditional_weights@np.array([sample['mean_X'] for sample in samples]))
    variance=float(conditional_weights@np.array([sample['variance_X']+
                        (sample['mean_X']-mean)**2 for sample in samples]))
    assert 0<p<1
    center=rows[3]
    assert abs(variance-center['variance_X'])>1e-3
    # The equal-time nonselective position law is unchanged: a finite accepted
    # interval plus its complement sum to the identity effect at every X.
    xgrid=np.linspace(-10,10,201)
    cdf=lambda v: .5*(1+math.erf(v/math.sqrt(2)))
    effects=np.array([cdf((high-x)/math.sqrt(tau2))-cdf((low-x)/math.sqrt(tau2))
                      for x in xgrid])
    assert np.min(effects)>=0 and np.max(effects)<=1
    complete_error=float(np.max(abs(effects+(1-effects)-1)))
    return dict(full_relative_quantum_DVR=True,grid_size=32,halfwidth=8.,b=1.,beta=1.,
        samples=rows,finite_bin=dict(bounds=[low,high],probability=p,
            mean_X=mean,variance_X=variance,bin_center_variance=center['variance_X']),
        accepted_and_failed_effect_completeness_error=complete_error,
        numerical_scope='one cell; not a 27-site nonlinear quantum simulation')


def exact_certificate():
    theta,g=THETA,GAIN
    y=F(13*g,81)
    bias=F(4,theta)+F(40,g)+8000*(y+1)/(theta*g)
    assert bias<F(1,512)
    assert F(13,81)+bias<F(1,4)
    assert F(2*theta,3)+F(207,4)<theta
    noise=F(1,256)-F(1,512)-F(1,2048)
    assert noise==F(3,2048)
    exponent=g*g*noise*noise/(2*theta)
    assert exponent==18
    overflow_exponent=F(g*g,theta)*F(9,64)
    assert overflow_exponent>28 and 162<2**8
    immediate_failure=F(2106,2**18)+F(1,2**20)
    delay_trace=F(1,2**12)
    assert immediate_failure+delay_trace<F(1,100)
    m=69+2*math.ceil(27*(y+21)**2/theta)
    upper=math.floor((y-1)**2/(4*theta))
    assert (m,upper)==(23337735,108037)
    assert F(1,81)>2*F(1,256)
    return dict(theta=theta,gradient_scale=g,anchor_center=fraction(y),
        source_bin_halfwidth=1,relative_gamma_upper=25,
        normalized_mean_difference_error_bound=fraction(bias),
        allowed_mean_difference_error=fraction(F(1,512)),
        final_pair_tolerance=fraction(F(1,256)),
        final_noise_allowance=fraction(noise),pair_tail_exponent=fraction(exponent),
        overflow_tail_exponent=fraction(overflow_exponent),
        immediate_failure_upper=fraction(immediate_failure),
        delayed_failure_upper=fraction(immediate_failure+delay_trace),
        source_probability_lower_log2=-m,source_probability_upper_log2=-upper,
        delay_maximum_log2=-(m+300),conditional_energy_upper_log2=m+274,
        source_reads=6,source_ternary_record_bits=12,
        final_reads=81,final_symbols_including_overflow=4098,
        final_bits_per_read=13,final_raw_payload_bits=1053,
        source_nonselective_injected_energy=12,final_injected_energy=162,
        finite_source_probability_not_rounded_to_zero=True,
        probability_bounds_not_claimed_sharp=True)


def finite_decoder_check():
    _,index,_,signature,_=graph()
    truth=signature/162.
    drift=np.where((np.arange(81).reshape(27,3)%3)==0,1.,-1.)/1024
    noise=np.where(np.arange(81).reshape(27,3)%2,1.,-1.)*3/4096
    # A correlated adversarial good-event diagnostic, not thermal Monte Carlo.
    analog=truth+drift+noise+.05
    assert np.max(abs(analog))<1
    codes=np.rint((analog+1)*2048).astype(int)
    decoded=codes/2048.-1
    maximum_error=0.
    for i in range(27):
        for j in range(i):
            maximum_error=max(maximum_error,float(np.max(abs(decoded[i]-decoded[j]-truth[i]+truth[j]))))
    assert maximum_error<F(1,256)
    origin=index[(0,0,0)]
    recovered=[]
    for i in range(27):
        relative=decoded[i]-decoded[origin]
        distances=np.max(abs(truth-relative),axis=1)
        recovered.append(int(np.argmin(distances)))
    assert recovered==list(range(27))
    assert max(codes.ravel())<=4096 and min(codes.ravel())>=0
    assert (4098-1).bit_length()==13
    return dict(all_27_given_dictionary_locations_recovered=True,
        maximum_pair_error=maximum_error,final_payload_bits=int(codes.size*13),
        declared_overflow_code=4097,bounded_test_records_not_quantum_samples=True)


def energy_and_time_check():
    # A trace comparison proves finite accepted energy, without simulating the
    # 54-coordinate relative Gibbs state or treating it as a classical density.
    assert .5 < math.sqrt(47/4)<4
    assert 1/(2*math.sinh(.125))<4
    assert 1/(2*math.sinh(2))>1/8
    common_upper=78*THETA+120
    assert common_upper<2**27
    assert 3*2**271+common_upper+12<2**274
    # 2 E t <= 2**(-25), so half trace distance < 2**(-12).
    assert 1+274-300==-25 and -25 < -24
    return dict(relative_coordinates_per_copy=54,
        relative_energy_upper_log2=271,three_copy_common_energy_upper=common_upper,
        source_nonselective_total_energy_upper_log2=274,
        twice_energy_times_delay_upper_log2=-25,
        delayed_half_trace_distance_upper=fraction(F(1,4096)),
        no_nonselective_remote_signal_from_source_instrument=True,
        record_communication_and_repreparation_assumed_not_generated=True)


def run():
    geometry=exact_geometry()
    covariance=common_covariance_check()
    quantum=one_cell_quantum_posterior()
    certificate=exact_certificate()
    decoder=finite_decoder_check()
    energy=energy_and_time_check()
    dependencies=('research_note_529.md','nonlinear_thermal_reference.py',
        'research_note_530.md','thermal_reference_dimension_model.py',
        'local_reference_probe_review.md','research_note_287.md','research_note_520.md')
    return dict(date='2026-09-30',round=531,scientific_base_through_round=530,
        tests_run=6,failures=0,errors=0,exact_dipole_geometry=geometry,
        common_thermal_comparison=covariance,actual_quantum_posterior_diagnostic=quantum,
        finite_source_and_readout_certificate=certificate,decoder_diagnostics=decoder,
        energy_and_time_certificate=energy,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                           for name in dependencies},
        primary_tool_source='https://arxiv.org/pdf/1404.5886',
        scope=dict(same_fixed_nonlinear_gate=True,spatial_mean_not_prepared_pointwise=True,
            actual_finite_bin_CP_source=True,conditional_U_relative_correlations_retained=True,
            rare_source_probability_explicit=True,positive_but_tiny_delayed_window=True,
            supplied_geometry_and_marked_anchors=True,unequal_sector_temperatures_supplied=True,
            accepted_state_not_claimed_Gaussian_or_thermal=True,
            sustainable_reference_source_proved=False,original_graph_to_lattice_proved=False,
            unique_dimension_or_GR_derived=False,stage_complete=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert json.loads(TARGET.read_text('utf8'))==json.loads(json.dumps(result))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(round=531,tests=6,failures=0,
        source_probability_log2_bounds=[-23337735,-108037],
        delayed_success_lower=1-result['finite_source_and_readout_certificate']['delayed_failure_upper']['decimal'])))
