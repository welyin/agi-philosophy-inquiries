"""Round 1145; read-only by default. --write stores this round's results only.
Finite matrices and coherent Gram identities diagnose the proof, not a field simulation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
GAMMA, VELOCITY = F(5, 4), F(3, 5)
D, E = F(864000, 10**9), F(22, 10**11)


def atan_bounds(x, n):
    s = sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)), F(0))
    other = s + (-1)**n*x**(2*n+1)/F(2*n+1)
    return min(s, other), max(s, other)


def group_provenance_and_rational_bounds():
    manifest = json.loads((HERE/'source_hashes.json').read_text(encoding='utf-8'))
    for rel, digest in manifest.items():
        assert hashlib.sha256((HERE/rel).resolve().read_bytes()).hexdigest() == digest, rel
    old = json.loads((HERE/'../1144/results.json').read_text(encoding='utf-8'))
    assert old['status'] == 'PASS'
    assert F(old['continuous_source_isometry_upper']) == D+E
    a, b = atan_bounds(F(1,5),20)
    c, d = atan_bounds(F(1,239),20)
    pi_lo, pi_hi = 16*a-4*d, 16*b-4*c
    assert F(314159,100000) < pi_lo < pi_hi < F(22,7)
    assert pi_hi-pi_lo < F(1,10**28)
    x = pi_lo/16
    sin_lo = x-x**3/6+x**5/120-x**7/5040
    # g'(x) >= 1-x^2/2-x^6/720 >0 for 0<x<1/5.
    assert pi_hi/16 < F(1,5)
    assert 1-F(1,5)**2/2-F(1,5)**6/720 > 0
    gap_lo = (1+sin_lo)/2-(1+GAMMA**2)*D-(1+GAMMA)*E
    assert gap_lo > F(59533,100000) > F(59,100)
    return {
        'frozen_sources': manifest,
        'pi_interval_width_upper': str(F(1,10**28)),
        'phase_error_sum_exact': str((1+GAMMA**2)*D),
        'mediator_error_sum_exact': str((1+GAMMA)*E),
        'gap_lower_rational': str(F((gap_lo*10**15).numerator//(gap_lo*10**15).denominator,10**15)),
        'gap_lower_float_view': float(gap_lo),
        'isometry_correct_upper': str(D+E),
        'isometry_wrong_target_upper': str(GAMMA**2*D+GAMMA*E),
        'strict_operational_gap_gt': '0.59533'
    }


def boost_matrix():
    g, v = float(GAMMA), float(VELOCITY)
    L = np.eye(4)
    L[0,0] = L[3,3] = g
    L[0,3] = L[3,0] = g*v
    return L


def group_source_geometry():
    L = boost_matrix()
    eta = np.diag([1.,-1.,-1.,-1.])
    assert np.max(np.abs(L.T@eta@L-eta)) < 1e-14
    assert abs(np.linalg.det(L)-1) < 1e-14
    g, v = float(GAMMA), float(VELOCITY)
    # Independent finite polynomial quadrature for the coordinate identity only.
    # Not the physical bump, not a new field or the physical T=1e9 evolution.
    z, w = np.polynomial.legendre.leggauss(16)
    rho = 15/16*(1-z*z)**2
    Tdiag, tc = 7., 0.43*7.*g
    profile = lambda u: np.where((u>=0)&(u<=1),30*u*u*(1-u)**2,0.)
    xpz = v*tc+z/g
    rest_tau = g*(tc-v*xpz)
    tau_texture = tc/g-v*z
    assert np.max(np.abs(rest_tau-tau_texture)) < 1e-14
    lhs = np.sum(w/g*profile(rest_tau/Tdiag)*rho)
    rhs = np.sum(w*profile(tau_texture/Tdiag)*rho)/g
    wrong = np.sum(w*profile(tau_texture/Tdiag)*g*rho)/g
    frozen_center = profile(tc/(g*Tdiag))/g
    assert abs(lhs-rhs) < 1e-14
    assert abs(wrong-g*lhs) < 1e-14
    assert abs(frozen_center-lhs) > 1e-3
    assert abs(np.sum(w*rho)-1) < 1e-14
    # Two coordinates with spacelike separation retain its exact sign under boost.
    q = np.array([0.25,0.3,0.2,0.7])
    q2 = L@q
    assert q@eta@q < 0 and abs(q2@eta@q2-q@eta@q) < 1e-14
    return {'Lorentz_metric_error':float(np.max(np.abs(L.T@eta@L-eta))),
            'time_texture_error':float(np.max(np.abs(rest_tau-tau_texture))),
            'correct_spatial_integral':float(lhs),
            'wrong_spatial_integral':float(wrong),
            'center_time_freezing_difference':float(frozen_center-lhs),
            'diagnostic_scope':'coordinate Jacobian and simultaneity texture only'}


def coherent_gram(hs):
    return np.array([[np.exp(-np.vdot(a,a).real/2-np.vdot(b,b).real/2+np.vdot(a,b))
                      for b in hs] for a in hs], dtype=complex)


def branch_isometry(hA,hB,alphaA,alphaB,theta):
    eye = np.eye(2)
    plus = np.array([1,1])/math.sqrt(2)
    minus = np.array([1,-1])/math.sqrt(2)
    branches=[(a,b,np.kron(eye[:,a],bv)) for a in (0,1) for b,bv in ((1,plus),(-1,minus))]
    hs=[np.zeros_like(hA)]+[a*hA+b*hB for a,b,_ in branches]
    gram=coherent_gram(hs)
    ev,U=np.linalg.eigh(gram)
    assert ev.min()>-1e-12
    factor=(U*np.sqrt(np.maximum(ev,0)))@U.conj().T
    V=np.zeros((20,4),complex)
    V0=np.zeros_like(V)
    for j,(a,b,vec) in enumerate(branches,1):
        phase=np.exp(1j*(alphaA*a*a+alphaB*b*b+theta*a*b))
        proj=np.outer(vec,vec.conj())
        V += phase*np.kron(proj,factor[:,j,None])
        V0 += phase*np.kron(proj,factor[:,0,None])
    return V,V0,gram


def group_full_record_algebra():
    hA=np.array([.041+.01j,-.027j,.017])
    hB=np.array([-.009j,.033+.007j,-.012])
    alphaA,alphaB,theta=.7,-1.2,math.pi/2+.0003
    V,V0,gram=branch_isometry(hA,hB,alphaA,alphaB,theta)
    iso_error=float(np.linalg.norm(V.conj().T@V-np.eye(4),2))
    displacement_bound=float(np.linalg.norm(hA)+np.linalg.norm(hB))
    comparison=float(np.linalg.norm(V-V0,2))
    assert iso_error<2e-13 and comparison<=displacement_bound+2e-12
    # An arbitrary joint input and a 3-dimensional passive reference; no input selection.
    rng=np.random.default_rng(1145)
    psi=rng.normal(size=12)+1j*rng.normal(size=12)
    psi/=np.linalg.norm(psi)
    ref_error=float(np.linalg.norm(np.kron(V-V0,np.eye(3))@psi))
    assert ref_error<=comparison+1e-12
    # A generic one-particle unitary checks the Gram representation of the exact
    # intertwiner. It is NOT a finite-dimensional implementation of a boost.
    Z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    Q,_=np.linalg.qr(Z)
    _,_,boost_gram=branch_isometry(Q@hA,Q@hB,alphaA,alphaB,theta)
    assert np.max(np.abs(boost_gram-gram))<1e-13
    g=float(GAMMA)
    bad,_,_=branch_isometry(g*hA,g*hB,g*g*alphaA,g*g*alphaB,g*g*theta)
    input10=np.array([0,0,1,0],complex)
    out=(V@input10).reshape(2,2,5)
    bout=(bad@input10).reshape(2,2,5)
    p=float(np.sum(np.abs(out[:,1,:])**2))
    pb=float(np.sum(np.abs(bout[:,1,:])**2))
    hp,hm=hA+hB,hA-hB
    overlap=coherent_gram([hp,hm])[0,1]
    p_formula=(1-np.real(np.exp(-2j*theta)*overlap))/2
    overlap_bad=coherent_gram([g*hp,g*hm])[0,1]
    pb_formula=(1-np.real(np.exp(-2j*g*g*theta)*overlap_bad))/2
    assert abs(p-p_formula)<1e-13 and abs(pb-pb_formula)<1e-13
    assert abs(p-math.sin(theta)**2)<=displacement_bound
    assert abs(pb-math.sin(g*g*theta)**2)<=g*displacement_bound
    # General data coherence retains alphaA. Removing it changes the joint output.
    zero=np.zeros(3,complex)
    exact,_,_=branch_isometry(zero,zero,alphaA,alphaB,math.pi/2)
    removed,_,_=branch_isometry(zero,zero,0.,alphaB,math.pi/2)
    superposed=np.array([1,0,1,0],complex)/math.sqrt(2)
    overlap_unknown=np.vdot(exact@superposed,removed@superposed)
    trace_distance=math.sqrt(max(0.,1-abs(overlap_unknown)**2))
    assert trace_distance>.34
    return {'full_isometry_error':iso_error,'all_input_displacement_error':comparison,
            'displacement_bound':displacement_bound,'passive_reference_error':ref_error,
            'coherent_gram_transport_error':float(np.max(np.abs(boost_gram-gram))),
            'diagnostic_pointer_correct':p,'diagnostic_pointer_wrong':pb,
            'dropping_self_phase_joint_trace_distance':trace_distance,
            'diagnostic_scope':'exact finite coherent Gram algebra, not continuum simulation'}


def group_four_momentum_and_energy():
    L=boost_matrix()
    spatial=np.array([[.2,-.1,.4],[-.3,.7,-.2],[.6,.1,.8],[-.2,-.5,.3]])
    k=np.column_stack((np.sqrt(1+np.sum(spatial**2,axis=1)),spatial))
    kp=k@L.T
    # Invariant mass-shell measure transported independently of displacement norm.
    measure=np.array([.17,.31,.23,.29])
    amp=np.array([.3+.2j,-.4+.1j,.11-.3j,.2-.17j])
    jac=kp[:,0]/k[:,0]
    measurep=measure*jac
    ampp=amp/np.sqrt(jac)
    weight=measure*np.abs(amp)**2
    weightp=measurep*np.abs(ampp)**2
    mean=weight@k
    meanp=weightp@kp
    covariance=np.einsum('a,ai,aj->ij',weight,k,k)
    covariancep=np.einsum('a,ai,aj->ij',weightp,kp,kp)
    err=float(np.max(np.abs(meanp-L@mean)))
    cerr=float(np.max(np.abs(covariancep-L@covariance@L.T)))
    assert err<1e-14 and cerr<1e-14
    assert np.max(np.abs(kp[:,0]**2-np.sum(kp[:,1:]**2,axis=1)-1))<1e-14
    assert np.all(kp[:,0]<=2*k[:,0])
    assert np.max(np.abs((float(GAMMA)*ampp)**2-float(GAMMA)**2*ampp**2))<1e-14
    rest=E**2
    boosted=2*rest
    mismatch=(GAMMA**2-1)*boosted
    assert rest==F(484,10**22)
    assert boosted==F(968,10**22)
    assert mismatch==F(5445,10**23)
    return {'four_momentum_transport_error':err,'connected_covariance_transport_error':cerr,
            'normal_rest_field_energy_upper':str(rest),
            'normal_boosted_field_energy_upper':str(boosted),
            'wrong_source_energy_difference_upper':str(mismatch),
            'diagnostic_scope':'mass-shell measure and coherent-state moments only',
            'whole_apparatus_four_momentum_closed':False}


def run():
    return {'round':1145,'status':'PASS','check_groups':4,
            'groups':{'frozen_sources_and_exact_budget':group_provenance_and_rational_bounds(),
                      'source_geometry':group_source_geometry(),
                      'full_record_algebra':group_full_record_algebra(),
                      'field_four_momentum':group_four_momentum_and_energy()},
            'scalar_source_type_is_model_input':True,
            'all_source_self_phases_retained':True,
            'proper_time_texture_retained':True,
            'common_source_family_transport':True,
            'autonomous_local_material_realization':False,
            'whole_apparatus_four_momentum_closed':False,
            'new_axioms':0,'scientific_calibration_increment':0,
            'finite_matrices_are_not_continuum_simulations':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=run()
    encoded=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
    if args.write:
        (HERE/'results.json').write_text(encoded,encoding='utf-8',newline='\n')
    print(encoded)
