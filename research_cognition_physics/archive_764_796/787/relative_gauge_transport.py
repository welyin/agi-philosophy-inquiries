"""787: general free gauge kernels and anchored formal relative BRST repair.

Finite matrices and exact polynomial complexes diagnose analytic arguments in
research_note_787.md. No continuum endpoint or interacting positivity claim.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import itertools
import json
import numpy as np
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('bv777_for787',HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv=importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)


def free_gauge_family():
    n=5
    eye=np.eye(n)
    d=np.diag(np.linspace(.7,1.3,n))@(np.roll(eye,1,axis=1)-eye)
    k=np.vstack((d,eye))
    f=np.hstack((eye,-d))
    l=np.hstack((np.zeros((n,n)),eye))+.11*f
    p=f.T@(eye+d.T@d)@f
    pi=np.eye(2*n)-k@l
    gp=pi@np.linalg.inv(p-k@k.T)@pi.T
    skew=np.roll(eye,1,axis=1)-np.roll(eye,-1,axis=1)
    param=eye+.17*skew
    ns=(-eye,eye+.3*d@d.T)
    cases=[]
    wrong_errors=[]
    for which,N in enumerate(ns):
        for t in (0.,.25,.5,.75,1.):
            F=(1-t)*param@k.T+t*l
            A=F@k
            ia=np.linalg.inv(A)
            T=np.eye(2*n)-k@ia@F
            top=T@gp@T.T-k@ia@N@ia.T@k.T
            G=np.block([[top,k@ia],[ia.T@k.T,np.zeros((n,n))]])
            D=np.block([[p,F.T],[F,N]])
            residual=max(np.max(np.abs(D@G-np.eye(3*n))),np.max(np.abs(G@D-np.eye(3*n))))
            physical=np.max(np.abs(pi@top@pi.T-gp))
            assert residual<5e-12 and physical<5e-12
            # A is not assumed self-adjoint; using A^{-1} for both legs fails.
            wrong=np.block([[T@gp@(np.eye(2*n)-F.T@ia@k.T)-k@ia@N@ia@k.T,k@ia],
                            [ia@k.T,np.zeros((n,n))]])
            wrong_errors.append(float(np.max(np.abs(D@wrong-np.eye(3*n)))))
            cases.append(dict(auxiliary_pairing_case=which,t=t,
                              inverse_residual=float(residual),physical_projection_residual=float(physical)))
    assert max(wrong_errors)>.1
    return dict(block_dimension=15,cases=cases,
                nonsymmetric_ghost_inverse_negative_control=max(wrong_errors),
                scope='Matrix calibration of both causal legs and general auxiliary pairing; not a continuum UV bound.')


v={name:bv.var(name) for name in bv.NAMES}
S=bv.add(bv.mul(v['p'],v['c']),bv.scale(bv.mul(v['a'],v['b']),-1))


def clean(p):
    return {i:x for i,x in p.items() if x}


def add(*polys):
    result={}
    for poly in polys:
        for i,value in poly.items():
            result[i]=bv.add(result.get(i,{}),value)
    return clean(result)


def scale(poly,c):
    return clean({i:bv.scale(value,c) for i,value in poly.items()})


def s(poly):
    return clean({i:bv.bracket(S,value) for i,value in poly.items()})


def deriv(poly):
    return clean({i-1:bv.scale(value,i) for i,value in poly.items() if i})


def integral(poly):
    return clean({i+1:bv.scale(value,Q(1,i+1)) for i,value in poly.items()})


def D(pair):
    a,b=pair  # chi is written on the LEFT, hence the minus sign in -s b.
    return s(a),add(deriv(a),scale(s(b),-1))


def H(pair):
    return integral(pair[1]),{}


def pair_add(*pairs):
    return add(*(p[0] for p in pairs)),add(*(p[1] for p in pairs))


def relative_homotopy():
    assert not bv.bracket(S,S)
    monomials=[bv.ONE]+list(v.values())
    for a,b in itertools.combinations_with_replacement(bv.NAMES,2):
        value=bv.mul(v[a],v[b])
        if value:
            monomials.append(value)
    checked=0
    for monomial in monomials:
        for degree in range(5):
            for slot in (0,1):
                pair=({degree:monomial},{}) if slot==0 else ({},{degree:monomial})
                expected=pair
                if slot==0 and degree==0:
                    expected=({}, {})
                assert D(D(pair))==({}, {})
                assert pair_add(D(H(pair)),H(D(pair)))==expected
                checked+=1
    r4=bv.prod(v['r'],v['r'],v['r'],v['r'])
    beta=Q(2,3)
    # The physical term survives the ORIGINAL complex's quartet projection.
    physical_anomaly=({}, {0:bv.scale(r4,beta)})
    assert D(physical_anomaly)==({}, {})
    repair=tuple(scale(x,-1) for x in H(physical_anomaly))
    assert pair_add(D(repair),physical_anomaly)==({}, {})
    assert repair[0].get(0,{})=={}
    assert repair[0][1]==bv.scale(r4,-beta)
    x=bv.prod(v['q'],v['r'],v['r'])
    y=bv.prod(v['h'],v['c'],v['r'])
    trial=({1:x,2:y},{})
    anomaly=D(trial)
    assert anomaly[0] and anomaly[1]
    assert D(anomaly)==({}, {})
    assert H(anomaly)==trial
    # H cannot remove an unrepaired anomaly at the anchor while preserving it.
    base=({0:bv.prod(v['c'],v['r'],v['r'])},{})
    assert D(base)==({}, {}) and H(base)==({}, {})
    return dict(coefficient_arithmetic='fractions.Fraction',external_parameter_degree_through=4,
                monomials=len(monomials),cochain_cases=checked,
                D_squared_zero=True,DH_plus_HD_equals_identity_minus_anchor=True,
                physical_ghost_zero_term='(2/3)*r^4',anchored_counterterm='-(2/3)*t*r^4',
                no_assumption_of_original_H0_vanishing=True,
                mixed_anomaly_components_repaired=True,
                unrepaired_anchor_anomaly_not_removed=True,
                scope='Exact external-parameter relative complex diagnostic, not an original anomaly coefficient calculation.')


def nonuniform_endpoint():
    rows=[]
    for delta in (Q(1,2),Q(1,10),Q(1,100),Q(1,1000)):
        t=1-delta
        spectral_value=delta**-2
        inverse=1/(delta*spectral_value+t)
        error=1-inverse
        assert error>Q(0)
        rows.append(dict(t=str(t),spectral_value=str(spectral_value),
                         inverse=str(inverse),difference_from_endpoint=str(error)))
    assert Q(rows[-1]['difference_from_endpoint'])>Q(999,1000)
    # Fixed spectral values DO converge: no claim of failure of strong convergence.
    fixed=Q(7)
    fixed_errors=[abs(1/(delta*fixed+1-delta)-1) for delta in (Q(1,10),Q(1,100),Q(1,1000))]
    assert fixed_errors[0]>fixed_errors[1]>fixed_errors[2]
    return dict(family='A_t=(1-t) A_0+t I, spectrum(A_0) subset [1,infinity)',
                fixed_spectrum_limit_is_identity=True,operator_norm_convergence=False,
                strong_or_distributional_convergence_disproved=False,
                increasing_spectral_sequence=rows,
                scope='Unbounded positive spectral counterexample to uniform endpoint inference; not a Lorentzian propagator model.')


def run():
    return dict(round=787,fresh_test_groups=3,free_gauge_family=free_gauge_family(),
                anchored_relative_complex=relative_homotopy(),endpoint=nonuniform_endpoint(),
                all_checks_passed=True,actual_auxiliary_endpoint_quantum_matching_proven=False,
                original_interacting_positive_state_proven=False)


if __name__=='__main__':
    result=run()
    HERE.joinpath('relative_gauge_transport_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
