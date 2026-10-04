"""577: original matter transduces a local phase to an adjacent bounded readout.
Full-graph statements are analytic; the numerical quantum grid is a labelled radial diagnostic.
No new pointer field, interaction term, quantum geometry or continuum limit is introduced.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_original_matter_record_results.json'
spec=importlib.util.spec_from_file_location('record_entry',HERE/'round577_drafts/radial_operator_entry_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
M=original.M
_,u,_=original.lattice.scalar.parameters()
HSTAR,SSTAR=np.sqrt(u)
CENTER=np.array([0.,HSTAR,0.,0.,SSTAR])


def mixed_edge_hessian(x,y,R=None,k=1.):
    if R is None:R=np.eye(5)
    v=R@y;fx=original.F(x);fy=original.F(v);den=np.sqrt(fx*fy)
    z=(np.dot(x-v,x-v)/12+(np.sqrt(fx)-np.sqrt(fy))**2/2)/den
    C=1+z;r,dr=entry.ratio(z)
    CA=-v/(6*den)+C*x/(6*fx)
    CB=-x/(6*den)+C*v/(6*fy)
    CAB=-np.eye(5)/(6*den)-np.outer(v,v)/(36*fy*den)+np.outer(x,CB)/(6*fx)
    return 6*k*(dr*np.outer(CA,CB)+r*CAB)@R


def response_coefficient(x,y,R=None,k=1.,wA=1.,wB=1.):
    dT=x.copy();dT[4]=0
    db=np.array([0.,0.,0.,0.,np.cos(y[4])])
    return -float((original.inverse(x)@dT)@mixed_edge_hessian(x,y,R,k)@(original.inverse(y)@db))/(wA*wB)


def representation(vector,phase):
    W=original.lattice.old.group_exp(np.array(vector),2)*np.exp(3j*phase)
    R=np.eye(5);R[:4,:4]=np.block([[W.real,-W.imag],[W.imag,W.real]])
    return R


def finite_mixed_hessian(x,y,R,step):
    out=np.zeros((5,5));unit=np.eye(5)*step
    def energy(a,b):return float(original.distance_squared(a,R@b)/2)
    for i in range(5):
        for j in range(5):
            out[i,j]=(energy(x+unit[i],y+unit[j])-energy(x+unit[i],y-unit[j])
                -energy(x-unit[i],y+unit[j])+energy(x-unit[i],y-unit[j]))/(4*step*step)
    return out


def endpoint_and_gauge_check():
    kin=original.metric(CENTER);hess=mixed_edge_hessian(CENTER,CENTER)
    zero_error=float(np.max(abs(hess+kin)))
    value=response_coefficient(CENTER,CENTER)
    expected=-original.F(CENTER)*SSTAR*HSTAR**2*np.cos(SSTAR)/(6*M)
    assert zero_error<1e-14 and abs(value-expected)<1e-14 and value<0
    x=CENTER+np.array([.014,-.021,.013,.004,.018])
    y=CENTER+np.array([-.012,.017,-.009,.011,-.019])
    R=representation((.12,-.21,.09),.04)
    exact=mixed_edge_hessian(x,y,R)
    errors=[float(np.max(abs(finite_mixed_hessian(x,y,R,s)-exact))) for s in (1e-3,3e-4)]
    assert errors[-1]<errors[0]/5 and errors[-1]<2e-7
    GA=representation((.3,.7,-.2),-.13);GB=representation((-.5,.1,.4),.21)
    gauge_error=abs(response_coefficient(GA@x,GB@y,GA@R@GB.T)-response_coefficient(x,y,R))
    assert gauge_error<1e-14
    assert np.max(abs(GA.T@GA-np.eye(5)))<1e-14
    rng=np.random.default_rng(577)
    nearby=[response_coefficient(CENTER+.004*rng.normal(size=5),
              CENTER+.004*rng.normal(size=5),representation(.01*rng.normal(size=3),.003*rng.normal())) for _ in range(25)]
    assert max(nearby)<0
    dT=CENTER.copy();dT[4]=0
    cost=float(dT@original.inverse(CENTER)@dT)
    cost_formula=float(original.F(CENTER)*(HSTAR**2-HSTAR**4/(6*M)))
    assert abs(cost-cost_formula)<1e-14 and cost>0
    return dict(vacuum_response_coefficient=value,formula=expected,
        mixed_diagonal_Hessian_error=zero_error,off_diagonal_finite_difference_errors=errors,
        gauge_covariance_error=gauge_error,nearby_sample_range=[min(nearby),max(nearby)],
        preparation_gradient_norm=cost,strict_sign_neighbourhood_proved_by_continuity_not_sampling=True)


def radial_operator(q,hbar):
    dx=q['dx']
    def derivative(f,axis):return (np.roll(f,-1,axis=axis)-np.roll(f,1,axis=axis))/(2*dx)
    def hamiltonian(psi):
        answer=q['V']*psi
        for start,kin,density in ((0,q['kx'],q['mux']),(2,q['ky'],q['muy'])):
            gradient=[derivative(psi,start+j) for j in range(2)]
            divergence=sum(derivative(density*sum(kin[...,i,j]*gradient[j] for j in range(2)),
                                      start+i) for i in range(2))
            answer=answer-hbar*hbar*divergence/(2*density)
        return answer
    def inner(a,b):return np.sum(q['mu']*a.conj()*b)*dx**4
    return hamiltonian,inner


def radial_quantum_comparison_check():
    q=entry.arrays(19);hbar=.07;theta=.2;H,inner=radial_operator(q,hbar)
    rows=[]
    for angle in (0.,theta,-theta):
        psi=q['amp']*np.exp(1j*angle*q['T']/hbar)
        h1=H(psi);h2=H(h1);h3=H(h2);b=q['b']
        first=-2*inner(h1,b*psi).imag/hbar
        second=-2*(inner(h2,b*psi).real-inner(h1,b*h1).real)/(hbar*hbar)
        third=(2*inner(h3,b*psi).imag-6*inner(h2,b*h1).imag)/(hbar**3)
        rows.append(dict(theta=angle,norm=float(inner(psi,psi).real),
            probability_zero=float(.5+inner(psi,b*psi).real/4),
            first=float(first),second=float(second),third=float(third),energy=float(inner(psi,h1).real)))
    baseline,plus,minus=rows
    first_two=max(abs(row[key]-baseline[key]) for row in rows[1:] for key in ('first','second'))
    same_initial=abs(plus['probability_zero']-minus['probability_zero'])
    equal_energy=abs(plus['energy']-minus['energy'])
    assert first_two<1e-9 and same_initial<1e-14 and equal_energy<1e-13
    assert plus['energy']>baseline['energy'] and minus['third']>0>plus['third']
    assert abs(plus['third']+minus['third'])<1e-12
    cbar=float(np.sum(q['mu']*q['amp']**2*q['coefficient'])*q['dx']**4)
    normalized_response=(plus['third']-minus['third'])/4
    expected_response=theta*cbar/2
    support=q['coefficient'][q['amp']>0]
    assert support.max()<0
    return dict(rows=rows,first_two_derivative_difference=first_two,
        initial_reader_probability_difference=same_initial,opposite_phase_energy_difference=equal_energy,
        common_extra_energy=plus['energy']-baseline['energy'],
        third_derivative_of_probability_difference=float(normalized_response),
        continuum_prediction_for_same_radial_quadrature=expected_response,
        sampled_support_coefficient_range=[float(support.min()),float(support.max())],
        effect='E_plus=1/2+sin(s_B)/4',
        radial_diagnostic_not_full_gauge_state_or_time_simulation=True)


def inherited_operator_check():
    saved=json.loads((HERE/'round577_drafts/radial_operator_entry_probe_results.json').read_text('utf8'))
    current=entry.run();assert saved==current
    rows=current['rows']
    assert rows[-1]['relative_third_derivative_error']<.08
    assert rows[-1]['relative_third_derivative_error']<rows[0]['relative_third_derivative_error']/3
    return dict(original_entry_results_reproduced=True,rows=rows,
                coarse_grid_not_exact_continuum_or_full_Gauss_proof=True)


def run():
    checks=('endpoint_and_gauge_check','radial_quantum_comparison_check','inherited_operator_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_curved_quantum_source.py','joint_pointer_matter_compatibility.py',
          'round577_drafts/radial_operator_entry_probe.py',
          'round577_drafts/radial_operator_entry_probe_results.json','research_round_576_checks.json')
    return dict(round=577,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='fixed finite graph and fixed geometry, original 574 matter Hamiltonian; analytic exact Gauss compact sources, equal initial reader statistics and equal preparation energy, different bounded short-time record probabilities at any fixed hbar; numerical radial diagnostic only; no added pointer interaction, no autonomous preparation/terminal readout or quantum geometry')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
