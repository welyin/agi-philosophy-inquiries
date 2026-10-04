"""685 entry: actual unsubtracted K determinant in the677 joint regulator limit."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_local_source_lift as body
base=body.base;prior=body.prior;err=body.err
TARGET=HERE/'body_weight_limit_probe_results.json'
CHARGES=np.array([1]*6+[-4]*3+[2]*3+[-3]*2+[6,0])

def background(theta):
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    r=base.prior.rep(np.eye(3),np.eye(2),np.exp(1j*theta))
    links[1]=r
    _,_,_,h,_=base.kernel(links)
    g5=np.kron(np.eye(4),np.kron(base.internal.spin.G5,np.eye(16)))
    return links,h,g5

def determinant(h,g5,a,layers):
    n=len(h);X=g5@h;B=2*np.eye(n)+a*X
    phase,lb=np.linalg.slogdet(B)
    assert abs(phase-1)<2e-12
    ha=g5@(a*X)@np.linalg.inv(B)
    assert err(ha-ha.conj().T)<2e-13
    ev=np.linalg.eigvalsh((ha+ha.conj().T)/2)
    assert np.max(abs(ev))<1
    return float((layers+1)*lb+np.sum(np.logaddexp(layers*np.log1p(ev),layers*np.log1p(-ev))))

def run():
    theta=np.pi/10;free,h0,g5=background(0.);links,h1,_=background(theta)
    r=base.prior.rep(np.eye(3),np.eye(2),np.exp(1j*1e-5))
    charges=np.sort(np.rint(np.angle(np.linalg.eigvals(r))/1e-5).astype(int))
    assert np.array_equal(charges,np.sort(CHARGES))
    s=np.sin(CHARGES*theta)
    analytic_delta=float(4*np.sum(np.sqrt(5+4*s)+np.sqrt(5-4*s)-2*np.sqrt(5)))
    observed_delta=float(np.sum(abs(np.linalg.eigvalsh(h1)))-np.sum(abs(np.linalg.eigvalsh(h0))))
    assert abs(analytic_delta-observed_delta)<3e-12 and analytic_delta<-.1
    assert abs(np.trace(g5@h1)-np.trace(g5@h0))<2e-12
    # All space-time plaquettes are identity, but closed temporal holonomy differs.
    for i in range(4):
        loop=links[0,i]@links[1,i]@links[0,i].conj().T@links[1,i].conj().T
        assert err(loop-np.eye(16))<1e-13
    holonomy=links[1,0]@links[1,2]
    assert err(holonomy-np.eye(16))>.1
    direct=[]
    for label,h in [('free',h0),('flat_holonomy',h1)]:
        for layers in (1,3):
            k,_,_=body.blocks(g5@h,g5,.2,layers)
            sk,lk=np.linalg.slogdet(k);formula=determinant(h,g5,.2,layers)
            assert abs(sk-1)<2e-12 and abs(lk-formula)<1e-9
            direct.append(dict(background=label,L=layers,direct_logdet=float(lk),
                formula_logdet=formula,absolute_error=float(abs(lk-formula))))
    rows=[];target=analytic_delta/2
    for a in (.1,.05,.025,.0125):
        layers=int(np.ceil(4/a**2));scale=a*layers
        difference=determinant(h1,g5,a,layers)-determinant(h0,g5,a,layers)
        rows.append(dict(a=a,L=layers,a_times_L=scale,
            log_unsubtracted_weight_ratio=difference,
            log_ratio_divided_by_aL=difference/scale,
            analytic_limit=target,error=float(abs(difference/scale-target))))
    assert rows[-1]['error']<rows[0]['error'] and rows[-1]['error']<.5
    return dict(entry_round=685,latest_formal_round=684,not_formal_round=True,
        theta=theta,original_hypercharges=CHARGES.tolist(),all_plaquettes_identity=True,
        nontrivial_temporal_holonomy=True,analytic_trace_absolute_H_difference=analytic_delta,
        numerical_trace_absolute_H_difference=observed_delta,direct_formula_checks=direct,rows=rows,
        joint_limit_not_physical_time_continuum=True,
        unsubtracted_body_not_claimed_positive=True,
        conditional_weight_obstruction_not_full_Gauss_RP_verdict=True,
        other_counterterms_or_other_models_not_excluded=True)

if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=685,all_checks_passed=True,
        analytic_limit=result['analytic_trace_absolute_H_difference']/2,
        normalized_log_ratios=[r['log_ratio_divided_by_aL'] for r in result['rows']])))
