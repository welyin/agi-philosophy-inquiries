"""682 executed entry: conditional mean versus full quantum bulk.
Two physical transverse oscillator modes, not original interacting Q0.
Retaining the bulk noise works only with the corresponding boundary Schur action.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'conditional_bulk_probe_results.json'
def covariance(k,times):
    vals,vec=np.linalg.eigh(k)
    return [(vec*(np.exp(-np.sqrt(vals)*t)/(2*np.sqrt(vals))))@vec.T for t in times]

def run():
    omega=.7;kappa=.8;d=omega**2+kappa;a=np.sqrt(d)
    big=np.sqrt(omega**2+2*kappa)
    k=np.array([[d,-kappa],[-kappa,d]])
    # Source identity at every frequency, with an actual nonzero source contact.
    rows=[]
    for p in (0.,.3,1.,2.):
        precision=p*p*np.eye(2)+k
        inv=np.linalg.inv(precision)
        D=p*p+d;schur=D-kappa*kappa/D;f=kappa/D
        j=np.array([.23,-.47])
        full=.5*j@inv@j
        mapped=.5*(j[0]+f*j[1])**2/schur+.5*j[1]**2/D
        assert abs(full-mapped)<2e-15
        assert abs(np.linalg.det(precision)-D*schur)<3e-14
        mean=f*f*inv[0,0];noise=1/D
        assert abs(inv[1,1]-mean-noise)<2e-15
        rows.append(dict(frequency=p,boundary_covariance=float(inv[0,0]),
            conditional_mean_covariance=float(mean),conditional_noise=float(noise),
            complete_bulk_covariance=float(inv[1,1]),full_source=float(full),
            mapped_source=float(mapped),omitted_contact=.5*j[1]**2/D))
    tau=.8;times=tau*np.arange(1,4)
    r1=np.exp(-omega*tau);r2=np.exp(-big*tau);ra=np.exp(-a*tau)
    c=np.array([r1*r2,-r1-r2,1.])
    tt=(times[:,None]+times[None,:]).ravel()
    qfull=np.array([v[1,1] for v in covariance(k,tt)]).reshape(3,3)
    noise=np.exp(-a*(times[:,None]+times[None,:]))/(2*a)
    mean=qfull-noise
    exact=-(ra*(ra-r1)*(ra-r2))**2/(2*a)
    qmean=float(c@mean@c);qnoise=float(c@noise@c);qtotal=float(c@qfull@c)
    assert exact<0 and abs(qmean-exact)<3e-16 and abs(qtotal)<3e-16
    assert abs(qmean+qnoise)<3e-16
    # Fixing boundary covariance to the single original oscillator requires a
    # time-nonlocal counterterm. Euclidean positive precision alone is not RP.
    pinned=[]
    for p in (0.,.3,1.,2.):
        D=p*p+d
        precision=np.array([[p*p+omega**2+kappa*kappa/D,-kappa],[-kappa,D]])
        inv=np.linalg.inv(precision)
        assert np.linalg.eigvalsh(precision)[0]>0
        assert abs(inv[0,0]-1/(p*p+omega**2))<2e-15
        scalar=1/(p*p+omega**2)-kappa/(p*p+d)**2
        assert abs(inv[1,1]-scalar)<2e-15
        pinned.append(dict(frequency=p,positive_precision_minimum=float(np.linalg.eigvalsh(precision)[0]),
            unchanged_boundary=float(inv[0,0]),new_bulk=float(inv[1,1])))
    tau2=2*np.log(2)/(a-omega);times2=tau2*np.arange(1,3)
    c2=np.array([np.exp(-omega*tau2),-1.])
    t=times2[:,None]+times2[None,:]
    qpin=np.exp(-omega*t)/(2*omega)-kappa*np.exp(-a*t)*(1+a*t)/(4*a**3)
    S=float(c2@np.exp(-a*times2))
    T=float(c2@(times2*np.exp(-a*times2)))
    exact2=-kappa*(S*S+2*a*S*T)/(4*a**3)
    actual2=float(c2@qpin@c2)
    assert S>0 and T>0 and exact2<0 and abs(actual2-exact2)<2e-18
    return dict(date='2026-10-02',entry_round=682,not_formal_round=True,
        omega=omega,kappa=kappa,quantum_precision_spatial_part=k.tolist(),
        original_complete_model_not_replaced=True,frequency_source_checks=rows,
        conditional_mean_witness=dict(times=times.tolist(),coefficients=c.tolist(),
            exact_negative_norm=float(exact),computed_negative_norm=qmean,
            noise_contribution=qnoise,full_quantum_norm=qtotal,
            full_RP_eigenvalues=np.linalg.eigvalsh(qfull).tolist()),
        boundary_preserving_counterterm=dict(frequency_checks=pinned,times=times2.tolist(),
            coefficients=c2.tolist(),analytic_bulk_negative_norm=exact2,
            computed_bulk_negative_norm=actual2,
            full_dynamical_RP_not_restored_by_noise_if_boundary_is_pinned=True),
        dependency_hashes={p:hashlib.sha256((HERE.parent/p).read_bytes()).hexdigest()
            for p in ('research_note_643.md','research_note_678.md','research_note_681.md')},
        scope='Declared two-mode free quantum bulk: conditional noise and full source contact are necessary in this mapping; boundary action changes. A separately defined boundary-preserving nonlocal counterterm is Euclidean positive but not reflection positive. Not original Gauss/S9, not chiral or gravity theorem.')

if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=682,not_formal_round=True,
        mean_RP=result['conditional_mean_witness']['exact_negative_norm'],
        pinned_bulk_RP=result['boundary_preserving_counterterm']['analytic_bulk_negative_norm'],
        all_checks_passed=True)))

