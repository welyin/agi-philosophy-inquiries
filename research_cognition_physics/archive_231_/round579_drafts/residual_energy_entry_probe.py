"""579 entry: residual full-potential lower bound after the target spectral shift.
Candidate evidence only; no completed numbered round or continuum solution.
"""
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'residual_energy_entry_results.json'
M=original.M
L,U0,_=original.lattice.scalar.parameters()
LAMBDA0=float(np.linalg.det(L)/np.trace(L))
USUM=float(np.sum(U0))
BETA=4/np.sqrt(6.)
RHO0=float(np.sqrt(6)*np.arctanh(np.sqrt((6*M+USUM)/(12*M))))
A=float(LAMBDA0*(6*M-USUM)**2/(512*M*M))
HBAR=.07


def threshold_radius():
    return max(RHO0,2/BETA)+1


def matching_radius(w):
    target=np.log(HBAR*HBAR/(8*A))-2*np.log(w)
    lo=0.;hi=max(1.,target/BETA+2.)
    def lhs(x):return 2*np.log(x)+BETA*x
    while lhs(hi)<target:hi*=2
    for _ in range(100):
        mid=(lo+hi)/2
        if lhs(mid)<target:lo=mid
        else:hi=mid
    return (lo+hi)/2


def lower_bound(w):
    r=matching_radius(w)
    assert r>=RHO0
    return HBAR*HBAR/(8*w*r*r)


def run():
    z,quad=np.polynomial.legendre.leggauss(240)
    rho=2.2+1.8*z;weights=1.8*quad
    amp=np.exp(-1/(1-z*z));amp_d=amp*(-2*z)/(1-z*z)**2/1.8
    u=amp*np.exp(.3j*rho);du=(amp_d+.3j*amp)*np.exp(.3j*rho)
    left=float(np.dot(weights,abs(du)**2-abs(u)**2/(4*rho*rho)))
    right=float(np.dot(weights,abs(du-u/(2*rho))**2))
    assert left>0 and abs(left-right)<1e-13
    assert LAMBDA0>0 and np.min(np.linalg.eigvalsh(L))>=LAMBDA0
    ratios=[]
    for r in np.linspace(RHO0,10.,35):
        rad=np.sqrt(6*M)*np.tanh(r/np.sqrt(6.))
        for fraction in np.linspace(0.,1.,19):
            phi=np.array([rad*np.sqrt(fraction),0.,0.,0.,rad*np.sqrt(1-fraction)])
            actual=float(original.node_potential(phi))
            bound=A*np.exp(BETA*r)
            ratios.append(actual/bound)
    assert min(ratios)>1
    rcrit=threshold_radius()
    wcrit=float(HBAR/(np.sqrt(8*A)*rcrit*np.exp(BETA*rcrit/2)))
    node_rows=[]
    for w in wcrit*np.array([.5,.1,.01,.001]):
        r=matching_radius(w);b=lower_bound(w)
        other=w*A*np.exp(BETA*r)
        log_derivative=-1+4/(BETA*r+2)
        assert abs(other/b-1)<2e-13 and log_derivative<0
        node_rows.append(dict(w=float(w),radius=r,bound=b,
            matching_relative_error=float(abs(other/b-1)),logarithmic_derivative=log_derivative))
    assert all(node_rows[i+1]['bound']>node_rows[i]['bound'] for i in range(3))
    # Full graph bound using the half of nodes with small cell volumes.
    total_volume=2.3
    start=max(64,int(np.ceil(2*total_volume/wcrit))*2)
    graph=[]
    for N in (start,2*start,4*start,8*start):
        wbar=2*total_volume/N
        robust=N/2*lower_bound(wbar)
        # A genuinely unequal volume allocation, below the monotonic threshold.
        raw=1+.35*np.sin(np.arange(N)*.37)
        w=raw/raw.sum()*total_volume
        assert np.max(w)<wcrit
        actual=sum(lower_bound(float(cell)) for cell in w)
        assert actual>=robust
        graph.append(dict(nodes=N,volume=float(w.sum()),robust_total_bound=robust,
                          summed_node_bound=actual,scaled=robust*np.log(N)**2/N**2))
    return dict(status='candidate entry only',checks_passed=4,
        constants=dict(M=M,lambda_lower=LAMBDA0,u_sum=USUM,beta=BETA,
                       rho0=RHO0,potential_prefactor=A,hbar=HBAR,small_volume_threshold=wcrit),
        hardy_identity=dict(left=left,right=right,error=abs(left-right)),
        potential_bound_min_sample_ratio=min(ratios),node_rows=node_rows,graph_rows=graph,
        scope='analytic residual lower-bound candidate for original full positive potential after xi=1/5 shift; no actual ground-state asymptotic, renormalized continuum or completed round')


if __name__=='__main__':
    data=run();payload=json.dumps(data,ensure_ascii=False,indent=2)+'\n'
    if '--write' in sys.argv:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==data
    print(json.dumps(data,ensure_ascii=False))
