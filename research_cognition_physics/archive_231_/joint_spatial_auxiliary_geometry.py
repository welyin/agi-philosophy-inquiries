"""656: the original free spatial auxiliary Pfaffian and its joint response.

No S9 integration, interacting gauge positivity or physical CAR reconstruction.
All lattices and internal matrices are inherited, not inferred dimensions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_subgroup_measure_source as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_spatial_auxiliary_geometry_results.json'


def frame(kappa=1.,nx=3,nt=2):
    sites=[(t,x) for t in range(nt) for x in range(nx)];n=len(sites)
    st=np.zeros((n,n));sx=st.copy()
    for i,(t,x) in enumerate(sites):
        st[i,sites.index(((t+1)%nt,x))]=-1 if t==nt-1 else 1
        sx[i,sites.index((t,(x+1)%nx))]=1
    def wilson(s,g):
        return np.kron(np.eye(n)-(s+s.T)/2,np.eye(4))+np.kron((s-s.T)/2,g)
    g5=np.kron(np.eye(n),old.spin.G5)
    dx=wilson(sx,old.spin.GAMMA[0])
    h=g5@(-np.eye(4*n)+wilson(st,old.spin.GAMMA[3])+kappa*dx)
    ev,vec=np.linalg.eigh(h);v=vec[:,ev<0];p=v@v.conj().T
    assert np.max(abs(h-h.conj().T))<1e-13 and min(abs(ev))>.4
    b=np.kron(np.eye(n),old.B)
    assert np.max(abs(b.conj().T@p.conj()@b-p))<2e-14
    return dict(sites=sites,H=h,dH=g5@dx,ev=ev,vec=vec,V=v,P=p,B=b)


def local_pairings(f):
    v=f['V']
    return [v[4*i:4*i+4].T@old.B@v[4*i:4*i+4] for i in range(len(f['sites']))]


def auxiliary(f,e):
    b=local_pairings(f)
    return sum(np.kron(bi,sum(ea*ta for ea,ta in zip(ei,old.T))) for bi,ei in zip(b,e))


def ratio(f,e):
    base=np.zeros_like(e);base[:,0]=1
    return old.pfaffian(auxiliary(f,e))/old.pfaffian(auxiliary(f,base))


def plane(f,theta):
    e=np.zeros((len(theta),10));e[:,0]=np.cos(theta);e[:,1]=np.sin(theta)
    v=f['V'];d=np.repeat(np.exp(1j*theta),4)
    m=v.conj().T@(d[:,None]*v)
    return e,m


def weights(f):
    p=f['P'];n=len(f['sites']);w=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i!=j:w[i,j]=8*np.linalg.norm(p[4*i:4*i+4,4*j:4*j+4])**2
    return w


def hessian_trace(f):
    b=local_pairings(f);n=len(b)
    a0=np.kron(sum(b),old.T[0]);inv=np.linalg.inv(a0)
    # First tangent a=1. Independent full antisymmetric-matrix trace evaluation.
    g=[inv@np.kron(bi,old.T[1]) for bi in b]
    h=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            h[i,j]=-.5*np.trace(g[i]@g[j]).real
            if i==j:h[i,j]-=.5*np.trace(inv@np.kron(b[i],old.T[0])).real
    return h


def planar_check():
    rng=np.random.default_rng(65601);rows=[];maxerr=0.
    for kappa,nx,nt in ((0.,3,2),(.5,3,2),(1.,3,2),(1.,2,3)):
        f=frame(kappa,nx,nt)
        for width in (.2,1.,2.7):
            theta=rng.uniform(-width,width,len(f['sites']));e,m=plane(f,theta)
            actual=ratio(f,e);predicted=abs(np.linalg.det(m))**8
            error=abs(actual-predicted);maxerr=max(maxerr,error)
            assert error<3e-13 and -1e-13<=actual.real<=1+1e-13
            rows.append(dict(kappa=kappa,nx=nx,nt=nt,angle_width=width,
                normalized_Pfaffian=[float(actual.real),float(actual.imag)],
                compressed_determinant_power=float(predicted),error=float(error)))
    f=frame(0);theta=np.zeros(6);theta[3]=np.pi
    e,m=plane(f,theta);zero=abs(ratio(f,e));assert zero<1e-24
    theta=np.full(6,.71);e,m=plane(frame(1),theta)
    constant=ratio(frame(1),e);assert abs(constant-1)<2e-13
    return dict(rows=rows,maximum_identity_error=float(maxerr),zero_example=float(zero),
                uniform_direction_error=float(abs(constant-1)),sign_checked_by_Pfaffian_not_square_root=True)


def tangent_check():
    ks=[old.T[0].conj().T@t for t in old.T[1:]];err=0.
    for a,ka in enumerate(ks):
        err=max(err,float(np.max(abs(ka+ka.conj().T))),float(abs(np.trace(ka))))
        for b,kb in enumerate(ks):
            err=max(err,float(np.max(abs(ka@kb+kb@ka+(2 if a==b else 0)*np.eye(16)))))
    assert err<1e-13
    rows=[]
    for kappa in (0.,.5,1.):
        f=frame(kappa);w=weights(f);lap=np.diag(w.sum(axis=1))-w
        h=hessian_trace(f);error=float(np.max(abs(h+lap)))
        assert error<3e-13 and min(np.linalg.eigvalsh(lap))>-1e-13
        rows.append(dict(kappa=kappa,first_Hessian_row=h[0].tolist(),trace_identity_error=error,
            Laplacian_spectrum=np.linalg.eigvalsh(lap).tolist()))
    # Test all nine components with the original full Pfaffian in one arbitrary tangent.
    rng=np.random.default_rng(65602);y=rng.normal(size=(6,9));y*=.2
    f=frame(1);w=weights(f);lap=np.diag(w.sum(axis=1))-w
    predicted=-float(np.einsum('ia,ij,ja->',y,lap,y));step=.002
    def logged(t):
        e=np.column_stack((np.sqrt(1-t*t*np.sum(y*y,axis=1)),t*y))
        return np.log(abs(ratio(f,e)))
    fd=(logged(step)+logged(-step)-2*logged(0))/step**2
    assert abs(fd-predicted)<2e-5
    return dict(internal_Clifford_error=err,rows=rows,nine_tangent_predicted_second=predicted,
                nine_tangent_Pfaffian_second=float(fd),finite_difference_error=float(abs(fd-predicted)))


def geometric_check():
    def analytic(k):
        r=np.sqrt(1+3*k*k)
        return np.array([4*k*k/r**2,8*(1+2/r)**2/9,8*(1-1/r)**2/9])
    rows=[]
    for k in (0.,.25,.5,1.,1.2):
        f=frame(k);w=weights(f);expected=analytic(k)
        actual=w[0,[1,3,4]];err=float(max(abs(expected-actual)))
        assert err<3e-13
        rows.append(dict(kappa=k,space_time_diagonal_weights=actual.tolist(),formula_error=err))
    f=frame(1);ev=f['ev'];vec=f['vec'];indicator=(ev<0).astype(float)
    factor=np.zeros((len(ev),len(ev)))
    for i in range(len(ev)):
        for j in range(len(ev)):
            if indicator[i]!=indicator[j]:factor[i,j]=(indicator[i]-indicator[j])/(ev[i]-ev[j])
    dp=vec@(factor*(vec.conj().T@f['dH']@vec))@vec.conj().T
    p=f['P'];deriv=[]
    for j in (1,3,4):
        deriv.append(float(16*np.vdot(p[:4,4*j:4*j+4],dp[:4,4*j:4*j+4]).real))
    expected=np.array([.5,-8/3,1/3]);assert max(abs(expected-deriv))<3e-13
    delta=1e-4;fd=(weights(frame(1+delta))[0,[1,3,4]]-weights(frame(1-delta))[0,[1,3,4]])/(2*delta)
    assert max(abs(fd-expected))<4e-8
    w=weights(f);lap=np.diag(w.sum(axis=1))-w
    exact_spectrum=np.array([0,11/3,11/3,8,31/3,31/3])
    assert max(abs(np.linalg.eigvalsh(lap)-exact_spectrum))<5e-13
    # No pure equal-time spatial multiplier can fix the copied chain on this slice.
    temporal=np.array([0.,0.,0.,1.,1.,1.]);angle=.01
    costs=[]
    for k in (0.,1.):
        ff=frame(k);e,_=plane(ff,angle*temporal)
        cost=-np.log(abs(ratio(ff,e)));ww=weights(ff);ll=np.diag(ww.sum(axis=1))-ww
        coefficient=.5*float(temporal@ll@temporal)
        assert abs(cost/angle**2-coefficient)<2e-4
        costs.append(dict(kappa=k,quadratic_coefficient=coefficient,actual_cost_over_angle_squared=float(cost/angle**2)))
    assert abs(costs[0]['quadratic_coefficient']-12)<1e-12
    assert abs(costs[1]['quadratic_coefficient']-6)<1e-12
    return dict(weight_rows=rows,original_weight_derivative=deriv,finite_difference_derivative=fd.tolist(),
        derivative_error=float(max(abs(fd-expected))),original_Laplacian_spectrum=exact_spectrum.tolist(),
        spatially_uniform_temporal_counterexample=costs,
        excludes_only_equal_time_multiplier_normalized_to_one_on_spatially_uniform_fields=True)


def leakage_locality_check():
    f=frame(1);theta=np.array([.12,-.17,.08,.21,-.09,.03]);e,m=plane(f,theta)
    v=f['V'];d=np.repeat(np.exp(1j*theta),4)
    r=(np.eye(len(v))-f['P'])@(d[:,None]*v);leak=r.conj().T@r
    residual=float(np.max(abs(m.conj().T@m+leak-np.eye(m.shape[0]))))
    assert residual<2e-14
    eig=np.linalg.eigvalsh(leak);norm=float(max(eig));cost=-np.log(abs(ratio(f,e)))
    lower=float(4*np.trace(leak).real);upper=lower/(1-norm)
    assert lower-1e-12<=cost<=upper+1e-12
    # Full S9, not restricted to a common circle: exact magnitude and leakage.
    rng=np.random.default_rng(65604);general_rows=[]
    vv=np.kron(v,np.eye(16));pp=vv@vv.conj().T;w=weights(f)
    for width in (.07,.25,.8):
        ee=np.eye(10)[0]+width*rng.normal(size=(6,10));ee/=np.linalg.norm(ee,axis=1)[:,None]
        uu=np.zeros((384,384),complex)
        for i,ei in enumerate(ee):
            ki=old.T[0].conj().T@sum(a*t for a,t in zip(ei,old.T))
            uu[64*i:64*(i+1),64*i:64*(i+1)]=np.kron(np.eye(4),ki)
        compressed=vv.conj().T@uu@vv;rr=(np.eye(384)-pp)@uu@vv
        ll=rr.conj().T@rr;leig=np.linalg.eigvalsh(ll);lmax=float(max(leig))
        actual=ratio(f,ee);logcost=-np.log(abs(actual))
        predicted=-.25*np.linalg.slogdet(compressed.conj().T@compressed)[1]
        assert abs(predicted-logcost)<2e-11
        pair=.5*sum(w[i,j]*np.linalg.norm(ee[i]-ee[j])**2 for i in range(6) for j in range(i+1,6))
        trace_lower=.25*np.trace(ll).real
        assert abs(pair-trace_lower)<1e-12 and pair-1e-12<=logcost<=pair/(1-lmax)+1e-11
        general_rows.append(dict(width=width,actual_negative_log_absolute=float(logcost),
            compressed_determinant_cost=float(predicted),pair_distance_lower_bound=float(pair),
            leakage_upper_bound=float(pair/(1-lmax)),magnitude_identity_error=float(abs(predicted-logcost)),
            phase_observed_only=float(np.angle(actual))))
    # Full four-dimensional free symbol: algebraic gap and norm certificate.
    rng=np.random.default_rng(65603);maximum=0.;minimum=100.;error=0.
    for k in rng.uniform(-np.pi,np.pi,(32,4)):
        a=1-np.cos(k);x=(sum(a)-1)*np.eye(4,dtype=complex)
        for mu in range(4):x+=1j*np.sin(k[mu])*old.spin.GAMMA[mu]
        h=old.spin.G5@x;r2=1+2*sum(a[i]*a[j] for i in range(4) for j in range(i+1,4))
        error=max(error,float(np.max(abs(h@h-r2*np.eye(4)))))
        minimum=min(minimum,np.sqrt(r2));maximum=max(maximum,np.sqrt(r2))
    assert error<3e-14 and minimum>=1-1e-13 and maximum<=7+1e-13
    # Consequence of the analytic binomial-series tail, not a fit to these samples.
    return dict(leakage_identity_error=residual,exact_planar_cost=float(cost),lower_bound=lower,
        upper_bound=upper,maximum_leakage=norm,full_S9_magnitude_rows=general_rows,
        full_4D_symbol_identity_error=error,
        sampled_gap_minimum=float(minimum),sampled_norm_maximum=float(maximum),
        analytic_gap_lower=1,analytic_norm_upper=7,locality_q=24/25,
        weight_tail_prefactor=9800,weight_tail_power='2*floor(graph_distance/2), distance>=2',
        no_hard_locality_or_emergent_dimension_claim=True)


def run():
    deps=('joint_subgroup_measure_source.py','joint_chiral_fibre_source.py','joint_spinor_subgroup_mass.py',
          'research_note_615.md','research_note_653.md','research_note_655.md')
    return dict(date='2026-10-02',round=656,tests_run=4,failures=0,errors=0,
        planar=planar_check(),tangent=tangent_check(),geometry=geometric_check(),
        leakage_and_locality=leakage_locality_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original free spatial auxiliary integrand: exact full-S9 magnitude/leakage identity, common-plane signed Pfaffian and all-nine tangent Hessian, shared spatial/time coupling and source derivative. No full S9 integral or general positivity, positive physical transfer, original CAR-state identification, interacting gauge field or GR reconstruction.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as stream:json.dump(result,stream,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
