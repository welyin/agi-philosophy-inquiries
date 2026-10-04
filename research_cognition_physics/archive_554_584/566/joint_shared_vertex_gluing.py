"""566: exact shared-matter gluing and collective records on a two-edge star.

Graph, SU(2), source preparations and apparatus couplings are explicit inputs.
Angular matrices below are exact invariant polynomial sectors, not time truncations.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
import joint_gauge_link_reference as link
import joint_radial_prediction_closure as old
import joint_finite_time_gauge_probe as inherited

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_shared_vertex_gluing_results.json'
A=.5; B=1.; K=.8; ETA=.5; NU=.5


def angular():
    ident=np.eye(5)
    ls=[];rs=[]
    for i,j in enumerate(link.J,1):
        right=np.column_stack([link.left(np.eye(4)[m])@np.eye(4)[i]/2 for m in range(4)])
        l=np.zeros((5,5));l[1:,1:]=j.T
        r=np.zeros((5,5));r[1:,1:]=right.T
        ls.append(l);rs.append(r)
    l1=[np.kron(l,ident) for l in ls];l2=[np.kron(ident,l) for l in ls]
    r2=[np.kron(ident,r) for r in rs]
    kc=-4*sum((x+y)@(x+y) for x,y in zip(l1,l2))
    k1=-4*sum(x@x for x in l1);k2=-4*sum(x@x for x in l2)
    reverse=np.kron(ident,np.diag([1,1,-1,-1,-1]))
    kr=-4*sum((-x+y)@(-x+y) for x,y in zip(l1,r2))
    assert np.max(abs(kr-reverse@kc@reverse))<1e-14
    f=[]
    for m in [np.eye(4)]+[2*j for j in link.J]:
        c=np.zeros((5,5));c[1:,1:]=m/2;f.append(c.ravel())
    vacuum=np.zeros(25);vacuum[0]=1.
    pairs=[(vacuum+sign*ETA*f[0])/math.sqrt(1+ETA*ETA) for sign in (1,-1)]
    for i,psi in enumerate(f):
        assert abs(psi@psi-1)<1e-14
        assert np.linalg.norm(kc@psi-(0 if i==0 else 8)*psi)<1e-14
        assert np.linalg.norm(k1@psi-3*psi)+np.linalg.norm(k2@psi-3*psi)<1e-14
        c=psi.reshape(5,5)
        target=np.diag([0]+[.25]*4)
        assert np.max(abs(c@c.T-target))<1e-14
        assert np.max(abs(c.T@c-target))<1e-14
    marginal=[]
    for psi in pairs:
        c=psi.reshape(5,5);marginal.append(c@c.T)
        assert np.linalg.norm(kc@psi)<1e-14
    assert np.max(abs(marginal[0]-marginal[1]))<1e-14
    return dict(kc=kc,k1=k1,k2=k2,functions=f,pairs=pairs,
                left1=l1,left2=l2,reverse_residual=float(np.max(abs(kr-reverse@kc@reverse))))


def cubature():
    points=[sign*np.eye(4)[i] for i in range(4) for sign in (-1,1)]
    points += [np.array(v)/2 for v in itertools.product((-1,1),repeat=4)]
    points=np.array(points);weights=np.full(24,1/24)
    assert np.max(abs((points.T*weights)@points-np.eye(4)/4))<1e-15
    assert abs(weights@(points[:,0]**4)-1/8)<1e-15
    assert abs(weights@(points[:,0]**2*points[:,1]**2)-1/24)<1e-15
    return points,weights


def angular_moments():
    q,w=cubature();pairw=w[:,None]*w[None,:];z=q[:,0]
    S=q@q.T;f0=2*S
    out=[]
    for sign in (1,-1):
        density=(1+sign*ETA*f0)**2/(1+ETA*ETA)
        measure=pairw*density
        result=dict(norm=float(measure.sum()),
            S=float(np.sum(measure*S)),
            z1=float(np.sum(measure*z[:,None])),
            z2=float(np.sum(measure*z[None,:])),
            z1z2=float(np.sum(measure*z[:,None]*z[None,:])),
            z1_squared=float(np.sum(measure*z[:,None]**2)))
        assert abs(result['norm']-1)<1e-14
        assert abs(result['S']-sign*ETA/(1+ETA*ETA))<1e-14
        assert abs(result['z1z2']-sign*ETA/(4*(1+ETA*ETA)))<1e-14
        assert abs(result['z1_squared']-.25)<1e-14
        out.append(result)
    return out


def full_H_check(ops):
    rng=np.random.default_rng(566)
    coef=np.zeros((5,5),dtype=complex);coef[0,0]=1
    coef[1,0]=.12;coef[0,2]=.09j
    coef[3,4]=.11;coef[1,1]=.07j
    coef=coef.astype(complex) if np.iscomplexobj(coef) else coef
    # Complex non-class functions are used; no angular averaging in this check.
    def evaluate(c,w):
        v=[np.r_[1.,2*t] for t in w]
        return v[0]@c.reshape(5,5)@v[1]
    def relative(q,U):
        n=q[:,:4]/np.linalg.norm(q[:,:4],axis=1)[:,None]
        nc=n[0].copy();nc[1:]*=-1
        return [link.left(nc)@U[e]@n[e+1] for e in range(2)]
    def bump(r):
        t=(r-1)/.4
        return math.exp(-1/(1-t*t)) if abs(t)<1 else 0.
    def lap(r):
        t=(r-1)/.4
        lp=-2*t/(.4*(1-t*t)**2)
        lpp=-2/(.4**2*(1-t*t)**2)-8*t*t/(.4**2*(1-t*t)**3)
        return lpp+lp*lp+3*lp/r
    def prefactor(q):
        return np.prod([bump(np.linalg.norm(v[:4])) for v in q])*math.exp(-np.sum(q[:,4]**2)/2)
    def wave(q,U):
        return prefactor(q)*evaluate(coef,relative(q,U))
    L,u,_=inherited.constants()
    residuals=[[],[],[]];gauge=0.;bad=[]
    for _ in range(8):
        q=rng.normal(size=(3,5))
        for row in q:row[:4]*=rng.uniform(.92,1.08)/np.linalg.norm(row[:4])
        uq=rng.normal(size=(2,4));uq/=np.linalg.norm(uq,axis=1)[:,None]
        U=[link.left(v) for v in uq]
        gq=rng.normal(size=(3,4));gq/=np.linalg.norm(gq,axis=1)[:,None]
        G=[link.left(v) for v in gq];moved=q.copy()
        for i in range(3):moved[i,:4]=G[i]@q[i,:4]
        transformed=[G[0]@U[e]@G[e+1].T for e in range(2)]
        ww=relative(q,U)
        gauge=max(gauge,float(np.max(abs(np.array(ww)-np.array(relative(moved,transformed))))))
        radii=np.linalg.norm(q[:,:4],axis=1)
        W=0.
        for r,sing in zip(radii,q[:,4]):
            d=np.array([r*r,sing*sing])-u;W+=d@L@d/4
        F=sum(K/2*np.sum((q[0,:4]-U[e]@q[e+1,:4])**2) for e in range(2))
        scalar=-A*sum(lap(r) for r in radii)+A*np.sum(1-q[:,4]**2)+W+F
        c=coef.ravel()
        ang=A*(ops['kc']@c/radii[0]**2+ops['k1']@c/radii[1]**2+ops['k2']@c/radii[2]**2)
        ang+=B/4*(ops['k1']+ops['k2'])@c
        target=prefactor(q)*(scalar*evaluate(c,ww)+evaluate(ang,ww))
        naive=target+prefactor(q)*A/radii[0]**2*evaluate((ops['k1']+ops['k2']-ops['kc'])@c,ww)
        bad.append(float(abs(naive-target)))
        psi=wave(q,U)
        for index,h in enumerate((.008,.004,.002)):
            laplacian=0j;glap=0j
            for i in range(3):
                for j in range(5):
                    delta=np.zeros((3,5));delta[i,j]=h
                    laplacian+=(-wave(q+2*delta,U)+16*wave(q+delta,U)-30*psi+16*wave(q-delta,U)-wave(q-2*delta,U))/(12*h*h)
            for e in range(2):
                for j in link.J:
                    vals=[]
                    for m in (-2,-1,1,2):
                        rot=math.cos(m*h/2)*np.eye(4)+2*math.sin(m*h/2)*j
                        new=U.copy();new[e]=rot@U[e];vals.append(wave(q,new))
                    glap+=(-vals[0]+16*vals[1]-30*psi+16*vals[2]-vals[3])/(12*h*h)
            actual=-A*laplacian-B*glap+(W+F)*psi
            residuals[index].append(float(abs(actual-target)/max(1.,abs(target))))
    errors=list(map(max,residuals))
    assert gauge<2e-14 and errors[-1]<2e-7
    assert errors[1]<errors[0]/8 and errors[2]<errors[1]/8
    assert max(bad)>1e-4
    return dict(gauge_residual=gauge,steps=[.008,.004,.002],
        full_fifteen_matter_six_group_direction_errors=errors,
        missed_cross_term_maximum_pointwise_residual=max(bad))


def radial_stats(n):
    r,w,lp,_=old.radial_grid(n)
    return r,w,dict(mean=float(w@r),second=float(w@(r*r)),
                   fourth=float(w@(r**4)),inverse_second=float(w@(r**-2)),
                   radial_kinetic=A*float(w@(lp*lp)))


def source_and_injection(n=160):
    _,_,m=radial_stats(n);L,u,c= inherited.constants()
    m1,m2,m4,mi=m['mean'],m['second'],m['fourth'],m['inverse_second']
    Wcell=(L[0,0]*(m4-2*u[0]*m2+u[0]**2)
       +2*L[0,1]*(m2-u[0])*(.5-u[1])
       +L[1,1]*(.75-u[1]+u[1]**2))/4
    base=3*(m['radial_kinetic']+A/2+Wcell)+2*K*m2
    occupation=ETA**2/(1+ETA**2)
    energy=base+occupation*(6*A*mi+1.5*B)
    angular=angular_moments();rows=[]
    for row in angular:
        # Moment evaluation through independent positive Haar cubature.
        Y2=K*K/4*(6*m4+10*m2*m2)+K*K/2*m2*m2+2*K*K*m2*m1*m1*row['z1z2']+NU**2
        centre_cross=K*K*(m2+m1*m1*row['S'])
        separate=8*A*K*K*m2+3*B*K*K*m2*m2/8
        shared=separate+2*A*centre_cross  # g^2 Var(P)=1 for NU=.5, g=1.
        rows.append(dict(source_mean_H=energy,mean_actual_Y=2*K*m2,
                         second_actual_Y=Y2,shared_pointer_source_injection=shared,
                         independent_pointer_source_injection=separate))
    delta_y2=rows[0]['second_actual_Y']-rows[1]['second_actual_Y']
    delta_energy=rows[0]['shared_pointer_source_injection']-rows[1]['shared_pointer_source_injection']
    assert abs(delta_y2-ETA*K*K*m2*m1*m1/(1+ETA*ETA))<1e-12
    assert abs(delta_energy-4*A*K*K*m1*m1*ETA/(1+ETA*ETA))<1e-12
    D0=36*K*K*u[0]**2;D1=48*K*K/c['ell']
    tau=.2;mass=tau*tau/D1
    return dict(radial_moments=m,cell_mean_W=float(Wcell),source_rows=rows,
        singlet_vs_triplet_source_energy_difference=8*A*mi,
        second_record_difference=delta_y2,collective_injection_difference=delta_energy,
        finite_pulse_existence_budget=dict(D0=float(D0),D1=float(D1),tau=tau,M=float(mass),
            theta=.5,lower_bound=-float(D0/(2*D1)),pointer_initial_energy=float(1/(2*mass)),
            no_uniform_distinguishing_tau_certificate_claim=True))


def record(n=96,nangle=120):
    r,w,m=radial_stats(n);z,wz=old.class_grid(nangle)
    frequency=.7;centre=2*K*m['second']
    x=frequency*K*r[:,None]*r[None,:]
    c=np.cos(x[:,:,None]*z);s=np.sin(x[:,:,None]*z)
    c0=c@wz;c2=c@(wz*z*z);sin1=s@(wz*z)
    total=np.zeros(2);cross_value=0.
    for i,r0 in enumerate(r):
        mean_phase=frequency*(K/2*(2*r0*r0+r[:,None]**2+r[None,:]**2)-centre)
        phase=np.cos(mean_phase);weights=w[i]*w[:,None]*w[None,:]
        common=np.outer(c0[i],c0[i])+4*ETA**2*(np.outer(c2[i],c2[i])
           +np.outer(c0[i]-c2[i],c0[i]-c2[i])/3)
        corr=4*ETA*np.outer(sin1[i],sin1[i])
        for index,sign in enumerate((1,-1)):
            total[index]+=np.sum(weights*phase*(common-sign*corr))/(1+ETA*ETA)
        cross_value+=np.sum(weights*phase*np.outer(sin1[i],sin1[i]))
    out=.5+.5*math.exp(-NU*NU*frequency**2/2)*total
    predicted=-4*ETA/(1+ETA*ETA)*math.exp(-NU*NU*frequency**2/2)*cross_value
    assert np.max(out)<=1 and np.min(out)>=0
    assert abs(out[0]-out[1]-predicted)<1e-14 and abs(predicted)>1e-4
    # Independent rigorous small-frequency Taylor criterion, unshifted cos statistic.
    delta=source_and_injection(n)['second_record_difference']
    bmax=4*K*1.4**2
    t_cert=math.sqrt(3*delta/bmax**4) # half the allowed t^2 ceiling, from ΔB².
    lower=math.exp(-NU*NU*t_cert*t_cert/2)*t_cert*t_cert*delta/8
    assert lower>0
    return dict(radial_nodes=n,angular_nodes=nangle,frequency=frequency,shift=centre,
        bounded_compensated_record_expectations=out.tolist(),difference=float(out[0]-out[1]),
        analytic_small_frequency_witness=dict(B_upper=bmax,frequency=t_cert,
            sufficient_effect_gap_lower=lower,
            formula_is_analytic_numerical_value_is_not_interval_certificate=True))


def gradient_checks():
    rng=np.random.default_rng(8566);err=0.;fd_error=0.
    for _ in range(40):
        x=rng.normal(size=(3,4));uq=rng.normal(size=(2,4));uq/=np.linalg.norm(uq,axis=1)[:,None]
        U=[link.left(q) for q in uq]
        d=[x[0]-U[e]@x[e+1] for e in range(2)]
        cross=K*K*np.dot(d[0],d[1])
        r=np.linalg.norm(x,axis=1);n=x/r[:,None]
        cn=n[0].copy();cn[1:]*=-1
        ww=[link.left(cn)@U[e]@n[e+1] for e in range(2)]
        target=K*K*(r[0]**2-r[0]*r[1]*ww[0][0]-r[0]*r[2]*ww[1][0]+r[1]*r[2]*np.dot(ww[0],ww[1]))
        err=max(err,abs(cross-target))
        def obs(q):
            return sum(K/2*np.sum((q[0]-U[e]@q[e+1])**2) for e in range(2))
        h=1e-5;fd=np.zeros((3,4))
        for i in range(3):
            for j in range(4):
                dx=np.zeros((3,4));dx[i,j]=h
                fd[i,j]=(obs(x+dx)-obs(x-dx))/(2*h)
        grad=np.array([K*(d[0]+d[1]),-K*U[0].T@d[0],-K*U[1].T@d[1]])
        fd_error=max(fd_error,float(np.max(abs(fd-grad))))
        squared=A*np.sum(grad*grad)
        edge_sum=2*A*K*K*sum(np.dot(v,v) for v in d)
        assert abs(squared-edge_sum-2*A*cross)<2e-13
    assert err<5e-14 and fd_error<1e-8
    return dict(shared_gradient_identity_residual=err,Cartesian_gradient_difference_residual=fd_error,
                electric_cross_terms_between_distinct_links_zero=True)


def run():
    op=angular()
    hcheck=full_H_check(op);am=angular_moments();source=source_and_injection()
    source2=source_and_injection(128)
    r1=record(96,100);r2=record(144,160)
    refinement=max(abs(r1['difference']-r2['difference']),
        abs(source['source_rows'][0]['source_mean_H']-source2['source_rows'][0]['source_mean_H']))
    assert refinement<2e-8
    gradients=gradient_checks()
    tests=['full_shared_vertex_H_and_orientation_reversal',
       'exact_angular_sector_central_Casimir_and_identical_link_marginals',
       'positive_Haar_cubature_and_complete_equal_energy_sources',
       'actual_collective_pointer_second_moment_and_bounded_record',
       'shared_Cartesian_gradient_and_collective_backreaction',
       'finite_duration_apparatus_semiboundedness_and_refinement']
    deps=['joint_radial_prediction_closure.py','joint_gauge_link_reference.py','joint_finite_time_gauge_probe.py']
    return dict(round=566,tests_run=len(tests),failures=0,errors=0,checks=tests,
       dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
       exact_angular=dict(link_Casimir=3,central_singlet_Casimir=0,central_triplet_Casimir=8,
           every_single_link_reduced_density_equal=True,reverse_residual=op['reverse_residual']),
       full_H=hcheck,positive_angular_moments=am,source_and_injection=source,
       bounded_record=r2,refinement_difference=refinement,gradients=gradients,
       scope=dict(exact_finite_graph_Gauss_form_reduction=True,
           same_local_marginals_not_same_joint_state=True,not_a_failure_of_local_tomography=True,
           independent_readers_can_use_joint_statistics=True,
           source_injection_identity_is_nonselective_instantaneous_only=True,
           finite_positive_duration_distinguishability_is_per_source_existence=True,
           graph_SU2_source_pointer_control_inputs_remain=True,
           no_dimension_continuum_gravity_or_unified_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','full_H','source_and_injection','bounded_record','refinement_difference')},ensure_ascii=False))
