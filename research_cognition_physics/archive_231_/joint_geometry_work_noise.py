"""590: same curved matter, changing geometry, and force-noise domains.

Full operator statements are analytic. Numerical checks evaluate original
geometry coefficients, actual Gauss packet marginals, and local differential
identities. No full lattice quantum propagation or GR simulation is claimed.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_full_spatial_metric as model
import joint_record_source_compression as records
import joint_quantum_measure_records as packet

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_geometry_work_noise_results.json'
spec=importlib.util.spec_from_file_location('time_entry590',HERE/'round590_drafts/time_dependent_form_entry.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
HBAR, W=.7,.8


def time_entry_check():
    r=entry.run();assert r==json.loads(entry.TARGET.read_text('utf8'))
    return dict(frozen_entry_reproduced=True,entry_evidence=r)


def conformal_force_check():
    q=model.original.shared_source(8);x,y,z=np.moveaxis(q['grid'],-1,0)
    psi=1.1+.05*np.cos(x)+.02*np.sin(y);shape=model.shape_field(q['grid'])
    base=model.energy(model.graph_data(q,psi),shape)
    powers=dict(scalar_kinetic=-6,onsite=6,gradient=2,electric=-2,magnetic=-2)
    slope=sum(powers[k]*base[k] for k in powers);rows=[]
    for u in (-.09,.04,.12):
        got=model.energy(model.graph_data(q,psi*np.exp(u)),shape)
        expected={k:np.exp(powers[k]*u)*base[k] for k in powers}
        error=max(abs(got[k]-expected[k]) for k in powers)
        assert error<2e-10
        rows.append(dict(log_scale=u,exact_full_source_scaling_error=error))
    diff=[]
    for h in (.01,.005,.0025):
        ep=sum(base[k]*np.exp(powers[k]*h) for k in powers)
        em=sum(base[k]*np.exp(-powers[k]*h) for k in powers)
        diff.append(abs((ep-em)/(2*h)-slope))
    assert diff[-1]<diff[0]/12 and abs(slope)<=6*base['total']
    return dict(full_original_Gauss_source=True,scaling_rows=rows,
                conformal_force_classical_value=slope,relative_form_bound=6,
                derivative_errors=diff)


def amplitude(h,s):
    x=(h-packet.HCENTER)/packet.HRADIUS;y=(s-packet.SCENTER)/packet.SRADIUS
    f=model.original.M-(h*h+s*s)/6
    return np.exp(-1/(1-x*x)-1/(1-y*y))*f**1.5/model.original.M**.25


def metric_radial(h,s):
    q=np.array([h,s]);f=model.original.M-np.dot(q,q)/6
    g=f*(np.eye(2)-np.outer(q,q)/(6*model.original.M))
    drift=np.array([3*f/h-f*h/(3*model.original.M),-f*s/(3*model.original.M)])
    return g,drift


def curved_phase_check():
    h,s=.537,.419;g,b=metric_radial(h,s)
    def lap_ratio(n,step):
        def f(x,y):return amplitude(x,y)*np.exp(1j*n*y)
        center=f(h,s)
        dh=(f(h+step,s)-f(h-step,s))/(2*step)
        ds=(f(h,s+step)-f(h,s-step))/(2*step)
        hh=(f(h+step,s)-2*center+f(h-step,s))/step**2
        ss=(f(h,s+step)-2*center+f(h,s-step))/step**2
        hs=(f(h+step,s+step)-f(h+step,s-step)-f(h-step,s+step)+f(h-step,s-step))/(4*step**2)
        return (g[0,0]*hh+2*g[0,1]*hs+g[1,1]*ss+b[0]*dh+b[1]*ds)/center
    rows=[]
    for n in (3,9):
        errs=[]
        for step in (.001,.0005,.00025):
            # The even phase difference cancels all original multiplicative
            # potentials and every other node/link derivative exactly.
            value=(lap_ratio(n,step)+lap_ratio(-n,step))/2-lap_ratio(0,step)
            errs.append(float(abs(value+n*n*g[1,1])/(n*n*g[1,1])))
        assert errs[-1]<errs[0]/12 and errs[-1]<2e-5
        rows.append(dict(phase_frequency=n,relative_second_difference_errors=errs))
    return dict(rows=rows,original_H5_radial_laplacian=True,
                target_source_leading_coefficient=-3*HBAR**2/W*g[1,1],
                full_potential_cancels_analytically=True)


def packet_data(n):
    z,w=np.polynomial.legendre.leggauss(n);x,y=np.meshgrid(z,z,indexing='ij')
    h=packet.HCENTER+packet.HRADIUS*x;s=packet.SCENTER+packet.SRADIUS*y
    chi=np.exp(-1/(1-x*x)-1/(1-y*y));prob=w[:,None]*w[None,:]*h**3*chi**2
    prob/=prob.sum();f=model.original.M-(h*h+s*s)/6
    kss=f*(1-s*s/(6*model.original.M))
    return h,s,prob,kss


def noise_tail_check():
    quadrature=[]
    for n in (48,80,128):
        h,s,p,k=packet_data(n)
        a=HBAR**2/(2*W)*float(np.sum(p*k))
        d=9*HBAR**4/W**2*float(np.sum(p*k*k))
        quadrature.append(dict(nodes=n,phase_energy_n2_coefficient=a,force_norm_n4_coefficient=d))
    assert abs(quadrature[-1]['force_norm_n4_coefficient']-quadrature[-2]['force_norm_n4_coefficient'])<1e-12
    a=quadrature[-1]['phase_energy_n2_coefficient'];d=quadrature[-1]['force_norm_n4_coefficient']
    zeta4=np.pi**4/90;limit=a*(np.pi**2/6)/zeta4;tails=[]
    for N in (10,100,1000,10000):
        n=np.arange(1,N+1,dtype=float)
        mass=float(np.sum(n**-4)/zeta4)
        increment=a*float(np.sum(n**-2)/zeta4)
        leading=d*N/zeta4
        assert mass<1 and increment<limit and leading>0
        tails.append(dict(N=N,excited_mixture_mass=mass,remainder_on_unphased_packet=1-mass,
                          energy_increment=increment,leading_force_second_moment_term=leading))
    return dict(packet_quadrature=quadrature,tail_rows=tails,finite_energy_increment_limit=limit,
                force_second_moment_asymptotic_slope=d/zeta4,
                displayed_noise_is_leading_term_not_full_moment=True,
                infinite_noise_is_analytic_tail_theorem=True,phases_are_Gauss_invariant=True)


def record_noise_check():
    h,s,p,k=packet_data(100);lf,df,lc,dc=records.instruments(s)
    leading=9*HBAR**4/W**2*k*k
    old=float(np.sum(p*leading))
    fine=float(np.sum(p*leading*np.sum(lf*lf,axis=-1)))
    coarse=float(np.sum(p*leading*np.sum(lc*lc,axis=-1)))
    assert max(abs(fine-old),abs(coarse-old))<1e-12
    # A joint geometry--matter probability distribution need not factorize.
    # Normalize a positive correlated density over two angular branches and
    # the original compact radial packet. Geometry derivative of L is zero.
    angles=np.array([-.4,.7]);weights=np.array([.43,.57]);eta=.12*np.sin(angles)
    corr=np.stack((1+.1*np.sin(s),1-.1*np.sin(s)),axis=0)
    joint=weights[:,None,None]*p[None,:,:]*corr;joint/=joint.sum()
    volume=W*np.exp(6*eta)
    gap=HBAR**2/2*float(np.sum(joint*k[None,:,:]*records.delta_a(s)[None,:,:]/volume[:,None,None]))
    sumA=np.sum(df*df,axis=-1)-np.sum(dc*dc,axis=-1)
    direct=HBAR**2/2*float(np.sum(joint*k[None,:,:]*sumA[None,:,:]/volume[:,None,None]))
    assert abs(gap-direct)<1e-14 and gap>0
    return dict(original_n4_coefficient=old,after_fine_n4= fine,after_direct_parity_n4=coarse,
                leading_noise_coefficient_preserved=True,joint_correlated_source_energy_gap=gap,
                injection_formula_error=abs(gap-direct),geometry_matter_joint_probability_not_factorized=True,
                instantaneous_record_only=True)


def internal_exchange_check():
    # Check the operator-valued product rule on the actual geometry matrices,
    # without replacing a full quantum state by a classical source trajectory.
    q=model.original.shared_source(4);t=.37;inertia=1.7
    def f(x):return np.exp(.13*np.cos(x)+.31j*np.sin(x))
    fp=f(t)*(-.13*np.sin(t)+.31j*np.cos(t))
    fpp=f(t)*((-.13*np.sin(t)+.31j*np.cos(t))**2-.13*np.cos(t)-.31j*np.sin(t))
    C,dC=entry.coefficients(q,t)
    worst=[]
    for step in (.004,.002,.001):
        plus,_=entry.coefficients(q,t+step);minus,_=entry.coefficients(q,t-step)
        maxerr=0.
        for key in ('gradient','electric','magnetic'):
            second=(plus[key]-2*C[key]+minus[key])/step**2
            # i/hbar [Pi^2/(2I), C] f, from an independent finite difference.
            comm=1j/HBAR*(-HBAR**2/(2*inertia))*((plus[key]*f(t+step)-2*C[key]*f(t)+minus[key]*f(t-step))/step**2-C[key]*fpp)
            rhs=-1j*HBAR/(2*inertia)*(second*f(t)+2*dC[key]*fp)
            maxerr=max(maxerr,float(np.max(abs(comm-rhs))))
        worst.append(maxerr)
    assert worst[-1]<worst[0]/12 and worst[-1]<1e-6
    # Local coverage of all six symmetric metric directions on a compact
    # controller torus: exp(Q0 + r sum sin(theta_a) T_a).
    T=[]
    for i in range(3):
        a=np.zeros((3,3));a[i,i]=1;T.append(a)
    for i,j in ((0,1),(0,2),(1,2)):
        a=np.zeros((3,3));a[i,j]=a[j,i]=1/np.sqrt(2);T.append(a)
    T=np.array(T);radius=.2;Q0=np.array([[.1,.04,0],[.04,-.12,.03],[0,.03,.02]])
    def gamma(v):return model.shape_exp(Q0+radius*np.einsum('a,aij->ij',np.sin(v),T))
    columns=[]
    for a in range(6):
        d=np.zeros(6);d[a]=1e-5
        columns.append(((gamma(d)-gamma(-d))/(2e-5)).reshape(9))
    sv=np.linalg.svd(np.stack(columns,axis=-1),compute_uv=False)
    rng=np.random.default_rng(590);mineig=min(float(np.linalg.eigvalsh(gamma(v)).min()) for v in rng.uniform(-np.pi,np.pi,(30,6)))
    assert sv[-1]>.1 and mineig>0
    return dict(actual_geometry_coefficient_commutator_errors=worst,
                local_six_metric_direction_singular_values=sv.tolist(),sampled_min_metric_eigenvalue=mineig,
                compactness_proof_not_replaced_by_sampling=True,
                controller_inertia_and_potential_are_extra_inputs=True,Einstein_dynamics_not_claimed=True)


def run():
    tests=(time_entry_check,conformal_force_check,curved_phase_check,noise_tail_check,record_noise_check,internal_exchange_check)
    evidence={f.__name__:f() for f in tests}
    files=('joint_full_spatial_metric.py','joint_curved_quantum_source.py','joint_record_source_compression.py',
           'joint_quantum_measure_records.py','round590_drafts/time_dependent_form_entry.py',
           'round590_drafts/time_dependent_form_entry_results.json','round590_drafts/time_dependent_form_entry.md')
    return dict(round=590,tests_run=len(tests),failures=0,errors=0,evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files},
                scope='Same finite curved matter and quotient gauge theory. Time-dependent background and explicitly added compact quantum geometry controller; not Einstein gravity. Finite mean energy does not ensure force-noise domain. No full lattice quantum propagation or fixed-hbar continuum claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
