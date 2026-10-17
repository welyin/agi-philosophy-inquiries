"""1241: explicit states and fixed-constant scale compatibility, not a detector.

NumPy only. Default recomputes and reads results.json; --write creates it once.
Analytic inequalities in proof.md carry the full-parameter claims.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]

def need(ok,label):
    if not bool(ok): raise AssertionError(label)

def quadrature(a,b,n):
    x,w=np.polynomial.legendre.leggauss(n)
    return (a+b)/2+(b-a)*x/2,(b-a)*w/2

def group_packet():
    x,w=quadrature(0,1,96)
    f=np.sqrt(2)*np.sin(np.pi*x)
    fp=np.sqrt(2)*np.pi*np.cos(np.pi*x)
    moments=[float(np.sum(w*f*f*x**j)) for j in range(3)]
    time_var=float(np.sum(w*fp*fp))
    exact=[1.,.5,1/3-1/(2*np.pi**2)]
    need(max(abs(a-b) for a,b in zip(moments,exact))<3e-14,'spectral moments')
    need(abs(time_var-np.pi**2)<1e-12,'time variance by Plancherel')
    points=np.array([0.,.7,np.pi,2*np.pi,7.,13.])
    integral=np.sum(w[None,:]*np.sin(np.pi*x)[None,:]*np.exp(-1j*points[:,None]*x),axis=1)
    stable=.5*(np.sinc((points-np.pi)/(2*np.pi))+np.sinc((points+np.pi)/(2*np.pi)))
    err=float(np.max(np.abs(integral-np.exp(-.5j*points)*stable)))
    need(err<3e-14,'independent Fourier integration including removable poles')
    # Integrate p(u), then give a rigorous analytic bound on the omitted tail.
    U=2048.
    gx,gw=np.polynomial.legendre.leggauss(16)
    centres=np.arange(-U+2,U,4.)
    u=centres[:,None]+2*gx[None,:]
    weights=2*gw[None,:]
    amp=.5*(np.sinc((u-np.pi)/(2*np.pi))+np.sinc((u+np.pi)/(2*np.pi)))
    dens=amp*amp/np.pi
    norm=float(np.sum(weights*dens))
    second=float(np.sum(weights*dens*u*u))
    tail0=128*np.pi/(27*U**3)
    tail2=128*np.pi/(9*U)
    need(norm<1+1e-12 and norm+tail0>1-1e-12,'normalization with analytic tail')
    need(second<np.pi**2 and second+tail2>np.pi**2,'time second moment with analytic tail')
    return dict(status='PASS',spectral_moments= moments,time_variance=time_var,
                Fourier_agreement=err,time_integral_U=U,integrated_probability=norm,
                omitted_probability_bound=float(tail0),integrated_time_second=second,
                omitted_second_moment_bound=float(tail2),mathematical_POVM_not_detector=True)

def group_endpoint():
    # Non-negligible chi is used only to diagnose the analytic relativistic bound.
    chi=.02; T=.82
    x,w=np.polynomial.hermite.hermgauss(28)
    v=chi/(2*T)
    p2=2*v*(x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2)
    weights=w[:,None,None]*w[None,:,None]*w[None,None,:]/np.pi**1.5
    p8=float(np.sum(weights*p2**4))
    need(abs(p8-945*v**4)<1e-14,'3D Gaussian eighth momentum moment')
    gap=(p2**2)/(4*(1+p2/2+np.sqrt(1+p2)))
    need(np.all(gap<=p2*p2/8+1e-15),'global positive-energy dispersion remainder')
    phase=T*gap/chi
    norm=float(np.sqrt(np.sum(weights*4*np.sin(phase/2)**2)))
    bound=math.sqrt(945)*chi/(32*T)
    need(0<norm<=bound,'actual free relativistic norm below all-spectrum bound')
    times=np.linspace(0,T,17)
    variance=chi*T/2*(1+(times/T)**2)
    need(abs(variance[-1]-chi*T)<1e-15,'no-covariance maximum variance')
    # An independent energy integral checks finite mean kinetic energy.
    kinetic=float(np.sum(weights*(p2/(np.sqrt(1+p2)+1))))
    kinetic_bound=3*chi/(4*T)
    need(0<kinetic<kinetic_bound,'complete COM kinetic energy retained')
    return dict(status='PASS',diagnostic_chi=chi,T=T,eighth_moment=p8,
                free_relativistic_norm=norm,analytic_norm_bound=bound,
                maximum_axis_variance=float(variance[-1]),
                positive_relativistic_mean_kinetic=kinetic,nonrelativistic_kinetic_upper=kinetic_bound,
                fixed_time_marginals_only=True)

def group_window():
    r=1e-8; T=.82; b=1e-3; zeta=1e-8
    dt=1e-8; dx=1e-9; pt=px=1e-4
    CR=math.sqrt(945)/(32*T)
    chi_time=r*dt*math.sqrt(pt)/math.pi
    chi_pos=px/(3*T/dx**2+CR)
    chi_star=min(chi_time,chi_pos); mu_star=b*zeta/2
    threshold=1/math.sqrt(chi_star*mu_star)
    ell=1e35; mass=1e5; B=1e12
    chi=1/(mass*ell); mu=mass/ell
    mlo=1/(ell*chi_star); mhi=mu_star*ell
    timing=math.pi**2/(B*dt)**2
    localisation=(3*T/dx**2+CR)*chi
    compactness=2*mu/b; signal_cap=chi*B
    need(2*b<.03 and chi<chi_star and mu<mu_star,'same physical scale interior')
    need(mlo<mass<mhi,'nonempty finite mass interval')
    need(signal_cap<r and timing<pt and localisation<px and compactness<zeta,'all screening criteria')
    need(abs(chi*mu*ell*ell-1)<1e-14,'fixed physical constants identity')
    rlow=1e-36; lowB=rlow/chi
    success_upper=min(1.,8*lowB*dt/math.pi**3)
    need(success_upper<1-pt,'low energy cannot preserve this timing preparation')
    increased=2*(mass*1e20/ell)/b
    need(increased>zeta,'heavier endpoint fails fixed-size compactness screen')
    small_ell=1e10; small_lo=1/(small_ell*chi_star); small_hi=mu_star*small_ell
    need(small_lo>small_hi,'selected sufficient-certificate interval empty, not universal no-go')
    return dict(status='PASS',inputs=dict(r=r,T=T,b=b,zeta=zeta,delta_t=dt,delta_x=dx,p_t=pt,p_x=px),
                C_R=CR,chi_time=chi_time,chi_position=chi_pos,chi_star=chi_star,mu_star=mu_star,
                certificate_L_over_lP_threshold=threshold,
                witness=dict(L_over_lP=ell,M_over_mP=mass,B=B,chi=chi,mu=mu,
                             mass_interval=[mlo,mhi],signal_energy_cap_ratio=signal_cap,
                             time_window_failure_bound=timing,position_failure_bound=localisation,
                             free_relativistic_norm_bound=CR*chi,body_compactness=compactness,
                             mean_COM_kinetic_ratio_upper=3*chi/(4*T)),
                deletions=dict(low_energy_ratio=rlow,low_energy_success_upper=success_upper,
                               heavy_body_compactness=increased,small_L_over_lP=small_ell,
                               small_L_certificate_mass_interval=[small_lo,small_hi]),
                actual_instrument=False,metric_error_remainder_certified=False)

def compute():
    return dict(schema='gr_response_scale_screen_1241_v1',round=1241,status='PASS',science_groups=3,
                group_1=group_packet(),group_2=group_endpoint(),group_3=group_window(),
                source_sha256={str((HERE/name).relative_to(RESEARCH)).replace('\\','/'):
                               hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                               for name in ['proof.md','admission.md','check.py']})

def compare(a,b):
    if isinstance(a,dict):
        need(isinstance(b,dict) and a.keys()==b.keys(),'result keys')
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        need(len(a)==len(b),'list length')
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):need(math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-300),'numeric replay')
    else:need(a==b,'result replay')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=compute();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        compare(data,json.loads(path.read_text('utf-8')));mode='READ_ONLY_MATCH'
    print(json.dumps(dict(round=1241,status='PASS',mode=mode,science_groups=3,
                         chi_star=data['group_3']['chi_star'],L_threshold=data['group_3']['certificate_L_over_lP_threshold']),ensure_ascii=False))
