"""One finite thermoelastic compatibility screen, not microscopic matter generation."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'thermoelastic_common_account_results.json'


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def n(x):return dict(exact=str(x),value=float(x))
def dot(x,y):return sum((a*b for a,b in zip(x,y)),F(0))
def multiply(a,x):return [dot(row,x) for row in a]


def operators(size):
    D=[[F(0) for _ in range(size)] for _ in range(size)]
    L=[[F(0) for _ in range(size)] for _ in range(size)]
    for j in range(size):
        D[j][(j+1)%size]+=F(1,2);D[j][(j-1)%size]-=F(1,2)
        L[j][j]-=2;L[j][(j+1)%size]+=1;L[j][(j-1)%size]+=1
    assert all(sum(row)==0 for row in D+L)
    assert all(D[i][j]==-D[j][i] and L[i][j]==L[j][i] for i in range(size) for j in range(size))
    assert all(L[i][j]>=0 for i in range(size) for j in range(size) if i!=j)
    return D,L


def rates(eps,v,T,D,L,rho,K,b,c,T0,kappa,thermal_b):
    theta=[t-T0 for t in T]
    sigma=[K*e-b*d for e,d in zip(eps,theta)]
    de=multiply(D,v)
    dv=[z/rho for z in multiply(D,sigma)]
    dT=[(kappa*l-thermal_b*t*x)/c for l,t,x in zip(multiply(L,T),T,de)]
    dE=dot([rho*x for x in v],dv)+dot([K*x+b*T0 for x in eps],de)+c*sum(dT)
    dS=b*sum(de)+c*sum(dt/t for dt,t in zip(dT,T))
    dA=dot([rho*x for x in v],dv)+dot([K*x for x in eps],de)+dot([c*(1-T0/t) for t in T],dT)
    return dict(strain_rate=[n(x) for x in de],velocity_rate=[n(x) for x in dv],
        temperature_rate=[n(x) for x in dT],energy_rate=n(dE),entropy_rate=n(dS),availability_rate=n(dA))


def run():
    rho,K,b,c,T0,kappa=F(1),F(4),F(1),F(2),F(1),F(1,2)
    eps=[F(1,100),F(0),F(-1,100),F(0)]
    v=[F(0),F(1,10),F(0),F(-1,10)]
    T=[F(6,5),F(1),F(4,5),F(1)]
    D,L=operators(4)
    good=rates(eps,v,T,D,L,rho,K,b,c,T0,kappa,b)
    bad=rates(eps,v,T,D,L,rho,K,b,c,T0,kappa,F(0))
    e=lambda row,key:F(row[key]['exact'])
    edge_entropy=kappa*sum((T[j]-T[(j+1)%4])**2/(T[j]*T[(j+1)%4]) for j in range(4))
    assert e(good,'energy_rate')==0
    assert e(good,'entropy_rate')==edge_entropy==F(1,12)>0
    assert e(good,'availability_rate')==-T0*edge_entropy
    expected_defect=b*dot([t-T0 for t in T],multiply(D,v))
    assert e(bad,'energy_rate')==expected_defect==F(1,25)
    assert e(bad,'entropy_rate')==edge_entropy
    # Independent floating Legendre and availability identities of the full potentials.
    E=float(dot([rho*x/2 for x in v],v)+dot([K*x/2 for x in eps],eps)+b*T0*sum(eps)+c*sum(t-T0 for t in T))
    S=float(b*sum(eps))+float(c)*sum(math.log(float(t/T0)) for t in T)
    availability=sum(float(rho*vel*vel/2+K*ep*ep/2)+float(c)*(float(t-T0)-float(T0)*math.log(float(t/T0)))
                     for ep,vel,t in zip(eps,v,T))
    assert abs(E-float(T0)*S-availability)<1e-14 and availability>0
    for ep,t in zip(eps,T):
        f=float(K*ep*ep/2-b*ep*(t-T0)+c*(t-T0))-float(c*t)*math.log(float(t/T0))
        s=float(b*ep)+float(c)*math.log(float(t/T0))
        u=float(K*ep*ep/2+b*T0*ep+c*(t-T0))
        assert abs(f+float(t)*s-u)<2e-15
    # Continuous symbol at one declared k, not the four-cell lattice dispersion.
    k=F(1)
    coefficients=[rho*c,rho*kappa*k*k,(K*c+b*b*T0)*k*k,K*kappa*k**4]
    routh=coefficients[1]*coefficients[2]-coefficients[0]*coefficients[3]
    assert all(x>0 for x in coefficients) and routh==rho*kappa*b*b*T0*k**4==F(1,2)>0
    roots=np.roots([float(x) for x in coefficients])
    symbol=np.array([[0,1j*float(k),0],
        [1j*float(k*K/rho),0,-1j*float(k*b/rho)],
        [0,-1j*float(k*b*T0/c),-float(kappa*k*k/c)]],dtype=complex)
    eigen=np.linalg.eigvals(symbol)
    assert max(z.real for z in roots)<0
    assert max(min(abs(z-w) for w in eigen) for z in roots)<1e-12
    # The centered four-cell D also has a checkerboard null vector: no false all-mode claim.
    checker=[F(1),F(-1),F(1),F(-1)]
    assert multiply(D,checker)==[F(0)]*4
    sources=[HERE/'drafts/adoption_decision.md',HERE/'drafts/STATUS.md',
        STAGE/'999/research_round_999_checks.json',STAGE/'999/overall_operation_hypothesis_v1.md',
        STAGE.parent/'archive_629_652/research_note_645.md',
        STAGE.parent/'archive_342_369/research_note_354.md',
        STAGE.parent/'archive_342_369/research_note_355.md',Path(__file__)]
    return dict(round=1000,date='2026-10-07',kind='classical_effective_thermoelastic_dynamic_screen',
        all_scientific_checks_passed=True,parameters={key:n(value) for key,value in
            zip(('rho','K','b','c','T0','kappa'),(rho,K,b,c,T0,kappa))},
        finite_state=dict(strain=[n(x) for x in eps],velocity=[n(x) for x in v],temperature=[n(x) for x in T]),
        matched=good,omitted_thermal_feedback=bad,edge_entropy_production=n(edge_entropy),
        state_function_diagnostics=dict(energy_relative_reference=E,entropy_relative_reference=S,availability=availability),
        continuous_linear_symbol=dict(k=n(k),polynomial_coefficients=[n(x) for x in coefficients],
            routh_margin=n(routh),roots=[dict(real=float(z.real),imag=float(z.imag)) for z in roots],
            isothermal_speed_squared=n(K/rho),adiabatic_speed_squared=n((K+b*b*T0/c)/rho)),
        discrete_checkerboard_null_mode_preserved=True,long_time_ode_simulated=False,
        nonlinear_temperature_constitutive_branch=True,microscopic_phase_derived=False,
        transport_coefficients_derived=False,noise_or_memory_lifetime_certified=False,
        all_discrete_nonuniform_modes_damped=False,einstein_geometry_derived=False,
        quantum_channel_constructed=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sources})


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14),(a,b)
    else:assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps(dict(round=1000,all_scientific_checks_passed=True,
        matched_energy_rate=out['matched']['energy_rate'],
        entropy_rate=out['matched']['entropy_rate'],
        incomplete_model_energy_rate=out['omitted_thermal_feedback']['energy_rate'],
        roots=out['continuous_linear_symbol']['roots'],full_goal_completed=False),ensure_ascii=False,indent=2))
