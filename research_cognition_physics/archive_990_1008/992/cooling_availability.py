"""992: differential redshift and athermality in an inherited finite material.
No cosmological trajectory, heat engine, reset device or infinite bath is built.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'cooling_availability_results.json'

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def entropy(p):return float(-p@np.log(p))

def run():
    base,h,q,f,_,native_residuals=load('native968_for992',STAGE/'968/internal_relay.py').material()
    h=h+.05*f
    e,v=np.linalg.eigh(h); eps=e-e[0]
    beta0=2*math.log(2);temp0=1/beta0
    def gibbs(beta):
        logz=float(np.logaddexp.reduce(-beta*eps))
        logp=-beta*eps-logz
        return np.exp(logp),logp,logz
    p0,logp0,_=gibbs(beta0)
    old=read(STAGE/'980/finite_thermal_records_results.json')['parameters']
    assert np.max(abs(e-np.array(old['energies'])))<1e-13
    assert np.max(abs(p0-np.array(old['thermal_populations'])))<1e-13
    s0=entropy(p0);energy_exc=float(p0@eps);energy=float(p0@e)
    rho=(v*p0)@v.T;logrho=(v*logp0)@v.T
    assert np.linalg.norm(h@rho-rho@h)<1e-14
    mass_offset=200.;M=mass_offset+energy
    frequency=base['J'];modes=12;V0=1.
    nbar=1/math.expm1(beta0*frequency)
    R=modes*frequency*nbar
    sg=modes*((1+nbar)*math.log1p(nbar)-nbar*math.log(nbar))
    def evaluate(a):
        temp=temp0/a
        pa,logpa,logz=gibbs(beta0*a)
        d=float(p0@(logp0-logpa))
        availability=temp*d
        direct=energy_exc-temp*s0+temp*logz
        logsigma=(v*logpa)@v.T
        matrix_d=float(np.trace(rho@(logrho-logsigma)))
        sref=entropy(pa)
        derivative=temp0/(a*a)*(s0-sref)
        # Compare free energies with the actual full mass baseline kept in each.
        full_f_actual=mass_offset+energy-temp*s0
        full_f_reference=mass_offset+float(pa@e)-temp*sref
        assert max(abs(availability-direct),abs(d-matrix_d),
                   abs(availability-(full_f_actual-full_f_reference)))<1e-12
        scaled_p=np.exp(-(beta0*a)*(eps/a));scaled_p/=sum(scaled_p)
        scaled_control=float(temp*np.sum(p0*np.log(p0/scaled_p)))
        assert abs(scaled_control)<1e-14
        pressure=R/(3*V0*a**4)
        balance=-R/a**2+pressure*3*V0*a*a
        assert abs(balance)<1e-12
        return dict(a=a,temperature=temp,actual_populations=p0.tolist(),
            reference_populations=pa.tolist(),actual_material_entropy=s0,
            reference_entropy=sref,relative_entropy=d,availability=availability,
            analytic_da=derivative,uniform_rescaling_control=scaled_control,
            equilibrium_state_trace_distance=float(np.sum(abs(p0-pa))/2),
            radiation_energy=R/a,total_cell_energy=M+R/a,
            pressure=pressure,pressure_work_balance_residual=balance,
            entropy_of_material_plus_radiation=s0+sg)
    rows=[evaluate(a) for a in (.5,1.,2.,4.,8.,32.)]
    assert rows[0]['availability']>0
    assert abs(rows[1]['availability'])<1e-14
    cold=[r for r in rows if r['a']>1]
    assert all(0<r['availability']<energy_exc and r['analytic_da']>0 for r in cold)
    assert all(x['availability']<y['availability'] for x,y in zip(cold,cold[1:]))
    derivative_errors=[]
    for row in rows:
        a=row['a'];step=2e-5*a
        numerical=(evaluate(a+step)['availability']-evaluate(a-step)['availability'])/(2*step)
        derivative_errors.append(abs(numerical-row['analytic_da']))
    assert max(derivative_errors)<1e-8
    variance=float(p0@eps**2-energy_exc**2)
    source_files=[Path(__file__),HERE/'drafts/selection.md',HERE/'drafts/STATUS.md',
        STAGE/'956/native_material_interface.py',STAGE/'968/internal_relay.py',
        STAGE/'980/finite_thermal_records_results.json',STAGE/'research_note_964.md',
        STAGE/'research_note_969.md',STAGE/'research_note_980.md',
        STAGE.parent/'archive_301_341/research_note_305.md',
        STAGE.parent/'archive_301_341/research_note_325.md']
    return dict(round=992,all_scientific_checks_passed=True,
        inputs=dict(eta=.05,beta0=beta0,temp0=temp0,mass_offset=mass_offset,
            photon_modes=modes,photon_initial_frequency=frequency,initial_mean_occupation=nbar,
            V0=V0,initial_R=R,material_mean_mass=M,material_energies=e.tolist()),
        initial_material_entropy=s0,radiation_entropy=sg,
        expansion_availability_limit=energy_exc,
        near_equilibrium_quadratic_coefficient=beta0*variance/2,
        rows=rows,derivative_max_residual=max(derivative_errors),
        native_material_residuals=native_residuals,
        decision='differential scaling generates a bounded athermality resource; not an automatic arrow or memory reset',
        physical_reset_implemented=False,work_extraction_implemented=False,
        initial_cosmological_boundary_explained=False,full_quantum_geometry_verified=False,
        additional_blank_registers_created=0,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-11),(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps(dict(round=992,all_scientific_checks_passed=True,
        expansion_availability_limit=out['expansion_availability_limit'],
        derivative_max_residual=out['derivative_max_residual'],
        rows=[{k:r[k] for k in ('a','temperature','availability','relative_entropy',
            'equilibrium_state_trace_distance')} for r in out['rows']]),ensure_ascii=False,indent=2))
