"""931: original charged mass data constrain a partial EM response and scalar source.
Canonical continuum free-fermion vacuum / one fermion loop only.
No prediction of bare gauge couplings, QCD response or spacetime generation.
"""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib,argparse
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
STORE=Layout()
with ResearchRuntime(STORE).installed():
    import joint_continuum_source_spectrum as inherited
TARGET=HERE/'charged_matter_response_results.json'
QSTAR=inherited.QSTAR

def f(q):return np.sqrt(6)*np.sinh(q/np.sqrt(6))
def kappa(q):return 1/(np.sqrt(6)*np.tanh(q/np.sqrt(6)))
def gauss(n=128):
    x,w=np.polynomial.legendre.leggauss(n)
    return (x+1)/2,w/2

def original():
    M,masses,weights,phi,direction=inherited.original_data()
    species=[]
    for label,start,count,charge in [('u',0,3,Fraction(2,3)),('d',6,3,Fraction(-1,3)),('e',12,1,Fraction(-1))]:
        m=abs(M[start,start+1]);raw=inherited.finite.Y[label]
        expected=abs(raw)*phi[0]/np.sqrt(inherited.old.original.F(phi))
        assert abs(m-expected)<1e-15
        species.append(dict(label=label,mass=float(m),colors=count,charge=float(charge),weight=float(count*charge**2)))
    electric=[]
    for label,start,count,charge in [('u',0,3,2/3),('d',6,3,-1/3),('e',12,1,-1)]:
        for _ in range(count):electric.extend([charge,-charge])
    electric.extend([0.,0.]);Q=np.diag(electric)
    assert np.max(abs(Q@M+M@Q))<1e-15
    charge_trace=float(np.trace(Q@Q).real/2)
    assert abs(charge_trace-8/3)<1e-15
    old=json.loads((ROOT/'research_cognition_physics/archive_629_652/630/joint_continuum_source_spectrum_results.json').read_text('utf-8'))
    olderr=float(np.max(abs(np.sort(masses)-old['original_mass_and_metric']['masses_sorted'])))
    assert olderr<1e-15
    return species,olderr,charge_trace

def kernel(z,q,species,n=128):
    """Pi(z)-Pi(0) real below pair threshold; z=-Q^2 for spacelike."""
    x,w=gauss(n);u=x*(1-x);a=f(q)/f(QSTAR)
    return float(sum(-s['weight']/(2*np.pi**2)*np.dot(w,u*np.log1p(-z*u/(s['mass']*a)**2)) for s in species))

def dispersion(z,q,species):
    v,w=gauss();a=f(q)/f(QSTAR)
    return float(sum(z*s['weight']/(12*np.pi**2)*np.dot(w,v*v*(3-v*v)/(4*(s['mass']*a)**2-z*(1-v*v))) for s in species))

def spectral(sval,q,species):
    a=f(q)/f(QSTAR);total=0.
    for row in species:
        m=row['mass']*a
        if sval>4*m*m:
            b=np.sqrt(1-4*m*m/sval)
            total+=row['weight']/(12*np.pi**2)*(1+2*m*m/sval)*b
    return float(total)

def radial_derivative(z,q,species):
    x,w=gauss();u=x*(1-x);a=f(q)/f(QSTAR)
    return float(sum(-s['weight']*kappa(q)/np.pi**2*np.dot(w,u*(z*u)/(s['mass']**2*a*a-z*u)) for s in species))

def run():
    species,olderr,charge_trace=original()
    # Independent current-projector evaluation of the cut, including spin.
    I=np.eye(2);Z=np.zeros((2,2));sigmas=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    aa=[np.block([[Z,s],[s,Z]]) for s in sigmas];beta=np.diag([1.,1.,-1.,-1.])
    cuterr=0.;paramcuterr=0.
    for row in species:
        m=row['mass']
        for r in (.2,.8,2.):
            p=r*m;E=np.hypot(p,m);H=aa[2]*p+beta*m
            pp=(np.eye(4)+H/E)/2;pm=np.eye(4)-pp
            avg=sum(np.trace(pm@A@pp@A).real for A in aa)/3
            rho_trace=row['weight']*(p/E)*avg/(16*np.pi**2)
            rho=row['weight']/(12*np.pi**2)*(1+m*m/(2*E*E))*(p/E)
            cuterr=max(cuterr,abs(rho-rho_trace))
            b=p/E;lo=(1-b)/2;hi=(1+b)/2
            primitive=lambda x:x*x/2-x*x*x/3
            rho_log=row['weight']/(2*np.pi**2)*(primitive(hi)-primitive(lo))
            paramcuterr=max(paramcuterr,abs(rho-rho_log))
    gap=2*min(r['mass'] for r in species)
    C1=sum(r['weight']/r['mass']**2 for r in species)/(60*np.pi**2)
    C2=sum(r['weight']/r['mass']**4 for r in species)/(560*np.pi**2)
    rows=[];disperr=0.;deriverr=0.;scalingerr=0.
    for ratio in (.1,.35,.65):
        scale=ratio*gap
        for sign in (-1,1):
            z=sign*scale**2;val=kernel(z,QSTAR,species);disp=dispersion(z,QSTAR,species)
            disperr=max(disperr,abs(val-disp))
            rem=abs(val-C1*z)
            bound=C2*z*z/(1-ratio**2) if sign==1 else C2*z*z
            assert rem<=bound*(1+1e-11)
            h=1e-6;fd=(kernel(z,QSTAR+h,species)-kernel(z,QSTAR-h,species))/(2*h)
            analytic=radial_derivative(z,QSTAR,species);deriverr=max(deriverr,abs(fd-analytic))
            for q in (.2,.4):
                a=f(q)/f(QSTAR)
                if z<0 or z<gap*gap*a*a:
                    scalingerr=max(scalingerr,abs(kernel(z,q,species)-kernel(z/a**2,QSTAR,species)))
            rows.append(dict(momentum_squared=z,type='spacelike' if sign<0 else 'timelike_below_threshold',threshold_fraction=ratio,kernel=val,leading=C1*z,remainder=rem,remainder_bound=bound,radial_source_derivative=analytic,finite_difference_error=abs(fd-analytic)))
    # Same finite off-shell source at two original scalar backgrounds.
    zprobe=-(.35*gap)**2
    backgrounds=[dict(q=q,subtracted_kernel=kernel(zprobe,q,species),charged_threshold=gap*f(q)/f(QSTAR)) for q in (.27,.4)]
    delta=abs(backgrounds[0]['subtracted_kernel']-backgrounds[1]['subtracted_kernel'])
    # The local determinant part; an independent finite local matching term remains legal.
    coefficient=charge_trace/(6*np.pi**2)
    qp=.4;local=-coefficient*np.log(f(qp)/f(QSTAR))
    slope=-coefficient*kappa(QSTAR)
    c=.07;local_added=c*(qp-QSTAR)
    spectral_rows=[dict(energy=w,density=spectral(w*w,QSTAR,species)) for w in (.9*gap,1.1*gap,.3,.7)]
    assert spectral_rows[0]['density']==0 and spectral_rows[1]['density']>0
    assert max(cuterr,paramcuterr,disperr,scalingerr)<1e-13 and deriverr<1e-9 and delta>1e-5
    deps={str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    for key in STORE.read_paths:
        actual=STORE.entries[key][0] if key in STORE.entries else ROOT/key
        if actual.is_file():deps[str(actual.relative_to(ROOT))]=hashlib.sha256(actual.read_bytes()).hexdigest()
    frozen=ROOT/'research_cognition_physics/archive_629_652/630/joint_continuum_source_spectrum_results.json'
    deps[str(frozen.relative_to(ROOT))]=hashlib.sha256(frozen.read_bytes()).hexdigest()
    return dict(round=931,date='2026-10-06',all_scientific_checks_passed=True,q_star=QSTAR,species=species,
        physical_Dirac_charge_squared_sum=charge_trace,frozen_full_mass_spectrum_error=olderr,
        current_projector_spectral_error=cuterr,Feynman_cut_spectral_error=paramcuterr,
        independent_dispersion_error=disperr,radial_scaling_error=scalingerr,radial_source_derivative_error=deriverr,
        charged_first_pair_threshold=gap,C1=C1,C2=C2,finite_momentum_rows=rows,spectral_rows=spectral_rows,
        shared_background_comparison=backgrounds,subtracted_response_difference=delta,
        local_determinant_log_coefficient=coefficient,local_determinant_change=local,local_determinant_radial_slope=slope,
        finite_matching_slope_shift=c,finite_matching_change=local_added,
        same_local_matching_at_q_star_but_different_radial_source=True,
        finite_local_polynomials_cannot_cancel_pair_cut=True,
        original_matter_used_no_extra_dilaton_species=True,
        scalar_source_and_EM_kernel_from_same_determinant=True,
        spacelike_subtraction_removes_only_zero_momentum_local_coefficient=True,
        no_radial_axion_generated_with_fixed_mass_phases=True,
        bare_gauge_normalization_selected=False,all_higher_derivative_matching_coefficients_fixed=False,
        full_interacting_SM_or_hadron_response_computed=False,
        graph_continuum_mapping_or_actual_reference_instrument_constructed=False,
        free_on_shell_photon_speed_modified_claimed=False,
        all_six_protocols_realized_in_continuum_candidate=False,
        full_backreaction_or_GR_generated=False,whole_stage_completed=False,full_goal_completed=False,source_hashes=deps)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    d=run()
    if a.write:
        with TARGET.open('x',encoding='utf-8') as fobj:json.dump(d,fobj,ensure_ascii=False,indent=2);fobj.write('\n')
    print(json.dumps({k:v for k,v in d.items() if k not in ('source_hashes','finite_momentum_rows','spectral_rows')},ensure_ascii=False,indent=2))
