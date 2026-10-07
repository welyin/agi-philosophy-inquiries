"""995: Z2-even internal Weinberg modulation, CP invariant and cyclic control.
No physical thermal collision kernel, generated charge or relic abundance is computed.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'internal_cp_bridge_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def exact(x):return dict(exact=str(x),value=float(x))
def inv(a,b):return float(np.trace(a.conj().T@b).imag)
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)

def run():
    c0=np.diag([.1,.2,.3]).astype(complex)
    c2=np.array([[.2j,.05,0],[.05,.1+.1j,.04j],[0,.04j,-.05+.15j]])
    assert np.max(abs(c2-c2.T))==0
    I=inv(c0,c2);assert abs(I-float(F(17,200)))<1e-15
    rng=np.random.default_rng(995)
    residuals=[]
    for _ in range(8):
        z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        u,_=np.linalg.qr(z)
        residuals.append(abs(inv(u.T@c0@u,u.T@c2@u)-I))
    assert max(residuals)<1e-14
    assert abs(inv(c0.conj(),c2.conj())+I)<1e-15
    # A real subtraction shift F_R -> F_R+c must also change C0.
    baseline=[]
    for shift in (-.4,.6):
        cp=c0-shift*c2
        for f in (.0,.07,.15):
            assert np.max(abs(cp+(f+shift)*c2-(c0+f*c2)))<1e-15
        baseline.append(abs(inv(cp,c2)-I))
    assert abs(inv(c0,.3*c0))<1e-15
    assert abs(inv(np.zeros((3,3)),c2))<1e-15

    old=read(STAGE/'991/dark_mode_screen_results.json')
    omega=F(1,2);assert omega**2==F(old['masses_squared']['dark']['exact'])
    vol=F(old['parameters']['volume']['exact']);lam=F(old['parameters']['lambda_dark']['exact'])
    scale=F(10) # New illustrative matching scale, not fitted to real neutrino masses.
    # Pure even squeezed Gaussian with exp(2r)=4; exact canonical moments.
    q2=1/(8*omega);p2=2*omega;q4=3*q2*q2
    energy=(p2+omega**2*q2)/2+lam*q4/(4*vol)
    f0=q2/(vol*scale**2)
    fdd=(2*p2-2*omega**2*q2-2*lam*q4/vol)/(vol*scale**2)
    assert q2*p2==F(1,4) and fdd>0
    # Quadratic control only. The quartic contribution above is NOT dropped
    # from the certified initial acceleration; long evolution below is H2 only.
    covariance=np.diag([float(q2),float(p2)])
    w=float(omega);period=math.pi/w
    harmonic=[];det_errors=[];energy_errors=[]
    for t in (0.,period/4,period/2,3*period/4,period):
        co=math.cos(w*t);si=math.sin(w*t)
        rot=np.array([[co,si/w],[-w*si,co]])
        vt=rot@covariance@rot.T
        ff=vt[0,0]/float(vol*scale**2)
        k=inv(c0+float(f0)*c2,c0+ff*c2)
        assert abs(k-I*(ff-float(f0)))<1e-14
        harmonic.append(dict(t=t,f=float(ff),two_time_cp_factor=k))
        det_errors.append(abs(np.linalg.det(vt)-.25))
        energy_errors.append(abs((vt[1,1]+w*w*vt[0,0])/2-float((p2+omega**2*q2)/2)))
    assert max(det_errors)<1e-13 and max(energy_errors)<1e-13

    # Unit-normalized causal kernel exp(-tau); algebraic control, not a thermal rate.
    A=(float(q2)+float(p2)/(w*w))/(2*float(vol*scale**2))
    B=(float(p2)/(w*w)-float(q2))/(2*float(vol*scale**2))
    Om=2*w;gamma=1.
    ts=np.arange(4096)*period/4096
    signal=I*B*(Om**2*np.cos(Om*ts)-gamma*Om*np.sin(Om*ts))/(gamma**2+Om**2)
    cycle=float(period*np.mean(signal))
    assert abs(cycle)<1e-16 and float(np.max(abs(signal)))>0
    # Explicit quadrature of f(t-tau)-f(t), compared with the convolution above.
    tau=np.linspace(0,32,8193);conv_errors=[]
    for t in (0.,.7,1.9,3.2):
        vals=np.exp(-tau)*(A-B*np.cos(Om*(t-tau))-(A-B*math.cos(Om*t)))
        numeric=I*float(np.trapezoid(vals,tau))
        analytic=I*B*(Om**2*math.cos(Om*t)-gamma*Om*math.sin(Om*t))/(gamma**2+Om**2)
        conv_errors.append(abs(numeric-analytic))
    assert max(conv_errors)<5e-12
    # Single-flavor local rephasing: theta_dot = Im(C* C_dot)/|C|^2.
    theta_errors=[]
    for f,fd in ((.01,.02),(.2,-.1),(.6,.03)):
        eps=1e-5;theta_plus=np.angle(1+1j*(f+eps*fd));theta_minus=np.angle(1+1j*(f-eps*fd))
        theta_errors.append(abs((theta_plus-theta_minus)/(2*eps)-fd/(1+f*f)))
    assert max(theta_errors)<1e-10
    inputs=[Path(__file__).resolve(),HERE/'drafts/STATUS.md',HERE/'drafts/selection.md',
        STAGE/'993/common_candidate_v1.md',STAGE/'994/boundary_adoption_v1.md',
        STAGE/'research_note_994.md',STAGE/'991/dark_mode_screen_results.json',
        STAGE.parent/'archive_629_652/research_note_629.md',
        STAGE.parent/'archive_702_741/research_note_710.md',
        STAGE.parent/'archive_301_341/research_note_332.md']
    return dict(round=995,all_scientific_checks_passed=True,
        kind='internal_CP_mechanism_adoption_not_charge_generation',
        invariant=exact(F(17,200)),flavor_basis_max_residual=max(residuals),
        baseline_shift_max_residual=max(baseline),CP_conjugate_invariant=inv(c0.conj(),c2.conj()),
        parity_even_initial_state=dict(omega=exact(omega),volume=exact(vol),scale=exact(scale),
            canonical_Q2=exact(q2),canonical_P2=exact(p2),canonical_Q4=exact(q4),
            energy_including_quartic=exact(energy),f_initial=exact(f0),
            f_second_derivative_including_quartic=exact(fdd),mean_sigma=0),
        quadratic_control=dict(rows=harmonic,covariance_det_max_residual=max(det_errors),
            energy_max_residual=max(energy_errors)),
        causal_periodic_control=dict(period=period,cycle_integral=cycle,
            pointwise_max=float(np.max(abs(signal))),quadrature_max_residual=max(conv_errors)),
        local_rephasing_derivative_max_residual=max(theta_errors),
        added_fundamental_species=0,added_effective_operator_dimension=7,
        actual_thermal_kernel_computed=False,net_lepton_charge_generated=False,
        expanding_cosmological_history_verified=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','quadratic_control')},ensure_ascii=False,indent=2))
