"""Read-only reproducible qualification, not a numbered research round."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math, cmath, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
FROZEN={'research_cognition_physics/archive_1120_/research_note_1125.md': '98b7bcd6cf8aca6854051f0b4b6277f63567009d55a34056f7457c7b194ac634', 'research_cognition_physics/archive_819_853/research_note_821.md': 'b3b87bb1a75a1ab7a19beb3caf126a7e2770515094a94898ed11e1c998c493df', 'research_cognition_physics/archive_819_853/research_note_853.md': '16f7b4f5aa31ec7f1d664d81f707e3e7fd4db98cc5317c14c31a580a16e8be35', 'research_cognition_physics/archive_742_763/research_note_745.md': '6c0173877324a30e82ef8b29a04cc0c0f96e6048b2392061ff7fddd99ef39c8f', 'research_cognition_physics/archive_742_763/research_note_757.md': '97f0cb9ace6204b2fa5aacc9a88fa7805f651e27e39dede964ab97401ce46972', 'research_cognition_physics/archive_854_872/research_note_872.md': 'fd0c16968a98e81e56c6055dc4da292ab8e65b89e6a98696b83021cd2bb3d892', 'research_cognition_physics/archive_873_893/research_note_876.md': 'dc0db4fb188d44968806f96acd430f506b80ec66cf104e02513f661651bcda99', 'research_cognition_physics/archive_873_893/research_note_881.md': '27089c6a594c938e725027143cf758093d3e21584ac155048a5da07bc5ce3134'}

def exp_upper(x,n=20):
    term=F(1); total=term
    for k in range(1,n+1):
        term*=x/k; total+=term
    first=term*x/(n+1)
    return total+first/(1-x/(n+2))

def sin_integral(a,T):
    return T if abs(a)<1e-12 else math.sin(a*T)/a

def exp_integral(a,T):
    return T if abs(a)<1e-12 else (cmath.exp(1j*a*T)-1)/(1j*a)

def oscillator(T,w):
    nu=2*math.pi/T; a=1/(2*math.sqrt(T))
    assert abs(w*w-nu*nu)>1e-6
    def D(v):
        return sin_integral(v,T)-(sin_integral(v+nu,T)+sin_integral(v-nu,T))/2
    phase=a*a*((D(0)-D(w))/(w*w)-(D(nu)-D(w))/(w*w-nu*nu))
    disp=a*(exp_integral(w,T)-(exp_integral(w+nu,T)+exp_integral(w-nu,T))/2)/math.sqrt(2*w)
    grid=np.linspace(0,T,16385)
    pulse=a*(1-np.cos(nu*grid))
    q=a*((1-np.cos(w*grid))/(w*w)-(np.cos(nu*grid)-np.cos(w*grid))/(w*w-nu*nu))
    qdd=a*(np.cos(w*grid)+(nu*nu*np.cos(nu*grid)-w*w*np.cos(w*grid))/(w*w-nu*nu))
    residual=float(np.max(np.abs(qdd+w*w*q-pulse)))
    integral=float(np.trapezoid(pulse*q,grid))
    phase_bound=2*math.pi/(T*w**3)
    vector_bound=4*math.pi/(math.sqrt(2)*T**1.5*w**2.5)
    assert abs(phase-(3/8)/w**2)<=phase_bound+1e-13
    assert abs(disp)<=vector_bound+1e-13
    assert residual<1e-12 and abs(integral-phase)<1e-8
    return dict(T=T,omega=w,phase=phase,static_phase=(3/8)/w**2,
                phase_error=abs(phase-(3/8)/w**2),phase_bound=phase_bound,
                residual_vector=abs(disp),vector_bound=vector_bound,
                ODE_residual=residual,quadrature_crosscheck=abs(integral-phase))

def main():
    for rel,digest in FROZEN.items():
        assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
    assert exp_upper(F(1,3))<F(3,2) and exp_upper(F(5),40)<150
    assert F(7,66000)>F(1,10000)
    rho2=F(3); lam2=F(12000); n2=F(24); T=10**9
    phase_bound=lam2*n2*rho2/T
    # sqrt(3)<7/4, sqrt(12000)<110, sqrt(2)>7/5; T^(3/2)>3*10^13.
    vector_bound=F(6600,3*10**13)
    assert phase_bound+vector_bound<F(865,10**6)
    rows=[oscillator(T,w) for T in [4.,8.,16.] for w in [1.,1.7,2.3]]
    # The cross term has no extra 1/2. A projector and B=X have eigenvalues 0,1 and -1,+1.
    labels=[(a,b) for a in [0,1] for b in [-1,1]]
    alpha,beta,theta=.31,.17,math.pi/2
    direct=np.diag([cmath.exp(1j*(alpha*a*a+beta*b*b+theta*a*b)) for a,b in labels])
    factor=np.diag([cmath.exp(1j*alpha*a) for a,b in labels])*cmath.exp(1j*beta)
    factor=factor@np.diag([cmath.exp(1j*theta*a*b) for a,b in labels])
    factor_error=float(np.linalg.norm(direct-factor))
    assert factor_error<1e-14
    # Vacuum coherent-state Gram retains the original quantum environment.
    hA=np.array([.002+.001j,.003-.001j]); hB=np.array([-.001+.0015j,.002+.002j])
    hs=[a*hA+b*hB for a,b in labels]
    gram=np.array([[cmath.exp(-np.vdot(x,x).real/2-np.vdot(y,y).real/2+np.vdot(x,y)) for y in hs] for x in hs])
    mineig=float(np.linalg.eigvalsh(gram).min())
    assert mineig>-1e-13 and np.max(np.abs(np.diag(gram)-1))<1e-14
    vector_exact=max(math.sqrt(2*(1-math.exp(-np.vdot(h,h).real/2))) for h in hs)
    vector_triangle=float(np.linalg.norm(hA)+np.linalg.norm(hB))
    assert vector_exact<=vector_triangle
    data={'status':'PASS_QUALIFICATION_ONLY','numbered_round':None,'fresh_numbered_groups':0,
          'checked_groups':4,'frozen_sources':FROZEN,
          'rational_finite_T_certificate':{'mass':1,'ball_radius':1,'separation':3,'T':T,
          'rho_L2_squared_upper':'3','K_lower':'1/10000','C_f_exact':'10/7','lambda_squared_upper':'12000',
          'phase_error_upper':str(phase_bound),'residual_vector_upper':str(vector_bound),
          'complete_isometry_upper':str(phase_bound+vector_bound),'target_isometry_tolerance':'865/1000000'},
          'frequency_diagnostics':rows,'full_QP_algebra':{'self_phase_factorization_error':factor_error,
          'coherent_Gram_min_eigenvalue':mineig,'vacuum_vector_distance':vector_exact,
          'triangle_bound':vector_triangle},
          'scope':'Analytic continuum source criterion plus independent per-frequency/algebra diagnostics. No discretized oscillator is identified with the full KG field.',
          'fixed_profile_exact_closure_no_go_is_analytic':True,
          'worldtube_cross_spacelike_phase_zero_is_analytic':True,
          'actual_material_Y_clock_autonomy_or_point_microcausality_proved':False,
          'same_1143_task_deadline_matched':False,'full_stage_completed':False}
    target=Path(__file__).with_name('after1143_native_record_reuse_results.json')
    if '--write' in sys.argv:
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        previous=json.loads(target.read_text(encoding='utf-8'))
        assert previous==data,'saved result differs'
    print(json.dumps({'status':data['status'],'checked_groups':4,'result_matches':True,
                      'complete_isometry_upper':float(phase_bound+vector_bound)},ensure_ascii=False))
if __name__=='__main__':
    main()
