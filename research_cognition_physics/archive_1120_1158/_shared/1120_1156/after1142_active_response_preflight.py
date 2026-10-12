"""Unnumbered active F3 response preflight; default read-only.
Reuses the frozen 1139 original-g certificate, without storing large arrays.
Analytic full-CAR response is in the sibling Markdown; finite matrices diagnose it.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCES={
 'local1139_tail_certificate.py':'fbc400a14c92788f00efd4579f4197e96ff579bacc2f42e1ffb92323ee478dd3',
 'local1139_tail_certificate.md':'a502e98a6bcf634e6976c394a57bcf226b402447cfad1b61b1fe5d5507427a90',
 'after1142_joint_contract_review.md':'dcb3316bcedcb86a4e0a7416652445b43d63e382d5e7c488a24929d950891b59',
}
def frozen():
    for name,want in SOURCES.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==want,name
    return SOURCES

def tail_check():
    p=HERE/'local1139_tail_certificate.py'
    spec=importlib.util.spec_from_file_location('original1139_tail',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    ns,k,c,e,z,dz=m.coefficients()
    S=float(c@c);Slo=float(np.maximum(abs(c)-e,0)@np.maximum(abs(c)-e,0))
    grad=math.sqrt(3*float((k*k*c)@c)*S*S);h=m.L/m.GRID
    dproj=math.sqrt(max(0,1-Slo**3));dcoeff=3.01*float(np.linalg.norm(e))
    dcell=h/math.sqrt(12)*grad;delta=dproj+dcoeff+dcell+1e-6
    pc,pd,roots=m.phase_tables()
    n2=ns[:,None,None]**2+ns[None,:,None]**2+ns[None,None,:]**2
    base=c[:,None,None]*c[None,:,None]*c[None,None,:];ds=base*pd[n2]
    small=[base*pc[n2],ds,1j*k[:,None,None]*ds,
           1j*k[None,:,None]*ds,1j*k[None,None,:]*ds]
    ix=np.ix_(ns%m.GRID,ns%m.GRID,ns%m.GRID)
    density=np.zeros((m.GRID,)*3);fft_error=0.
    for a in small:
        full=np.zeros((m.GRID,)*3,dtype=np.complex128);full[ix]=a
        val=m.inverse_fft(full,roots)*m.GRID**3/m.L**1.5
        diagnostic=np.fft.ifftn(full)*m.GRID**3/m.L**1.5
        fft_error=max(fft_error,float(np.max(abs(val-diagnostic))))
        density+=abs(val)**2
    x=np.fft.fftfreq(m.GRID)*m.L
    r=np.sqrt(x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2)
    rows={}
    for label,rin,rout in [('old_1139',1.9,2.),('active_F3',1.1,1.25)]:
        v=np.clip((np.maximum(r-math.sqrt(3)*h/2,0)-rin)/(rout-rin),0,1)
        low=35*v**4-84*v**5+70*v**6-20*v**7
        a_cell=float(np.sum(low*density)*h**3)
        b_cell=float(np.sum(low*low*density)*h**3)
        mu=max(0,math.sqrt(max(0,a_cell-1e-6))-delta)**2
        nu=max(0,math.sqrt(max(0,b_cell-1e-6))-delta)**2
        rows[label]=dict(rin=rin,rout=rout,mu_lower=mu,nu_lower=nu,
                         q_lower=mu,cell_mu=a_cell,cell_nu=b_cell)
    saved=json.loads((HERE/'local1139_tail_certificate.json').read_text(encoding='utf-8'))
    assert math.isclose(rows['old_1139']['mu_lower'],saved['q_lower'],
                        rel_tol=1e-10,abs_tol=1e-12)
    assert rows['old_1139']['mu_lower']>.035
    assert rows['active_F3']['mu_lower']>.4
    assert grad<3.1 and fft_error<1e-10
    return dict(projection_error=dproj,coefficient_error=dcoeff,cell_error=dcell,
                arithmetic_norm_budget=1e-6,complete_L2_error=delta,
                norm_diagnostic=float(np.sum(density)*h**3),
                independent_fft_max_error=fft_error,masks=rows)

def car_check():
    # All 16 states of four CAR modes g1,g2,e1,e2; not one-particle truncation.
    n=4;dim=1<<n;cs=[]
    for j in range(n):
        a=np.zeros((dim,dim),complex)
        for b in range(dim):
            if (b>>j)&1:
                a[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
        cs.append(a)
    def q(A):
        return sum(A[i,j]*(cs[i].conj().T@cs[j]) for i in range(n) for j in range(n))
    mu=.4;nu=.35;b=math.sqrt(nu-mu*mu);kap=math.pi/4
    sy=np.array([[0,-1j],[1j,0]])
    spatial=np.array([[mu,b/2],[b/2,0.]])
    A=np.kron(spatial,sy);Q=q(A);Y=q(np.kron(np.diag([1.,0.]),sy))
    effect=cs[0].conj().T@cs[0]
    outside=np.zeros((4,4),complex);outside[2:,2:]=[[.7,.2],[.2,1.1]]
    H0=q(outside)+kap*Y
    norm=float(np.max(abs(np.linalg.eigvalsh(Q))))
    assert abs(norm-math.sqrt(nu))<1e-12
    def expectation(theta,state):
        ev,V=np.linalg.eigh(H0+theta*kap*Q)
        out=V@(np.exp(-1j*ev)*(V.conj().T@state))
        return float(np.vdot(out,effect@out).real)
    errs=[];finite_shifts=[];eps=1e-5
    for outside_bits in range(4):
        state=np.zeros(dim,complex);state[1+(outside_bits<<2)]=1
        derivative=(expectation(eps,state)-expectation(-eps,state))/(2*eps)
        errs.append(abs(derivative+kap*mu))
        finite_shifts.append(.5-expectation(1/math.pi,state))
    assert max(errs)<1e-8
    lower=mu/4-nu/8
    assert min(finite_shifts)>lower-1e-12
    return dict(covariance='arbitrary original sea in analytic proof; all four outside occupations in diagnostic',
                fock_dimension=16,quadratic_norm=norm,expected_norm=math.sqrt(nu),
                derivative_max_error=max(errs),finite_shifts=finite_shifts,
                fixed_mu_nu_lower=lower,
                scope='finite full-CAR diagnostic, not actual BT evolution')

def budget():
    # Normal Z gain h=1+s, E h=1, |h|<=11/10, x=theta*kappa=1/4.
    mu=F(2,5);x=F(1,4);h2=F(121,100)
    ideal=(x-2*x*x*h2)*mu
    actual_if_each_e_le_005=ideal-4*F(1,200)
    assert ideal==F(79,2000) and ideal>F(3,100)
    assert actual_if_each_e_le_005==F(39,2000)>F(1,100)
    return dict(mu_safe=str(mu),x=str(x),normal_gain_second_moment_upper=str(h2),
                ideal_probability_gap_lower=str(ideal),ideal_probability_gap_float=float(ideal),
                per_program_isometry_error_new_requirement='1/200 (not proved here)',
                conditional_actual_gap_lower=str(actual_if_each_e_le_005),
                actual_BT_connection='UNPROVED; do not inherit old 1137 error for new program')

def run():
    return dict(status='PASS_IDEAL_PREFLIGHT_ONLY',numbered_round=None,
                source_hashes=frozen(),tail=tail_check(),full_CAR=car_check(),budget=budget(),
                new_inputs=['relational radial mask 1.1 to 1.25','explicit off/source/material program',
                            'normal Z gain 1+s(Z)','program mass and actual BT output comparison still to qualify'],
                scientific_calibration_increment=0)
def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
    else: assert a==b,(path,a,b)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run();target=Path(__file__).with_suffix('.json')
    if args.write_results:
        if target.exists():raise FileExistsError(target)
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif target.exists():compare(result,json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps(result,ensure_ascii=False,indent=2))
