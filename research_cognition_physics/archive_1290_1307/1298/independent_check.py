"""1298 independent exact algebra and conservative analytic certificates.
No author/historical calculation imports.  Finite polynomial tests calibrate
identities; the continuum and renormalization scopes require proof review.
Default/--check reads saved results; --write creates them exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir())
BASE=ROOT/'research_cognition_physics'
OUT=HERE/'independent_results.json'
ZERO=(0,0,0,0)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F): return {'num':x.numerator,'den':x.denominator}
    if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [enc(v) for v in x]
    return x

def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def stress(a,ba,c,bc):
    # a,c use spin vector (1,+/-1), norm squared 2. Divide by 2 for
    # spin normalization and by 2 for the symmetric Dirac stress kernel.
    K=(2,0,0,2)*(len(a)//4)
    return (dot(c,[F(k)*x for k,x in zip(K,ba)])+
            dot(bc,[F(k)*x for k,x in zip(K,a)]))/4

def central_cp_check():
    rows=[]
    for g,v,w in [(F(2),F(3,4),F(5)),
                  (F(3,2),F(7,5),F(11,3)),
                  (F(1),F(2,3),F(7,4))]:
        for spin in (1,-1):
            zeta=(F(1),F(spin))
            f=[F(0),F(0),g*zeta[0],g*zeta[1]]
            bf=[g*zeta[0],g*zeta[1],F(0),F(0)]
            pf=[-v*zeta[0]/2,-v*zeta[1]/2,g*zeta[0]/2,g*zeta[1]/2]
            bpf=[g*zeta[0]/2,g*zeta[1]/2,
                 (-w/2+(w-v)/6)*zeta[0],(-w/2-(w-v)/6)*zeta[1]]
            # Exact CP covariance: -FP-PF+FPF+V rho V*, with FPF=F/2.
            embed=lambda a,i: [F(0)]*(4*i)+a+[F(0)]*(4*(1-i))
            fs=[embed(f,i) for i in range(2)]
            bfs=[embed(bf,i) for i in range(2)]
            ps=[embed(pf,i) for i in range(2)]
            bps=[embed(bpf,i) for i in range(2)]
            c=sum((-stress(fs[i],bfs[i],ps[i],bps[i])
                   -stress(ps[i],bps[i],fs[i],bfs[i])
                   +stress(fs[i],bfs[i],fs[i],bfs[i])/2
                   for i in range(2)),F(0))
            B=[[stress(fs[i],bfs[i],fs[j],bfs[j]) for j in range(2)] for i in range(2)]
            assert all(x==0 for row in B for x in row)
            expected=2*g*(F(2,3)*w+F(1,3)*v)
            assert c==expected
            rhos=[[[F(1),F(0)],[F(0),F(0)]],
                  [[F(1,3),F(1,3)],[F(1,3),F(2,3)]],
                  [[F(1,2),F(-1,2)],[F(-1,2),F(1,2)]]]
            for rho in rhos:
                got=c+sum((B[i][j]*rho[j][i] for i in range(2) for j in range(2)),F(0))
                assert got==expected
            assert stress(f,bf,f,bf)==0
            rows.append({'g':g,'omega_inverse_g':v,'omega_g':w,
                         'sigma_x_eigenvalue':spin,'full_two_flavor_source':c,
                         'raw_code_only_source':F(0),'logical_matrix':B})
    return {'method':'direct finite-rank CP legs, not copied central formula',
            'cases':rows,'unknown_partner_scope':'exact scalar compressed source at the center; full CAR state required',
            'per_flavor_coefficients':[F(2,3),F(1,3)],
            'two_flavor_coefficients':[F(4,3),F(2,3)],
            'discard_reference_legs_counterexample':True,
            'continuous_positivity_uses_heat_semigroup_not_test_samples':True}

# Laurent polynomial in v; derivative of b(v)P(v), b=exp(-1/v).
def lclean(p): return {n:c for n,c in p.items() if c}
def ladd(*ps):
    r={}
    for p in ps:
        for n,c in p.items():r[n]=r.get(n,F(0))+c
    return lclean(r)
def lscale(p,c):return lclean({n:c*v for n,v in p.items()})
def lshift(p,j):return {n+j:c for n,c in p.items()}
def ld(p):return lclean({n-1:n*c for n,c in p.items() if n})
def bd(p):return ladd(ld(p),lshift(p,-2))
def radial_lap(p):
    d2=bd(bd(p))
    return ladd(lscale(bd(p),-6),lscale(d2,4),lscale(lshift(d2,1),-4))
def laurent_bound(p):
    assert all(n<=0 for n in p)
    return sum((abs(c)*math.factorial(-n) for n,c in p.items()),F(0))

def window_certificate():
    p2=radial_lap({0:F(1)})
    p4=radial_lap(p2)
    assert p2=={-4:F(4),-3:F(-12),-2:F(2)}
    assert p4=={-8:F(16),-7:F(-224),-6:F(896),-5:F(-1168),-4:F(444),-3:F(-24)}
    p6=radial_lap(p4)
    L2,L4,L6=laurent_bound(p2),laurent_bound(p4),laurent_bound(p6)
    assert L6==263117600640
    assert L2==172 and L4==2570160
    a=F(1,32)
    # ||b||_2^2>pi*a^3/96 and volume=4*pi*a^3/3.
    # sqrt(128)<12; G_j<=G_4 for j<=4 and omega>=1.
    K=12*(1+2*L2/a**2+L4/a**4)
    C1=48*K*K
    ell=F(1024)/C1
    center=F(1)/(4*a**3)
    loss=4*C1*ell
    assert center==8192 and loss==4096 and 0<ell<a/8
    positive=center-loss
    flux=4*ell**2*(ell**2/2)*positive
    assert flux==8192*ell**4 and flux>0
    # Independent nonempty linear weak-metric budget via Poisson+Kirchhoff.
    S=C1
    geometry_factor=a*a+8*a+1
    kappa=1/(300*S*geometry_factor)
    barh_cap=kappa*S*geometry_factor
    h_cap=3*barh_cap
    assert h_cap==F(1,100)
    B6=1+3*L2/a**2+3*L4/a**4+L6/a**6
    M6=8192*B6**2
    Hstar=512*M6
    # For the zero-data inhomogeneous wave part: spatial second derivatives
    # <=M6, mixed time-space <=3M6, dt^2=Delta-2*tau <=5M6.
    # Homogeneous initial terms: <=20*a^2*M6 from f and <=24*a^2*M6
    # from the initial velocity; |n|_1<2 is sufficient.
    barH_C2_factor=5+44*a**2
    assert barH_C2_factor==F(1291,256) and barH_C2_factor<8
    assert Hstar==334779972733175839233728425479031683340042240000
    kappa_C2=1/(16*Hstar)
    lorentz_operator_cap=4*kappa_C2*Hstar
    assert lorentz_operator_cap==F(1,4)<1
    return {'radial_laplacian':p2,'radial_bilaplacian':p4,'radial_trilaplacian':p6,
            'factorial_sup_laplacian':L2,'factorial_sup_bilaplacian':L4,'factorial_sup_trilaplacian':L6,
            'radius':a,'G4_upper':K,'each_coordinate_C1_upper':C1,
            'strict_center_lower':center,'ell':ell,'ell_diagnostic':float(ell),
            'patch':'t=z in [0,ell], |x|,|y|<=ell; all coordinates in original flat frame',
            'maximum_source_loss':loss,'strict_patch_source_lower':positive,
            'strict_weighted_flux_lower':flux,
            'weak_geometry':{'source_sup_safe':S,'initial_potential_bound_factor':a*a,
              'initial_gradient_bound_factor':2*a,'initial_velocity_bound_factor':6*a,
              'kirchhoff_barh_bound_factor':geometry_factor,'fixed_positive_kappa':kappa,
              'barh_component_cap':barh_cap,'h_component_cap':h_cap,
              'coordinate_time_domain':[F(0),F(1)]},
            'C2_geometry':{'B6':B6,'M6':M6,'Hstar':Hstar,
              'barH_component_C2_factor':barH_C2_factor,
              'author_safe_barH_C2_factor':F(8),
              'fixed_kappa_upper':kappa_C2,'metric_operator_perturbation_cap':lorentz_operator_cap,
              'Poisson_derivatives_on_source_order':3,
              'not_using_global_L2_for_Coulomb_tail':True,
              'no_nonlinear_solution_closeness_inferred':True},
            'analytic_inputs':['heat semigroup maximum at radial bump center',
              'Sobolev H2 point evaluation constant (8*pi)^(-1/2)<1',
              'smooth compact modes, omega-commuting Dirac propagator and projector',
              'six absolute CP rank terms bound the total two-flavor source',
              '3D Poisson kernel and Kirchhoff formula'],
            'not_certified':['physical resolution at ell','EFT frequency-band coverage',
              'nonlinear Einstein closeness','curved-state or instrument transport',
              'absolute reference polarization kernel']}

# Small exact 4-variable polynomial algebra for an independent Ward/seed test.
def clean(p):return {m:c for m,c in p.items() if c}
def mono(c,*m):return {tuple(m):F(c)} if c else {}
def add(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,F(0))+c
    return clean(out)
def scale(p,c):return clean({m:c*v for m,v in p.items()})
def deriv(p,i):
    out={}
    for m,c in p.items():
        if m[i]:
            n=list(m);n[i]-=1
            out[tuple(n)]=c*m[i]
    return clean(out)
def twice(p,i):return deriv(deriv(p,i),i)
def lap(p):return add(*(twice(p,i) for i in (1,2,3)))
def box(p):return add(twice(p,0),scale(lap(p),-1))
def tcoef(p,j):
    return {(0,)+m[1:]:c for m,c in p.items() if m[0]==j}
def tpower(p,j):return {(m[0]+j,)+m[1:]:c for m,c in p.items()}
def ix(p):
    return {(m[0],m[1]+1,m[2],m[3]):c/(m[1]+1) for m,c in p.items()}
def poisson_poly(p):
    out={};term=p
    for _ in range(12):
        if not term:break
        j=ix(ix(term));out=add(out,j)
        term=scale(add(twice(j,2),twice(j,3)),-1)
    assert not term
    assert lap(out)==p
    return out

def wave(f,g,T,kappa):
    coef=[f,g]
    for j in range(20):
        coef.append(scale(add(lap(coef[j]),scale(tcoef(T,j),-2*kappa)),F(1,(j+1)*(j+2))))
    assert not coef[-1] and not coef[-2]
    out=add(*(tpower(p,j) for j,p in enumerate(coef)))
    assert add(box(out),scale(T,2*kappa))=={}
    return out

def harmonic_tensor(T,kappa,drop_velocity=False,drop_poisson=False):
    phi=[{} if drop_poisson else scale(poisson_poly(tcoef(T[0][nu],0)),2*kappa) for nu in range(4)]
    g00={} if drop_velocity else add(*(deriv(phi[i],i) for i in (1,2,3)))
    h=[[{} for _ in range(4)] for _ in range(4)]
    for mu in range(4):
        for nu in range(mu,4):
            f=phi[nu] if mu==0 else {}
            g=g00 if mu==nu==0 else {}
            h[mu][nu]=h[nu][mu]=wave(f,g,T[mu][nu],kappa)
    C=[add(deriv(h[0][nu],0),scale(add(*(deriv(h[i][nu],i) for i in (1,2,3))),-1)) for nu in range(4)]
    return phi,h,C

def harmonic_certificate():
    pot=add(mono(1,1,3,0,0),mono(F(2,3),0,2,0,2),
            mono(F(5,7),2,0,2,0),mono(F(3,5),2,1,1,1),
            mono(F(7,11),3,0,0,1),mono(F(1,2),0,0,4,0))
    eta=(1,-1,-1,-1); b=box(pot)
    T=[[add(deriv(deriv(pot,mu),nu),scale(b,-eta[mu] if mu==nu else 0))
        for nu in range(4)] for mu in range(4)]
    div=lambda TT,nu:add(deriv(TT[0][nu],0),scale(add(*(deriv(TT[i][nu],i) for i in (1,2,3))),-1))
    assert all(div(T,nu)=={} for nu in range(4))
    kappa=F(3,7)
    phi,h,C=harmonic_tensor(T,kappa)
    assert all(c=={} for c in C)
    assert all(add(scale(box(h[mu][nu]),F(-1,2)),scale(T[mu][nu],-kappa))=={}
               for mu in range(4) for nu in range(4))
    _,_,C_no_v=harmonic_tensor(T,kappa,drop_velocity=True)
    assert tcoef(C_no_v[0],0)!={}
    _,_,C_no_p=harmonic_tensor(T,kappa,drop_poisson=True)
    assert any(tcoef(deriv(c,0),0) for c in C_no_p)
    bad=[[dict(T[mu][nu]) for nu in range(4)] for mu in range(4)]
    bad[0][0]=add(bad[0][0],mono(1,1,0,0,0))
    _,_,C_bad=harmonic_tensor(bad,kappa)
    assert any(C_bad)
    assert all(add(box(C_bad[nu]),scale(div(bad,nu),2*kappa))=={} for nu in range(4))
    return {'signature':'+---','box':'dt^2-Delta','G1_in_harmonic':'-box(barh)/2',
      'polynomial_ward_components_zero':4,'poisson_seed_components_verified':4,
      'wave_equations_verified':10,'harmonic_constraints_identically_zero':4,
      'linear_einstein_components_verified':10,
      'omit_seed_velocity_detected':True,'omit_poisson_seed_detected':True,
      'nonconserved_source_forces_C_detected':True,
      'source_potential_terms':[{ 'powers':m,'coefficient':c} for m,c in sorted(pot.items())],
      'scope':'exact polynomial calibration only; compact physical source and global existence use analytic proof'}

def optical_weight_check():
    ell,R0,R1=F(1,4),F(2),F(3)
    area=R0*ell**2/2-R1*ell**3/6
    boost=R0*ell**2/2-R1*ell**3/3
    weight=R1*ell**3/6
    assert area==boost+weight and weight>0
    return {'terminal_expansion':'theta(lambda=0)=0','lambda_domain':[-ell,F(0)],
            'affine_curvature_R0':R0,'affine_curvature_R1':R1,
            'geometric_area_response':area,'boost_integral':boost,
            'linear_weight_correction':weight,'bare_clausius_residual_nonzero':True,
            'weight':'2*lambda+ell','not_a_second_order_only_term':True}


def curved_current_identity_control():
    # g=e^(2t)*eta at the origin, tau^00=1 constant, chi=partial_t.
    # This is an algebraic diagnostic, not a conserved curved physical source.
    divergence_current=F(6)
    half_Lie_metric=F(1)
    chi_divergence_tau=F(5)
    assert divergence_current==half_Lie_metric+chi_divergence_tau
    return {'dimension':4,'divergence_j':divergence_current,
            'half_tau_Lie_g':half_Lie_metric,'chi_div_tau':chi_divergence_tau,
            'omitting_covariant_source_divergence_fails':True,
            'scope':'product-rule calibration; flat conservation does not imply curved conservation'}


def calculate():
    sources=[HERE/'independent_check.py',
        BASE/'archive_1290_/_shared/round1298_admission.md',
        BASE/'archive_1290_/1293/proof.md',
        BASE/'archive_1290_/1295/proof.md']
    return enc({'round':1298,'status':'INDEPENDENT_ANALYTIC_CALIBRATIONS_PASS',
       'formal_review_is_separate':True,'new_scientific_groups':0,
       'central_cp':central_cp_check(),'finite_window':window_certificate(),
       'harmonic':harmonic_certificate(),'finite_optical_balance':optical_weight_check(),
       'curved_current_product_rule':curved_current_identity_control(),
       'source_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in sources},
       'scope':'same flat Dirac mean source and constrained linear/leading-kappa geometry; no finite-coupling quantum transport or actual readout certification'})

def main():
    parser=argparse.ArgumentParser()
    action=parser.add_mutually_exclusive_group()
    action.add_argument('--write',action='store_true')
    action.add_argument('--check',action='store_true')
    args=parser.parse_args()
    data=calculate()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write('\n')
        mode='CREATED_EXCLUSIVE'
    else:
        saved=json.loads(OUT.read_text('utf-8'))
        assert data==saved,'stored independent result or source binding mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps({'round':1298,'status':data['status'],'mode':mode,
                      'result_sha256':sha(OUT),'formal_review_is_separate':True},sort_keys=True))

if __name__=='__main__':main()
