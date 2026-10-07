"""897: absolute-source jet residual certificate, original128 matrix calibration.
The continuum statement is a conditional Cauchy-kernel/diagonal-trace theorem.
This finite symbol check does NOT evaluate a curved renormalized vacuum source.
"""
from pathlib import Path
import json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'896'))
import compatible_instrument_source_bridge as bridge
old=bridge.old
TARGET=HERE/'absolute_source_residual_certificate_results.json'
INDICES=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))

def within(b,a):return b[0]<=a[0] and b[1]<=a[1]
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def fro(x):return float(np.linalg.norm(x,'fro'))
def value(poly,t):return sum(v*t**n for n,v in poly.items())
def add_int(dst,src,factor):
    for n,c in src.items():dst[n+1]=dst.get(n+1,0.)+factor*c/(n+1)

def jet_h(K,M):
    out={a:np.zeros_like(K) for a in INDICES}
    out[(0,0)]=K+M
    for a in INDICES[1:]:
        if a[1]==0:out[a]=(-1.)**a[0]*K/math.factorial(a[0])
        if a[0]==0:out[a]=M/math.factorial(a[1])
    return out

def generator(H,A):
    ans={a:np.zeros_like(H[(0,0)]) for a in INDICES}
    for a in INDICES:
        for b in INDICES:
            if within(b,a):
                X=H[b];Y=A[sub(a,b)];ans[a]+=-1j*(X@Y-Y@X)
    return ans

def propagate_coefficients(H,N,order):
    v={a:np.zeros_like(N) for a in INDICES};v[(0,0)]=N
    powers=[v]
    for _ in range(order+1):powers.append(generator(H,powers[-1]))
    return powers

def evaluate(powers,order,t):
    return {a:sum(t**j/math.factorial(j)*powers[j][a] for j in range(order+1)) for a in INDICES}

def residual_bounds(H,powers,order,t):
    # Normalized parameter Taylor coefficients, so no binomial coefficients.
    norms={b:2*float(np.linalg.norm(H[b],2)) for b in INDICES[1:]}
    polys={}
    for a in INDICES:
        polys[a]={order+1:fro(powers[order+1][a])/math.factorial(order+1)}
        for b in INDICES[1:]:
            if within(b,a):add_int(polys[a],polys[sub(a,b)],norms[b])
    return {a:value(polys[a],t) for a in INDICES}

def vertex_jets(K,M,H):
    z=np.zeros_like(K);out={'energy':H}
    out['conformal']={a:(-1.)**(a[0]+1)*K/math.factorial(a[0]) if a[1]==0 else z for a in INDICES}
    out['mass_scale']={a:M/math.factorial(a[1]) if a[0]==0 else z for a in INDICES}
    return out

def source_jet(V,N,a):
    return sum(old.source(N[sub(a,b)],V[b]) for b in INDICES if within(b,a))
def source_bound(V,bounds,a):
    return .5*sum(fro(V[b])*bounds[sub(a,b)] for b in INDICES if within(b,a))

def exact_covariance(K,M,N,gamma,xi,t):
    h=np.exp(-gamma)*K+np.exp(xi)*M
    vals,v=np.linalg.eigh(h);U=(v*np.exp(-1j*t*vals))@v.conj().T
    return U@N@U.conj().T

def run():
    bg=old.background(.04);zero=[np.zeros_like(x) for x in bg]
    B,vertices=bridge.source_vertices(bg)
    mass64=old.mass_jet(bg[1],np.zeros(5))[0];M=old.assemble([mass64,mass64]);K=B-M
    H=jet_h(K,M)
    # Verify the full nonlinear original mass-coordinate/conformal family.
    checks=[]
    phi=bg[1];F=old.old.old.matter.original.F(phi);x=phi/np.sqrt(F)
    for gamma,xi in ((.03,.02),(-.02,.04)):
        xx=np.exp(xi)*x;pp=xx*np.sqrt(2/(1+np.dot(xx,xx)/6))
        family=(np.exp(2*gamma)*bg[0],pp,bg[2],bg[3])
        original=old.matrix_and_source(*family,*zero)[0]
        checks.append(float(np.max(abs(original-(np.exp(-gamma)*K+np.exp(xi)*M)))))
    assert max(checks)<3e-13
    N0,dN,_,_,C=old.initial_modes();N0=N0+dN
    # N denotes occupation, in accordance with730/804. Never replace by 1-N.
    val,q=np.linalg.eigh(B);S=q[:,val<0]@q[:,val<0].conj().T
    T=.35;orders=(1,2,4,6);reference_order=20
    powers=propagate_coefficients(H,N0,reference_order)
    reference=evaluate(powers,reference_order,T)
    bref=residual_bounds(H,powers,reference_order,T)
    Nexact=exact_covariance(K,M,N0,0.,0.,T)
    assert fro(reference[(0,0)]-Nexact)<2e-12
    exact_spectrum=np.linalg.eigvalsh(Nexact)
    assert exact_spectrum.min()>-3e-12 and exact_spectrum.max()<1+3e-12
    V=vertex_jets(K,M,H)
    # Independent finite difference checks of first, pure and mixed second jets.
    h=.0008;fd=[]
    for a in ((1,0),(0,1),(2,0),(0,2),(1,1)):
        if a==(1,0):estimate=(exact_covariance(K,M,N0,h,0,T)-exact_covariance(K,M,N0,-h,0,T))/(2*h)
        elif a==(0,1):estimate=(exact_covariance(K,M,N0,0,h,T)-exact_covariance(K,M,N0,0,-h,T))/(2*h)
        elif a==(2,0):estimate=(exact_covariance(K,M,N0,h,0,T)-2*Nexact+exact_covariance(K,M,N0,-h,0,T))/(2*h*h)
        elif a==(0,2):estimate=(exact_covariance(K,M,N0,0,h,T)-2*Nexact+exact_covariance(K,M,N0,0,-h,T))/(2*h*h)
        else:estimate=(exact_covariance(K,M,N0,h,h,T)-exact_covariance(K,M,N0,h,-h,T)-exact_covariance(K,M,N0,-h,h,T)+exact_covariance(K,M,N0,-h,-h,T))/(4*h*h)
        error=fro(estimate-reference[a]);assert error<2e-6
        fd.append(dict(normalized_jet=list(a),finite_difference_error=error))
    rows=[]
    pairs=bridge.gauge_pairs(bg,B)
    for n in orders:
        approx=evaluate(powers,n,T);bounds=residual_bounds(H,powers,n,T)
        errors={a:fro(approx[a]-reference[a]) for a in INDICES}
        assert all(errors[a]<=bounds[a]+bref[a]+3e-12 for a in INDICES)
        records=[]
        for name,v in V.items():
            for a in INDICES:
                err=abs(source_jet(v,approx,a)-source_jet(v,reference,a))
                bd=source_bound(v,bounds,a)+source_bound(v,bref,a)
                assert err<=bd+2e-11
                records.append(dict(source=name,normalized_jet=list(a),error=err,bound=bd))
        difference=approx[(0,0)]-reference[(0,0)]
        all_errors=[abs(old.source(difference,v)) for _,v in vertices]
        assert all(err<=.5*fro(v)*(bounds[(0,0)]+bref[(0,0)])+1e-12 for (_,v),err in zip(vertices,all_errors))
        dot=sum(T**(j-1)/math.factorial(j-1)*powers[j][(0,0)] for j in range(1,n+1))
        residual=dot+1j*(B@approx[(0,0)]-approx[(0,0)]@B)
        explicit=-T**n/math.factorial(n)*powers[n+1][(0,0)]
        assert fro(residual-explicit)<3e-13
        ward_identity=max(abs(old.source(dot,Q)+old.source(approx[(0,0)],D)-old.source(residual,Q)) for Q,D in pairs)
        defect=max(abs(old.source(dot,Q)+old.source(approx[(0,0)],D)) for Q,D in pairs)
        assert ward_identity<3e-13
        assert defect<=max(.5*fro(Q)*fro(residual) for Q,_ in pairs)+1e-13
        ev=np.linalg.eigvalsh(approx[(0,0)])
        rows.append(dict(time_polynomial_order=n,
            kernel_jet_errors={str(a):errors[a] for a in INDICES},
            kernel_jet_residual_bounds={str(a):bounds[a] for a in INDICES},
            joint_source_jet_errors=records,
            maximum_other_source_error=max(all_errors),
            instantaneous_kernel_residual=fro(residual),
            actual_gauge_exchange_defect=defect,off_shell_balance_identity_error=ward_identity,
            numerical_covariance_min=float(ev.min()),numerical_covariance_max=float(ev.max())))
    # These are two forbidden substitutions, not a second numbered experiment.
    contact=abs(source_jet(V['conformal'],reference,(1,0))-old.source(reference[(1,0)],V['conformal'][(0,0)]))
    vacuum_omission=abs(old.source(S,B))
    assert contact>.01 and vacuum_omission>1
    assert rows[-1]['kernel_jet_errors']['(0, 0)']<rows[0]['kernel_jet_errors']['(0, 0)']
    return dict(round=897,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3682,
        argument_scope='Conditional full-kernel residual certificate for the fixed renormalized one-fermion-loop absolute source and its first two background jets. Calibration retains original128 local symbol, full occupation reference and source contacts; it is not a numerical continuum vacuum or original-Q/E matching proof.',
        numerical_Nambu_dimension=128,numerical_physical_modes=64,time=T,
        actual_original_parameter_family_errors=checks,independent_parameter_jet_checks=fd,
        certified_reference_order=reference_order,
        reference_jet_bounds={str(a):bref[a] for a in INDICES},rows=rows,
        omitted_conformal_Hessian_contact_error=contact,
        omitted_reference_energy_error=vacuum_omission,
        finite_subtraction_reference_is_not_continuum_parametrix=True,
        numerical_covariance_not_declared_a_physical_state=True,
        physical_target_state_is_original_positive_instrument_state=True,
        actual_curved_absolute_source_computed=False,
        full_common_model_Q_E_error_certificate_completed=False,full_goal_completed=False,
        floating_point_not_interval_arithmetic=True)

if __name__=='__main__':
    result=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},ensure_ascii=False,indent=2))
    for r in result['rows']:print(r['time_polynomial_order'],r['kernel_jet_errors'],r['kernel_jet_residual_bounds'],r['actual_gauge_exchange_defect'])
