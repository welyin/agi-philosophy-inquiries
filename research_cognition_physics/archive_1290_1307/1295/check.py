"""1295 author group: exact finite constants and joint matrix/CAR diagnostics.
Default/--check only read; --write exclusively creates results.json.
No continuum discretization, inherited checker import or apparatus claim.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'
ALPHA=F(31,2**27)
EPS=ALPHA/128
TOL=3e-11


def encode(v):
    if isinstance(v,F): return f'{v.numerator}/{v.denominator}'
    if isinstance(v,dict): return {str(k):encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [encode(x) for x in v]
    if isinstance(v,np.generic): return v.item()
    return v


def residual(a,b=0):
    x=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    assert x<TOL, x
    return x


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def rank(a):
    a=[[F(x) for x in row] for row in a]
    r=0
    for c in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        z=a[r][c]; a[r]=[x/z for x in a[r]]
        for i in range(r+1,len(a)):
            if a[i][c]:
                z=a[i][c]; a[i]=[x-z*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r


def commutator_map(a,b):
    n=len(a); rows=[]
    for h in (a,b):
        for i in range(n):
            for j in range(n):
                rows.append([(h[q][j] if p==i else F(0))-
                             (h[i][p] if q==j else F(0))
                             for p in range(n) for q in range(n)])
    return rows


def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def poly_rem(a,b):
    a=[F(x) for x in a]; b=trim([F(x) for x in b])
    while len(trim(a))>=len(b) and any(a):
        q=a[-1]/b[-1]; k=len(a)-len(b)
        for i in range(len(b)):a[i+k]-=q*b[i]
    return trim(a)


def jacobi_polynomial(n):
    prev=[F(1)]; cur=[F(-1),F(1)]
    for j in range(2,n+1):
        nxt=[F(0)]*(len(cur)+1)
        for i,c in enumerate(cur):nxt[i]-=j*c; nxt[i+1]+=c
        for i,c in enumerate(prev):nxt[i]-=c
        prev,cur=cur,trim(nxt)
    return cur


def exact_matrices(n):
    a=[[F(i+1) if i==j else F(int(abs(i-j)==1)) for j in range(n)] for i in range(n)]
    a2=mm(a,a)
    hn=[[a[i][j]+EPS*a2[i][j]/2 for j in range(n)] for i in range(n)]
    he=[[F(i+1)+EPS*F((i+1)**2,4) if i==j else F(0) for j in range(n)] for i in range(n)]
    r=rank(commutator_map(hn,he))
    assert r==n*n-1
    hd=[[F(i+1)+EPS*F((i+1)**2,2) if i==j else F(0) for j in range(n)] for i in range(n)]
    rd=rank(commutator_map(hd,he)); assert rd==n*n-n
    edge=[2*EPS+EPS**2*(2*i+1) for i in range(1,n)]
    assert all(x>0 for x in edge)
    # Exact distinctness of A's spectrum is separately certified by charpoly gcd.
    poly=jacobi_polynomial(n); deriv=[i*poly[i] for i in range(1,len(poly))]
    x,y=poly,deriv
    while any(y):x,y=y,poly_rem(x,y)
    assert len(trim(x))==1
    assert 1-EPS>0 and 1+(n+2)*EPS<2
    row=max(sum(abs(v) for v in rr) for rr in a)
    assert row<=n+1<n+2
    numeric=np.array(hn,float)
    ev=np.linalg.eigvalsh(numeric)
    assert min(np.diff(ev))>0
    return {'n':n,'commutator_rank':r,'commutant_complex_dimension':n*n-r,
            'delete_path_commutant_dimension':n*n-rd,'scaled_nu_matrix':hn,
            'scaled_e_matrix':he,'adjacent_Hnu_entries':edge,
            'A_characteristic_polynomial_ascending':poly,'charpoly_derivative_gcd_degree':0,
            'A_row_sum_bound':row,'nu_mass_norm_upper':1+(n+2)*EPS,
            'scaled_Hnu_min_eigen_gap_diagnostic':float(min(np.diff(ev))),
            'exact_rank_uses_rescaled_H_not_tiny_float_threshold':True}


def exact_budget():
    d=3*EPS
    assert ALPHA**2-2*ALPHA*d==F(61,64)*ALPHA**2
    assert ALPHA-d==F(125,128)*ALPHA
    assert 4*(F(512)*32**2+4)<1450**2
    return {'alpha_star':ALPHA,'epsilon':EPS,'center_d_upper':d,
            'center_probability_strict_lower':F(61,64)*ALPHA**2,
            'center_amplitude_smin_strict_lower':F(125,128)*ALPHA,
            'open_probability_strict_lower':ALPHA**2/2,
            'open_amplitude_smin_strict_lower':3*ALPHA/4,
            'time_window':[F(7,16),F(1,2)],'source_radius':F(1,32),
            'center_distance':F(1,4),'gradient_squared_strict_upper':524288,
            'energy_square_strict_upper':4*524292,'energy_strict_upper':1450,
            'open_d_bound_only_for_neutrino_block':True,
            'each_model_requires_own_Hadamard_reference':True}


SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)


def matrices(n):
    d=np.diag(np.arange(1,n+1,dtype=float))
    a=d+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
    return np.eye(n)+float(EPS)*a,2*np.eye(n)+float(EPS)*d,a


def dirac(m,p):
    n=len(m); sp=p[0]*SX+p[1]*SY+p[2]*SZ
    k=np.kron(sp,np.eye(n)); mass=np.kron(np.eye(2),m)
    return np.block([[-k,mass],[mass.conj().T,k]])


def evolve(h,t):
    vals,vec=np.linalg.eigh(h)
    return (vec*np.exp(-1j*t*vals))@vec.conj().T


def spectral_energy_and_evolution(n):
    m,_,a=matrices(n); maxerr=0.; samples=[]
    for p in ((0,0,0),(F(1,3),F(-1,5),F(2,7))):
        p=np.asarray(p,float); h=dirac(m,p); h0=dirac(np.eye(n),p)
        for t in (float(F(7,16)),.5):
            du=float(np.linalg.norm(evolve(h,t)-evolve(h0,t),2))
            bound=t*float(EPS)*float(np.linalg.norm(a,2))
            assert du<=bound+2e-14
            samples.append({'p':p.tolist(),'t':t,'norm_difference':du,'Duhamel_bound':bound})
        hzero=dirac(np.zeros((n,n)),p)
        maxerr=max(maxerr,residual(evolve(hzero,.5)[:2*n,2*n:]))
    # Independent finite +-momentum analogue of the radial-mode identity.
    p=np.array([.3,-.2,.4]); hp=dirac(m,p); hm=dirac(m,-p)
    z=np.zeros_like(hp); h=np.block([[hp,z],[z,hm]])
    vals,vec=np.linalg.eigh(h)
    pref=(vec*(vals<0))@vec.conj().T
    absolute=(vec*np.abs(vals))@vec.conj().T
    v=np.zeros((len(h),2),complex)
    for j in range(2):v[2*n+j,j]=1/math.sqrt(2);v[len(hp)+2*n+j,j]=1/math.sqrt(2)
    ff=v@v.conj().T
    maxerr=max(maxerr,residual(v.conj().T@h@v))
    expected=float(np.trace(ff@absolute).real)
    energies=[]
    for rho in (np.diag([1.,0.]),np.array([[.5,.3j],[-.3j,.5]])):
        dc=-ff@pref-pref@ff+ff@pref@ff+v@rho@v.conj().T
        energy=np.trace(h@dc)
        maxerr=max(maxerr,residual(energy,expected));energies.append(float(energy.real))
    assert expected>0
    return {'n':n,'finite_symbol_Duhamel_diagnostics':samples,'max_identity_residual':maxerr,
            'relative_energies_finite_momentum_analogue':energies,
            'trace_F_absolute_h':expected,'not_continuum_propagation_or_energy_estimate':True}


def ward_sources(n):
    yn,ye,_=matrices(n); zz=np.zeros((n,n),complex)
    m=np.block([[yn,zz],[zz,ye]]); zero=np.zeros_like(m)
    hh=np.block([[zero,m],[m.conj().T,zero]])
    def scalar_variation(dhp,dh0):
        dm=np.block([[yn*np.conj(dh0),ye*dhp],[-yn*np.conj(dhp),ye*dh0]])
        return np.block([[zero,dm],[dm.conj().T,zero]])
    data=[]; maxerr=0.; omitted=[]
    generators=[('W1',SX/2),('W2',SY/2),('W3',SZ/2),('Y',np.eye(2)/2)]
    for name,th in generators:
        dh=1j*th@np.array([0,1.],complex)
        if name=='Y':
            tl=-.5*np.eye(2*n);tr=np.kron(np.diag([0.,-1.]),np.eye(n))
        else:tl=np.kron(th,np.eye(n));tr=np.zeros_like(tl)
        tg=np.block([[tl,zero],[zero,tr]])
        dj=scalar_variation(dh[0],dh[1])
        comm=1j*(hh@tg-tg@hh)
        er=residual(comm+dj)
        maxerr=max(maxerr,er); omitted.append(float(np.linalg.norm(comm,2)))
        data.append({'generator':name,'Ward_residual':er,'drop_Higgs_defect':omitted[-1]})
    assert min(omitted)>.4
    # Each of the four real Higgs source directions is retained as Hermitian.
    for hp,h0 in ((1,0),(1j,0),(0,1),(0,1j)):
        j=scalar_variation(hp,h0);maxerr=max(maxerr,residual(j,j.conj().T))
    return {'n':n,'all_four_real_Higgs_sources_retained':True,'Ward_checks':data,
            'max_identity_residual':maxerr,'spurion_Y_derivatives_in_menu':False,
            'diagnostic_is_mass_current_Higgs_Ward_not_spacetime_PDE':True}


def annihilation(i,n):
    a=np.zeros((2**n,2**n),complex);bit=1<<(n-1-i)
    for c in range(2**n):
        if c&bit:a[c^bit,c]=(-1)**((c>>(n-i)).bit_count())
    return a


def trace_out(rho,dims,keep):
    a=rho.reshape(tuple(dims)*2);cur=list(dims)
    for i in reversed(range(len(cur))):
        if i not in keep:a=np.trace(a,axis1=i,axis2=i+len(cur));cur.pop(i)
    d=math.prod(cur);return a.reshape(d,d)


def covariance(rho,ops):
    return np.array([[np.trace(rho@ops[j].conj().T@ops[i]) for j in range(len(ops))]
                     for i in range(len(ops))])


def dg(j,ops):
    return sum((j[a,b]*ops[a].conj().T@ops[b] for a in range(len(ops)) for b in range(len(ops))),
               np.zeros_like(ops[0]))


def car_and_matrix_sources():
    ops=[annihilation(i,4) for i in range(4)]
    # Correlated reference, not a product between the code modes and complement.
    href=np.array([[-2,.2,.7,.1j],[.2,-1,.3j,.6],[.7,-.3j,1,.4],[-.1j,.6,.4,2]],complex)
    ev,u=np.linalg.eigh(href); pv=u[:,ev<0]@u[:,ev<0].conj().T
    psi=np.zeros(16,complex);psi[0]=1
    for col in u[:,ev<0].T:
        creation=sum((col[i]*ops[i].conj().T for i in range(4)),np.zeros((16,16),complex))
        psi=creation@psi
    psi/=np.linalg.norm(psi);ref=np.outer(psi,psi.conj())
    er=residual(covariance(ref,ops),pv)
    complement=trace_out(ref,[4,4],[1])
    code_ref=trace_out(ref,[4,4],[0])
    defect=float(np.linalg.norm(ref-np.kron(code_ref,complement),2));assert defect>.05
    w=np.zeros((4,2),complex);w[2,0]=1;w[1,1]=1
    rho=np.array([[.5,.3j],[-.3j,.5]],complex)
    state=np.kron(w@rho@w.conj().T,complement)
    ff=np.diag([1.,1.,0.,0.]);qq=np.eye(4)-ff
    v=np.eye(4)[:,:2]
    expected=qq@pv@qq+v@rho@v.conj().T
    er=max(er,residual(covariance(state,ops),expected))
    n01=ops[0].conj().T@ops[0]@ops[1].conj().T@ops[1]
    double=float(np.trace(state@n01).real);assert abs(double)<TOL and np.linalg.det(rho).real>.1
    # Entangled logical input; final tensor ordering is code, complement, R.
    bell=np.array([1,0,0,1],complex)/math.sqrt(2)
    qr=np.outer(bell,bell.conj());wr=np.kron(w,np.eye(2))
    fr=wr@qr@wr.conj().T
    full=np.kron(fr,complement).reshape(4,2,4,4,2,4).transpose(0,2,1,3,5,4).reshape(32,32)
    er=max(er,residual(trace_out(full,[4,4,2],[1,2]),np.kron(complement,np.eye(2)/2)))
    er=max(er,residual(trace_out(full,[4,4,2],[0,2]),fr))
    onebody=np.array([[1,.2j,.7,0],[-.2j,-.4,.3,.1j],[.7,.3,.2,-.5],[0,-.1j,-.5,.8]],complex)
    src=dg(onebody,ops);bm=v.conj().T@onebody@v
    c=np.trace(onebody@qq@pv@qq)
    er=max(er,residual(np.trace(state@src),c+np.trace(rho@bm)))
    bR=SZ
    actual=np.trace(full@np.kron(src,bR))
    rhs=c*np.trace(qr@np.kron(np.eye(2),bR))+np.trace(qr@np.kron(bm,bR))
    er=max(er,residual(actual,rhs));assert abs(actual)>.1
    # Same identity with the full reference baseline explicitly retained.
    cdelta=np.trace(onebody@(-ff@pv-pv@ff+ff@pv@ff))
    er=max(er,residual(np.trace((state-ref)@src),cdelta+np.trace(rho@bm)))
    scalar_wrong=float(np.linalg.norm(bm-np.trace(bm)*np.eye(2)/2,2));assert scalar_wrong>.1
    return {'modes':4,'reference_factorization_defect':defect,'max_identity_residual':er,
            'nonquasifree_actual_double_occupation':double,
            'wrong_Wick_double_occupation':float(np.linalg.det(rho).real),
            'full_remote_R_product_verified':True,'full_code_R_state_preserved':True,
            'source_R_correlation_real':float(actual.real),'wrong_scalar_source_defect':scalar_wrong,
            'source_baseline_retained':True,'finite_algebra_diagnostic_not_continuum_model':True}


def bindings():
    items={'proof':'proof.md','check':'check.py','admission':'../_shared/round1295_admission.md',
           '1011':'../../archive_1009_1043/research_note_1011.md',
           '730':'../../archive_702_741/research_note_730.md',
           '803':'../../archive_797_818/research_note_803.md',
           '805':'../../archive_797_818/research_note_805.md',
           '813':'../../archive_797_818/research_note_813.md',
           '818':'../../archive_797_818/research_note_818.md',
           '1293proof':'../1293/proof.md','1293results':'../1293/results.json',
           'background_dictionary':'../_admission/background_dressing_screen_after1293.md'}
    return {key:{'path':rel,'sha256':hashlib.sha256((HERE/rel).read_bytes()).hexdigest()}
            for key,rel in items.items()}


def compare(a,b):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isfinite(a) and math.isfinite(b)
        assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-14),(a,b)
    else:assert a==b,(a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.check and args.write:parser.error('choose --check or --write')
    payload=encode({'round':1295,'author_check_groups':1,'verdict':'PASS_AUTHOR_SCOPED',
                    'exact_budget':exact_budget(),'Yukawa_families':[exact_matrices(n) for n in (3,4)],
                    'Dirac_and_energy':[spectral_energy_and_evolution(n) for n in (3,4)],
                    'same_action_Ward':[ward_sources(n) for n in (3,4)],
                    'CAR_and_source_matrix':car_and_matrix_sources(),
                    'scope':{'all_window_continuum_bound_by_proof':True,
                             'own_reference_per_model':True,'source_matrices_not_assumed_scalar':True,
                             'open_family_exists_without_optimized_radius':True,
                             'finite_effect_error_controls_unbounded_source':False,
                             'full_interacting_SM_certified':False,'apparatus_certified':False,
                             'generation_number_selected':False},'sources':bindings()})
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:json.dump(payload,f,ensure_ascii=False,indent=2,sort_keys=True);f.write('\n')
        print('PASS_AUTHOR_SCOPED EXCLUSIVE_WRITE')
    else:
        saved=json.loads(OUT.read_text(encoding='utf-8'));compare(payload,saved)
        print('PASS_AUTHOR_SCOPED READ_ONLY_MATCH')
    print('n=3,4; exact commutants=1; center gap >=61 alpha^2/64; Erel<1450')
    print('finite CAR/Dirac/source diagnostics only; continuum and open-family claims are analytical')


if __name__=='__main__':main()
