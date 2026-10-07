"""773: local gauge extraction from existing material jets, and jet homotopy.

Numerical checks test identities, not a full spacetime anomaly or an instrument.
The background-patch coefficient class is explicitly wider than a polynomial EFT.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_reference_constraint_strata as original

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_material_local_brst_results.json'


def maxabs(x):
    return float(np.max(np.abs(x)))


def realvec(x):
    a=np.asarray(x).reshape(-1)
    return np.r_[a.real,a.imag]


def internal_extraction():
    old=original.geo.old
    point=np.array([0.,np.pi/2,np.pi/4])
    data=old.fields(point)
    H=data['X']; B=sum(data['a'][0,i]*old.T[i] for i in range(3))+3*data['a0'][0]*np.eye(2)
    dH=np.array([0.,.05*np.sqrt(old.PAR['h2'])*np.cos(point[0]+point[1])],complex)
    DH=dH+1j*B@H
    generators=[*old.T,3*np.eye(2)]
    EW=np.column_stack([realvec([1j*t@H,1j*t@DH]) for t in generators])
    Lw=np.linalg.solve(EW.T@EW,EW.T)
    alpha=np.array([.19,-.23,.31,-.07]); dalpha=np.array([.43,.11,-.29,.17])
    t=sum(x*y for x,y in zip(alpha,generators)); dt=sum(x*y for x,y in zip(dalpha,generators))
    deltaH=1j*t@H; deltadH=1j*dt@H+1j*t@dH
    deltaB=1j*(t@B-B@t)-dt
    deltaDH=deltadH+1j*deltaB@H+1j*B@deltaH
    ewerr=max(maxabs(deltaDH-1j*t@DH),maxabs(Lw@realvec([deltaH,deltaDH])-alpha))
    omitted=maxabs(deltadH+1j*B@deltaH-1j*t@DH)
    col=original.color(1.)
    magpairs=((0,1),(0,2),(1,2))
    mag=[col['F'][i,j] for i,j in magpairs]
    magC=np.column_stack([realvec([1j*(t@f-f@t) for f in mag]) for t in original.T])
    # Original canonical electric density fixes the temporal curvature at lapse 1.
    _,psi,_,_,_=original.completed(16,1.)
    psi_point=float(psi[0,4,2]); kc=original.geo.old.PAR['K'][0]
    electric=col['E']/(kc*psi_point**2)
    Avec=np.concatenate((np.zeros((1,3,3),complex),col['A']))
    dAvec=np.zeros((4,4,3,3),complex);dAvec[0,1:]=electric
    pairs=tuple(itertools.combinations(range(4),2))
    curv=[dAvec[i,j]-dAvec[j,i]+1j*(Avec[i]@Avec[j]-Avec[j]@Avec[i]) for i,j in pairs]
    C=np.column_stack([realvec([1j*(t@f-f@t) for f in curv]) for t in original.T])
    Lc=np.linalg.solve(C.T@C,C.T)
    coeff=np.array([.21,-.17,.13,.29,-.31,.11,-.07,.19])
    ac=sum(x*y for x,y in zip(coeff,original.T))
    dac=[sum((.03*(i+1)*(j+1)-.07)*t for j,t in enumerate(original.T)) for i in range(4)]
    # A symmetric second gauge-parameter jet cancels from d(delta A).
    d2ac=[[sum(.004*(i+j+1)*(k+1)*t for k,t in enumerate(original.T)) for j in range(4)] for i in range(4)]
    da=[1j*(ac@a-a@ac)-d for a,d in zip(Avec,dac)]
    d_da=[[1j*(dac[i]@Avec[j]-Avec[j]@dac[i]+ac@dAvec[i,j]-dAvec[i,j]@ac)-d2ac[i][j] for j in range(4)] for i in range(4)]
    df=[d_da[i][j]-d_da[j][i]+1j*(da[i]@Avec[j]-Avec[j]@da[i]+Avec[i]@da[j]-da[j]@Avec[i]) for i,j in pairs]
    cerr=max(maxabs(Lc@realvec(df)-coeff),maxabs(np.array(df)-np.array([1j*(ac@f-f@ac) for f in curv])))
    ew_eigs=np.linalg.eigvalsh(EW.T@EW); color_eigs=np.linalg.eigvalsh(C.T@C)
    assert np.linalg.matrix_rank(EW)==4 and np.linalg.matrix_rank(C)==8
    assert max(ewerr,cerr,maxabs(Lw@EW-np.eye(4)),maxabs(Lc@C-np.eye(8)))<2e-13
    assert abs(np.linalg.det(np.column_stack([H,DH])))>.1 and omitted>.01
    assert np.linalg.matrix_rank(magC)==7
    eta=np.diag([-1.,1.,1.,1.]); e=np.diag([1.,psi_point**2,psi_point**2,psi_point**2])
    raw=np.arange(16).reshape(4,4)*.019-.08
    lorentz=(raw-eta@raw.T@eta)/2
    de=lorentz@e
    measured=de@np.linalg.inv(e)
    lorerr=max(maxabs((measured-eta@measured.T@eta)/2-lorentz),maxabs(de.T@eta@e+e.T@eta@de))
    assert lorerr<1e-14
    return dict(point=point.tolist(), electroweak_rank=4,color_curvature_rank=8,
                magnetic_only_rank=7,original_psi_point=psi_point,local_Lorentz_extraction_error=lorerr,
                electroweak_Gram_min_eigenvalue=float(ew_eigs[0]),
                color_Gram_min_eigenvalue=float(color_eigs[0]),
                electroweak_identity_error=ewerr,color_identity_error=cerr,
                omitted_connection_variation_error=omitted,
                Higgs_DH_determinant_abs=float(abs(np.linalg.det(np.column_stack([H,DH])))),
                scope='Actual original initial fields, electric curvature and internal Lie directions; arbitrary nonconstant gauge jets. Old background collocation is reused, but gauge extraction uses no Green inverse. No physical gauge breaking or full anomaly calculation.')


def spacetime_extraction():
    # Off-shell jet calibration. Original on-shell chart rank is inherited from 573/753.
    g=np.diag([-1.,1.,1.,1.]); gi=np.linalg.inv(g)
    dh=np.array([1.,.1,0.,0.]); ds=np.array([.2,1.,.1,0.])
    hh=np.array([[.1,.2,.3,.4],[.2,-.2,.1,.2],[.3,.1,.4,-.1],[.4,.2,-.1,.3]])
    hs=np.array([[.2,-.1,.1,.2],[-.1,.3,.2,-.3],[.1,.2,-.2,.2],[.2,-.3,.2,.1]])
    J=np.stack([dh,ds,2*hh@gi@dh,2*hs@gi@ds])
    xi=np.array([.23,-.19,.31,.17]); dxi=np.arange(16).reshape(4,4)*.013-.07
    # dxi[mu,nu]=partial_mu xi^nu.
    dg=dxi@g+g@dxi.T; dgi=-gi@dg@gi
    vdh=dxi@dh+hh@xi; vds=dxi@ds+hs@xi
    dX=np.array([dh@xi,ds@xi,dh@dgi@dh+2*dh@gi@vdh,ds@dgi@ds+2*ds@gi@vds])
    recovered=np.linalg.solve(J,dX)
    missing=dX.copy(); missing[2:]-=[dh@dgi@dh,ds@dgi@ds]
    bad=float(np.linalg.norm(np.linalg.solve(J,missing)-xi))
    assert abs(np.linalg.det(J))>.01 and maxabs(recovered-xi)<1e-14 and bad>.01
    return dict(reference_jet_determinant=float(np.linalg.det(J)),
                scalar_covariance_error=maxabs(dX-J@xi),left_inverse_error=maxabs(recovered-xi),
                omitted_metric_variation_error=bad,
                scope='General finite-jet covariance check of the material extractor. These diagnostic jets are not claimed to solve the original field equations.')


# Exact super-polynomials. Four local jet levels test prolongation without
# confusing this algebra with the Fock-state contraction of round 772.
EVEN=('q0','q1','q2','q3','r','z')  # z: ghost antifield, ghost number -2
ODD=('c0','c1','c2','c3','p','a')   # p: gauge antifield; a: physical antifield
NE=len(EVEN); ZERO=(0,)*NE
def var(name):
    if name in EVEN:
        v=list(ZERO);v[EVEN.index(name)]=1
        return {(tuple(v),0):Q(1)}
    return {(ZERO,1<<ODD.index(name)):Q(1)}
def add(*args):
    out={}
    for p in args:
        for k,v in p.items():out[k]=out.get(k,Q(0))+v
    return {k:v for k,v in out.items() if v}
def scale(p,x):return {k:v*x for k,v in p.items() if v*x}
def mul(p,q):
    out={}
    for (a,m),v in p.items():
        for (b,n),w in q.items():
            if m&n:continue
            inv=sum((n&((1<<i)-1)).bit_count() for i in range(len(ODD)) if m&(1<<i))
            key=(tuple(x+y for x,y in zip(a,b)),m|n)
            out[key]=out.get(key,Q(0))+v*w*(-1)**inv
    return {k:v for k,v in out.items() if v}
ONE={(ZERO,0):Q(1)}
def product(names):
    out=ONE
    for n in names:out=mul(out,var(n))
    return out
def derive(poly,rules,parity):
    out={}
    for (powers,mask),value in poly.items():
        names=[n for n,k in zip(EVEN,powers) for _ in range(k)]+[n for i,n in enumerate(ODD) if mask&(1<<i)]
        prefix=0
        for i,name in enumerate(names):
            if name in rules:
                term=mul(mul(product(names[:i]),rules[name]),product(names[i+1:]))
                out=add(out,scale(term,value*(-1)**(parity*prefix)))
            prefix+=int(name in ODD)
    return out
SR={f'q{i}':var(f'c{i}') for i in range(4)}|{'z':var('p'),'a':var('r')}
HR={f'c{i}':var(f'q{i}') for i in range(4)}|{'p':var('z')}
DR={f'q{i}':var(f'q{i+1}') for i in range(3)}|{f'c{i}':var(f'c{i+1}') for i in range(3)}
def s(p):return derive(p,SR,1)
def hraw(p):return derive(p,HR,1)
def dx(p):return derive(p,DR,0)
def count_gauge(key):
    powers,mask=key
    return sum(powers[:4])+powers[5]+(mask&31).bit_count()
def homotopy(p):
    return hraw({k:v/Q(count_gauge(k)) for k,v in p.items() if count_gauge(k)})
def reduced(p):return {k:v for k,v in p.items() if not count_gauge(k)}
def gh(key):
    powers,mask=key
    return -2*powers[5]+(mask&15).bit_count()-int(bool(mask&16))-int(bool(mask&32))


def local_jet_homotopy():
    generators=EVEN+ODD
    checked=0; positive=0
    for size in range(4):
        for names in itertools.combinations_with_replacement(generators,size):
            p=product(names)
            if not p:continue
            assert not s(s(p))
            assert add(s(homotopy(p)),homotopy(s(p)))==add(p,scale(reduced(p),-1))
            assert dx(homotopy(p))==homotopy(dx(p))
            if gh(next(iter(p)))>0:
                assert not reduced(p)
                positive+=1
            checked+=1
    seed=add(product(['q0','r','r']),product(['q0','q1','r']),product(['z','c0','c1','r']))
    anomaly=s(seed); primitive=homotopy(anomaly)
    assert anomaly and not s(anomaly) and s(primitive)==anomaly
    assert dx(primitive)==homotopy(dx(anomaly))
    return dict(exact_monomials_checked=checked,positive_ghost_number_monomials=positive,
                closed_diagnostic_terms=len(anomaly),primitive_terms=len(primitive),
                nilpotency_and_contraction_exact=True,spacetime_derivative_commutation_exact=True,
                original_quantum_anomaly_computed=False,
                scope='Exact local super-polynomial/jet signs including antifields and horizontal derivative. General local differential splitting and formal recursion are analytic; this is not a computed original loop anomaly.')


def run():
    dependencies=['research_note_573.md','research_note_753.md','research_note_765.md',
                  'research_note_771.md','research_note_772.md','joint_reference_constraint_strata.py',
                  'joint_gauss_continuum_sampling.py','joint_gauss_einstein_initial_data.py',
                  'round773_drafts/research_note_773_working.md']
    return dict(round=773,tests_run=3,failures=0,errors=0,
                internal=internal_extraction(),spacetime=spacetime_extraction(),
                local_jet=local_jet_homotopy(),
                scope='Original regular material patch admits a finite differential left inverse of infinitesimal gauge generators; enlarged smooth background-jet class has positive-ghost local formal contraction. No global polynomial counterterm theorem, complete QME normalization, interacting state, instrument, continuum equivalence or generated GR.',
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
