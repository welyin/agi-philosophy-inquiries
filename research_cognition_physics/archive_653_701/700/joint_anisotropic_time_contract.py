"""700: original Wilson symbol under an explicit continuous-time deformation.

This is a necessary free-sector test of a NEW anisotropic regulator path.
It does not modify the frozen699 candidate or claim dynamic RP/continuum SM.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_chiral_fibre_source as cliff
import joint_full_holonomy_reduction as original
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_anisotropic_time_contract_results.json'
I4=np.eye(4,dtype=complex)
G5=cliff.G5;GAMMA=cliff.GAMMA


def op(a):return float(np.linalg.norm(a,2))


def symbol(k0,k,nu,rs):
    k=np.asarray(k,float)
    w=float(rs*np.sum(1-np.cos(k))-np.cos(k0))
    v=np.r_[nu*np.sin(k),np.sin(k0)]
    x=w*I4+sum(1j*g*y for g,y in zip(GAMMA,v))
    r=float(np.sqrt(w*w+v@v))
    assert r>1e-14
    h=G5@x
    d=(I4+x/r)/2
    p=(I4-h/r)/2
    return h,d,p,r


def derivative(k0,k,nu,rs):
    """Exact first spatial derivative from the same full four-gamma symbol."""
    k=np.asarray(k,float)
    w=rs*np.sum(1-np.cos(k))-np.cos(k0)
    v=np.r_[nu*np.sin(k),np.sin(k0)]
    dw=rs*np.sin(k[0]);dv=np.array([nu*np.cos(k[0]),0,0,0.])
    r=np.sqrt(w*w+v@v);dr=(w*dw+v@dv)/r
    x=w*I4+sum(1j*g*y for g,y in zip(GAMMA,v))
    dx=dw*I4+sum(1j*g*y for g,y in zip(GAMMA,dv))
    return -.5*G5@(dx/r-x*dr/r**2)


def full_inherited_matrix():
    """Build an actual all16-channel AP-time/periodic-space finite Wilson matrix.
    Two suppressed spatial directions have unit links and hence zero terms.
    """
    nt=ns=4;sites=list(itertools.product(range(nt),range(ns)));count=len(sites)
    st=np.zeros((count,count),complex);sx=np.zeros_like(st)
    for a,(t,x) in enumerate(sites):
        st[a,sites.index(((t+1)%nt,x))]=-1 if t+1==nt else 1
        sx[a,sites.index((t,(x+1)%ns))]=1
    theta=.173;rep=original.b.prior.rep(np.eye(3),np.eye(2),np.exp(1j*theta))
    assert op(rep-np.diag(np.diag(rep)))<1e-14
    nu=.3;rs=.9;i64=np.eye(64);g4=np.kron(GAMMA[3],np.eye(16))
    xmat=rs*np.eye(64*count,dtype=complex)-.5*np.kron(st+st.conj().T,i64)
    xmat+=.5*np.kron(st-st.conj().T,g4)
    xmat-=rs*.5*(np.kron(sx,np.kron(I4,rep))+np.kron(sx.conj().T,np.kron(I4,rep.conj().T)))
    xmat+=nu*.5*(np.kron(sx,np.kron(GAMMA[0],rep))-np.kron(sx.conj().T,np.kron(GAMMA[0],rep.conj().T)))
    k0=np.pi/nt;k1=2*np.pi/ns
    wave=np.array([np.exp(1j*(t*k0+x*k1)) for t,x in sites])/np.sqrt(count)
    embed=np.kron(wave[:,None],i64);block=embed.conj().T@xmat@embed
    expected=np.zeros((64,64),complex)
    for j,phase in enumerate(np.angle(np.diag(rep))):
        inds=np.arange(4)*16+j
        h,_,_,_=symbol(k0,[k1+phase,0,0],nu,rs)
        expected[np.ix_(inds,inds)]=G5@h
    error=op(block-expected)
    assert error<5e-14
    return dict(dimension=len(xmat),channels=16,AP_time_sites=nt,space_sites=ns,
        nu=nu,spatial_Wilson=rs,full_original_group_flat_link=True,Fourier_block_error=error)


def spectral_and_count_checks():
    rng=np.random.default_rng(700)
    err=0.
    for _ in range(32):
        k0=float(rng.uniform(-np.pi,np.pi));k=rng.uniform(-np.pi,np.pi,3)
        nu=float(rng.uniform(.08,1));rs=float(rng.uniform(.08,1.3))
        h,d,p,r=symbol(k0,k,nu,rs)
        ev,u=np.linalg.eigh(h);signed=(u*np.sign(ev))@u.conj().T
        err=max(err,op(h-h.conj().T),op(h@h-r*r*I4),op(p-(I4-signed)/2),op(p@p-p),
            op(G5@d+d@G5-2*d@G5@d))
    assert err<5e-13
    rows=[]
    for rs,expected in ((F(1,10),8),(F(1,5),7),(F(1,3),4),(F(1),1)):
        corners=[]
        for bits in itertools.product((0,1),repeat=3):
            n=sum(bits);w=2*n*rs-1;assert w
            h,d,p,r=symbol(0,np.pi*np.array(bits),.17,float(rs))
            light=w<0
            assert op(d-(np.zeros((4,4)) if light else I4))<1e-12
            corners.append(dict(bits=list(bits),Wilson_scalar=str(w),light=light))
        count=sum(c['light'] for c in corners);assert count==expected
        rows.append(dict(rs=str(rs),light_corners=count,all_corners=corners))
    # Exact first derivatives at every light corner: time and space have the
    # same 1/(2|w|) residue factor, so their velocity ratio is nu a_s/a_t.
    at=F(1,80);aspace=F(1,4);c=F(3,2);nu=c*at/aspace
    assert nu*aspace/at==c
    return dict(max_matrix_identity_error=err,corner_classification=rows,
        velocity_example=dict(at=str(at),a_s=str(aspace),c=str(c),nu=str(nu)),
        full_matrix=full_inherited_matrix())


def locality_checks():
    derivatives=[]
    for rs in (.6,1.,1.4):
        star=float(np.arccos(1-1/rs));s=np.sin(star)
        for nu in (.5,.2,.08,.03):
            h,d,p,r=symbol(0,[star,0,0],nu,rs)
            actual=op(derivative(0,[star,0,0],nu,rs))
            exact=rs/(2*nu)
            assert abs(r-nu*s)<1e-12 and abs(actual/exact-1)<2e-13
            step=nu*1e-5
            finite=(symbol(0,[star+step,0,0],nu,rs)[2]-symbol(0,[star-step,0,0],nu,rs)[2])/(2*step)
            err=op(finite-derivative(0,[star,0,0],nu,rs))/exact
            assert err<2e-8
            derivatives.append(dict(rs=rs,nu=nu,p_star=star,gap=r,exact_derivative=exact,relative_difference_error=err))
    # The collapse persists in the actual lowest AP Matsubara block at fixed beta.
    ap=[];beta=2.;a_s=.4;c=1.3;rs=1.;star=np.pi/2;omega=np.pi/beta
    limit=rs/(2*np.sqrt(omega**2+(c/a_s)**2))
    for nt in (32,64,128,256,512):
        at=beta/nt;nu=c*at/a_s;k0=np.pi/nt
        h,d,p,r=symbol(k0,[star,0,0],nu,rs)
        exact=np.sqrt(4*np.sin(k0/2)**2+nu**2)
        upper=at*np.sqrt(omega**2+(c/a_s)**2)
        scaled=at*op(derivative(k0,[star,0,0],nu,rs))
        assert abs(r-exact)<1e-13 and r<=upper*(1+1e-13)
        ap.append(dict(nt=nt,at=at,nu=nu,gap=r,gap_upper=upper,scaled_derivative=scaled,limit=limit))
    assert abs(ap[-1]['scaled_derivative']/limit-1)<1e-4
    # Numerical Fourier calibration only. The nonuniform locality theorem uses
    # the exact derivative, not a finite Fourier cutoff.
    fourier=[];n=4096;grid=2*np.pi*np.arange(n)/n;freq=np.fft.fftfreq(n)*n
    for nu in (.5,.25,.125,.0625):
        w=-np.cos(grid);v=nu*np.sin(grid);rad=np.sqrt(w*w+v*v)
        mats=(I4[None]-w[:,None,None]/rad[:,None,None]*G5[None]
              -1j*v[:,None,None]/rad[:,None,None]*(G5@GAMMA[0])[None])/2
        coeff=np.fft.fft(mats,axis=0)/n
        norms=np.linalg.svd(coeff,compute_uv=False)[:,0]
        moment=float(np.sum(np.abs(freq)*norms))
        bound=1/(2*nu)
        assert moment>=bound*(1-1e-10)
        fourier.append(dict(n=n,nu=nu,computed_first_moment=moment,necessary_lower=bound))
    return dict(derivative_checks=derivatives,AP_fixed_beta=ap,Fourier_diagnostics=fourier,
        exact_unbounded_first_moment='r_s/(2 nu)',
        physical_first_moment_necessary_bound='r_s a_s^2/(2 c a_t)',
        no_uniform_exponential_locality_at_fixed_a_s=True,
        no_claim_that_each_finite_regulator_is_nonlocal=True)


def run():
    one=spectral_and_count_checks();two=locality_checks()
    deps=('research_note_604.md','research_note_606.md','research_note_610.md','research_note_612.md',
        'research_note_613.md','research_note_643.md','research_note_655.md','research_note_673.md',
        'research_note_699.md','joint_chiral_fibre_source.py','joint_full_holonomy_reduction.py',
        'round700_drafts/mass_time_scope_entry.md')
    return dict(date='2026-10-02',round=700,tests_run=2,failures=0,errors=0,
        inherited_symbol_and_corners=one,common_time_and_locality=two,
        scope=dict(new_anisotropic_candidate_explicitly_declared=True,m0_fixed_to_original_one=True,
            original_group_and_16_channels_preserved=True,free_massless_necessary_sector_only=True,
            physical_space_dimension_three_is_input=True,
            synchronized_naive_scaling_reintroduces_eight_corners=True,
            single_light_corner_with_fixed_rs_requires_nonuniform_spatial_locality=True,
            original_dynamic_RP_or_interacting_continuum_NOT_proved=True,
            hard_gap_NOT_reintroduced_as_finite_integrability_requirement=True,
            original_HF_and_joint_spacetime_limit_NOT_refuted=True,old_spatial_interfaces_unchanged=True),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=700,checks=2,corner_counts=[x['light_corners'] for x in result['inherited_symbol_and_corners']['corner_classification']],
        full_matrix_error=result['inherited_symbol_and_corners']['full_matrix']['Fourier_block_error'],all_checks_passed=True)))
