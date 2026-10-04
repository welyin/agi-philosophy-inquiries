"""657: original free auxiliary measure as a half-time BCS Gram kernel.

Exact finite reflection factorization, not a physical one-step transfer claim.
Strict integrated normalization is proved here for two time slices; general
even-time unnormalized reflection positivity is kept separate.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_spatial_auxiliary_geometry as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_auxiliary_reflection_gluing_results.json'


def fun(a,fn):
    ev,v=np.linalg.eigh((a+a.conj().T)/2)
    return (v*fn(ev))@v.conj().T


def reflection(f):
    sites=f['sites'];nt=max(t for t,x in sites)+1;nx=max(x for t,x in sites)+1
    assert nt%2==0
    hs=nx*nt//2;n=4*hs
    reflected=[4*((nt-1-t)*nx+x)+s for t in range(nt//2) for x in range(nx) for s in range(4)]
    p=f['P'];g=np.kron(np.eye(hs),old.old.spin.G5@old.old.spin.GAMMA[3])
    d=p[:n,:n];c=g.conj().T@p[np.ix_(reflected,range(n))]
    lower=g.conj().T@p[np.ix_(reflected,reflected)]@g
    a=fun(d,lambda x:np.sqrt(np.clip(x,0,1)))
    b=fun(d,lambda x:np.sqrt(np.clip(1-x,0,1)))
    return dict(D=d,C=c,lower=lower,A=a,B=b,G=g,reflected=reflected,half_sites=hs,
        C_expected=fun(d,lambda x:np.sqrt(np.clip(x*(1-x),0,None))))


def field_pair(e):
    n=len(e);m=np.zeros((64*n,64*n),complex)
    for i,ei in enumerate(e):
        m[64*i:64*(i+1),64*i:64*(i+1)]=np.kron(old.old.B,sum(x*t for x,t in zip(ei,old.old.T)))
    return m


def full_configuration(e,f,nt,nx):
    return np.concatenate((e.reshape(nt//2,nx,10),f.reshape(nt//2,nx,10)[::-1]),axis=0).reshape(nt*nx,10)


def block_check():
    rows=[]
    for nx,nt in ((2,2),(3,2),(3,4),(3,6),(4,8)):
        f=old.frame(1,nx,nt);r=reflection(f);n=len(r['D']);identity=np.eye(n)
        error=max(np.max(abs(r['C']-r['C'].conj().T)),np.max(abs(r['lower']-(identity-r['D']))),
            np.max(abs(r['D']@r['C']-r['C']@r['D'])),np.max(abs(r['C']@r['C']-(r['D']-r['D']@r['D']))))
        ceig=np.linalg.eigvalsh((r['C']+r['C'].conj().T)/2)
        sqrt_error=float(np.max(abs(r['C']-r['C_expected'])))
        assert error<1e-13 and min(ceig)>-2e-14 and sqrt_error<4e-8
        d=r['D'];bil=np.kron(np.eye(r['half_sites']),old.old.B)
        pairing_error=float(np.max(abs(d.T@bil-bil@d)))
        assert pairing_error<2e-14
        rows.append(dict(nx=nx,nt=nt,block_identity_error=float(error),C_minimum=float(min(ceig)),
            C_maximum=float(max(ceig)),square_root_error=sqrt_error,charge_pairing_error=pairing_error,
            D_minimum=float(min(np.linalg.eigvalsh(d))),D_maximum=float(max(np.linalg.eigvalsh(d)))))
    return dict(rows=rows,square_root_roundoff_near_zero_eigenvalues_explicit=True,
                arbitrary_even_time_positivity_uses_primary_free_overlap_spectral_lemma=True)


def physical_pairing_check():
    rng=np.random.default_rng(65701);rows=[]
    for nx,nt in ((3,2),(3,4),(2,6)):
        f=old.frame(1,nx,nt);r=reflection(f);hs=r['half_sites']
        aa=np.kron(r['A'],np.eye(16));bb=np.kron(r['B'],np.eye(16))
        zero=np.zeros((hs,10));zero[:,0]=1;baseline=field_pair(zero)
        norm=old.old.pfaffian(baseline)
        base_error=float(np.max(abs(aa.T@baseline@aa+bb.T@baseline@bb-baseline)))
        assert base_error<3e-13
        for width in (.08,.35):
            e=np.eye(10)[0]+width*rng.normal(size=(hs,10));e/=np.linalg.norm(e,axis=1)[:,None]
            ff=np.eye(10)[0]+width*rng.normal(size=(hs,10));ff/=np.linalg.norm(ff,axis=1)[:,None]
            me=field_pair(e);mf=field_pair(ff)
            phase_error=abs(old.old.pfaffian(me)/norm-1)
            assert phase_error<2e-13
            from_half=old.old.pfaffian(aa.T@me@aa+bb.T@mf@bb)/norm
            original=old.ratio(f,full_configuration(e,ff,nt,nx))
            error=abs(from_half-original)
            assert error<4e-8*max(abs(original),1e-3)
            rows.append(dict(nx=nx,nt=nt,width=width,original_Pfaffian=[float(original.real),float(original.imag)],
                half_pairing_Pfaffian=[float(from_half.real),float(from_half.imag)],
                difference=float(error),reference_matrix_error=base_error,constant_pairing_phase_error=float(phase_error)))
    return dict(rows=rows,full_S9_fields_not_plane_restricted=True,
                original_Pfaffian_not_replaced_by_absolute_value=True)


def subsets(n):return [[i for i in range(n) if mask>>i&1] for mask in range(1<<n)]


def bcs(z):
    out=np.zeros(1<<len(z),complex);out[0]=1
    for mask,s in enumerate(subsets(len(z))):
        if s and len(s)%2==0:out[mask]=old.old.pfaffian(z[np.ix_(s,s)])
    return out


def exterior(g):
    subs=subsets(len(g));out=np.zeros((len(subs),len(subs)),complex);out[0,0]=1
    for i,s in enumerate(subs):
        if not s:continue
        for j,t in enumerate(subs):
            if len(s)==len(t):out[i,j]=np.linalg.det(g[np.ix_(s,t)])
    return out


def unitary(rng,n):
    raw=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));h=(raw+raw.conj().T)/2
    h-=np.trace(h)*np.eye(n)/n
    return fun(h,lambda x:np.exp(1j*x))


def fock_check():
    # An independent explicit Fock check of signs and singular D endpoints.
    rng=np.random.default_rng(65702);n=6;j=np.kron(np.eye(n//2),np.array([[0.,1.],[-1.,0.]]))
    fields=[]
    for _ in range(5):
        q=unitary(rng,n);fields.append(q@j@q.T)
    u=unitary(rng,n);gamma_u=exterior(u.conj());rows=[]
    for ds in (np.array([.1,.2,.35,.65,.8,.9]),np.array([0.,.2,.35,.65,.8,1.])):
        d=(u*ds)@u.conj().T;a=(u*np.sqrt(ds))@u.conj().T;b=(u*np.sqrt(1-ds))@u.conj().T
        diagonal=[]
        for s in subsets(n):
            diagonal.append(np.prod([(1-ds[i])**.25 if i in s else ds[i]**.25 for i in range(n)]))
        w=(gamma_u*np.array(diagonal))@gamma_u.conj().T
        states=np.array([w@bcs(z) for z in fields]);gram=states.conj()@states.T
        pf=np.zeros_like(gram)
        for i,e in enumerate(fields):
            for k,f in enumerate(fields):pf[i,k]=old.old.pfaffian(a.T@e@a+b.T@f@b)/old.old.pfaffian(e)
        error=float(np.max(abs(pf-gram)));assert error<5e-13
        assert min(np.linalg.eigvalsh(w))>-2e-14 and min(np.linalg.eigvalsh(gram))>-2e-13
        rows.append(dict(D_eigenvalues=ds.tolist(),Pfaffian_Fock_error=error,
            Gram_minimum=float(min(np.linalg.eigvalsh(gram))),W_minimum=float(min(np.linalg.eigvalsh(w))),
            W_maximum=float(max(np.linalg.eigvalsh(w)))))
    return dict(rows=rows,Fock_modes=n,explicit_Fock_dimension=1<<n,
        diagnostic_pairings_not_claimed_as_full_original_64_mode_Fock_enumeration=True)


def create_pair(vector,i,j):
    out=np.zeros_like(vector)
    for mask,z in enumerate(vector):
        if z==0 or mask>>i&1 or mask>>j&1:continue
        # a_i^dagger a_j^dagger: the rightmost creator acts first.
        sign=(-1)**((mask&((1<<j)-1)).bit_count()+((mask|(1<<j))&((1<<i)-1)).bit_count())
        out[mask|(1<<i)|(1<<j)]+=sign*z
    return out


def sphere_check():
    # Exact degree-four sphere rule, testing the nilpotent moment formula.
    rng=np.random.default_rng(65703);n=8;d=10;mats=[]
    for _ in range(d):
        x=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));mats.append((x-x.T)*.07)
    def qa(z,v):
        return sum(z[i,k]*create_pair(v,i,k) for i in range(n) for k in range(i+1,n))
    def square_sum(v):return sum(qa(z,qa(z,v)) for z in mats)
    vac=np.zeros(1<<n,complex);vac[0]=1;s1=square_sum(vac);s2=square_sum(s1)
    formula=vac+s1/(2*d)+s2/(8*d*(d+2))
    selected=([0,1,2,3],[4,5,6,7],list(range(8)));actual=np.zeros(3,complex)
    def add(e,weight):
        z=sum(x*a for x,a in zip(e,mats))
        for k,s in enumerate(selected):actual[k]+=weight*old.old.pfaffian(z[np.ix_(s,s)])
    for i in range(d):
        for sign in (-1,1):
            e=np.zeros(d);e[i]=sign;add(e,1/(d*(d+2)))
    for signs in itertools.product((-1,1),repeat=d):add(np.array(signs)/np.sqrt(d),d/((d+2)*2**d))
    predicted=np.array([formula[sum(1<<i for i in s)] for s in selected]);error=float(max(abs(predicted-actual)))
    assert error<3e-14
    f=old.frame(1,3,2);r=reflection(f)
    log_lower=8*np.linalg.slogdet(r['D'])[1]
    assert abs(log_lower+160*np.log(2))<1e-11
    coefficients=[]
    rising=1;factorial=1
    for k in range(17):
        if k:rising*=k+4;factorial*=k
        coefficients.append(str(Fraction(1,4**k*factorial*rising)))
    return dict(moment_rule_error=error,exact_S9_average_coefficients=coefficients,
        nilpotent_Q_powers_per_original_site=16,original_two_slice_log_partition_lower_bound=float(log_lower),
        original_two_slice_partition_lower_bound_exact='2^(-160)',
        sphere_moments_are_exact_integration_not_planar_measure=True,
        full_original_partition_numerical_value_not_computed=True,
        general_even_time_strict_normalization_not_inferred_from_semidefinite_limit=True)


def run():
    deps=('joint_spatial_auxiliary_geometry.py','joint_subgroup_measure_source.py','joint_chiral_fibre_source.py',
          'research_note_612.md','research_note_653.md','research_note_656.md')
    return dict(date='2026-10-02',round=657,tests_run=4,failures=0,errors=0,
        reflection_blocks=block_check(),original_pairing=physical_pairing_check(),
        explicit_Fock=fock_check(),sphere_and_normalization=sphere_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original free finite even-time auxiliary Pfaffian admits an unnormalized reflection Gram representation using the free overlap spectral lemma. Full S9 integration gives a squared norm; strict normalization proved for two slices. No arbitrary-time strictly positive normalization proof, common one-step transfer, physical original CAR identification, interacting gauge measure, continuum or GR reconstruction.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as stream:json.dump(result,stream,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
