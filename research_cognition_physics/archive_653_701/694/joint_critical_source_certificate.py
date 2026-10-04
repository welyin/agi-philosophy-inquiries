"""694: exact rational certificate for a one-sided original-source jump.

Floating eigenvectors and inverses propose dyadic witnesses only. Every decisive
residual, inertia and final inequality is checked with Python integers/Fractions.
Fixed original two-direction box and E=e0; not an S9/Haar average or RP certificate.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as base

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_critical_source_certificate_results.json'
CHARGES=(1,-4,2,-3,6,0);MULT=(6,3,3,2,1,1)
LEFT=F(100967,10**6);RIGHT=F(100968,10**6)
BITS=38;SCALE=1<<BITS


def cp(a):
    assert np.array_equal(a.real,np.rint(a.real)) and np.array_equal(a.imag,np.rint(a.imag))
    return np.array(a.real,dtype=object).astype(object),np.array(a.imag,dtype=object).astype(object)


def integer_pair(a):
    r,i=cp(a)
    return np.vectorize(int,otypes=[object])(r),np.vectorize(int,otypes=[object])(i)


def mm(a,b):
    r,i=a;s,t=b
    return r@s-i@t,r@t+i@s


def adj(a):return a[0].T,-a[1].T
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def times(a,k):return a[0]*k,a[1]*k
def eye(n):return np.eye(n,dtype=object),np.zeros((n,n),dtype=object)


def norm_upper(a,den=1):
    # |Re|+|Im| dominates each complex modulus; max(1-norm,infinity-norm)
    # dominates spectral norm, with no floating square root.
    major=abs(a[0])+abs(a[1])
    return F(int(max(np.max(np.sum(major,axis=0)),np.max(np.sum(major,axis=1)))),int(den))


def proposal(a,den):
    return np.array(a[0]/den,dtype=float)+1j*np.array(a[1]/den,dtype=float)


def dyadic(a):
    return (np.vectorize(lambda x:int(round(float(x)*SCALE)),otypes=[object])(a.real),
            np.vectorize(lambda x:int(round(float(x)*SCALE)),otypes=[object])(a.imag))


def gaussian_power(z,n):
    out=(1,0)
    for _ in range(n):out=(out[0]*z[0]-out[1]*z[1],out[0]*z[1]+out[1]*z[0])
    return out


def sector_coefficients(q,sign):
    """Exact integer Laurent coefficients of four times the original8x8 H."""
    sp=base.internal.spin;sites=base.prior.SITES
    b=np.kron(np.eye(4),np.array([[1,0],[-sign,0],[0,1],[0,-sign]]))
    g5=np.kron(np.eye(4),sp.G5)
    out={0:2*b.T@g5@b}
    for mu in (0,1):
        gamma=np.kron(np.eye(4),sp.GAMMA[3 if mu else 0])
        for idx,(t,x) in enumerate(sites):
            target=((1-t,x) if mu else (t,1-x));j=sites.index(target)
            exponent=(q if mu else 2*q*t) if x else 0
            shift=np.zeros((16,16),complex)
            shift[4*idx:4*idx+4,4*j:4*j+4]=(-1 if mu and t else 1)*np.eye(4)
            for k,v in ((exponent,-b.T@g5@(np.eye(16)-gamma)@shift@b),
                        (-exponent,-b.T@g5@(np.eye(16)+gamma)@shift.conj().T@b)):
                out[k]=out.get(k,np.zeros((8,8),complex))+v
    return {k:integer_pair(v) for k,v in out.items()}


def rational_sector(t,q,sign):
    coeff=sector_coefficients(q,sign);power=max(abs(k) for k in coeff)
    p,r=t.numerator,t.denominator;d=r*r+p*p;z=(r*r-p*p,2*r*p)
    numerator=(np.zeros((8,8),dtype=object),np.zeros((8,8),dtype=object))
    for k,c in coeff.items():
        x,y=gaussian_power((z[0],z[1] if k>=0 else -z[1]),abs(k))
        term=((c[0]*x-c[1]*y)*d**(power-abs(k)),
              (c[0]*y+c[1]*x)*d**(power-abs(k)))
        numerator=add(numerator,term)
    assert np.array_equal(numerator[0],numerator[0].T)
    assert np.array_equal(numerator[1],-numerator[1].T)
    return numerator,4*d**power


def inertia(a):
    # Reuse676 rational symmetric-congruence proof, generalized to large integers.
    a=[[F(int(x)) for x in row] for row in a]
    counts=[0,0,0]
    while a:
        n=len(a);i=next((i for i in range(n) if a[i][i]),None)
        if i is not None:
            order=[i]+[j for j in range(n) if j!=i]
            a=[[a[i][j] for j in order] for i in order];p=a[0][0]
            counts[0 if p>0 else 1]+=1
            a=[[a[i][j]-a[i][0]*a[0][j]/p for j in range(1,n)] for i in range(1,n)]
        else:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if a[i][j]),None)
            if pair is None:counts[2]+=n;break
            i,j=pair;order=[i,j]+[k for k in range(n) if k not in (i,j)]
            a=[[a[i][j] for j in order] for i in order];p=a[0][1]
            counts[0]+=1;counts[1]+=1
            a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/p
                for j in range(2,n)] for i in range(2,n)]
    return counts


def exact_counts(t,q,sign):
    (r,i),_=rational_sector(t,q,sign)
    counts=inertia(np.block([[r,-i],[i,r]]))
    assert all(n%2==0 for n in counts)
    return [n//2 for n in counts]


def spectral_witness(q,sign):
    h,den=rational_sector(LEFT,q,sign)
    values,vectors=np.linalg.eigh(proposal(h,den))
    lam=[int(round(float(v)*SCALE)) for v in values];v=dyadic(vectors)
    orth=norm_upper(sub(mm(adj(v),v),times(eye(8),SCALE**2)),SCALE**2)
    assert orth<F(1,10**8)
    vl=(v[0]*np.array(lam,dtype=object)[None,:],v[1]*np.array(lam,dtype=object)[None,:])
    approx=mm(vl,adj(v))
    residual=norm_upper(sub(times(h,SCALE**3),times(approx,den)),den*SCALE**3)
    lmax=F(max(abs(x) for x in lam),SCALE)
    spectral_error=residual+lmax*(2+orth)*orth
    half_gap=F(lam[4]-lam[3],2*SCALE)
    assert half_gap>F(1,10)
    anchor_gap=half_gap-spectral_error
    projection_error=spectral_error/anchor_gap+(2+orth)*orth
    variation=6*abs(q)*(RIGHT-LEFT)
    assert anchor_gap>variation
    total_error=projection_error+variation/(anchor_gap-variation)
    small_v=(v[0][:,:4],v[1][:,:4]);proj=mm(small_v,adj(small_v))
    counts_left=exact_counts(LEFT,q,sign);counts_right=exact_counts(RIGHT,q,sign)
    assert counts_left==[4,4,0] and counts_right[2]==0
    # No other charge can cross zero in the interval; q6 has one isolated
    # possible zero per S sector. Eigenvalue proposals are validated above.
    distances=[F(abs(x),SCALE)-spectral_error-variation for x in lam]
    possible=[i for i,d in enumerate(distances) if d<=0]
    assert len(possible)==(1 if q==6 else 0)
    return proj,dict(charge=q,spin_sign=sign,inertia_left=counts_left,inertia_right=counts_right,
        orthogonality_error=str(orth),diagonalization_error=str(spectral_error),
        half_cluster_gap=str(half_gap),anchor_cut_gap_lower=str(anchor_gap),
        interval_H_variation_bound=str(variation),projector_error_bound=str(total_error),
        possible_zero_eigenvalue_indices=possible),total_error


def full_projector(projs):
    # Direct sum of original left charges, conjugated by the original signed J.
    j=base.mass.dictionary.dictionary();jr,ji=integer_pair(j)
    assert not np.any(ji) and np.array_equal(jr.T@jr,np.eye(16,dtype=object))
    out=(np.zeros((256,256),dtype=object),np.zeros((256,256),dtype=object))
    charge_list=[q for q,m in zip(CHARGES,MULT) for _ in range(m)]
    for idx,q in enumerate(charge_list):
        one=(np.zeros((16,16),dtype=object),np.zeros((16,16),dtype=object))
        for sign in (-1,1):
            b=np.kron(np.eye(4,dtype=object),np.array([[1,0],[-sign,0],[0,1],[0,-sign]],dtype=object))
            one=add(one,(b@projs[q,sign][0]@b.T,b@projs[q,sign][1]@b.T))
        target=int(np.flatnonzero(np.array(jr[:,idx],dtype=int))[0])
        inds=np.array([16*i+target for i in range(16)])
        out[0][np.ix_(inds,inds)]+=one[0];out[1][np.ix_(inds,inds)]+=one[1]
    return out,2*SCALE**2


def full_N(p,den):
    g5=np.kron(np.kron(np.eye(4),base.internal.spin.G5),np.eye(16))
    g=integer_pair(g5);qp=integer_pair((np.eye(256)+g5)/2);qm=integer_pair((np.eye(256)-g5)/2)
    m=integer_pair(np.kron(np.eye(4),np.kron(base.internal.B,base.internal.T[0])))
    mb=integer_pair(np.kron(np.eye(4),np.kron(base.internal.B,base.internal.T[0].conj().T)))
    d=sub(times(qp,den),mm(g,p))
    top=times(mm(m,qp),den);bottom=times(mm(mb,qm),den)
    n=(np.block([[top[0],-d[0].T],[d[0],bottom[0]]]),
       np.block([[top[1],-d[1].T],[d[1],bottom[1]]]))
    assert np.array_equal(n[0],-n[0].T) and np.array_equal(n[1],-n[1].T)
    return n


def connected_blocks(n):
    adjacent=(n[0]!=0)|(n[1]!=0);left=set(range(len(n[0])));blocks=[]
    while left:
        todo=[min(left)];found=set(todo)
        while todo:
            k=todo.pop()
            new=set(np.flatnonzero(adjacent[k]))-found
            found.update(new);todo.extend(new)
        left-=found;blocks.append(sorted(found))
    return blocks


def inverse_certificate(n,den):
    rows=[];bounds=[]
    for indices in connected_blocks(n):
        a=(n[0][np.ix_(indices,indices)],n[1][np.ix_(indices,indices)])
        b=dyadic(np.linalg.inv(proposal(a,den)))
        res=norm_upper(sub(times(eye(len(indices)),den*SCALE),mm(a,b)),den*SCALE)
        beta=norm_upper(b,SCALE)
        assert res<1
        bound=(1-res)/beta;bounds.append(bound)
        rows.append(dict(size=len(indices),inverse_residual_upper=str(res),
                         inverse_norm_upper=str(beta),singular_value_lower=str(bound)))
    return min(bounds),rows


def run():
    projs={};rows=[];errors=[]
    for q in CHARGES:
        for sign in (-1,1):
            p,row,error=spectral_witness(q,sign)
            projs[q,sign]=p;rows.append(row);errors.append(error)
    counts=[]
    for side in ('left','right'):
        counts.append([sum(m*next(row['inertia_'+side][1] for row in rows
            if row['charge']==q and row['spin_sign']==sign) for q,m in zip(CHARGES,MULT))
            for sign in (-1,1)])
    assert counts==[[64,64],[65,63]]
    p,den=full_projector(projs);n=full_N(p,den)
    lower,blocks=inverse_certificate(n,den);error=max(errors)
    assert lower>error
    # Cross-check construction against original kernel; this diagnostic is not
    # used for the rational residual proof or the strict inequality above.
    path=HERE/'round694_drafts/critical_source_entry.py'
    spec=importlib.util.spec_from_file_location('entry694',path)
    entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
    s=4*np.arctan(float(LEFT))/np.pi
    u,v,d,h,gap=base.kernel(entry.links_at(s))
    e=np.tile(np.eye(10)[0],(4,1));mat=base.fixed_matrices(e,base.mass.car.PHI)
    original=base.regular(u,v,d,mat,lam=0)
    difference=float(np.linalg.norm(proposal(n,den)-original['N'],2))
    assert difference<1e-8
    deps=('research_note_673.md','research_note_676.md','research_note_693.md',
          'round676_drafts/exact_flux_seed.py','round694_drafts/critical_source_entry.py',
          'round694_drafts/critical_source_entry_results.json','round694_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=694,tests_run=2,failures=0,errors=0,
        rational_path=dict(z='(1+i*t)/(1-i*t)',left=str(LEFT),right=str(RIGHT),
            inherited_spatial_z_power='2*time_index',temporal_z_power=1,
            full_negative_spin_counts=counts,sector_certificates=rows,
            crossings_exist_by_exact_inertia=True,root_uniqueness_not_required_or_claimed=True),
        nonzero_continued_source=dict(dyadic_bits=BITS,exact_inverse_blocks=blocks,
            exact_approximate_N_singular_lower=str(lower),
            uniform_continuation_N_error_upper=str(error),
            continued_N_singular_lower=str(lower-error),
            strict_positive_margin=float(lower-error),
            original_full_N_crosscheck_error=difference,
            original_anchor_weight=base.old.cpair(original['weight']),
            scalar_Pfaffian_one_sided_limit_nonzero_by_certificate=True,
            B_source_limit_not_certified=True),
        dependency_hashes={x:hashlib.sha256((HERE/x).read_bytes()).hexdigest() for x in deps},
        scope=dict(original_all16_channels_retained=True,original_full512_N=True,
            fixed_E_scalar_source_jump=True,no_new_transversality_or_hard_gap=True,
            full_auxiliary_or_Haar_average_jump_not_proved=True,
            no_RP_HF_continuum_GR_completion=True,old_space_interfaces_inherited=True),
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=694,all_checks_passed=True,
        counts=result['rational_path']['full_negative_spin_counts'],
        margin=result['nonzero_continued_source']['strict_positive_margin'],
        blocks=[r['size'] for r in result['nonzero_continued_source']['exact_inverse_blocks']])) )
