"""695: two exact original full-group auxiliary configurations of opposite sign.
Floating eigensolvers/inverses only propose dyadic witnesses; rational residuals,
Gaussian-integer determinants and strict rational inequalities decide the result.
"""
import argparse,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as base
import joint_critical_source_certificate as exact
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_auxiliary_sign_certificate_results.json'
FIXTURE=HERE/'round695_drafts/negative_auxiliary_fixture.json'
SCALE=exact.SCALE


class Q:
    def __init__(self,r=0,i=0):self.r=F(r);self.i=F(i)
    @staticmethod
    def cast(v):return v if isinstance(v,Q) else Q(v)
    def __add__(self,v):v=Q.cast(v);return Q(self.r+v.r,self.i+v.i)
    __radd__=__add__
    def __neg__(self):return Q(-self.r,-self.i)
    def __sub__(self,v):return self+-Q.cast(v)
    def __rsub__(self,v):return Q.cast(v)+-self
    def __mul__(self,v):v=Q.cast(v);return Q(self.r*v.r-self.i*v.i,self.r*v.i+self.i*v.r)
    __rmul__=__mul__
    def __truediv__(self,v):
        v=Q.cast(v);d=v.r*v.r+v.i*v.i
        return Q((self.r*v.r+self.i*v.i)/d,(self.i*v.r-self.r*v.i)/d)
    def conj(self):return Q(self.r,-self.i)
    def __bool__(self):return bool(self.r or self.i)
    def __complex__(self):return complex(float(self.r),float(self.i))


def qi(n):return [[Q(int(i==j)) for j in range(n)] for i in range(n)]
def qadj(a):return [[a[j][i].conj() for j in range(len(a))] for i in range(len(a[0]))]
def qmul(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))] for i in range(len(a))]
def qinverse(a):
    n=len(a);b=[list(row)+ident for row,ident in zip(a,qi(n))]
    for k in range(n):
        p=next(i for i in range(k,n) if b[i][k]);b[k],b[p]=b[p],b[k]
        z=b[k][k];b[k]=[x/z for x in b[k]]
        for i in range(n):
            if i!=k:
                z=b[i][k];b[i]=[x-z*y for x,y in zip(b[i],b[k])]
    return [r[n:] for r in b]
def qdet(a):
    b=[list(row) for row in a];out=Q(1)
    for k in range(len(b)):
        p=next(i for i in range(k,len(b)) if b[i][k])
        if p!=k:b[k],b[p]=b[p],b[k];out=-out
        d=b[k][k];out=out*d
        for i in range(k+1,len(b)):
            z=b[i][k]/d
            for j in range(k+1,len(b)):b[i][j]=b[i][j]-z*b[k][j]
    return out
def qpower(z,k):
    if k<0:z=z.conj();k=-k
    out=Q(1)
    for _ in range(k):out=out*z
    return out


def cayley(data,grid):
    n=len(data['real']);a=[[Q(F(data['real'][i][j],grid),F(data['imag'][i][j],grid)) for j in range(n)] for i in range(n)]
    assert all(a[i][j].r==a[j][i].r and a[i][j].i==-a[j][i].i for i in range(n) for j in range(n))
    eye=qi(n);plus=[[eye[i][j]+Q(0,1)*a[i][j] for j in range(n)] for i in range(n)]
    minus=[[eye[i][j]-Q(0,1)*a[i][j] for j in range(n)] for i in range(n)]
    u=qmul(plus,qinverse(minus));det=qdet(u)
    for i in range(n):u[i][0]=u[i][0]/det
    assert qdet(u).r==1 and qdet(u).i==0
    unit=qmul(qadj(u),u)
    assert all(unit[i][j].r==int(i==j) and unit[i][j].i==0 for i in range(n) for j in range(n))
    return u


def group_modules(data,grid):
    c=cayley(data['color'],grid);w=cayley(data['weak'],grid)
    t=F(data['u1_cayley'],grid);z=Q(1,t)/Q(1,-t)
    qw=[[c[i//2][j//2]*w[i%2][j%2]*z for j in range(6)] for i in range(6)]
    cc=[[x.conj() for x in row] for row in c]
    return [qw,[[x*qpower(z,-4) for x in row] for row in cc],
            [[x*qpower(z,2) for x in row] for row in cc],
            [[x*qpower(z,-3) for x in row] for row in w],[[qpower(z,6)]],[[Q(1)]]]


def encode(a):
    den=1
    for row in a:
        for x in row:den=math.lcm(den,x.r.denominator,x.i.denominator)
    r=np.array([[int(x.r*den) for x in row] for row in a],dtype=object)
    i=np.array([[int(x.i*den) for x in row] for row in a],dtype=object)
    return (r,i),den


def sector_h(links,irrep):
    d=len(links[0][0][irrep]);n=8*d
    small=np.array([[1,0],[1,0],[0,1],[0,1]],complex)
    sp=base.internal.spin;g=small.T@sp.G5@small/2
    gamma=[small.T@sp.GAMMA[k]@small/2 for k in (0,3)]
    def qentry(x):return Q(int(x.real),int(x.imag))
    h=[[Q() for _ in range(n)] for _ in range(n)]
    for site in range(4):
        for a in range(2):
            for b in range(2):
                for j in range(d):h[(2*site+a)*d+j][(2*site+b)*d+j]+=qentry(g[a,b])
    for mu in (0,1):
        forward=-g@(np.eye(2)-gamma[mu])/2
        reverse=-g@(np.eye(2)+gamma[mu])/2
        for site,(t,x) in enumerate(base.prior.SITES):
            target=base.prior.SITES.index((1-t,x) if mu else (t,1-x))
            sign=-1 if mu and t else 1;r=links[mu][site][irrep];ra=qadj(r)
            for aa in range(2):
                for bb in range(2):
                    cf=Q(F(float(forward[aa,bb].real)),F(float(forward[aa,bb].imag)))*sign
                    cr=Q(F(float(reverse[aa,bb].real)),F(float(reverse[aa,bb].imag)))*sign
                    for i in range(d):
                        for j in range(d):
                            h[(2*site+aa)*d+i][(2*target+bb)*d+j]+=cf*r[i][j]
                            h[(2*target+aa)*d+i][(2*site+bb)*d+j]+=cr*ra[i][j]
    return encode(h)


def certified_frame(h,den):
    assert np.array_equal(h[0],h[0].T) and np.array_equal(h[1],-h[1].T)
    ev,vec=np.linalg.eigh(exact.proposal(h,den));v=exact.dyadic(vec)
    lam=np.array([int(round(float(x)*SCALE)) for x in ev],dtype=object)
    orth=exact.norm_upper(exact.sub(exact.mm(exact.adj(v),v),exact.times(exact.eye(len(ev)),SCALE**2)),SCALE**2)
    approx=exact.mm((v[0]*lam[None,:],v[1]*lam[None,:]),exact.adj(v))
    res=exact.norm_upper(exact.sub(exact.times(h,SCALE**3),exact.times(approx,den)),den*SCALE**3)
    eps=res+F(max(abs(x) for x in lam),SCALE)*(2+orth)*orth
    gap=F(min(abs(x) for x in lam),SCALE)-eps
    assert orth<F(1,10**6) and gap>0
    # Polarize V, then rotate its negative projector to the exact negative
    # projector. The polar intertwiner differs from identity by<=2||P-Q||.
    nu=orth+2*eps/gap
    neg=np.array([x<0 for x in lam]);pos=~neg
    return (v[0][:,neg],v[1][:,neg]),(v[0][:,pos],v[1][:,pos]),nu,orth,dict(
        dimension=len(ev),negative_count=int(sum(neg)),positive_count=int(sum(pos)),
        exact_orthogonality_error=str(orth),spectral_residual_upper=str(eps),
        Wilson_gap_lower=str(gap),joint_frame_error_upper=str(nu))


def frames(links):
    pieces=[];reports=[];offset=0;negcount=poscount=0
    for k,d in enumerate((6,3,3,2,1,1)):
        h,den=sector_h(links,k);u,v,nu,eta,row=certified_frame(h,den)
        ix=np.array([(2*s+a)*16+offset+j for s in range(4) for a in range(2) for j in range(d)])
        pieces.append((ix,u,v,nu,eta));reports.append(row);offset+=d
        negcount+=u[0].shape[1];poscount+=v[0].shape[1]
    assert negcount==poscount==64
    u=(np.zeros((128,64),dtype=object),np.zeros((128,64),dtype=object))
    v=(np.zeros((128,64),dtype=object),np.zeros((128,64),dtype=object));cu=cv=0
    for ix,un,vn,nu,eta in pieces:
        m=un[0].shape[1];n=vn[0].shape[1]
        for j in (0,1):u[j][np.ix_(ix,np.arange(cu,cu+m))]=un[j];v[j][np.ix_(ix,np.arange(cv,cv+n))]=vn[j]
        cu+=m;cv+=n
    return u,v,max(x[3] for x in pieces),max(x[4] for x in pieces),reports


def E_fields(fixture):
    g=fixture['grid'];rows=[]
    for ints in fixture['negative_candidate_E_stereographic']:
        norm=sum(x*x for x in ints);den=g*g+norm
        row=[F(g*g-norm,den)]+[F(2*g*x,den) for x in ints]
        assert sum(x*x for x in row)==1;rows.append(row)
    return [[F(1)]+[F(0)]*9 for _ in range(4)],rows


def L_matrix(e):
    bspin=np.array([[1,0],[1,0],[0,1],[0,1]],complex)
    k=base.internal.spin.G5@base.internal.spin.GAMMA[1]
    l2=bspin.T@base.internal.B@k@bspin/2
    j=base.mass.dictionary.dictionary();ts=[j.T@t@j for t in base.internal.T]
    out=[[Q() for _ in range(128)] for _ in range(128)]
    for site,row in enumerate(e):
        for a in range(2):
            for b in range(2):
                spin=Q(int(l2[a,b].real),int(l2[a,b].imag))
                for i in range(16):
                    for j in range(16):
                        z=sum((x*Q(int(t[i,j].real),int(t[i,j].imag)) for x,t in zip(row,ts)),Q())
                        out[(2*site+a)*16+i][(2*site+b)*16+j]=spin*z
    value,den=encode(out)
    assert np.array_equal(value[0],value[0].T) and np.array_equal(value[1],value[1].T)
    return value,den


def inv_bound(a,den):
    b=exact.dyadic(np.linalg.inv(exact.proposal(a,den)))
    delta=exact.norm_upper(exact.sub(exact.times(exact.eye(len(a[0])),den*SCALE),exact.mm(a,b)),den*SCALE)
    beta=exact.norm_upper(b,SCALE);assert delta<1
    return beta/(1-delta),dict(inverse_residual=str(delta),inverse_norm_upper=str(beta/(1-delta)))


def gmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def gsub(a,b):return a[0]-b[0],a[1]-b[1]
def gdiv(a,b):
    den=b[0]*b[0]+b[1]*b[1];r=a[0]*b[0]+a[1]*b[1];i=a[1]*b[0]-a[0]*b[1]
    assert r%den==0 and i%den==0
    return r//den,i//den


def gaussian_determinant(a):
    n=len(a[0]);b=[[(int(a[0][i,j]),int(a[1][i,j])) for j in range(n)] for i in range(n)]
    previous=(1,0);sign=1
    for k in range(n-1):
        p=next(i for i in range(k,n) if b[i][k]!=(0,0))
        if p!=k:b[k],b[p]=b[p],b[k];sign=-sign
        pivot=b[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):b[i][j]=gdiv(gsub(gmul(pivot,b[i][j]),gmul(b[i][k],b[k][j])),previous)
        for i in range(k+1,n):b[i][k]=(0,0)
        previous=pivot
    return b[-1][-1][0]*sign,b[-1][-1][1]*sign


def original_crosscheck(links,fields,proposed_ratio):
    numeric=np.zeros((2,4,16,16),complex);j=base.mass.dictionary.dictionary()
    for mu in range(2):
        for site in range(4):
            left=np.zeros((16,16),complex);start=0
            for block in links[mu][site]:
                d=len(block);left[start:start+d,start:start+d]=np.array([[complex(x) for x in row] for row in block]);start+=d
            numeric[mu,site]=j@left@j.T
    u,v,d,h,gap=base.kernel(numeric);weights=[]
    g=np.kron(np.kron(np.eye(4),base.internal.spin.GAMMA[1]),np.eye(16))
    z=np.zeros_like(g);theta=np.block([[z,1j*g],[1j*g,z]])
    reality=[]
    for e in fields:
        mat=base.fixed_matrices(np.array(e,dtype=float),base.mass.car.PHI)
        data=base.regular(u,v,d,mat,lam=0)
        weights.append(data['weight'])
        reality.append(float(np.max(abs(theta.T@data['N'].conj()@theta-data['N']))))
    ratio=weights[1]/weights[0]
    assert abs(ratio/proposed_ratio-1)<2e-6 and max(reality)<2e-12
    return dict(original_512_weights=[base.old.cpair(w) for w in weights],
        original_weight_ratio=base.old.cpair(ratio),compressed_ratio_relative_error=float(abs(ratio/proposed_ratio-1)),
        original_full_Wilson_gap=gap,antiunitary_scalar_reality_errors=reality,
        numerical_values_not_used_for_certified_sign=True)


def run():
    f=json.loads(FIXTURE.read_text('utf8'));links=[[group_modules(x,f['grid']) for x in row] for row in f['links']]
    u,v,nu,eta,reports=frames(links);delta_c=(2+eta)*nu
    coeffs=[];dets=[];bounds=[];diagnostic=[]
    for e in E_fields(f):
        l,dl=L_matrix(e);c=exact.mm((u[0].T,u[1].T),exact.mm(l,v));dc=dl*SCALE**2
        inv,report=inv_bound(c,dc);r=inv*delta_c;assert r<1
        determinant=gaussian_determinant(c);dets.append(determinant);bounds.append(r)
        coeffs.append(report|dict(relative_matrix_error=str(r),
            determinant_real_integer=str(determinant[0]),determinant_imag_integer=str(determinant[1]),
            positive_denominator_power=str(dc),dimension=64))
        s,log=np.linalg.slogdet(exact.proposal(c,dc));diagnostic.append((complex(s),float(log)))
    product=gmul(dets[1],(dets[0][0],-dets[0][1]))
    assert product[0]<0 and 100*abs(product[1])<-product[0]
    phase_error=64*sum(r/(1-r) for r in bounds);assert phase_error<F(1,100)
    # Same physical prefactor nonzero: the actual chiral projection of the
    # positive full H frame is unitarily equivalent to these two64x64 blocks.
    sp=base.internal.spin;plus=0 if sp.G5[0,0].real>0 else 1;minus=1-plus
    inds=lambda a:np.array([(2*s+a)*16+j for s in range(4) for j in range(16)])
    physical=[]
    for frame,ind in ((v,inds(plus)),(u,inds(minus))):
        block=(frame[0][ind],frame[1][ind]);inv,row=inv_bound(block,SCALE)
        assert inv*nu<1
        physical.append(row|dict(frame_relative_error=str(inv*nu),true_chiral_block_invertible=True))
    g=sp.GAMMA[1]
    assert np.array_equal(g,g.T) and np.max(abs(g.imag))==0
    assert np.array_equal(g.T@base.internal.B@g,base.internal.B)
    proposed_ratio=(diagnostic[1][0]/diagnostic[0][0])*np.exp(diagnostic[1][1]-diagnostic[0][1])
    crosscheck=original_crosscheck(links,E_fields(f),proposed_ratio)
    deps=('research_note_673.md','research_note_674.md','research_note_694.md',
          'round695_drafts/self_complement_entry.py','round695_drafts/self_complement_entry_results.json',
          'round695_drafts/negative_auxiliary_fixture.json','joint_critical_source_certificate.py')
    return dict(date='2026-10-02',round=695,tests_run=2,failures=0,errors=0,
        rational_original_group=dict(grid=f['grid'],modules=[6,3,3,2,1,1],
            Sminus_negative_rank=64,Sminus_positive_rank=64,sector_certificates=reports,
            joint_frame_error=str(nu),compressed_matrix_error=str(delta_c)),
        opposite_auxiliary_signs=dict(compressed_certificates=coeffs,
            proposed_ratio_negative_real_cone_verified_exactly=True,
            determinant_ratio_phase_error_bound=str(phase_error),
            phase_error_diagnostic=float(phase_error),
            physical_prefactor_certificates=physical,
            diagnostic_relative_phase=[float((diagnostic[1][0]/diagnostic[0][0]).real),float((diagnostic[1][0]/diagnostic[0][0]).imag)],
            diagnostic_log_absolute_ratio=diagnostic[1][1]-diagnostic[0][1],
            exact_physical_scalar_weights_real=True,
            both_weights_nonzero_and_opposite_sign=True,at_least_one_strict_negative=True,
            original_matrix_crosscheck=crosscheck),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_full_group_and_16_channels_retained=True,
            universal_pointwise_auxiliary_nonnegativity_refuted=True,
            fixed_flux694_auxiliary_average_not_decided=True,
            original_S9_Haar_Hb_average_or_RP_not_refuted=True,
            no_HF_continuum_GR_completion=True,old_space_interfaces_inherited=True),
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=695,all_checks_passed=True,phase_error=r['opposite_auxiliary_signs']['phase_error_diagnostic'],
        log_ratio=r['opposite_auxiliary_signs']['diagnostic_log_absolute_ratio'])))
