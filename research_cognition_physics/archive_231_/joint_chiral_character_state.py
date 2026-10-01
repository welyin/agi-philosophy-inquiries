"""653: exact one-cell subgroup auxiliary character and original-CAR obstruction.

All character certificates are integer polynomial identities. Full matrices
test the actual 615 chiral frames/Pfaffian on general diagonal subgroup links.
No full interacting transfer matrix, continuum limit or new particles claimed.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np
import joint_subgroup_measure_source as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_chiral_character_state_results.json'
DEN=495*4**8
ZERO=(0,0,0,0)
RHO=(2,1,0,1,0)


def add(a,b,scale=1):
    for x,v in b.items():a[x]=a.get(x,0)+scale*v
    return a


def clean(a):return {k:v for k,v in a.items() if v}


def multiply(a,b):
    out=defaultdict(int)
    for x,u in a.items():
        for y,v in b.items():out[tuple(i+j for i,j in zip(x,y))]+=u*v
    return clean(out)


def canonical(v):return tuple(v[i]-v[4] for i in range(4))


def sign(perm):return (-1)**sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))


PERMS=[(a+tuple(3+j for j in b),sign(a)*sign(b))
       for a in itertools.permutations(range(3)) for b in itertools.permutations(range(2))]


def alternant(weight):
    out={}
    for perm,sg in PERMS:
        key=canonical(tuple(weight[i] for i in perm));out[key]=out.get(key,0)+sg
    return clean(out)


def numerator():
    unit=[tuple(int(i==j) for i in range(4)) for j in range(4)]+[(-1,)*4]
    h=[{ZERO:1}]+[{} for _ in range(8)]
    for e in unit:
        u={ZERO:2,e:1,tuple(-x for x in e):1}
        for r in range(1,9):add(h[r],multiply(u,h[r-1]))
    return h[8]


def complete_h8(values):
    h=np.zeros(9);h[0]=1.
    for a in values:
        for k in range(1,9):h[k]+=a*h[k-1]
    return float(h[8])


def auxiliary_from_v5(v):
    vals=np.linalg.eigvals(v)
    assert np.max(abs(abs(vals)-1))<1e-12
    return complete_h8((1+vals.real)/2)/495


def frames(alpha):
    internal=np.array([sum(alpha[list(s)]) if s else 0 for s in old.old.STATES])
    u=np.zeros((64,32),complex);v=np.zeros_like(u);d=np.zeros((64,64),complex)
    for i,b in enumerate(internal):
        x=np.pi+b
        rot=np.cos(x/2)*np.eye(4)+1j*np.sin(x/2)*old.spin.GAMMA[3]
        xx=-np.cos(x)*np.eye(4)+1j*np.sin(x)*old.spin.GAMMA[3]
        u[4*i:4*i+4,2*i:2*i+2]=rot@old.VP
        v[4*i:4*i+4,2*i:2*i+2]=rot@old.VM
        d[4*i:4*i+4,4*i:4*i+4]=(np.eye(4)+xx)/2
    return u,v,d,internal


def actual_measure_check():
    rng=np.random.default_rng(653)
    inputs=[np.zeros(5),np.array([-2,-2,-2,3,3])*.17]
    for _ in range(8):
        a=rng.normal(size=5);a-=a.mean();inputs.append(a)
    pairing_errors=[];pf_errors=[];physical_errors=[];rows=[]
    for a in inputs:
        u,v,d,internal=frames(a)
        cosine=np.repeat(np.cos(a/2),2)
        for _ in range(4):
            e=rng.normal(size=10);e/=np.linalg.norm(e)
            mat=u.T@np.kron(sum(x*t for x,t in zip(e,old.T)),old.B)@u
            expected=-np.kron(sum(x*t for x,t in zip(cosine*e,old.T)),old.EPS)
            pairing_errors.append(float(np.linalg.norm(mat-expected,2)))
            pf_errors.append(float(abs(old.pfaffian(mat)-np.dot(cosine*e,cosine*e)**8)))
        actual=np.linalg.det(np.kron(np.eye(16),old.VP.conj().T)@d@v)
        physical=float(np.prod(np.cos(internal/2)**2))
        physical_errors.append(float(abs(actual-physical)))
        rows.append(dict(alpha=a.tolist(),auxiliary=complete_h8(np.cos(a/2)**2)/495,physical=physical))
    assert max(pairing_errors+pf_errors+physical_errors)<2e-12
    old_errors=[]
    for theta in (.07,.17,.31,1.03,2.14):
        a=np.array([-2,-2,-2,3,3])*theta
        old_errors.append(abs(complete_h8(np.cos(a/2)**2)/495-old.measure(theta)[0]))
    assert max(old_errors)<1e-14
    return dict(actual_subgroup_frames=len(inputs),auxiliary_matrices=4*len(inputs),
                pairing_error=max(pairing_errors),pfaffian_error=max(pf_errors),
                physical_determinant_error=max(physical_errors),old615_restriction_error=max(old_errors),rows=rows)


def character_certificate():
    p=numerator();irreps=[];rhs={}
    arho=alternant(RHO)
    for l4 in sorted(p):
        l=l4+(0,)
        if not (l[0]>=l[1]>=l[2] and l[3]>=0):continue
        c=0
        for perm,sg in PERMS:
            shift=tuple(l[i]+RHO[i]-RHO[perm[i]] for i in range(5))
            c+=sg*p.get(canonical(shift),0)
        if not c:continue
        assert c>0
        a,b,j=l[0]-l[1],l[1]-l[2],l[3]
        charge=-2*sum(l[:3])+3*l[3]
        assert (2*a+4*b+3*j+charge)%6==0
        dimension=(a+1)*(b+1)*(a+b+2)*(j+1)//2
        irreps.append(dict(a=a,b=b,j=j,Q=charge,multiplicity=c,dimension=dimension,highest=list(l)))
        add(rhs,alternant(tuple(l[i]+RHO[i] for i in range(5))),c)
    lhs=multiply(p,arho)
    assert clean(rhs)==lhs
    dimension=sum(r['multiplicity']*r['dimension'] for r in irreps)
    assert sum(p.values())==dimension==DEN and len(irreps)==495 and len(p)==6661
    # Exact restriction of the full four-torus polynomial to original hypercharge.
    restricted={}
    for k,v in p.items():
        q=-2*sum(k[:3])+3*k[3];restricted[q]=restricted.get(q,0)+v
    spec=importlib.util.spec_from_file_location('probe653',HERE/'round653_drafts/holonomy_charge_probe.py')
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    old_result=probe.run()
    expected={int(k):Fraction(v) for k,v in old_result['coefficients']['auxiliary'].items()}
    assert {q:Fraction(v,DEN) for q,v in restricted.items()}==expected
    return dict(laurent_monomials=len(p),weyl_group_order=len(PERMS),irreducible_types=len(irreps),
        positive_integer_multiplicities=True,dimension=dimension,all_Z6_descend=True,
        minimum_multiplicity=min(r['multiplicity'] for r in irreps),
        maximum_multiplicity=max(r['multiplicity'] for r in irreps),
        exact_alternant_identity=True,alternant_monomials=len(lhs),exact_hypercharge_restriction=True,
        irreps=irreps,original_CAR_obstruction={k:old_result[k] for k in (
            'fixed_original_CAR_modes','original_charge_range','auxiliary_fourier_range','full_fourier_range',
            'full_highest_coefficient','outside_original_charge_weight','outside_original_charge_weight_float','charge_variances')})


def common_state_check():
    rng=np.random.default_rng(6531);group=[]
    for _ in range(7):
        c=old.old.old.gauge.group_exp(rng.normal(size=8)*.7,3)
        w=old.old.old.gauge.group_exp(rng.normal(size=3)*.7,2)
        z=np.exp(1j*rng.normal()*.3)
        group.append(old.old.carrier(c,w,z))
    gm=np.zeros((len(group),len(group)),complex);gz=gm.copy();errors=[]
    for i,x in enumerate(group):
        for j,y in enumerate(group):
            v=x.conj().T@y;m=auxiliary_from_v5(v)
            r=old.old.exterior(v)
            physical=np.linalg.det((np.eye(16)+r)/2)**2
            gm[i,j]=m;gz[i,j]=m*physical
            errors.append(abs(auxiliary_from_v5(group[0]@v@group[0].conj().T)-m))
    assert np.max(abs(gm-gm.conj().T))<1e-12 and np.max(abs(gz-gz.conj().T))<1e-12
    mineig_m=float(np.linalg.eigvalsh(gm).min());mineig_z=float(np.linalg.eigvalsh(gz).min())
    assert mineig_m>-1e-11 and mineig_z>-1e-11 and max(errors)<1e-12
    # Same conditional weight fixes the original response; no source refit.
    theta=.17;step=1e-5
    def full(t):
        v=np.diag(np.exp(1j*np.array([-2,-2,-2,3,3])*t))
        r=old.old.exterior(v)
        return float(np.real(np.linalg.det((np.eye(16)+r)/2)**2))*auxiliary_from_v5(v)
    fd=(np.log(full(theta+step))-np.log(full(theta-step)))/(2*step)
    m,dm=old.measure(theta)
    expected=float(-np.sum(old.Q*np.tan(old.Q*theta/2))+dm/m)
    assert abs(fd-expected)<3e-7
    return dict(noncommuting_group_elements=len(group),auxiliary_Gram_minimum=mineig_m,
        full_Gram_minimum=mineig_z,conjugacy_error=float(max(errors)),
        full_log_source=float(fd),old615_log_source=expected,source_error=abs(fd-expected),
        auxiliary_representation_dimension=DEN,formal_full_dimension=2**32*DEN,
        scope='Positive fixed conditional state and original subgroup; no Hamiltonian-time reconstruction.')


def rising(a,n):
    out=Fraction(1)
    for j in range(n):out*=a+j
    return out


def harmonic_spectrum():
    a=Fraction(9,2)
    c=[[Fraction(1)],[Fraction(0),Fraction(8)]]
    for n in range(1,8):
        nxt=[Fraction(0)]*(n+2)
        for j,v in enumerate(c[n]):nxt[j+1]+=2*(n+4)*v/(n+1)
        for j,v in enumerate(c[n-1]):nxt[j]-=Fraction(n+7,n+1)*v
        c.append(nxt)
    lam0=rising(a,8)/rising(2*a,8)
    rows=[]
    unit=[tuple(int(i==j) for i in range(4)) for j in range(4)]+[(-1,)*4]
    h=[{ZERO:1}]+[{} for _ in range(8)]
    for e in unit+[tuple(-x for x in e) for e in unit]:
        for r in range(1,9):add(h[r],multiply({e:1},h[r-1]))
    recovered={}
    for l in range(9):
        lam=lam0
        for j in range(l):lam*=Fraction(8-j,17+j)
        # Independent exact Funk-Hecke integral after t=2u-1, u~Beta(9/2,9/2).
        integral=Fraction(0)
        for j,co in enumerate(c[l]):
            moment=sum(Fraction(math.comb(j,k)*2**k*(-1)**(j-k))*rising(a,8+k)/rising(2*a,8+k) for k in range(j+1))
            integral+=co*moment
        integral/=sum(c[l])
        assert integral==lam and lam>0
        d=math.comb(l+9,9)-(math.comb(l+7,9) if l>=2 else 0)
        harmonic=dict(h[l])
        if l>=2:add(harmonic,h[l-2],-1)
        coeff=lam*DEN;assert coeff.denominator==1
        add(recovered,harmonic,int(coeff))
        rows.append(dict(l=l,eigenvalue=str(lam),dimension=d,scaled_integer=int(coeff)))
    assert clean(recovered)==numerator()
    assert sum(Fraction(r['eigenvalue'])*r['dimension'] for r in rows)==1
    assert sum(r['dimension'] for r in rows)==35750
    return rows


def temporal_gluing_check():
    spectrum=harmonic_spectrum();rng=np.random.default_rng(6532)
    values,vectors=np.linalg.eigh(old.spin.GAMMA[3])
    wp=vectors[:,values>.5];wm=old.spin.G5@wp;eps=wp.T@old.B@wp
    pp=wp@wp.conj().T;pm=wm@wm.conj().T
    rows=[]
    for n in (1,2,3):
        reps=[]
        for _ in range(n):
            cc=old.old.old.gauge.group_exp(rng.normal(size=8)*.2,3)
            ww=old.old.old.gauge.group_exp(rng.normal(size=3)*.2,2)
            zz=np.exp(1j*rng.normal()*.1)
            reps.append(old.old.exterior(old.old.carrier(cc,ww,zz)))
        shift=np.zeros((16*n,16*n),complex)
        for x,r in enumerate(reps):shift[16*x:16*(x+1),16*((x+1)%n):16*((x+1)%n+1)]=r*(-1 if x==n-1 else 1)
        u=(np.kron(np.eye(16*n),wp)+np.kron(shift.conj().T,wm))/np.sqrt(2)
        v=(np.kron(np.eye(16*n),wp)-np.kron(shift.conj().T,wm))/np.sqrt(2)
        xx=-np.kron(shift.conj().T,pp)-np.kron(shift,pm)
        d=(np.eye(64*n)+xx)/2
        g5=np.kron(np.eye(16*n),old.spin.G5)
        projection_error=float(np.linalg.norm(u@u.conj().T-(np.eye(64*n)-g5@xx)/2,2))
        scalar=rng.normal(size=(n,10));scalar/=np.linalg.norm(scalar,axis=1)[:,None]
        te=np.zeros((16*n,16*n),complex)
        for x,e in enumerate(scalar):te[16*x:16*(x+1),16*x:16*(x+1)]=sum(a*t for a,t in zip(e,old.T))
        pairing=u.T@np.kron(te,old.B)@u
        block=(te+shift.conj()@te@shift.conj().T)/2
        pairing_error=float(np.linalg.norm(pairing-np.kron(block,eps),2))
        expected=1.;coefficient_error=0.
        for x in range(n):
            b=block[16*x:16*(x+1),16*x:16*(x+1)]
            e=np.array([np.trace(t.conj().T@b)/16 for t in old.T])
            coefficient_error=max(coefficient_error,float(np.max(abs(e.imag))))
            expected*=float(np.real(e.conj()@e))**8
        actual=old.pfaffian(pairing)
        relative_pf=float(abs(actual-expected)/max(expected,1e-25))
        hol=np.eye(16,dtype=complex)
        for r in reps:hol=hol@r
        # bar-positive gamma5 basis in this gamma4-adapted frame.
        bar=np.kron(np.eye(16*n),(wp+wm).conj().T/np.sqrt(2))
        physical=np.linalg.det(bar@d@v)
        target=np.linalg.det(np.eye(16)+hol.conj().T)**2/2.**(32*n)
        physical_error=float(abs(physical-target)/max(abs(target),1e-25))
        assert projection_error<2e-12 and pairing_error<2e-12 and coefficient_error<1e-12
        assert relative_pf<1e-10 and physical_error<1e-10
        trace=sum(r['dimension']*Fraction(r['eigenvalue'])**n for r in spectrum)
        rows.append(dict(temporal_sites=n,projector_error=projection_error,pairing_error=pairing_error,
            Pfaffian_relative_error=relative_pf,physical_determinant_relative_error=physical_error,
            trivial_holonomy_auxiliary_partition=str(trace),partition_float=float(trace)))
    return dict(actual_overlap_chain_rows=rows,spherical_harmonic_spectrum=spectrum,
        exact_harmonic_character_identity=True,positive_transfer_support_dimension=35750,
        source='Same 615 m0=1, zero spatial hopping in Euclidean kernel, zero scalar background; AP temporal shift.',
        scope='Exact auxiliary transfer and physical determinant on this temporal-only branch, not original full graph H.')


def run():
    actual=actual_measure_check();certificate=character_certificate();state=common_state_check();gluing=temporal_gluing_check()
    names=('research_note_598.md','research_note_614.md','research_note_615.md','research_note_616.md',
           'research_note_628.md','research_note_629.md','research_note_646.md',
           'joint_subgroup_measure_source.py','joint_spinor_subgroup_mass.py',
           'round653_drafts/holonomy_charge_probe.py','round653_drafts/holonomy_charge_probe_results.json')
    return dict(date='2026-10-02',round=653,tests_run=4,failures=0,errors=0,
        actual_measure=actual,character_certificate=certificate,common_state=state,temporal_gluing=gluing,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
        scope='Actual spatially trivial zero-scalar AP-time overlap branch: original-CAR obstruction, full subgroup positive character and exact auxiliary temporal transfer. No original full graph H equivalence or interacting/continuum/GR completion.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
