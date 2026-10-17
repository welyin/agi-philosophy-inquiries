"""1262 independent finite checks: rational source algebra and original material.
No author code/results imported. Default read-only; --write-results is exclusive.
Finite examples diagnose the analytic theorem; they do not prove its quantifiers.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'independent_results.json'
OLD = ROOT / 'archive_1106_1119/1112/check.py'
OLD_SHA = '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f'
ETA = (-1,1,1,1)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(x,y):
    return sum(a*b for a,b in zip(x,y))


def dec(q):
    q=F(q)
    return Decimal(q.numerator)/Decimal(q.denominator)


def shell(m,v):
    """Rational parametrization, not numerical energy square roots."""
    m=F(m); v=tuple(F(x) for x in v); v2=dot(v,v)
    assert m>0 and v2<1
    e=m*(1+v2)/(1-v2)
    p=tuple(2*m*x/(1-v2) for x in v)
    assert e*e-dot(p,p)==m*m
    return e,p


def numerator(m,ep,p,eq,q):
    p4=(ep,)+tuple(p); q4=(eq,)+tuple(q)
    mink=sum(F(ETA[i])*p4[i]*q4[i] for i in range(4))
    return [[p4[i]*q4[j]+q4[i]*p4[j]-(ETA[i]*(mink+m*m) if i==j else 0)
             for j in range(4)] for i in range(4)]


def central(m,ep,p,eq,q):
    mid=tuple((p[i]+q[i])/2 for i in range(3))
    delta=tuple((p[i]-q[i])/2 for i in range(3))
    d2=dot(delta,delta); e2=m*m+dot(mid,mid)
    a=(ep+eq)/2; b=(ep-eq)/2; s2=ep*eq
    assert s2==e2+d2-2*b*b
    assert b*b<=d2 and a*a>=e2 and abs(s2-e2)<=d2
    n=[[F(0) for _ in range(4)] for _ in range(4)]
    n[0][0]=s2+e2-d2
    for i in range(3):
        n[0][i+1]=n[i+1][0]=2*a*mid[i]-2*b*delta[i]
        for j in range(3):
            n[i+1][j+1]=2*mid[i]*mid[j]-2*delta[i]*delta[j]
            if i==j: n[i+1][j+1]+=s2-e2+d2
    return n,mid,delta,e2,s2


def group_source():
    vs=[(0,0,0),(F(1,4),0,0),(-F(1,3),F(1,5),0),
        (F(1,7),-F(1,6),F(1,5)),(F(999,1000),0,0)]
    cases=[]
    for m in (F(1),F(26),F(31)):
        for i,v in enumerate(vs):
            for j,w in enumerate(vs):
                ep,p=shell(m,v); eq,q=shell(m,w)
                cases.append((m,ep,p,eq,q,m))
    nulls=[(F(5),(F(3),F(4),F(0))), (F(5),(F(-3),F(-4),F(0))),
           (F(13),(F(0),F(5),F(12))), (F(1000),(F(0),F(0),F(1000)))]
    for ep,p in nulls:
        for eq,q in nulls:
            cases.append((F(0),ep,p,eq,q,F(1)))
    rows=[]; maxratio=Decimal(0); maxlinear=Decimal(0); diagonal_count=0; zero_midpoint_count=0
    with localcontext() as ctx:
        ctx.prec=75
        for m,ep,p,eq,q,emin in cases:
            n=numerator(m,ep,p,eq,q)
            c,mid,delta,e2,s2=central(m,ep,p,eq,q)
            assert n==c
            transfer=(-(eq-ep),)+tuple(q[i]-p[i] for i in range(3))
            assert all(sum(transfer[i]*n[i][j] for i in range(4))==0 for j in range(4))
            if p==q:
                diagonal_count+=1
                for i in range(4):
                    for j in range(4):
                        p4=(ep,)+p
                        assert n[i][j]/(2*ep)==p4[i]*p4[j]/ep
            s=dec(s2).sqrt(); e=dec(e2).sqrt()
            if e==0: zero_midpoint_count+=1
            t=[[Decimal(0) for _ in range(4)] for _ in range(4)]
            t[0][0]=e
            for i in range(3):
                t[0][i+1]=t[i+1][0]=dec(mid[i])
                for j in range(3):
                    t[i+1][j+1]=dec(mid[i]*mid[j])/e if e else Decimal(0)
            bound=3*dec(dot(delta,delta))/dec(emin)
            err=max(abs(dec(n[i][j])/(2*s)-t[i][j]) for i in range(4) for j in range(4))
            assert err<=bound+Decimal('1e-60')
            # No mass or energy gap is needed for the weaker linear bound.
            radius=dec(dot(delta,delta)).sqrt()
            linear_bound=Decimal(7)*radius/2
            assert err<=linear_bound+Decimal('1e-60')
            assert e*e<=s*s+radius*radius+Decimal('1e-60')
            assert max(abs(dec(n[i][j])/(2*s)) for i in range(4) for j in range(4))<=2*s+Decimal('1e-60')
            assert max(abs(t[i][j]) for i in range(4) for j in range(4))<=e+Decimal('1e-60')
            if linear_bound: maxlinear=max(maxlinear,err/linear_bound)
            if bound: maxratio=max(maxratio,err/bound)
            # Independent stronger energy identity checks sign and normalization.
            energy=dec(n[0][0])/(2*s)-e
            assert abs(energy-((s-e)**2-dec(dot(delta,delta)))/(2*s))<Decimal('1e-60')
            assert -dec(dot(delta,delta))/(2*s)-Decimal('1e-60')<=energy<=Decimal('1e-60')
        for emin in (F(26),F(1),F(5,2)):
            for length in (F(1),F(10),F(100)):
                kernel_coefficient=F(3,4)/emin
                gaussian_second_moment=F(3)/length**2
                bound=kernel_coefficient*gaussian_second_moment
                assert bound==F(9,4)/emin/length**2
                rows.append({'emin':str(emin),'L':str(length),'Schur_bound':str(bound)})
    # Endpoint extension for a massless zero-energy leg: the numerator is zero,
    # the kernel extends continuously by zero, but its Weyl midpoint need not.
    endpoint_cases=0
    with localcontext() as ctx:
        ctx.prec=75
        for ep,p in [(F(0),(F(0),F(0),F(0)))]+nulls:
            for reverse in (False,True):
                eq=F(0); q=(F(0),F(0),F(0))
                if reverse: ep1,p1,eq1,q1=eq,q,ep,p
                else: ep1,p1,eq1,q1=ep,p,eq,q
                n=numerator(F(0),ep1,p1,eq1,q1)
                assert all(x==0 for row in n for x in row)
                _,mid,delta,e2,s2=central(F(0),ep1,p1,eq1,q1)
                assert s2==0
                e=dec(e2).sqrt(); radius=dec(dot(delta,delta)).sqrt()
                # Every Weyl component is bounded by e, hence K=0 suffices.
                assert e<=Decimal(7)*radius/2
                endpoint_cases+=1
    # Exact rational certificate for the split constant: sqrt(2)<3/2.
    assert F(2)<F(3,2)**2 and F(2)+F(3,2)==F(7,2)
    # For a Gaussian in three dimensions C1=2 sqrt(2/pi)/L < 2/L,
    # using pi>3>2; the continuum Gaussian integral is proved in the text.
    nogap=[]
    for length in (F(1),F(10),F(100),F(10)**20):
        upper=F(7,4)*F(2)/length
        assert upper==F(7,2)/length
        nogap.append({'L':str(length),'full_space_gap_free_Schur_upper':str(upper)})
    # Opposite null momenta exhibit failure of a uniform quadratic coefficient
    # down to zero energy. The full-massless linear estimate above remains valid.
    infrared=[]
    for lam in (F(1),F(1,16),F(1,256)):
        p=(lam,F(0),F(0)); q=(-lam,F(0),F(0))
        n=numerator(F(0),lam,p,lam,q)
        assert n[2][2]/(2*lam)==lam
        ratio=lam/(4*lam*lam)
        infrared.append({'energy':str(lam),'difference_over_transfer_squared':str(ratio)})
    # Band compression cannot certify leakage-free insertion: nonzero kernel
    # links p inside [1,3] to q outside. Gaussian Fourier factor is positive.
    ep,p=F(2),(F(0),F(0),F(2))
    eq,q=F(1,2),(F(0),F(0),F(1,2))
    assert ep*eq==1
    cross=numerator(F(0),ep,p,eq,q)[0][0]/2
    assert cross==1
    return {'method':'exact rational mass shells, all 16 tensor components; 75-digit bounds are diagnostics',
            'exact_pair_cases':len(cases),'exact_tensor_comparisons':16*len(cases),
            'exact_Ward_contractions':4*len(cases),'diagonal_cases':diagonal_count,
            'null_zero_midpoint_cases':zero_midpoint_count,
            'max_finite_diagnostic_bound_ratio':str(maxratio),
            'max_gap_free_linear_diagnostic_ratio':str(maxlinear),
            'massless_zero_leg_continuous_extension_cases':endpoint_cases,
            'Gaussian_gap_free_Schur_certificates':nogap,
            'positive_band_needed_only_for_finer_quadratic_bound':True,
            'analytic_continuous_bound_is_proof_not_sampling':True,
            'Gaussian_Schur_certificates':rows,'infrared_quadratic_deletion':infrared,
            'band_crossing_K00_exact':str(cross),
            'compression_not_dynamic_invariance':True}


def group_material():
    """Direct reconstruction of the old physical 864-by-36 accessible corner.
    Other spectator copies and full six-block source bound are analytic.
    """
    assert sha(OLD)==OLD_SHA
    x=np.array([[0,1],[1,0]],complex); y=np.array([[0,-1j],[1j,0]],complex)
    i2=np.eye(2); ur=np.kron((i2-y)/2,i2)-1j*np.kron((i2+y)/2,y)
    um=np.kron(np.diag([1,0]),i2)-1j*np.kron(np.diag([0,1]),x)
    aux0=np.kron(np.array([1,1j])/math.sqrt(2),np.kron([1,0],[1,0]))
    aux=[aux0,aux0,np.kron(ur,i2)@aux0,np.kron(i2,um)@np.kron(ur,i2)@aux0]
    clock=1.5*np.eye(4)
    for j in range(3): clock[j,j+1]=clock[j+1,j]=math.sqrt((j+1)*(3-j))/2
    swap=np.zeros((9,9))
    for q in range(3):
        for r in range(3): swap[3*r+q,3*q+r]=1
    emb=np.zeros((864,36),complex); view=emb.reshape(4,8,3,3,3,36)
    for h in range(4):
        for q in range(3):
            for r in range(3):
                qc,qr=(q,0) if h==0 else (0,q)
                view[h,:,qc,qr,r,9*h+3*q+r]=aux[h]
    assert np.linalg.norm(emb.conj().T@emb-np.eye(36))<1e-12
    physical_flag=np.broadcast_to((np.arange(8)%2)[None,:,None,None,None],(4,8,3,3,3)).reshape(864)
    flag=emb.conj().T@(physical_flag[:,None]*emb)
    assert np.linalg.norm(flag-np.kron(np.diag([0,0,0,1]),np.eye(9)))<1e-12
    baseline=17+2*math.sqrt(19)
    gates=[np.kron(np.eye(8),swap),np.kron(ur,np.eye(18)),np.kron(np.kron(i2,um),np.eye(9))]
    blocks=emb.reshape(4,72,-1); applied=(baseline+1.5)*blocks.copy()
    for j,gate in enumerate(gates):
        c=math.sqrt((j+1)*(3-j))/2
        applied[j+1]+=c*(gate@blocks[j]); applied[j]+=c*(gate.conj().T@blocks[j+1])
    assert np.linalg.norm(applied.reshape(864,36)-emb@np.kron(baseline*np.eye(4)+clock,np.eye(9)))<1e-11
    eig,vec=np.linalg.eigh(clock)
    assert np.max(abs(eig-np.arange(4)))<1e-13
    projectors=[np.zeros((36,36),complex) for _ in range(6)]
    ranks=[0]*6
    for j in range(4):
        for sign in (1,-1):
            k=j+1-sign
            projectors[k]+=np.kron(np.outer(vec[:,j],vec[:,j].conj()),(np.eye(9)+sign*swap)/2)
            ranks[k]+=24*(6 if sign==1 else 3)
    assert sum(ranks)==864 and ranks==[144,144,216,216,72,72]
    assert np.linalg.norm(sum(projectors)-np.eye(36))<1e-12
    masses=np.array([baseline+1+j for j in range(6)])
    assert min(masses)>26 and max(masses)<32
    age=vec@(np.exp(-1j*(math.pi/4)*eig)*(vec.conj().T@np.eye(4)[:,0]))
    load=np.zeros((9,2),complex); load[3,0]=load[6,1]=1
    prepared=np.kron(age[:,None],load)
    weights=[]
    for j,w in enumerate((1,3,4,4,3,1)):
        pull=prepared.conj().T@projectors[j]@prepared
        assert np.linalg.norm(pull-w*np.eye(2)/16)<1e-12
        weights.append(str(F(w,16)))
    mass=sum(m*pj for m,pj in zip(masses,projectors))
    assert np.linalg.norm(prepared.conj().T@np.kron(np.eye(4),np.eye(9)-swap)@prepared-np.eye(2))<1e-12
    assert np.linalg.norm(mass@flag-flag@mass)>3
    records=[]
    # Full effect compression, not two basis-state probability checks.
    for momentum in (0.,5.,50.):
        energies=np.sqrt(masses*masses+momentum*momentum)
        u=sum(np.exp(-1j*(math.pi/2)*(en-energies[0]))*pj for en,pj in zip(energies,projectors))
        evolved=u@prepared
        effect=evolved.conj().T@flag@evolved
        probability=float(np.trace(effect).real/2)
        assert np.linalg.norm(effect-probability*np.eye(2))<1e-12
        dephased=sum(prepared.conj().T@pj@flag@pj@prepared for pj in projectors)
        pd=float(np.trace(dephased).real/2)
        assert abs(pd-5/16)<1e-12
        # Arbitrary-reference claim follows from operator identity; this is one
        # non-product diagnostic, with the reference itself retained.
        bell=np.array([1,0,0,1j],complex)/math.sqrt(2)
        full=np.kron(evolved,np.eye(2))@bell
        pb=float(np.vdot(full,np.kron(flag,np.eye(2))@full).real)
        assert abs(pb-probability)<1e-12
        records.append({'momentum':momentum,'coherent_Mr':probability,'mass_dephased_Mr':pd})
    assert records[0]['coherent_Mr']-records[0]['mass_dephased_Mr']>.30
    # Independent classical kernel values on each actual mass block. These
    # scalars define a diagonal insertion without measuring the mass label.
    p=np.array([2.,-1.,3.]); q=np.array([-4.,2.,1.]); mid=(p+q)/2
    coeffs=[]; delta2=float(np.dot(p-q,p-q))/4
    for m in masses:
        ep=math.sqrt(m*m+np.dot(p,p)); eq=math.sqrt(m*m+np.dot(q,q))
        pp=np.r_[ep,p]; qq=np.r_[eq,q]; mink=-ep*eq+np.dot(p,q)
        kernel=(np.outer(pp,qq)+np.outer(qq,pp)-np.diag(ETA)*(mink+m*m))/(2*math.sqrt(ep*eq))
        e=math.sqrt(m*m+np.dot(mid,mid)); vv=np.r_[e,mid]
        target=np.outer(vv,vv)/e
        assert np.max(abs(kernel-target))<=3*delta2/26+1e-12
        coeffs.append(float(kernel[0,0]))
    insertion=sum(v*pj for v,pj in zip(coeffs,projectors))
    assert np.linalg.norm(insertion@mass-mass@insertion)<1e-10
    # A direct spectral lift preserves off-diagonal density matrices; dephasing
    # is a different CP channel. No locality of Mr follows from these matrices.
    return {'physical_dimension':864,'checked_accessible_corner_dimension':36,
            'full_mass_block_dimensions':ranks,'six_code_weights':weights,
            'masses':list(map(float,masses)),'actual_clock_records':records,
            'mass_dephasing_changes_original_Mr':True,
            'exact_material_Luders_energy_change_upper':'5',
            'query_total_apparatus_energy_or_locality_not_certified':True,
            'source_block_diagonality_is_not_state_dephasing':True,
            'unknown_reference_quantifier_by_tensor_operator_bound':True}


def calculate():
    return {'status':'PASS_INDEPENDENT_SOURCE_BRIDGE_SCOPE','round':1262,'independent_groups':2,
            'independent_code_sha256':sha(Path(__file__)), 'original_1112_definition_sha256':sha(OLD),
            'group1_full_tensor_and_finite_smear':group_source(),
            'group2_original_mass_and_coherence':group_material(),
            'scope':{'minimal_free_KG_one_particle_input':True,'normal_massive_momentum_tails_not_cut':True,
                     'massless_compression_not_evolution_invariance':True,
                     'gap_free_full_massless_insertion_bound':True,
                     'Gaussian_not_compact_local_instrument':True,'individual_sources_not_bounded':True,
                     'no_1261_AW_or_interacting_source_match':True,'no_new_CT_or_GR_completion':True}}


def compare(saved,now,path='root'):
    if isinstance(now,dict):
        assert isinstance(saved,dict) and saved.keys()==now.keys(),path
        for key in now: compare(saved[key],now[key],path+'.'+key)
    elif isinstance(now,list):
        assert isinstance(saved,list) and len(saved)==len(now),path
        for i,(s,n) in enumerate(zip(saved,now)): compare(s,n,f'{path}[{i}]')
    elif isinstance(now,float):
        assert math.isclose(float(saved),now,rel_tol=1e-10,abs_tol=1e-11),(path,saved,now)
    else: assert saved==now,(path,saved,now)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=calculate()
    if args.write_results:
        with OUT.open('x',encoding='utf-8') as f: json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        print('WRITE_ONCE_PASS_INDEPENDENT_1262: independent_groups=2')
    else:
        compare(json.loads(OUT.read_text(encoding='utf-8')),result)
        print('READ_ONLY_MATCH_INDEPENDENT_1262: independent_groups=2')


if __name__=='__main__': main()