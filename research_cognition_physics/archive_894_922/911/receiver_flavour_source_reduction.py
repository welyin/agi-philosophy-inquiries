"""911: contract the original870 preparation difference through all872 source terms.
Exact algebra calibrates a flavour-selection identity for the original PDE;
the finite boson/Dirac matrices below are NOT the original Green operators.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, random, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'846'));sys.path.insert(0,str(STAGE/'870'))
import constrained_six_leg_recursion as ex
import receiver_four_point_backreaction as car
TARGET=HERE/'receiver_flavour_source_reduction_results.json'

def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a):return sum(a[i][i] for i in range(len(a)))
def plus(a,b):return [[x+y for x,y in zip(row,col)] for row,col in zip(a,b)]
def times(a,t):return [[x*t for x in row] for row in a]
def trans(a):return list(map(list,zip(*a)))
def nv(m,v):return [sum(x*y for x,y in zip(row,v)) for row in m]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def tensor(a,b):return [[x*y for x in row for y in rr] for row in a for rr in b]
I=[[F(1),F(0)],[F(0),F(1)]];Y=[[F(1),F(0)],[F(0),F(-1)]]
# Y eigenbasis: Delta rho is invariant under common flavour rotations.
def bil(bar,m,psi):
    return ex.add(*(ex.scale(ex.mul(bar[i],psi[j]),m[i][j]) for i in range(len(bar)) for j in range(len(psi))))
def mixed(bar,m,v,barv,psi):return ex.add(bil(bar,m,v),bil(barv,m,psi))
def delta(p):return -p.get(15,F(0))/2  # unit Delta r; Delta(n1 n2)=-1/2

def finite_case(seed=911,old_scale=F(1),old_identity=False):
    rng=random.Random(seed);d=3
    def sym(den):
        a=F(rng.randint(-3,3),den);b=F(rng.randint(-3,3),den);c=F(rng.randint(-3,3),den)
        return [[a,b],[b,c]]
    mg=[sym(3) for _ in range(d)]
    mo=[times(sym(5),old_scale) for _ in range(d)]
    m1=[sym(7) for _ in range(d)]
    nu=[F(1,3),F(-1,4),F(2,5)]
    h1=[[plus(times(m1[j],nu[i]),times(m1[i],nu[j])) for j in range(d)] for i in range(d)]
    h0=[[times(I,F(i+j+1,19)) for j in range(d)] for i in range(d)]
    C=[[[F(i+j+k+1,11) for k in range(d)] for j in range(d)] for i in range(d)]
    G=[[F(2,3),F(1,7),F(0)],[F(1,7),F(3,4),F(1,9)],[F(0),F(1,9),F(4,5)]]
    K=[[F(2,3),F(1,4)],[F(-1,5),F(3,4)]]
    u=[F(2,3),F(-1,4)];w=[F(3,5),F(1,2)]
    cs=[{2:F(1)},{8:F(1)}];bars=[{1:F(1)},{4:F(1)}]
    psi=[ex.scale(cs[f],u[s]) for s in range(2) for f in range(2)]
    bar=[ex.scale(bars[f],w[s]) for s in range(2) for f in range(2)]
    mG=[tensor(x,I) for x in mg];mOld=[tensor(x,I if old_identity else Y) for x in mo]
    m0=[plus(x,y) for x,y in zip(mG,mOld)];mOne=[tensor(x,Y) for x in m1]
    hOne=[[tensor(x,Y) for x in row] for row in h1]
    hZero=[[tensor(x,I) for x in row] for row in h0]
    KK=tensor(K,I)
    def response(b,ms,pp=psi,bb=bar):
        right=ex.va(*(ex.pv(bi,ex.mv(M,pp)) for bi,M in zip(b,ms)))
        left=ex.va(*(ex.pv(bi,ex.mv(trans(M),bb)) for bi,M in zip(b,ms)))
        return ex.vs(ex.mv(KK,right),-1),ex.vs(ex.mv(trans(KK),left),-1)
    b20=ex.vs(ex.mv(G,[bil(bar,M,psi) for M in m0]),-1)
    b21=ex.vs(ex.mv(G,[bil(bar,M,psi) for M in mOne]),-1)
    A_poly=ex.vs(ex.mv(G,[bil(bar,M,psi) for M in mG]),-1)
    C_poly=ex.vs(ex.mv(G,[bil(bar,M,psi) for M in mOld]),-1)
    v32,z32=response(b21,mOne)
    vg,zg=response(b21,mG);vo,zo=response(b21,mOld)
    va,za=response(A_poly,mOne);vc,zc=response(C_poly,mOne)
    v31=ex.va(vg,vo,va,vc);z31=ex.va(zg,zo,za,zc)
    terms=[
        [ex.scale(ex.add(*(ex.scale(ex.mul(b21[j],b21[k]),C[i][j][k]) for j in range(d) for k in range(d))),F(1,2)) for i in range(d)],
        [mixed(bar,M,v32,z32,psi) for M in m0],
        [mixed(bar,M,v31,z31,psi) for M in mOne],
        [ex.add(*(ex.mul(b21[j],bil(bar,hOne[i][j],psi)) for j in range(d))) for i in range(d)]]
    full=ex.va(*terms)
    dropped=[[mixed(bar,M,v32,z32,psi) for M in mOld],
             [mixed(bar,M,vo,zo,psi) for M in mOne],
             [mixed(bar,M,vc,zc,psi) for M in mOne]]
    def numeric_response(b,ms):
        m=[[sum(b[i]*ms[i][j][k] for i in range(d)) for k in range(2)] for j in range(2)]
        return [-x for x in nv(K,nv(m,u))],[-x for x in nv(trans(K),nv(trans(m),w))]
    A=[-x for x in nv(G,[dot(w,nv(M,u)) for M in mg])]
    B=[-x for x in nv(G,[dot(w,nv(M,u)) for M in m1])]
    rB,sB=numeric_response(B,m1);rg,sg=numeric_response(B,mg);rA,sA=numeric_response(A,m1)
    def mix(m,v,z):return dot(w,nv(m,v))+dot(z,nv(m,u))
    reduced_terms=[
        [sum(C[i][j][k]*B[j]*B[k] for j in range(d) for k in range(d))/2 for i in range(d)],
        [mix(m,rB,sB) for m in mg],
        [mix(m,[x-y for x,y in zip(rg,rA)],[x-y for x,y in zip(sg,sA)]) for m in m1],
        [sum(B[j]*dot(w,nv(h1[i][j],u)) for j in range(d)) for i in range(d)]]
    reduced=[sum(t[i] for t in reduced_terms) for i in range(d)]
    if not old_identity:
        assert all(delta(p)==0 for ds in dropped for p in ds)
        assert [delta(p) for p in full]==reduced
        assert [[delta(p) for p in t] for t in terms]==reduced_terms
    # Independent complete simultaneous Euler substitution, then extract lambda^2.
    def solve(lam):
        ms=[plus(a,times(b,lam)) for a,b in zip(m0,mOne)]
        hs=[[plus(hZero[i][j],times(hOne[i][j],lam)) for j in range(d)] for i in range(d)]
        def forces(b,pp,bb):
            out=[]
            for i in range(d):
                cubic=ex.scale(ex.add(*(ex.scale(ex.mul(b[j],b[k]),C[i][j][k]) for j in range(d) for k in range(d))),F(1,2))
                out.append(ex.add(cubic,bil(bb,ms[i],pp),*(ex.mul(b[j],bil(bb,hs[i][j],pp)) for j in range(d))))
            return out
        def fermions(b,pp,bb):
            vr,zr=response(b,ms,pp,bb)
            for i in range(d):
                for j in range(d):
                    q=ex.scale(ex.mul(b[i],b[j]),F(1,2))
                    vr=ex.va(vr,ex.vs(ex.mv(KK,ex.pv(q,ex.mv(hs[i][j],pp))),-1))
                    zr=ex.va(zr,ex.vs(ex.mv(trans(KK),ex.pv(q,ex.mv(trans(hs[i][j]),bb))),-1))
            return ex.va(psi,vr),ex.va(bar,zr)
        b=[{} for _ in range(d)];pp=psi;bb=bar
        for iteration in range(12):
            bn=ex.vs(ex.mv(G,forces(b,pp,bb)),-1);pn,zn=fermions(b,pp,bb)
            if (bn,pn,zn)==(b,pp,bb):break
            b,pp,bb=bn,pn,zn
        else:raise AssertionError('Euler substitution did not stabilize')
        assert not any(ex.va(b,ex.mv(G,forces(b,pp,bb))))
        return [{15:p[15]} if 15 in p else {} for p in forces(b,pp,bb)]
    pm,p0,pp=solve(F(-1)),solve(F(0)),solve(F(1))
    independent=ex.vs(ex.va(pm,pp,ex.vs(p0,-2)),F(1,2))
    assert full==independent
    return dict(seed=seed,old_scale=str(old_scale),old_communication_is_identity=old_identity,
        source_parts=[[str(delta(p)) for p in t] for t in terms],
        full_preparation_difference_per_unit_delta_r=[str(delta(p)) for p in full],
        reduced_source_per_unit_delta_r=[str(x) for x in reduced],
        three_old_communication_paths=[[str(delta(p)) for p in t] for t in dropped],
        full_recursion_equals_independent_Euler_solution=True,
        reduced_matches_full=([delta(p) for p in full]==reduced),
        original_continuum_source_evaluated=False)

def car_check():
    rng=np.random.default_rng(911);cs=car.annihilators(4)
    ext=np.diag([.28,.52,.07,.13])
    drho=np.kron(car.receiver_state(1.)-car.receiver_state(.5),ext)
    dr=.5
    def op(M):return sum(M[i,j]*cs[i].conj().T@cs[j] for i in range(4) for j in range(4))
    err=0.;one=0.
    for _ in range(40):
        a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));a=(a+a.conj().T)/2
        b=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));b=(b+b.conj().T)/2
        want=dr/2*(np.trace(a[:2,:2]@b[:2,:2])-np.trace(a[:2,:2])*np.trace(b[:2,:2]))
        got=np.trace(drho@op(a)@op(b))
        err=max(err,float(abs(got-want)));one=max(one,float(abs(np.trace(drho@op(a)))))
    ni=np.zeros((4,4),complex);ni[:2,:2]=np.eye(2)
    yi=np.zeros((4,4),complex);yi[:2,:2]=car.Y
    n=op(ni);y=op(yi)
    vals=[float(np.trace(drho@o).real) for o in (n@n,y@y,n@y,y@n)]
    assert np.max(abs(np.array(vals)-np.array([-dr,dr,0,0])))<1e-14
    assert err<2e-14 and one<2e-14
    return dict(original870_states_used=True,CAR_modes=4,retained_modes=2,mode_exterior_modes=2,
        tested_general_full_bilinear_pairs=40,complete_pair_identity_max_error=err,
        one_body_difference_max_error=one,delta_N2_Y2_NY_YN=vals,
        arbitrary_mode_exterior_cross_blocks_included=True)

def run():
    exact=[]
    for seed in (911,1911,2911):
        cases=[finite_case(seed,F(s)) for s in (0,1,3)]
        assert all(c['full_preparation_difference_per_unit_delta_r']==cases[0]['full_preparation_difference_per_unit_delta_r'] for c in cases)
        assert any(F(x) for x in cases[0]['full_preparation_difference_per_unit_delta_r'])
        exact+=cases
    negative=finite_case(911,F(1),True)
    assert negative['reduced_matches_full'] is False
    case0=finite_case(911,F(0),True)
    assert negative['full_preparation_difference_per_unit_delta_r']!=case0['full_preparation_difference_per_unit_delta_r']
    return dict(round=911,date='2026-10-06',fresh_test_groups=1,cumulative_numbered_groups=3696,
        all_checks_passed=True,full_CAR_preparation_check=car_check(),exact_rational_cases=exact,
        negative_changed_old_flavour_to_identity=negative,
        actual_theorem_scope='Original870 same-two-point preparation difference, flavour-degenerate free receiver, old860 and new869 communication both sigma_y, pure lambda_alpha^2 hbar^2 four-receiver-leg coefficient of complete872 source. The old communication paths vanish after Delta4; identity-flavour stress, full mixed boson propagation, density vertices and compatible initial data remain.',
        original_continuum_total_source_numerically_evaluated=False,
        actual_finite_observable_error_certified=False,full872_numerical_feedback_completed=False,
        autonomous_instrument_completed=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__),STAGE/'846/constrained_six_leg_recursion.py',STAGE/'870/receiver_four_point_backreaction.py')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
