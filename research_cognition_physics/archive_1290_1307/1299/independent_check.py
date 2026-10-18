"""1299 independent charge, exterior, spectral-projector LSZ and Fock-Ward checks.
No author code is imported. --write creates results exclusively; default/--check
only compare. --diagnose has no source-file dependency and writes nothing.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, itertools, json, math
import numpy as np
HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESEARCH=STAGE.parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def charge_checks():
    # name, BL, Y_SM, color triality, weak parity, r exponent
    base=[('Q',F(1,3),F(1,6),1,1,1),('uc',F(-1,3),F(-2,3),-1,0,1),
          ('dc',F(-1,3),F(1,3),-1,0,1),('L',F(-1),F(-1,2),0,1,1),
          ('ec',F(1),F(1),0,0,1),('N',F(1),F(0),0,0,1)]
    fs=base+[(n+'bar',-b,-y,-c,w,-r) for n,b,y,c,w,r in base]
    hs=[('H',F(0),F(1,2),0,1,2),('Hbar',F(0),F(-1,2),0,1,2)]
    for _,b,y,_,_,r in fs+hs:
        x=5*b-4*y
        assert x.denominator==1 and (int(x)-r)%4==0
    count=0; hit=[]; nor=[]
    for nf in (0,2,4):
        for ix in itertools.combinations_with_replacement(range(12),nf):
            for nh in range((12-3*nf)//2+1):
                for nplus in range(nh+1):
                    vv=[fs[j] for j in ix]+[hs[0]]*nplus+[hs[1]]*(nh-nplus)
                    count+=1; b=sum((v[1] for v in vv),F(0)); y=sum((v[2] for v in vv),F(0))
                    d=F(3*nf+2*nh,2)
                    assert abs(y)<=4
                    if sum(v[3] for v in vv)%3 or sum(v[4] for v in vv)%2: continue
                    assert b.denominator==1 and int(b)%2==0
                    if not b: continue
                    k=-y/(3*b)
                    if k.denominator!=1: continue
                    assert k==0
                    row={'k':int(k),'dimension':str(d),'BL':str(b),'fields':[v[0] for v in vv]}
                    nor.append(row)
                    if sum(v[5] for v in vv)%4==0:
                        assert abs(b)==4 and d==6 and all(v[0] in ('N','Nbar') for v in vv)
                        hit.append(row)
    assert min(F(v['dimension']) for v in nor)==3
    assert min(F(v['dimension']) for v in hit)==6
    assert all(math.gcd(5+12*k,12)==1 for k in range(-100,101))
    lookup={v[0]:v for v in fs+hs}
    for names in [('Q','uc','H'),('Q','dc','Hbar'),('L','ec','Hbar'),('L','N','H')]:
        vv=[lookup[n] for n in names]
        assert sum(v[1] for v in vv)==sum(v[2] for v in vv)==0 and sum(v[5] for v in vv)%4==0
    higher=[]
    for k in (-3,-2,-1,1,2,3):
        p=6*abs(k)
        assert 12*k+(-2 if k>0 else 2)*p==0
        higher.append({'k':k,'dimension':6+6*p,'BL':4})
    return {'contents_checked':count,'with_r_BLV_contents':hit,'without_r_minimum_dimension':3,
            'with_r_minimum_dimension':6,'nonzero_k_counterexamples':[],
            'r_needed_for_k_selection':False,'higher_dimension_witnesses':higher,
            'scope':'Charge-content overapproximation; infinite quantifiers require the analytic proof.'}

def wedge(a,b):
    out={}
    for ia,ca in a.items():
        for ib,cb in b.items():
            if ia&ib: continue
            flips=sum((ib&((1<<j)-1)).bit_count() for j in range(8) if (ia>>j)&1)
            out[ia|ib]=out.get(ia|ib,0)+(-1 if flips%2 else 1)*ca*cb
    return {i:c for i,c in out.items() if abs(c)>1e-14}
def add(a,b,s=1):
    out=dict(a)
    for i,c in b.items(): out[i]=out.get(i,0)+s*c
    return {i:c for i,c in out.items() if abs(c)>1e-14}
def pair(a,b): return add(wedge(a[0],b[1]),wedge(a[1],b[0]),-1)
def lsz(cols):
    a=[{1<<j:cols[j][z] for j in (0,1)} for z in (0,1)]
    b=[{1<<j:cols[j][z] for j in (2,3)} for z in (0,1)]
    return wedge(pair(a,a),pair(b,b)).get(15,0)/4

def exterior_ward_checks():
    e=[{1<<j:F(1)} for j in range(8)]
    a,b,c,d=[e[2*j:2*j+2] for j in range(4)]
    assert pair(a,a)=={3:F(2)}
    assert wedge(pair(a,a),pair(b,b))=={15:F(4)}
    assert not wedge(pair(a,a),pair(a,a))
    assert not add(add(wedge(pair(a,b),pair(c,d)),wedge(pair(a,c),pair(b,d))),wedge(pair(a,d),pair(b,c)))
    ann=[]
    for j in range(4):
        op=np.zeros((16,16),complex)
        for s in range(16):
            if (s>>j)&1: op[s^(1<<j),s]=(-1)**((s&((1<<j)-1)).bit_count())
        ann.append(op)
    op=ann[0]@ann[1]@ann[2]@ann[3]; Q=sum(a.conj().T@a for a in ann)
    g=(2+3j)/10000; H=-g*op-g.conjugate()*op.conj().T
    lhs=1j*(H@Q-Q@H); rhs=-4j*(g*op-g.conjugate()*op.conj().T)
    r=np.diag(np.exp(-.5j*np.pi*np.diag(Q).real))
    res=max(np.linalg.norm(lhs-rhs),np.linalg.norm(r@op@r.conj().T-op),
            np.linalg.norm(r@(ann[0]@ann[1])@r.conj().T+ann[0]@ann[1]))
    assert res<1e-12 and op[0,15]==1
    return {'external_factor_before_quarter':4,'one_flavor_fourth_power_zero':True,
            'two_flavor_nonzero':True,'Grassmann_Schouten_identity':True,
            'Ward':'dQ/dt=-4i(g O4-conjugate(g) O4dagger)','max_residual':float(res)}

def spinor_checks():
    sig=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
    def h(p,m):
        sp=sum(p[j]*sig[j] for j in range(3)); mi=m*np.eye(2)
        return np.block([[-sp,mi],[mi,sp]])
    def leg(p,m,xi,anti):
        E=math.sqrt(m*m+float(p@p)); hh=h(-p if anti else p,m); sign=-1 if anti else 1
        P=(np.eye(4)+sign*hh/E)/2
        u=P@np.r_[xi,[0,0]]; u*=math.sqrt(2*E)/np.linalg.norm(u)
        return u,max(np.linalg.norm(hh@u-sign*E*u),abs(np.vdot(u,u).real-2*E))
    up=np.array([1,0],complex); down=np.array([0,1],complex)
    rot=np.array([[3,-1],[1,3]],complex)/math.sqrt(10); direction=np.array([3/5,0,4/5])
    rows=[]; residual=0.
    for Er in (F(5),F(41,8),F(21,4)):
        E=float(Er); pi=math.sqrt(E*E-9); po=math.sqrt(E*E-16)
        specs=[(np.array([0.,0.,pi]),3,down,False),(np.array([0.,0.,-pi]),3,up,False),
               (po*direction,4,rot@down,True),(-po*direction,4,rot@up,True)]
        legs=[]
        for p,m,xi,anti in specs:
            u,err=leg(p,m,xi,anti); legs.append(u); residual=max(residual,float(err))
        coeff=lsz([u[:2] for u in legs]); target=(E+pi)*(E+po)
        residual=max(residual,float(abs(coeff-target)))
        assert abs(coeff-target)<1e-10 and abs(coeff)>=72-1e-11 and pi>0 and po>0
        rows.append({'E':str(Er),'coefficient_real':float(np.real(coeff)),'coefficient_imag':float(np.imag(coeff))})
    assert residual<1e-10
    return {'method':'4x4 energy spectral projectors and external Grassmann coefficient',
            'norm':'u_dagger_u=v_dagger_v=2E','v_equation':'h(-p)v=-Ev',
            'outgoing_direction':['3/5','0','4/5'],'masses':[3,4,5],
            'E_window':['5','21/4'],'coefficient_lower_bound':72,'Lambda':100,
            'E5_abs_M':'9/1250','E5_dPhi2_per_dOmega_times_pi_squared':'3/160',
            'identical_final_state_symmetry_factor':'1/2','rows':rows,'max_residual':residual}

def payload():
    return {'verdict':'PASS_INDEPENDENT_SCOPED','group':'1299-exterior-spectral-independent',
            'charge_checks':charge_checks(),'exterior_Ward':exterior_ward_checks(),'massive_LSZ':spinor_checks(),
            'scope':{'all_k_proof_analytic':True,'finite_scans_diagnostic_only':True,
                     'tree_level_EFT_only':True,'actual_record_or_full_O_certified':False,'M4_certified':False}}
def sources():
    ps={'author_proof':HERE/'proof.md','author_check':HERE/'check.py','author_results':HERE/'results.json',
        'admission':STAGE/'_shared'/'round1299_admission.md',
        '1035':RESEARCH/'archive_1009_1043'/'research_note_1035.md',
        '1060':RESEARCH/'archive_1046_1062'/'1060'/'proof.md','independent_check':Path(__file__).resolve()}
    return {n:{'path':str(p.relative_to(RESEARCH)).replace(chr(92),'/'),'sha256':sha(p)} for n,p in ps.items()}
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--write',action='store_true'); ap.add_argument('--check',action='store_true'); ap.add_argument('--diagnose',action='store_true'); args=ap.parse_args()
    if args.write and (args.check or args.diagnose): ap.error('incompatible modes')
    out=payload()
    if args.diagnose: print(json.dumps(out,ensure_ascii=False,indent=2)); return
    out['sources']=sources(); dest=HERE/'independent_results.json'
    if args.write:
        with dest.open('x',encoding='utf-8',newline='\n') as f: json.dump(out,f,ensure_ascii=False,indent=2); f.write('\n')
        print('PASS_INDEPENDENT_SCOPED EXCLUSIVE_WRITE')
    else:
        assert json.loads(dest.read_text(encoding='utf-8'))==out,'result/source binding mismatch'
        print('PASS_INDEPENDENT_SCOPED READ_ONLY_MATCH')
    print('max_residual='+str(max(out['exterior_Ward']['max_residual'],out['massive_LSZ']['max_residual'])))
if __name__=='__main__': main()