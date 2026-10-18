"""Count-zero continuation after1299. Exact finite diagnostics; general m is proved in Markdown.
Default/--check read-only; --write exclusively creates its own result, never 1299/results.json.
"""
from fractions import Fraction as F
from math import comb, gcd
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=HERE/'global_form_corollary_results.json'
# q=q0+m*q1, color-center power, weak-center power, B-L.
FIELDS={
 'Q':(1,1,1,1,F(1,3)), 'uc':(-4,-1,-1,0,F(-1,3)),
 'dc':(2,-1,-1,0,F(-1,3)), 'L':(-3,-3,0,1,F(-1)),
 'ec':(6,3,0,0,F(1)), 'N':(0,3,0,0,F(1)),
 'H':(3,0,0,1,F(0))}


def enc(x):
    if isinstance(x,F):return {'numerator':x.numerator,'denominator':x.denominator}
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    return x


def charges(m):return {n:a+b*m for n,(a,b,_,_,_) in FIELDS.items()}


def phase(n,m,c,w,t):
    q0,q1,ct,wt,_=FIELDS[n]
    return ct*c+wt*w+(q0+q1*m)*t


def calculate():
    # Exact affine identities, not a scan over integer m.
    for n,(q0,q1,ct,wt,bl) in FIELDS.items():
        assert q1==3*bl
        X0=5*bl-F(2,3)*q0
        assert X0.denominator==1 and int(X0)%4==(2 if n=='H' else 1)
    assert all(gcd(5+2*j,4)==1 for j in range(4))
    types={'N':(0,0),'Q':(1,1),'uc':(2,0),'ec':(0,2),'dc':(2,2),'L':(3,1)}
    weight_rows=[]
    for n,(nc,nw) in types.items():
        affine=(-2*nc+3*nw,3-2*nc)
        assert affine==FIELDS[n][:2]
        weight_rows.append({'field':n,'occupancy':[nc,nw],
                            'q_affine':affine,'multiplicity':comb(3,nc)*comb(2,nw)})
    assert sum(r['multiplicity'] for r in weight_rows)==16
    assert 3*(-2)+2*3==0 and 3*(-2)==-6  # sum(v_m)=-6m, always even.
    # The complete finite four-lepton count follows dimension saturation analytically.
    contents=[]
    for e in range(5):
        for l in range(5-e):
            n=4-e-l
            m=F(-2*e-l,4)
            if m.denominator==1:
                assert int(m) in (0,-1,-2) and l%2==0
                contents.append({'ec':e,'Lbar':l,'N':n,'m':int(m)})
    assert {r['m'] for r in contents}=={-2,-1,0}
    kernel_rows=[]
    for m in (-2,-1,0):
        js=[]; representatives=[]
        for j in range(6):
            t=F(j,6); c=(2+2*m)*t; w=-3*t
            assert (3*c).denominator==1 and (2*w).denominator==1
            phases={n:phase(n,m,c,w,t) for n in FIELDS}
            assert all((v-(0 if n=='H' else F(m*j,2))).denominator==1
                       for n,v in phases.items())
            is_kernel=all(v.denominator==1 for v in phases.values())
            assert is_kernel==(m*j%2==0)
            if is_kernel:
                js.append(j)
                representatives.append({'j':j,'C_phase':c%1,'W_phase':w%1,'U1_phase':t})
        assert len(js)==(6 if m%2==0 else 3)
        assert all((a+b)%6 in js for a in js for b in js)
        if m==-1:
            assert all(r['C_phase']==0 and r['W_phase']==0 for r in representatives)
            assert [charges(m)[n]//3 for n in ['Q','uc','dc','L','ec','N','H']]==[0,-1,1,0,1,-1,1]
        # Fixed original z and the new m-dependent maximal common kernel are different.
        oldz={n:phase(n,m,F(1,3),F(1,2),F(1,6)) for n in FIELDS}
        assert (oldz['Q']-F(m,6)).denominator==1
        assert all(v.denominator==1 for v in oldz.values())==(m%6==0)
        kernel_rows.append({'m':m,'kernel_j':js,'representatives':representatives,
                            'descends_to_original_fixed_G6':m%6==0})
    # Odd-m gauge fermion parity: inspect retained m=-1 representative j=1, not in ordinary kernel.
    ff={n:phase(n,-1,F(0),F(-1,2),F(1,6))%1 for n in FIELDS}
    assert ff['H']==0 and all(v==F(1,2) for n,v in ff.items() if n!='H')
    # Conditional physical relabeling: hypercharge inversion plus singlet exchanges, H conjugation.
    qminus=charges(-2)
    rename={'Q':'Q','uc':'dc','dc':'uc','L':'L','ec':'N','N':'ec'}
    mapped={n:-qminus[old] for n,old in rename.items()}
    assert mapped=={n:charges(0)[n] for n in rename}
    assert -(-qminus['H'])==3
    names=[
      'research_cognition_physics/archive_1290_/_shared/global_form_corollary_after1299.md',
      'research_cognition_physics/archive_1290_/1299/global_form_corollary_check.py',
      'research_cognition_physics/archive_1290_/1299/proof.md',
      'research_cognition_physics/archive_1290_/1299/check.py',
      'research_cognition_physics/archive_1290_/1299/results.json',
      'research_cognition_physics/archive_1009_1043/1035/drafts/joint_selection_derivation.md',
      'research_cognition_physics/archive_1046_1062/1060/proof.md',
      'research_cognition_physics/archive_1046_1062/1060/independent_review.md',
      'research_cognition_physics/archive_1009_1043/research_note_1021.md']
    return enc({'status':'PASS_COUNT_ZERO_SCOPED','new_scientific_groups':0,
       'new_axioms':0,'charge_affine_weights':weight_rows,'four_lepton_contents':contents,
       'finite_kernel_representatives':kernel_rows,'odd_m_gauge_fermion_parity_example':ff,
       'conditional_m_minus2_relabel':mapped,
       'scope':{'full_integer_claim_is_analytic':True,
                'maximal_faithful_ordinary_group_not_unique_physical_global_form':True,
                'Spin_total_group_quotient_not_silently_selected':True,
                'm_minus1_survives':True,'new_instrument_or_amplitude':False},
       'source_hashes':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}})


def main():
    ap=argparse.ArgumentParser(); mode=ap.add_mutually_exclusive_group()
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=ap.parse_args();r=calculate()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            f.write(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
        print('CREATED PASS_COUNT_ZERO_SCOPED')
    else:
        assert json.loads(OUT.read_text(encoding='utf-8'))==r,'Result/source mismatch; no files modified.'
        print('READ_ONLY_MATCH PASS_COUNT_ZERO_SCOPED')
    print(hashlib.sha256(OUT.read_bytes()).hexdigest())


if __name__=='__main__':main()
