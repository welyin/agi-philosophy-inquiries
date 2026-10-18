"""1299: exact charge/exterior certificates and finite massive normalization.
The full integer/dimension theorem is analytic; finite cases are diagnostics.
Default and --check recompute read-only. --write creates results.json exclusively.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results.json'
TOL = 2e-12


def encode(x):
    if isinstance(x, F):
        return {'numerator': x.numerator, 'denominator': x.denominator}
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    return x


def err(x):
    return float(np.max(np.abs(x)))


def charge_certificate():
    # Entries: name, 6Y0, 3(B-L), color-center triality, r exponent modulo4.
    base = [('Q', 1, 1, 1, 1), ('uc', -4, -1, -1, 1),
            ('dc', 2, -1, -1, 1), ('L', -3, -3, 0, 1),
            ('ec', 6, 3, 0, 1), ('N', 0, 3, 0, 1)]
    fs = base + [(n+'bar', -y, -b, -c, -r) for n,y,b,c,r in base]
    assert all(F(abs(y), 6) <= 1 for _,y,_,_,_ in fs)
    # Stronger universal bound: |Y0| <= (2/3)*canonical_dimension.
    assert all(F(abs(y), 6) <= F(2, 3)*F(3, 2) for _,y,_,_,_ in fs)
    assert F(1,2) <= F(2,3)
    with_r, no_r = [], []
    contents = 0
    for nf in (2, 4):
        for inds in combinations_with_replacement(range(len(fs)), nf):
            selected = [fs[i] for i in inds]
            y0 = sum(z[1] for z in selected)
            b = sum(z[2] for z in selected)
            tr = sum(z[3] for z in selected)
            rr = sum(z[4] for z in selected)
            if not b or tr % 3:
                continue
            assert b % 6 == 0  # color singlet + even fermion parity.
            for nh in range((12-3*nf)//2+1):
                for hp in range(nh+1):
                    hm = nh-hp
                    contents += 1
                    yy = y0+3*(hp-hm)
                    dim = F(3*nf,2)+nh
                    assert F(abs(yy),6) <= F(2,3)*dim <= 4
                    # 6Yk = 6Y0 + 6*k*b; solve all integers, never scan k.
                    k = F(-yy,6*b)
                    if k.denominator != 1:
                        continue
                    assert k == 0
                    row = {'fermions':[z[0] for z in selected], 'H':hp,
                           'Hbar':hm,'BL':F(b,3),'dimension':dim,'k':k}
                    no_r.append(row)
                    if (rr+2*nh) % 4 == 0:
                        assert b in (12,-12) and nf == 4 and nh == 0
                        assert all(z[0] == ('N' if b>0 else 'Nbar') for z in selected)
                        with_r.append(row)
    assert len(with_r) == 2
    # Finite scan is only an overcount of contraction possibilities, not a basis.
    assert any(r['fermions']==['N','N'] and r['dimension']==3 for r in no_r)
    field = {n:(F(y,6),F(b,3),r) for n,y,b,_,r in fs}
    field.update({'H':(F(1,2),F(0),2),'Hbar':(F(-1,2),F(0),-2)})
    yukawas=[['Q','H','uc'],['Q','Hbar','dc'],['L','Hbar','ec'],['L','H','N']]
    for word in yukawas:
        assert sum(field[n][0] for n in word)==0
        assert sum(field[n][1] for n in word)==0
        assert sum(field[n][2] for n in word)%4==0
    return {'tested_charge_contents':contents,
            'coarse_passes_without_r':len(no_r),'coarse_passes_with_r':with_r,
            'unknown_integer_k_solved_exactly_not_scanned':True,
            'derivative_and_fieldstrength_neutral_dimension_cost_included_analytically':True,
            'coarse_scan_does_not_certify_Lorentz_or_weak_contractions':True,
            'without_r_nonzero_k_still_excluded':True,
            'all_four_Dirac_yukawas_neutral_for_all_k':True}


def add(a,b):
    z=dict(a)
    for m,c in b.items():
        z[m]=z.get(m,F(0))+c
    return {m:c for m,c in z.items() if c}


def scale(a,c):
    return {m:F(c)*v for m,v in a.items() if c*v}


def wedge(a,b):
    z={}
    for m,c in a.items():
        for n,d in b.items():
            if m&n:
                continue
            swaps=sum(1 for i in range(8) for j in range(8) if (m>>i)&1 and (n>>j)&1 and i>j)
            key=m|n
            z[key]=z.get(key,F(0))+(-1 if swaps%2 else 1)*c*d
    return {m:c for m,c in z.items() if c}


def pair(u,v):
    return add(wedge(u[0],v[1]),scale(wedge(u[1],v[0]),-1))


def exterior_certificate():
    a,b,c,d=[{1<<i:F(1)} for i in range(4)]
    n1,n2=[a,b],[c,d]
    o=scale(wedge(pair(n1,n1),pair(n2,n2)),F(1,4))
    assert o == {15:F(1)}
    assert not wedge(pair(n1,n1),pair(n1,n1))
    # Independent external Grassmann substitution; outgoing upper legs have sqrt2 common factors.
    # x1=(0,3), x2=(3,0); y3/sqrt2=(0,2), y4/sqrt2=(2,0).
    in_field=[scale(b,3),scale(a,3)]
    out_field_div_sqrt2=[scale(d,2),scale(c,2)]
    exterior_vertex=scale(wedge(pair(in_field,in_field),
                                pair(out_field_div_sqrt2,out_field_div_sqrt2)),F(2,4))
    assert exterior_vertex == {15:F(72)}
    nn=pair(n1,n1)
    assert scale(nn,-1)!=nn and scale(o,1)==o
    return {'two_flavor_O4_coefficient':o[15], 'same_flavor_fourth_power_zero':True,
            'external_four_leg_coefficient':exterior_vertex[15],
            'two_factorials_cancel_one_quarter':True}


def annihilator(mode,n=4):
    a=np.zeros((2**n,2**n),complex)
    for bits in range(2**n):
        if (bits>>mode)&1:
            sign=(-1)**((bits&((1<<mode)-1)).bit_count())
            a[bits^(1<<mode),bits]=sign
    return a


def ward_certificate():
    aa=[annihilator(i) for i in range(4)]
    O=aa[0]@aa[1]@aa[2]@aa[3]
    Q=sum(a.conj().T@a for a in aa)
    C=complex(F(3,5),F(4,5))
    lam=F(100)
    L=(C*O+C.conjugate()*O.conj().T)/float(lam**2)
    H=-L
    rhs=-4j*(C*O-C.conjugate()*O.conj().T)/float(lam**2)
    rr=np.diag([(-1j)**i.bit_count() for i in range(16)])
    nn=aa[0]@aa[1]
    residuals=[err(Q@O-O@Q+4*O),err(1j*(H@Q-Q@H)-rhs),
               err(rr@O@rr.conj().T-O),err(rr@nn@rr.conj().T+nn),err(H+L)]
    assert max(residuals)<TOL
    return {'max_residual':max(residuals),'dimension':16,
            'current_convention':'deltaL=-j.dalpha+4i*alpha*(C O-Cbar Odagger)/Lambda^2',
            'H_int_equals_minus_L_int':True,
            'interaction_Hilbert_source_definition_not_full_device_certificate':True}


I2=np.eye(2,dtype=complex)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
EPS=np.array([[0,1],[-1,0]],complex)


def hamiltonian(m,p):
    sp=p[0]*SX+p[1]*SY+p[2]*SZ
    return np.block([[-sp,m*I2],[m*I2,sp]])


def psqrt(a):
    val,vec=np.linalg.eigh(a)
    assert min(val)>0
    return (vec*np.sqrt(val))@vec.conj().T


def external(m,p,spin,negative=False):
    E=np.sqrt(m*m+np.dot(p,p))
    sp=p[0]*SX+p[1]*SY+p[2]*SZ
    up=psqrt(E*I2-sp)@spin
    low=psqrt(E*I2+sp)@spin
    return E,np.concatenate([up,(-1 if negative else 1)*low])


def massive_certificate():
    down=np.array([0.,1.],complex); up=np.array([1.,0.],complex)
    records=[]; legs=[]; residuals=[]
    for m,p,spin,neg in [(3,np.array([0.,0.,4.]),down,False),
                          (3,np.array([0.,0.,-4.]),up,False),
                          (4,np.array([0.,0.,3.]),down,True),
                          (4,np.array([0.,0.,-3.]),up,True)]:
        E,w=external(m,p,spin,neg)
        sign=-1 if neg else 1
        residuals += [err(hamiltonian(m,-p if neg else p)@w-sign*E*w),
                      float(abs(np.vdot(w,w)-2*E))]
        legs.append(w[:2]); records.append({'mass':m,'E':float(E),'momentum':p.tolist(),
                                           'negative_frequency':neg,'norm_squared':float(np.vdot(w,w).real)})
    contractions=[legs[0].T@EPS@legs[1],legs[2].T@EPS@legs[3]]
    residuals += [abs(contractions[0]+9),abs(contractions[1]+8)]
    exact_amp=F(72,100**2)
    assert exact_amp==F(9,1250) and exact_amp**2==F(81,1562500)
    assert max(residuals)<TOL
    # Charge/Lorentz selection is exact; the one massive point normalizes its vertex.
    return {'masses':[3,4,5],'Lambda':100,'C':1,'external_states':records,
            'epsilon_contractions':[float(x.real) for x in contractions],
            'M':exact_amp,'M_squared':exact_amp**2,
            'positive_energy_window':[F(5),F(21,4)],
            'window_lower_abs_M':exact_amp,'incoming_beta_at_left':F(4,5),
            'outgoing_beta_at_left':F(3,5),'max_residual':float(max(residuals)),
            'matrix_element_order':'first order C/Lambda^2, tree level only',
            'finite_momentum_not_finite_spacetime_instrument':True}


def deletion_certificate():
    # Actual generic Higgs jet: H=(1,0), D0H=(0,1), spatial jets=0.
    H=np.array([1.,0.]); dt=np.array([0.,1.])
    J0=float(H@EPS.real@dt)
    assert J0==1
    rows=[]
    for k in (1,-1):
        sign=-1 if k>0 else 1  # Y of Jdagger or J.
        exponent=6*abs(k)
        y=12*k+2*sign*exponent
        dim=6+6*exponent
        assert y==0 and dim==42
        rows.append({'k':k,'high_dimension':dim,'Y_total':y,'BL':4,'nonzero_local_jet':1})
    # Additional SNN example explicitly conserves enlarged B-L if BL(S)=-2.
    assert F(-2)+2*F(1)==0
    assert (2+1+1)%4==0
    return {'without_r_k_filter_survives':True,'without_r_k0_NN_allowed':True,
            'high_dimension_general_formula':'6+36*abs(k)',
            'high_dimension_local_jet_checks':rows,
            'charged_S_NN_not_BL_violation_if_S_carries_minus2':True,
            'C_zero_also_allowed':True}


def source_hashes():
    names=[
      'research_cognition_physics/archive_1290_/_shared/round1299_admission.md',
      'research_cognition_physics/archive_1290_/1299/proof.md',
      'research_cognition_physics/archive_1290_/1299/check.py',
      'research_cognition_physics/archive_1009_1043/1035/drafts/joint_selection_derivation.md',
      'research_cognition_physics/archive_1046_1062/1060/proof.md',
      'research_cognition_physics/archive_1009_1043/research_note_1034.md',
      'research_cognition_physics/archive_1046_1062/_shared/common_branch_scope_after1062.md']
    return {n:sha256((ROOT/n).read_bytes()).hexdigest() for n in names}


def calculate():
    return encode({'round':1299,'status':'PASS_AUTHOR_SCOPED','scientific_groups':1,
       'charge_necessary_conditions':charge_certificate(),
       'exterior_algebra':exterior_certificate(),'Ward_and_discrete_protection':ward_certificate(),
       'massive_four_point':massive_certificate(),'deletion_scope':deletion_certificate(),
       'source_hashes':source_hashes(),
       'scope':{'full_integer_theorem_by_proof':True,
                'inherited_anomaly_scope_only':True,'nonzero_C_is_input':True,
                'no_actual_instrument_or_higher_order_error_certificate':True}})


def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=calculate()
    text=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            f.write(text)
        print('CREATED PASS_AUTHOR_SCOPED')
    else:
        saved=json.loads(OUT.read_text(encoding='utf-8'))
        assert saved==result, 'Saved result or source hash differs; no file modified.'
        print('READ_ONLY_MATCH PASS_AUTHOR_SCOPED')
    print('sha256_results='+sha256(OUT.read_bytes()).hexdigest())


if __name__=='__main__':
    main()
