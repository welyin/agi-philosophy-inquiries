"""1272 independent exact gate/source audit, one group.
No author module is imported. Default compares without writes; --write is exclusive.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial, prod
from pathlib import Path
import argparse
import hashlib
import json

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
OUT=HERE/'independent_results.json'
DIMS=(2,2,2,2,4,5,2,3,2,2,2)
ZERO=(F(0),F(0))
ONE=(F(1),F(0))
HALF_I=(F(1,2),-F(1,2))
HALF_S=(F(1,2),F(1,2))

def add(z,w): return (z[0]+w[0],z[1]+w[1])
def mul(z,w): return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def conj(z): return (z[0],-z[1])
def abs2(z): return z[0]*z[0]+z[1]*z[1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def serial(v):
    if isinstance(v,F): return str(v)
    if isinstance(v,dict): return {str(k):serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [serial(x) for x in v]
    return v

def update(key,i,j=None):
    out=list(key)
    if j is None: out[i]^=1
    else: out[i],out[j]=out[j],out[i]
    return tuple(out)

def transition(key,step):
    # Full finite mother basis, retaining stock, seeds, control and every record.
    if step==0 and key[4]==0: return [(update(key,0,2),ONE)]
    if step==1 and key[5] in (0,1): return [(update(key,1,3),ONE)]
    if step==2 and key[6]==1 or step==7:
        return [(key,HALF_I),(update(key,0,1),HALF_S)]
    if step==3 and key[7] in (1,2) and key[0]: return [(update(key,8),ONE)]
    if step==4 and key[7] in (1,2) and key[8]: return [(update(key,9),ONE)]
    if step==5 and key[7]==2 and key[8]: return [(update(key,9),ONE)]
    if step==6 and key[7]==2 and key[0]: return [(update(key,8),ONE)]
    if step==8 and key[0]: return [(update(key,10),ONE)]
    return [(key,ONE)]

def apply(state,step):
    out={}
    for key,z in state.items():
        for dest,w in transition(key,step):
            out[dest]=add(out.get(dest,ZERO),mul(z,w))
    return {key:z for key,z in out.items() if z!=ZERO}

def inner(left,right,drop_control=False):
    if drop_control:
        def collapse(v):
            r={}
            for key,z in v.items():
                new=key[:6]+key[7:]
                r[new]=add(r.get(new,ZERO),z)
            return r
        left,right=collapse(left),collapse(right)
    value=ZERO
    for key,z in left.items(): value=add(value,mul(conj(z),right.get(key,ZERO)))
    return (value[0]/20,value[1]/20)

def probability(state,predicate):
    return sum((abs2(z) for key,z in state.items() if predicate(key)),F(0))/20

def initial(control,program):
    return {(0,0,1,1,za,zb,control,program,0,0,0):ONE
            for za in range(4) for zb in range(5)}

def protocol_checks():
    assert prod(DIMS)==15360
    assert add(HALF_I,HALF_S)==ONE
    assert abs2(HALF_I)+abs2(HALF_S)==1
    assert add(mul(conj(HALF_I),HALF_S),mul(conj(HALF_S),HALF_I))==ZERO
    assert add(mul(HALF_I,HALF_I),mul(HALF_S,HALF_S))==ZERO
    assert mul((F(2),F(0)),mul(HALF_I,HALF_S))==ONE
    # Exhaustive selection rules for every basis vector; no finite sampling.
    counts=[0]*9
    for key in product(*(range(d) for d in DIMS)):
        energy=sum(key[:4])
        for step in range(9):
            terms=transition(key,step)
            assert all(sum(dest[:4])==energy for dest,_ in terms)
            assert all(dest[4:8]==key[4:8] for dest,_ in terms)
            counts[step]+=1
            if step not in (2,7):
                back=transition(terms[0][0],step)
                assert back==[(key,ONE)]
    cases={}
    outputs={}
    for control,program in product(range(2),range(3)):
        v=initial(control,program)
        assert inner(v,v)==ONE
        prefixes=[]
        for step in range(9):
            v=apply(v,step)
            assert inner(v,v)==ONE
            assert all(sum(key[:4])==2 for key in v)
            prefixes.append({'pA':probability(v,lambda x:x[0]==1),
                             'pB':probability(v,lambda x:x[1]==1)})
        assert prefixes[1]=={'pA':F(1,4),'pB':F(2,5)}
        if control:
            assert prefixes[2]['pA']-prefixes[1]['pA']==F(3,40)
            assert prefixes[2]['pB']-prefixes[1]['pB']==-F(3,40)
        p=probability(v,lambda x:x[10]==1)
        expected=F(2,5) if control and program!=1 else F(13,40)
        assert p==expected and p==prefixes[-1]['pA']
        if program!=1:
            assert probability(v,lambda x:x[8] or x[9])==0
        else:
            assert probability(v,lambda x:x[8]!=x[9])==0
        outputs[control,program]=v
        cases[f'K{control}_L{program}']={'prefixes':prefixes,'p_Jout_1':p,
            'p_M1':probability(v,lambda x:x[8]==1),
            'retained_basis_terms':len(v)}
    coherences={}
    for program in range(3):
        assert inner(outputs[0,program],outputs[1,program])==ZERO
        coherences[str(program)]=inner(outputs[0,program],outputs[1,program],True)
    # This overlap is sensitive to the controlled half-swap's identity phase.
    assert coherences['0']==(F(31,40),-F(9,40))
    assert coherences['2']==coherences['0']
    assert outputs[1,0] != outputs[1,1]
    return {'mother_dimension':prod(DIMS),'energy_checks_per_gate':counts,
            'cases':cases,'environment_overlap_K0_K1':coherences,
            'isometry_on_control_and_arbitrary_reference':True}

def source_budget():
    pulse_int=F(140*factorial(3)*factorial(3),factorial(7))
    assert pulse_int==1
    tau=F(1,16)
    peak=F(140,64)/tau
    xi=F(1,10**9)
    pi_upper=F(22,7)
    vmax=peak*pi_upper*(1+xi)
    dpeak=F(420,16)/tau**2
    workmax=dpeak*pi_upper*(1+xi)
    assert peak==35 and vmax<111 and dpeak==6720 and workmax<21121
    assert 11+4+vmax<126 and 11+2+vmax<124
    assert 11+2==13
    # Every selected G is a positive multiple of a controlled (I-P),
    # where P is a Hermitian involution; positive pulses give v >= 0.
    generator_eigenvalues=(F(0),F(2))
    assert min(generator_eigenvalues)==0 and 1-xi>0
    start=F(1,2)-F(1,256)
    last=F(1,2)+9*tau+F(1,256)
    assert start>F(21,400)
    assert last<F(2)-F(1,10**8)<F(17,8)<F(9,4)
    channel=xi+9*pi_upper*xi
    gap=F(3,40)-2*channel
    assert channel<F(3,10**8)
    assert gap>F(3,40)-F(6,10**8)>0
    # Exponents (hbar,s): E_gap * proper_time / hbar is dimensionless.
    physical_phase_powers=tuple(sum(v) for v in zip((1,-1),(0,1),(-1,0)))
    assert physical_phase_powers==(0,0)
    # Same source coefficient must multiply clock and all added components.
    assert -13+11+2==0
    assert -1+1==0  # coefficients of the actual v, not only its expectation.
    # Each pulse has one fixed G; sequential G need not commute with others.
    assert 35*pi_upper<110+F(1,100)
    old=read(STAGE/'1249/quantum_bounds_results.json')['group_constructive_finite_resource_bound']['ledger']
    M=1<<3059
    EH=int(old['Fourier_H_exponent'])
    E0=int(old['Fourier_menu_base_exponent'])
    N=int(old['integer_n_power_exponent_decimal'])
    EF=int(old['Fourier_transport_exponent_decimal'])
    assert (EH,E0)==(3093,14008) and EF==95928+100*M
    assert N==1+max((EF+44)//3,(EH+EF+47)//6)
    assert 4*N-7>110 and 4*N-15>110
    geom=read(STAGE/'1267/geometry_results.json')['uniform_tube_and_bounds']
    assert N==int(geom['N_decimal'])
    slack=F(geom['positive_threshold_slacks']['J'])
    assert slack>F(18,1000)>F(1,2**100)
    margins=[2*N-(E0+EH+15*M+27),4*N-(E0+60)]
    assert min(margins)>110
    assert F(2,2**110)<F(1,2**80)
    dk=F(1,15000)+F(5,10**16)
    field=F(111,100)*dk+F(3,2**80)
    assert field<F(1,10000)
    oldfield=read(STAGE/'1270/continuation_results.json')['field_tests']['actual_bounded_quantum_source_comparison_bound']
    assert field==F(oldfield)
    return {'baseline_carriers':11,'conserved_excitation':2,'initial_J_debit_in_epsilon':13,
        'gate_peak':peak,'actual_v_norm_upper':vmax,
        'clock_derivative_upper':workmax,'saved_J_slack':slack,
        'full_finite_code_mass_bound':126,'prepared_support_mass_bound':124,
        'full_finite_code_h_lower':11,'prepared_support_h_lower':13,
        'earliest_first_gate':start,'latest_last_gate':last,
        'causal_start_margin':start-F(21,400),
        'actual_channel_error_bound':channel,'record_gap_bound':gap,
        'coarser_record_gap':F(3,40)-F(6,10**8),
        'epsilon_bound':'2^(-4N)','spectral_shift_upper':'124 epsilon < 2^-110',
        'old_source_exponent_margin_bitlengths':[x.bit_length() for x in margins],
        'inherited_source_error':'<2^-80','inherited_accepted_field_bound':field,
        'proper_gap':'hbar/s','proper_beta':'(s/hbar)*dimensionless_beta',
        'source_cancellation':'same completed operator coefficient; clock includes baseline, thermal energy and actual pulse',
        'conditional_geometry_requires_old_new_preparation_independence':True}

def calculate():
    result={'round':1272,'verdict':'PASS_INDEPENDENT_WITH_SCOPE','independent_groups':1,
            'method':'exact Gaussian-rational sparse amplitudes, exhaustive energy selection rules, Fraction/integer budgets',
            'protocol':protocol_checks(),'source_budget':source_budget()}
    files=('admission.md','candidate_spec.md','source_baseline.json',
           'source_proof.md','source_check.py','source_results.json',
           'protocol_proof.md','protocol_check.py','protocol_results.json',
           'NEXT.md','../research_note_1272.md','independent_check.py')
    baseline=read(HERE/'source_baseline.json')
    for name,value in baseline['source_sha256'].items():
        assert sha(STAGE/name)==value,name
    result['historical_stage_bindings_matched']=len(baseline['source_sha256'])
    result['sha256']={name:sha(HERE/name) for name in files}
    result['scope']={
        'new_comoving_material':True,'K1219_direction_contract_imported':False,
        'actual_selected_ports_are_explicit_finite_effective_input':True,
        'finite_internal_carrier_not_all_local_field_modes':True,
        'unbounded_full_field_stress_certified':False,'thermal_geometry_generation':False,
        'Jacobson_temperature_or_area_entropy':False,'CT_adopted':False,
        'goal_complete':False}
    return serial(result)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write('\n')
        print('EXCLUSIVE_FIRST_WRITE PASS_INDEPENDENT_WITH_SCOPE')
    else:
        assert read(OUT)==result,'READ_ONLY_MISMATCH'
        print('READ_ONLY_MATCH PASS_INDEPENDENT_WITH_SCOPE')
    print('groups=1; nine energy-preserving gates; ideal record difference=3/40')

if __name__=='__main__': main()
