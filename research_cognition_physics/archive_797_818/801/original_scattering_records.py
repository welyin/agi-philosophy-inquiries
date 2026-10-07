"""801: common old-source/output-record identity and retained joint-state blocks.

Finite exact algebra diagnostics. These matrices are not the original continuum
Einstein-matter coefficients and do not establish a physical readout apparatus.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import json
import sys
sys.dont_write_bytecode=True
import native_scattering_blocks_probe as blocks
a=blocks.a
C=a.C
HERE=Path(__file__).resolve().parent
LIMIT=(3,2,2)  # epsilon, old source z, new late source y
ZERO=(0,0,0)
ONE={ZERO:a.I}
def clean(s):return {k:v for k,v in s.items() if any(v.flat)}
def add(*ss):
    out={}
    for s in ss:
        for k,v in s.items():out[k]=out[k]+v if k in out else v.copy()
    return clean(out)
def sc(s,c):return clean({k:a.scale(v,c) for k,v in s.items()})
def mul(s,t):
    out={}
    for k,x in s.items():
        for l,y in t.items():
            v=tuple(i+j for i,j in zip(k,l))
            if all(i<=j for i,j in zip(v,LIMIT)):
                out[v]=out[v]+x@y if v in out else x@y
    return clean(out)
def adj(s):return clean({k:a.dag(v) for k,v in s.items()})
def same(s,t):return not add(s,sc(t,-1))
def exp(g):
    assert ZERO not in g
    term=out=ONE
    for n in range(1,sum(LIMIT)+1):
        term=mul(term,g)
        out=add(out,sc(term,Q(1,factorial(n))))
    return out
def imag(s):return sc(s,C(0,1))
def zero_source(s,axis):return {k:v for k,v in s.items() if not k[axis]}
def first(s):
    k=min(s)
    i,j=next((i,j) for i in range(2) for j in range(2) if s[k][i,j])
    return dict(epsilon=k[0],old_z=k[1],late_y=k[2],row=i,column=j,value=str(s[k][i,j]))
def joint_menu_probe():
    past=exp(imag({(1,0,0):a.scale(a.Y,Q(1,5))}))
    mid_generator={(1,0,0):a.scale(a.Z,Q(1,3))}
    future=exp(imag({(1,0,0):a.scale(a.X,Q(1,7))}))
    middle=exp(imag(mid_generator))
    total=mul(mul(future,middle),past)
    old_linear={(0,1,0):a.X,(1,1,0):a.scale(a.Y,Q(1,11))}
    old_contact={(2,2,0):a.scale(a.Z,Q(1,13))}
    mid_old=exp(imag(add(mid_generator,old_linear,old_contact)))
    total_old=mul(mul(future,mid_old),past)
    p0=a.proj(2,0)
    late=exp(imag({(0,0,1):p0}))
    old=mul(adj(total),total_old)
    output=mul(mul(adj(total),late),total)
    combined=mul(mul(adj(total),late),total_old)
    assert same(combined,mul(output,old))
    assert same(zero_source(combined,2),old)
    assert same(zero_source(combined,1),output)
    assert same(mul(adj(combined),combined),ONE)
    # The joint identity contains genuinely nonzero mixed coefficients.
    mixed=[k for k in combined if k[1] and k[2]]
    assert mixed
    wrong_old=mul(mul(future,exp(imag(add(mid_generator,old_linear)))),past)
    missing=add(combined,sc(mul(mul(adj(total),late),wrong_old),-1))
    assert missing
    # The same free mode labels at different interaction times are not the same observable.
    early_total=mul(middle,past)
    internal=mul(mul(adj(early_total),{ZERO:p0}),early_total)
    outgoing=mul(mul(adj(total),{ZERO:p0}),total)
    different=add(outgoing,sc(internal,-1))
    assert different
    # Old records/sources are retained by a common output-picture map, without equality to p0.
    to_out=lambda x:mul(mul(total,x),adj(total))
    from_out=lambda x:mul(mul(adj(total),x),total)
    for x in (old,internal):assert same(from_out(to_out(x)),x)
    assert same(to_out(outgoing),{ZERO:p0})
    assert not same(to_out(internal),{ZERO:p0})
    assert same(to_out(mul(outgoing,internal)),mul(to_out(outgoing),to_out(internal)))
    return dict(all_checks_passed=True,epsilon_order=3,old_source_order=2,late_source_order=2,
                full_joint_generating_family_preserved=True,
                setting_new_source_to_zero_recovers_old_family=True,
                joint_family_unitarity=True,nonzero_mixed_coefficients=len(mixed),
                omitted_old_quantum_contact_defect=first(missing),
                identifying_internal_and_outgoing_records_defect=first(different),
                original_records_and_ordered_products_recovered_in_output_picture=True)
def run():
    return dict(round=801,all_checks_passed=True,
                joint_menus=joint_menu_probe(),retained_state_blocks=blocks.run(),
                combined_test_groups=1,
                finite_diagnostic_is_original_continuum_action=False,
                native_continuum_distinguishing_signal_numerically_certified=False,
                autonomous_terminal_or_preparation_proven=False,
                finite_coupling_probability_proven=False)
if __name__=='__main__':
    result=run()
    (HERE/'original_scattering_records_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
