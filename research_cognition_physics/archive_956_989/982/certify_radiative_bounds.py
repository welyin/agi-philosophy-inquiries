"""Exact rational checks of the inequalities used in research_note_982.
These arithmetic checks supplement the analytic operator proof in the note.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
HERE=Path(__file__).resolve().parent;TARGET=HERE/'radiative_bounds_certificate.json'

def describe(x):return dict(exact=str(x),decimal=float(x))
def run():
    om_lo=F('1.077');om_hi=F('1.078')
    assert om_lo**2<F('1.16')<om_hi**2
    Jlo=(om_lo-1)/2;Jhi=(om_hi-1)/2
    assert Jlo==F('.0385') and Jhi==F('.039')
    # r^2=J/(1+J); the rational inequalities imply r<1/5 and 1/r<6.
    assert Jhi/(1+Jhi)<F(1,25) and Jlo/(1+Jlo)>F(1,36)
    # sqrt(2)<3/2, sqrt(3)<7/4: AP<68, V on n<=2 has norm<32.
    # For the latter use tighter sqrt(2)<1.415 and sqrt(3)<1.733.
    assert F('1.415')**2>2 and F('1.733')**2>3
    MA=1/(2*Jlo)+F(2,5)+F(1,2)+12+F(1,2)
    assert MA<27 and F(5,2)*27<68
    assert (F('1.415')+F('1.733'))*10<32
    assert 50*om_hi<54
    bracket=95+F(11,7)*(2312+54)
    bound=Jhi/F(40000)*bracket
    assert bracket==3813 and bound<F('.00372')
    # At gT=pi/2 the two sines have arguments a,b in (0,1).
    s5lo=F('2.236');s5hi=F('2.237')
    assert s5lo**2<5<s5hi**2
    pi_lo=F('3.14159');pi_hi=F(22,7)  # Known mathematical enclosures, not numerically inferred.
    a_lo=pi_lo*(3-s5hi)/4;a_hi=pi_hi*(3-s5lo)/4
    b_lo=pi_lo*(s5lo-1)/4;b_hi=pi_hi*(s5hi-1)/4
    assert 0<a_lo<a_hi<1 and 0<b_lo<b_hi<1
    lower=lambda x:x-x**3/6
    upper=lambda x:x-x**3/6+x**5/120
    pg_lo=(lower(a_lo)+lower(b_lo))**2/5
    pg_hi=(upper(a_hi)+upper(b_hi))**2/5
    assert pg_lo>F('.38') and pg_hi<F('.4')
    isolated_difference=1-F('.4')-F('.00372')
    assert isolated_difference==F('.59628')
    return dict(round=982,all_rational_checks_passed=True,
        assumptions_not_arithmetic_proofs=[
            'pi in (3.14159,22/7)',
            'operator and finite-support integration-by-parts argument in research_note_982',
            'selected matched three-mode Hamiltonian; no bound here on omitted physical modes'],
        sqrt_omega_interval=[str(om_lo),str(om_hi)],
        A_coefficient_sum_upper=describe(MA),
        full_Fock_isometry_upper_before_rounding=describe(bound),
        strict_isometry_upper=describe(F('.00372')),
        rwa_probability_rational_lower=describe(pg_lo),
        rwa_probability_rational_upper=describe(pg_hi),
        isolated_RWA_difference_lower=describe(isolated_difference))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
