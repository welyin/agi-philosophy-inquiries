"""Independent source-contract review. Default verifies, --write creates one result."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, argparse
import numpy as np
BASE=Path(__file__).resolve().parent

def run():
    t=F(22,7); g=F(1,10**15)
    h1=2*300*t; h2=h1*h1+60000*t
    da=t/F(10**9); db=(t*t/2+t)/F(10**9)
    # Separate the ket and effect errors; do not import author calculations.
    ket_error=2*((2*h1+da)*da+2*db)
    effect_error=2*((h1+da)**2+2*(h2+db))/F(10**12)
    coeff_error=ket_error+effect_error
    alpha=2*g*301*t; beta=g*g*60001*t
    second=alpha*alpha+beta
    ket_tail=alpha**3/(1-alpha)+alpha*beta/((1-alpha)*(1-beta))+beta*beta/(1-beta)
    prob_tail=2*(2*ket_tail*(1+alpha+second)+ket_tail**2+2*alpha*second+second**2)
    c=F(22700475,13958643712)
    radius=(coeff_error+prob_tail/(g*g))/c
    assert F(6325181164,10**11)<radius<F(6325181166,10**11)
    assert radius<F(8,125)
    saved=json.loads((BASE/'results.json').read_text(encoding='utf-8'))
    for key, value in [('coefficient',c),('A',h1),('B',h2),('delta1',da),('delta2',db),('coefficient_error',coeff_error),('actual_pair_tail_over_g2',prob_tail/(g*g)),('xi_radius',radius)]:
        assert F(saved[key]['exact'])==value, key
    # Gaussian moments: radial Gaussian in d=3 has E[(P^2)^k]=prod(3+2j)*v^k.
    variance=F(7500)
    assert 3*variance==150**2
    assert 15*variance**2<30000**2
    assert 945*variance**4<1800000000**2
    mass=F(10**9)
    free_error=t*(F(3*30000)/mass**2+F(1800000000)/(8*mass**3))
    rms=F(1,100)+t*150/mass
    assert F(saved['com_free_only']['evolution_upper']['exact'])==free_error
    assert F(saved['com_free_only']['rms_upper']['exact'])==rms
    # Independent matrix checks of the algebraic probability-coefficient inequality.
    rng=np.random.default_rng(1185)
    ratios=[]
    def random_vector(): return rng.normal(size=4)+1j*rng.normal(size=4)
    def effect():
        a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
        a=a.conj().T@a
        return a/(1+np.linalg.norm(a,2))
    for unused in range(64):
        e=random_vector(); e/=np.linalg.norm(e)
        a=random_vector(); b=random_vector()
        ap=a+0.03*random_vector(); bp=b+0.02*random_vector()
        f=effect(); fp=0.98*f+0.02*effect()
        q=lambda aa,bb,ff: float(np.real(np.vdot(aa,ff@aa)+2*np.vdot(e,ff@bb)))
        actual=abs(q(ap,bp,fp)-q(a,b,f))
        aa=np.linalg.norm(a); bb=np.linalg.norm(b)
        deltaa=np.linalg.norm(ap-a); deltab=np.linalg.norm(bp-b)
        epsilon=np.linalg.norm(fp-f,2)
        upper=(2*aa+deltaa)*deltaa+2*deltab+epsilon*((aa+deltaa)**2+2*(bb+deltab))
        assert actual<=upper+1e-12
        ratios.append(actual/upper)
    # Exact integer coefficient of omitted preparation derivative: E2=|1>, F=|+><+|.
    e0=np.array([1,0]); e2=np.array([0,1]); twice_f=np.array([[1,1],[1,1]])
    assert int(e0@twice_f@e2)==1
    # Effect-only variant is sin(g)^2; its coefficient is one, despite F_g=F at g=0.
    x=np.array([[0,1],[1,0]]); f=np.diag([0,1])
    d1=-1j*x@e0; d2=-0.5*e0
    readout_coeff=float(np.real(np.vdot(d1,f@d1)+2*np.vdot(e0,f@d2)))
    assert readout_coeff==1
    boundfiles=['scope_decision.md','coefficient_bound.md','com_scope.md','check.py','results.json','reviewer_check.py']
    hashes={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in boundfiles}
    parent=BASE.parent/'after1185_dynamic_quadrupole'/'response_proof.md'
    parenthash=hashlib.sha256(parent.read_bytes()).hexdigest()
    assert parenthash=='50b7c1c527217a8714b241255b36caac5fcfd82f2e3762acef27f7c7a40d7686'
    return {'schema':'independent_conditional_source_review_v1','status':'PASS_CONDITIONAL_ONLY',
        'bindings':hashes,'frozen_response_proof_sha256':parenthash,
        'coefficient_error_exact':str(coeff_error),'actual_pair_tail_over_g2_exact':str(prob_tail/(g*g)),
        'xi_radius_exact':str(radius),'certified_interval':['0.06325181164','0.06325181166'],
        'matrix_inequality_samples':64,'matrix_max_ratio':max(ratios),
        'omitted_preparation_p2':1,'omitted_effect_p2':1,
        'com_free_error_exact':str(free_error),'com_rms_exact':str(rms),
        'scope':{'full_insert_contract_added':True,'actual_word_bounds_independent':True,
                 'local_source_qualification_proved':False,'actual_source_exists_proved':False,
                 'strict_COM_support_proved':False,'new_round_approved':False,
                 'new_axiom_approved':False,'GR_completion_approved':False}}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run();dest=BASE/'review_checks.json'
    if args.write:
        with dest.open('x',encoding='utf-8') as out: json.dump(result,out,ensure_ascii=False,indent=2);out.write('\n')
        print('SAVED_INDEPENDENT_CONDITIONAL_REVIEW')
    else:
        old=json.loads(dest.read_text(encoding='utf-8'));assert old==result
        print('READ_ONLY_MATCH: conditional source coefficient and COM scope reviewed')
