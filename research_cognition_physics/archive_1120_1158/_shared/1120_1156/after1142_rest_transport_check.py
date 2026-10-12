"""Exact constants for original w; readonly by default. No continuum simulation."""
from pathlib import Path
from fractions import Fraction as F
import json, hashlib, argparse
BASE=Path(__file__).resolve().parent

def run():
    sources={n:hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in [
        'after1142_common_account_screen.md','after1142_rest_transport_qualification.md',
        '../1133/proof.md','../1142/proof.md']}
    r3=2048*(4**4+2*3**3+2*2**2)
    r4=2048*(6**6+6*5**5+12*4**4+12*3**3)
    assert r3<10**6 and r4<2*10**8
    eta2=128+F(8,4)
    eta3=10**6+F(3*128,4)+F(8*6,16)
    eta4=2*10**8+F(6*10**6,4)+F(3*128,16)+F(4*128*6,16)+F(8*36,64)
    assert eta2==130 and eta3<1000100 and eta4<202000000
    w2=130+2*8*8+128
    w4=202000000+4*1000100*8+6*130*128+4*8*10**6+2*10**8
    assert w2<400 and w4<5*10**8
    I2=23*(3*400+9*5*10**8)
    I2decl=104000000000
    rest=F(3,4)*I2decl
    integral=F(5,4)*rest
    assert I2<I2decl and rest==78000000000 and integral<10**11
    return dict(status='PASS',sources=sources,groups={
        'original_smooth_weight':dict(r3=r3,r4=r4,eta2=str(eta2),eta3=str(eta3),eta4=str(eta4),w2=w2,w4=w4),
        'full_fourier_and_transport':dict(I2_upper=I2,I2_declared=I2decl,rest_defect_upper=str(rest),integrated_upper=str(integral))},
        continuum_simulated=False,analytic_derivative_proof_checked_by_code=False,
        microcausality_proved=False,complete_common_balance_proved=False,new_axioms=0,science_increment=0)

def main():
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    got=run();out=BASE/'after1142_rest_transport_results.json'
    if a.write_results:
        assert not out.exists(),'preserve existing result'
        out.write_text(json.dumps(got,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert json.loads(out.read_text(encoding='utf-8'))==got
    print(json.dumps(dict(status=got['status'],groups=len(got['groups']),saved_results_match=not a.write_results),ensure_ascii=False))
if __name__=='__main__': main()
