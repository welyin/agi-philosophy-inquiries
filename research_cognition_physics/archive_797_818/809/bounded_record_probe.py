"""809 working: binary Weyl-effect identity and Gaussian conversion factor.

No completion of the original perturbative observable domain is claimed here.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'bounded_record_probe_results.json'
def adj(a):return {-k:v.conjugate() for k,v in a.items()}
def mul(a,b):
    out={}
    for k,v in a.items():
        for j,w in b.items():out[k+j]=out.get(k+j,0)+v*w
    return {k:v for k,v in out.items() if v!=0}
def scale(a,c):return {k:v*c for k,v in a.items()}
def run():
    # V is unitary exp(i alpha X/2); exact Gaussian-integer coefficients.
    hp={1:1-1j,-1:1+1j};hm={1:1+1j,-1:1-1j}
    ep=scale(mul(adj(hp),hp),1/8);em=scale(mul(adj(hm),hm),1/8)
    assert ep=={2:-.25j,0:.5,-2:.25j}
    assert em=={2:.25j,0:.5,-2:-.25j}
    assert adj(ep)==ep and adj(em)==em
    assert {k:ep.get(k,0)+em.get(k,0) for k in (0,2,-2)}=={0:1.,2:0.,-2:0.}
    nodes,weights=np.polynomial.hermite.hermgauss(80);rows=[]
    for variance in (.2,1.,3.):
        x=np.sqrt(2*variance)*nodes
        for alpha in (.1,.5,1.,-1.):
            direct=float(np.dot(weights,np.sin(alpha*x)*x)/np.sqrt(np.pi))
            expected=float(variance*alpha*np.exp(-alpha*alpha*variance/2))
            assert abs(direct-expected)<2e-13
            rows.append(dict(variance=variance,alpha=alpha,
                binary_connected_response_factor=float(alpha*np.exp(-alpha*alpha*variance/2)/2),
                gaussian_identity_residual=abs(direct-expected)))
    return dict(working_round=809,all_checks_passed=True,
        exact_Laurent_square_identity=True,
        effects='E+ = H+*H+/8 = (1+sin(alpha X))/2; E- = 1-E+',
        gaussian_checks=rows,original_Wick_Weyl_domain_extension_proven=False,
        original_continuous_response_value_computed=False,
        finite_coupling_measurement_implemented=False,formal_round_completed=False,
        new_numbered_test_groups=0)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
