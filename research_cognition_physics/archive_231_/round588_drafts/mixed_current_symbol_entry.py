"""588 entry: actual curved-edge mixed principal symbol, not a completed round."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_matter_energy_current as prior
TARGET=HERE/'mixed_current_symbol_entry_results.json'

def run():
    phi=np.array([.14,.53,-.17,.1,.42])
    w,k=.8,.9
    def current(x,y,p,q):
        gx,gy=prior.edge_gradients(x,y)
        return float(k/(2*w)*(gy@prior.original.inverse(y)@q-gx@prior.original.inverse(x)@p))
    def kinetic(x,y,p,q,L):
        return float((L[0]*p@prior.original.inverse(x)@p+L[1]*q@prior.original.inverse(y)@q)/(2*w))
    def poisson(p,q,L,step):
        # At the equal-configuration point V_B=V_C=0, so the d_p J d_q H
        # terms vanish exactly. Only derivatives in the singlet p-directions
        # are needed; the other components of K^{-1}p must still be included.
        vx=L[0]/w*(prior.original.inverse(phi)@p)
        vy=L[1]/w*(prior.original.inverse(phi)@q)
        ans=0.
        for a in range(5):
            d=np.eye(5)[a]*step
            ans+=(current(phi+d,phi,p,q)-current(phi-d,phi,p,q))/(2*step)*vx[a]
            ans+=(current(phi,phi+d,p,q)-current(phi,phi-d,p,q))/(2*step)*vy[a]
        return ans
    rows=[]
    for L in ((.7,1.3),(1.,1.)):
        expected=k*(L[1]-L[0])/(2*w*w)*prior.original.inverse(phi)[-1,-1]
        convergence=[]
        for step in (.02,.01,.005,.0025):
            values=[]
            for a,b in ((1,1),(1,-1),(-1,1),(-1,-1)):
                p=np.zeros(5);q=np.zeros(5);p[-1]=a;q[-1]=b
                values.append(poisson(p,q,L,step))
            got=(values[0]-values[1]-values[2]+values[3])/4
            convergence.append(dict(step=step,mixed_coefficient=float(got),error=float(abs(got-expected))))
        assert convergence[-1]['error']<1e-6
        if L[0]!=L[1]:
            assert abs(expected)>.1
            assert convergence[-1]['error']<convergence[0]['error']/40
        else:assert max(x['error'] for x in convergence)<1e-12
        rows.append(dict(lapse=L,analytic_mixed_coefficient=float(expected),rows=convergence))
    # The original lapse Hamiltonian has no cross-node momentum term.
    e=np.eye(5)[-1]
    hs=[kinetic(phi,phi,a*e,b*e,(.7,1.3)) for a,b in ((1,1),(1,-1),(-1,1),(-1,-1))]
    assert abs((hs[0]-hs[1]-hs[2]+hs[3])/4)<1e-14
    return dict(status='588 analytic entry only; not a completed round',checks_passed=True,
                original_target_inverse_ss=float(prior.original.inverse(phi)[-1,-1]),
                cases=rows,original_lapse_H_mixed_coefficient=0.,
                scope='fixed-geometry original finite-generator family; not a general GR or cognition no-go')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    out=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==out
    print(json.dumps(out,ensure_ascii=False,indent=2))
