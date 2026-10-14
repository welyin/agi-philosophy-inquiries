"""Unnumbered scope diagnostic. GR and classical scalar matter are inputs."""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent

def run():
    m=mp=a0=A=B=1.0
    H0=np.sqrt((B*B+m*m*A*A)/(6*mp*mp))
    def flow(x):
        a,h,p,v=x
        return np.array([a*h,-v*v/(2*mp*mp),v,-3*h*v-m*m*p])
    def pressure(x): return .5*(x[3]**2-m*m*x[2]**2)
    def density(x): return .5*(x[3]**2+m*m*x[2]**2)
    def hdot(x): return flow(x)[1]
    def addot(x): return x[0]*(x[1]**2+hdot(x))
    def directional(fn,x,step):
        v=flow(x)
        return (fn(x+step*v)-fn(x-step*v))/(2*step)
    pair=[np.array([a0,H0,s*A,B]) for s in (1.,-1.)]
    constraints=[float(3*mp*mp*x[1]**2-density(x)) for x in pair]
    stresses=[np.diag([density(x)]+[x[0]**2*pressure(x)]*3) for x in pair]
    predicted=np.array([-4*m*m*A*B,2*m*m*A*B/(mp*mp),2*a0*m*m*A*B/(mp*mp)])
    differences=[];conservation=[]
    for step in [1e-5,1e-6]:
        vals=np.array([[directional(fn,x,step) for fn in (pressure,hdot,addot)] for x in pair])
        diff=vals[0]-vals[1]
        differences.append(dict(step=step,delta_pdot_Hddot_a3=diff.tolist(),formula_residual=float(np.max(abs(diff-predicted)))))
        conservation.append(float(max(abs(directional(density,x,step)+3*x[1]*(density(x)+pressure(x))) for x in pair)))
    assert max(abs(np.array(constraints)))<1e-14
    assert np.max(abs(stresses[0]-stresses[1]))==0
    assert max(row['formula_residual'] for row in differences)<1e-8
    assert max(conservation)<1e-8
    return dict(status='PASS_SCOPE_ONLY',kind='classical_given_GR_initial_source_not_full_material_state',parameters=dict(a0=a0,A=A,B=B,m=m,Mp=mp),initial_H=float(H0),initial_stress_common=stresses[0].tolist(),constraint_residual=constraints,analytic_differences=predicted.tolist(),directional_checks=differences,source_conservation_residual=conservation,counts_as_numbered_round=False,quantum_parent_or_finite_instrument_certificate=False)

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float): assert abs(a-b)<1e-8
    else: assert a==b

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    data=run();target=HERE/'results.json'
    if args.write:
        with target.open('x',encoding='utf-8') as f: json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    else: compare(data,json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps({'status':data['status'],'mode':'first_save' if args.write else 'read_only_compare','analytic_differences':data['analytic_differences']},ensure_ascii=False))
