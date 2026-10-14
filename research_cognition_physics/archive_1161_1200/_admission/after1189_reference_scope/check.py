"""Unnumbered inner-action diagnostic; not a GR or full-cognitive countermodel."""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent

def run():
    I=np.eye(2);X=np.array([[0.,1.],[1.,0.]]);Z=np.diag([1.,-1.]);eye=np.eye(4)
    Hs=[3*eye+np.kron(Z,I)+np.kron(I,Z)+g*np.kron(X,X) for g in (0.,1.)]
    def expi(H,t):
        w,v=np.linalg.eigh(H)
        return (v*np.exp(1j*t*w))@v.conj().T
    t=np.pi/(2*np.sqrt(5));rows=[];cov=0.;dress=0.;shift=0.
    matrices=[]
    for i in range(4):
        for j in range(4):
            a=np.zeros((4,4),complex);a[i,j]=1;matrices.append(a)
    for g,H in enumerate(Hs):
        U=expi(H,-t)
        unit=float(np.linalg.norm(U.conj().T@U-eye,2))
        energy=float(np.linalg.norm(U.conj().T@H@U-H,2))
        p=float(abs(U[3,0])**2)
        exact=float(g*g/(4+g*g)*np.sin(np.sqrt(4+g*g)*t)**2)
        assert abs(p-exact)<1e-13
        assert unit<1e-13 and energy<1e-12
        rows.append(dict(g=g,eigenvalues=np.linalg.eigvalsh(H).tolist(),probability=p,analytic_probability=exact,unitary_residual=unit,energy_residual=energy))
        for s in (-.6,0.,.9):
            for tau in (-.3,.4):
                v=expi(H,s);minus=expi(H,-s)
                shift=max(shift,float(np.linalg.norm(v@expi(H,-(s-tau))-expi(H,tau),2)))
                for a in matrices:
                    pi=minus@a@v
                    dress=max(dress,float(np.linalg.norm(v@pi@minus-a,2)))
                    left=expi(H,-(s-tau))@a@expi(H,s-tau)
                    right=minus@expi(H,tau)@a@expi(H,-tau)@v
                    cov=max(cov,float(np.linalg.norm(left-right,2)))
    assert min(rows[0]['eigenvalues'])>.99 and min(rows[1]['eigenvalues'])>.76
    assert abs(rows[0]['probability'])<1e-13 and abs(rows[1]['probability']-.2)<1e-13
    assert max(cov,dress,shift)<1e-12
    return dict(status='PASS_REFERENCE_SCOPE_ONLY',kind='inner_action_crossed_product_does_not_select_marked_dynamics',fixed_time=float(t),models=rows,all_matrix_units_checked=16,covariance_residual=cov,dressing_residual=dress,translation_dressing_residual=shift,probability_gap=float(rows[1]['probability']-rows[0]['probability']),numbered_round=False,physical_clock_or_gravity_included=False,exact_scope='All parameters and arbitrary references follow from the analytic unitary identities; floating checks do not certify an infinite algebra type.')

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-11
    else:assert a==b

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();data=run()
    target=HERE/'results.json'
    if args.write:
        with target.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(data,json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps(dict(status=data['status'],mode='first_save' if args.write else 'read_only_compare',probability_gap=data['probability_gap']),ensure_ascii=False))
