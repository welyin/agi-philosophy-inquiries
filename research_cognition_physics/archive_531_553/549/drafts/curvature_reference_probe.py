"""Preliminary exact source mapping; NOT round 549 completion or a GR solution."""
from fractions import Fraction as F
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def run():
    # f0*T=2*pi² and T2=2*T*(h²+s²) imply Delta K=-(h²+s²)/12.
    xi=F(1,6)
    coefficient_K=-F(2*2,48)
    assert coefficient_K==-F(1,12) and 2*coefficient_K==-xi
    # Keep n=Tr Z EXPLICIT. Do not mix one-generation and three-generation traces.
    r,T,C0=F(7,3),F(1),F(1,4)
    q=F(1,4); S=(T-3*q)/r
    h2=C0/q; s2=C0*r/T*(1-S/q)
    assert h2+s2==C0*(3+r)/T
    rows=[]
    for n in (F(8), F(64), F(128)):
        M0sq=n*C0/(24*T)
        actual=M0sq-(h2+s2)/6
        expected=C0*(n-4*(3+r))/(24*T)
        assert actual==expected
        rows.append(dict(identity_trace=str(n),F_at_bare_stationary_point=str(actual)))
    # An arbitrary positive stationary h=s=1 witness for stress algebra only.
    # It is NOT declared to equal any frozen RG or spectral vacuum.
    eps=F(1,4); h=s=F(1)
    eta=[-1,1,1,1]
    p=[eps,0,0,0]; qv=[0,eps,0,0]
    H=[[F(0) for _ in range(4)] for _ in range(4)]
    Sjet=[[F(0) for _ in range(4)] for _ in range(4)]
    H[0][2]=H[2][0]=eps; Sjet[1][3]=Sjet[3][1]=eps
    qq=[[2*(p[i]*p[j]+h*H[i][j]+qv[i]*qv[j]+s*Sjet[i][j])
         for j in range(4)] for i in range(4)]
    box=sum(eta[i]*qq[i][i] for i in range(4))
    assert box==0
    improved=[[p[i]*p[j]+qv[i]*qv[j]+xi*((eta[i]*box if i==j else 0)-qq[i][j])
               for j in range(4)] for i in range(4)]
    assert improved[0][0]==improved[1][1]==(1-2*xi)*eps**2
    assert improved[0][2]==-h*eps/3 and improved[1][3]==-s*eps/3
    assert improved[0][2]!=0  # a constant potential cannot cancel off-diagonal terms
    return dict(status='preliminary_not_a_completed_round',K_quadratic_coefficient=str(coefficient_K),
        same_scale_bare_F_examples=rows,identity_trace_not_selected=True,
        improved_flat_stress=[[str(v) for v in row] for row in improved],
        checks=dict(same_scale_normalization=True,retained_identity_trace=True,nonminimal_flat_residual=True),
        scope=dict(no_low_energy_curvature_matching=True,no_dynamic_constraint_solution=True,
                   no_completed_unification=True))


if __name__=='__main__':
    result=run(); target=HERE/'curvature_reference_probe_results.json'
    if target.exists():
        assert json.loads(target.read_text('utf8'))==result
    else:
        with target.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'checks':result['checks']}))
