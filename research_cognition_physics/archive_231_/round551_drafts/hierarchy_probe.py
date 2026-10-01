"""Preliminary 551 exact constant-background compatibility after vacuum matching."""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_reference_gravity_constraints as previous


def run():
    p=previous.bare_parameters(); rows=[]
    for beta in (Q(-1,10),Q(0),Q(1,10),Q(9,10),Q(99,100)):
        t=1-beta; c_eff=t*p['C0']; r=6*p['C0']*beta
        h2=c_eff/p['q']; s2=c_eff*p['r']/p['T']*(1-p['S']/p['q'])
        f=p['M02']-(h2+s2)/6
        assert f>0 and h2>0 and s2>0
        v0_new=(p['M02']*r+p['C0']*(h2+s2))/4
        v=v0_new-p['C0']*(h2+s2)/2+(p['lh']*h2*h2+2*p['p']*h2*s2+p['ls']*s2*s2)/4
        assert -p['C0']+p['lh']*h2+p['p']*s2+r/6==0
        assert -p['C0']+p['ls']*s2+p['p']*h2+r/6==0
        assert f*r==4*v
        bound=(3+p['r'])*(2*p['ng']/t-1)/18
        assert f/h2==p['q']*(3+p['r'])*(2*p['ng']/t-1)/(6*p['T'])
        assert f/h2<=bound
        rows.append(dict(R_over_moment_scale_squared=str(r/(p['C0']/4)),
            field_scale_factor=str(t),V0_new=str(v0_new),delta_V0=str(v0_new-p['V0']),
            F_over_h_squared=str(f/h2),upper_bound=str(bound),
            scalar_and_metric_equations_exact=True))
    return dict(status='preliminary_not_completed_round',examples=rows,
                stability_or_quantum_matching_proved=False,completed_unification=False)


if __name__=='__main__':
    result=run()
    with (HERE/'hierarchy_probe_results.json').open('x',encoding='utf8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps(result,ensure_ascii=False))
