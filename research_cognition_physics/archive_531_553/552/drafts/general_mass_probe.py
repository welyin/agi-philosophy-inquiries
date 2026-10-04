"""Preliminary exact unequal-mass curvature/Hessian compatibility, not a round."""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from joint_vacuum_hierarchy_matching import parameters


def run():
    p=parameters()
    a,b,c=p['lh'],p['p'],p['ls']
    det=a*c-b*b
    v=[(c-b)/det,(a-b)/det]
    A=sum(v)
    rows=[]
    examples=[('healthy_unequal_mass',(Q(3,10),Q(1,5)),Q(1,10)),
              ('positive_F_but_radial_saddle',(Q(9,5),Q(17,10)),Q(3)),
              ('degenerate_curvature',(Q(3,2),Q(3,2)),Q(3)),
              ('degenerate_curvature',(Q(3,2),Q(3,2)),Q(6))]
    for label,m,R in examples:
        u=[(c*m[0]-b*m[1])/det,(a*m[1]-b*m[0])/det]
        x=[u[i]-R*v[i]/6 for i in range(2)]
        F0=p['M02']-sum(u)/6
        F=p['M02']-sum(x)/6
        assert min(x)>0 and F>0
        assert F==F0+R*A/36
        V0=(F0*R+sum(m[i]*u[i] for i in range(2)))/4
        V=V0-(m[0]*x[0]+m[1]*x[1])/2+(a*x[0]**2+2*b*x[0]*x[1]+c*x[1]**2)/4
        assert a*x[0]+b*x[1]==m[0]-R/6
        assert b*x[0]+c*x[1]==m[1]-R/6
        assert F*R==4*V
        z=R/(36*F)
        ratio=1-z*A
        assert ratio==F0/F
        assert (a-z)*(c-z)-(b-z)**2==det*ratio
        rows.append(dict(label=label,mass_coefficients=list(map(str,m)),R=str(R),
                         h2=str(x[0]),s2=str(x[1]),F=str(F),flat_F0=str(F0),V0=str(V0),
                         rank_one_positivity_factor=str(ratio)))
    assert Q(rows[0]['rank_one_positivity_factor'])>0
    assert Q(rows[1]['rank_one_positivity_factor'])<0
    assert rows[2]['rank_one_positivity_factor']==rows[3]['rank_one_positivity_factor']=='0'
    assert rows[2]['V0']==rows[3]['V0']=='3'
    return dict(status='preliminary_552_not_independently_reviewed_or_completed',examples=rows,
                running_curvature_couplings_or_vacuum_not_included=True)


if __name__=='__main__':
    result=run()
    with (HERE/'general_mass_probe_results.json').open('x',encoding='utf8',newline='\n') as stream:
        stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
