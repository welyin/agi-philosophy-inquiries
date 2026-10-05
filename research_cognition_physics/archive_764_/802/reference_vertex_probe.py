"""Working 802: a raw singlet variation is not a physical slice direction.

Exact rational prolongation of the off-shell diagnostic jet from round 773.
This is not an original on-shell background or a continuous response simulation.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
def mat(rows):return np.array([[Q(str(x)) for x in row] for row in rows],dtype=object)
def vec(xs):return np.array([Q(str(x)) for x in xs],dtype=object)
def solve(a,b):
    n=len(b);rows=[[Q(x) for x in a[i]]+[Q(b[i])] for i in range(n)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if rows[i][k])
        rows[k],rows[pivot]=rows[pivot],rows[k]
        d=rows[k][k];rows[k]=[x/d for x in rows[k]]
        for i in range(n):
            if i!=k:
                d=rows[i][k]
                rows[i]=[x-d*y for x,y in zip(rows[i],rows[k])]
    return np.array([r[-1] for r in rows],dtype=object)
def run():
    g=mat([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    gi=g.copy()
    dh=vec([1,.1,0,0]);ds=vec([.2,1,.1,0])
    hh=mat([[.1,.2,.3,.4],[.2,-.2,.1,.2],[.3,.1,.4,-.1],[.4,.2,-.1,.3]])
    hs=mat([[.2,-.1,.1,.2],[-.1,.3,.2,-.3],[.1,.2,-.2,.2],[.2,-.3,.2,.1]])
    jac=np.stack([dh,ds,2*hh@gi@dh,2*hs@gi@ds])
    amplitude=Q(2,7)
    # Raw delta s is constant near the point, with raw delta h=delta g=0.
    dx=vec([0,amplitude,0,0])
    xi=solve(jac,dx)
    # Quadratic scalar background jets and constant metric. d(raw X variation)=0.
    dxi=[]
    for mu in range(4):
        dj=np.stack([hh[:,mu],hs[:,mu],2*hh@gi@hh[:,mu],2*hs@gi@hs[:,mu]])
        dxi.append(solve(jac,-dj@xi))
    dxi=np.stack(dxi)
    projected_g=-(dxi@g+g@dxi.T)
    projected_h=-dh@xi;projected_s=amplitude-ds@xi
    projected_dh=-(dxi@dh+hh@xi)
    projected_ds=-(dxi@ds+hs@xi)
    inverse_metric_variation=-gi@projected_g@gi
    projected_x=np.array([projected_h,projected_s,
        dh@inverse_metric_variation@dh+2*dh@gi@projected_dh,
        ds@inverse_metric_variation@ds+2*ds@gi@projected_ds],dtype=object)
    assert not any(projected_x) and not any(projected_dh) and not any(projected_ds)
    assert any(projected_g.flat)
    # Same original mass-coordinate function x5=s/sqrt(F), tested via its square.
    h=s=Q(1);f=Q(2)-(h*h+s*s)/6
    raw_mass_square_change=(2*s/f+s**3/(3*f*f))*amplitude
    physical_mass_square_change=(2*s/f+s**3/(3*f*f))*projected_s+s*s*h/(3*f*f)*projected_h
    volume_change=sum((gi@projected_g)[i,i] for i in range(4))/2
    assert raw_mass_square_change==Q(66,175)
    assert physical_mass_square_change==0 and volume_change!=0
    return dict(working_round=802,all_checks_passed=True,
        diagnostic_jet_source='773 spacetime_extraction, with exactly rational coefficients',
        raw_singlet_amplitude=str(amplitude),
        compensating_vector=[str(x) for x in xi],
        projected_reference_variation=[str(x) for x in projected_x],
        projected_scalar_first_jets_zero=True,
        nonzero_projected_metric_components=sum(bool(x) for x in projected_g.flat),
        relative_volume_density_variation=str(volume_change),
        raw_original_mass_coordinate_square_variation=str(raw_mass_square_change),
        projected_mass_coordinate_square_variation=str(physical_mass_square_change),
        conclusion='The reference slice fixes h and s, while the compensated metric varies. A bare singlet vertex alone is not the reduced physical response.',
        original_on_shell_background=False,
        total_reduced_vertex_or_response_computed=False,
        physical_interaction_disappears=False,
        formal_round_completed=False,new_numbered_test_groups=0)
if __name__=='__main__':
    result=run()
    (HERE/'reference_vertex_probe_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
