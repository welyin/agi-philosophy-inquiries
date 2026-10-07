"""Working 795: cutoff response versus fixed Cauchy data.

Exact finite symplectic response diagnostic. Does not prove original continuum
star-comparison, quantum state positivity or all-order initial-state matching.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
M=((Q(3,5),Q(4,5)),(Q(-4,5),Q(3,5)))


def apply(v):
    return tuple(sum((M[i][j]*v[j] for j in range(2)),Q(0)) for i in range(2))


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def run():
    old_past=[Q(1,3),Q(2,5),Q(1,7)]
    new_past=[Q(3,11),Q(-1,9),Q(2,13)]
    old=new=(Q(0),Q(0))
    for a,b in zip(old_past,new_past):
        old=add(apply(old),(Q(0),a))
        new=add(apply(new),(Q(0),b))
    correction=sub(old,new)
    assert any(correction)
    matched=add(new,correction)
    assert matched==old
    cov=((Q(1,2),Q(0)),(Q(0),Q(1,2)))
    rotated=tuple(tuple(sum((M[i][a]*cov[a][b]*M[j][b] for a in range(2) for b in range(2)),Q(0))
                        for j in range(2)) for i in range(2))
    assert rotated==cov
    assert M[0][0]*M[1][1]-M[0][1]*M[1][0]==1
    rows=[]
    for step in range(7):
        assert matched==old
        defect=sub(old,new)
        assert defect==correction
        # x + beta*x^2 at first fluctuation order: mu_x + beta*C_xx.
        beta=Q(2,3)
        assert old[0]+beta*cov[0][0]==matched[0]+beta*cov[0][0]
        rows.append(dict(step=step, old_mean=[str(x) for x in old],
                         unmatched_difference=[str(x) for x in defect],
                         matched_difference=['0','0']))
        drive=(Q(0),Q(5,7)+Q(step,17))
        old=add(apply(old),drive)
        new=add(apply(new),drive)
        matched=add(apply(matched),drive)
        correction=apply(correction)
    return dict(working_round=795, same_post_initial_surface_source=True,
                initial_mean_difference=[str(x) for x in rows[0]['unmatched_difference']],
                covariance_unchanged_by_mean_displacement=True,
                symplectic_transition_exact=True, samples=rows,
                first_order_composite_mean_preserved_after_initial_data_match=True,
                all_checks_passed=True,
                scope='Finite linear response/centered-covariance diagnostic; continuum constrained displacement and full interacting state remain to be checked.',
                original_interacting_state_constructed=False, formal_round_completed=False)


if __name__=='__main__':
    result=run()
    HERE.joinpath('cauchy_state_probe_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
