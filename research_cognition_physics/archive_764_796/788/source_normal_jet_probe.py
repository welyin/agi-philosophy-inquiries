"""Unpublished 788: on-shell determinant equality does not fix a tadpole.

Finite Euclidean gauge model only, with an explicit relational-field correction.
No original Lorentzian in-in source or endpoint theorem is asserted.
"""
from fractions import Fraction as Q
from pathlib import Path
import json


def run():
    mass=Q(3,2)
    rows=[]
    for a,b in ((Q(1),Q(0)),(Q(2),Q(1,3)),(Q(3,2),Q(-2,5))):
        # At q=0,r=x: Hessian of m(r-q^2/2)^2/2 plus the linear gauge square.
        # H(x) = [[a^2-m*x,a*b],[a*b,m+b^2]], ghost A=a at q=0.
        det0=a*a*mass
        det_derivative=-mass*(mass+b*b)
        source_from_logdet=det_derivative/(2*det0)
        covariance_qq=(mass+b*b)/det0
        covariance_rr=a*a/det0
        source_from_cubic=-mass*covariance_qq/2
        assert source_from_cubic==source_from_logdet
        mean_r=-source_from_cubic/mass
        invariant_mean=mean_r-covariance_qq/2
        assert invariant_mean==0 and covariance_rr==1/mass
        rows.append(dict(a=str(a),b=str(b),
                         gauge_cancelled_on_shell_determinant=str(det0/(a*a)),
                         one_loop_source_r=str(source_from_logdet),
                         mean_r_coefficient=str(mean_r),
                         relational_quadratic_correction=str(-covariance_qq/2),
                         relational_mean_coefficient=str(invariant_mean),
                         physical_free_covariance=str(covariance_rr)))
    assert len({row['one_loop_source_r'] for row in rows})==len(rows)
    return dict(round=788,status='working_probe_not_signed_off',
                mass=str(mass),same_on_shell_effective_action=True,
                same_untransported_tadpole=False,
                invariant_relational_mean_after_joint_transport=True,cases=rows,
                original_continuum_source_matching_proven=False)


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('source_normal_jet_probe_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
