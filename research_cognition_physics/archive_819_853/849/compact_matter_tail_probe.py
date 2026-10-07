"""849: compact range obstruction in a declared oscillator diagnostic.

Zero polynomial moments do not imply zero radiation for a massive equation.
The actual original 53-component range criterion is proved in the note, not
inferred from these finite oscillator channels.
"""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'compact_matter_tail_probe_results.json'

def run():
    def weights(n):
        s,w=np.polynomial.legendre.leggauss(n);u=1-s*s
        base=np.exp(-1/u);norm=w@base
        second=base*(4*s*s/u**4-2/u**2-8*s*s/u**3)
        return s,w*base/norm,w*second/norm
    s,base,d2=weights(256);width=.2
    moments=[float(d2.sum()),float(d2@s)]
    assert max(map(abs,moments))<1e-11
    rows=[];max_formula_error=0.;max_compact_repair_error=0.;max_quad_error=0.
    for omega in (.7,1.3):
        for kappa in (2.,8.,32.,128.):
            frequency=omega*width/kappa
            co=np.cos(frequency*s);si=np.sin(frequency*s)
            cosine=float(d2@co);sine=float(d2@si)
            expected=-frequency**2*float(base@co)
            max_formula_error=max(max_formula_error,abs(cosine-expected))
            tail=float(np.hypot(cosine,sine)/omega)
            bound=frequency**2*np.cos(frequency)/omega
            assert tail>bound-1e-11 and tail>0
            # For f=k/delta*chi(k t/delta), Df has zero on-shell moments.
            repaired=(kappa/width)**2*d2+omega**2*base
            error=float(np.hypot(repaired@co,repaired@si))
            max_compact_repair_error=max(max_compact_repair_error,error)
            if kappa==8.:
                ss,_,dd=weights(192)
                max_quad_error=max(max_quad_error,abs(cosine-float(dd@np.cos(frequency*ss))))
            rows.append(dict(omega=omega,kappa=kappa,source_cosine_moment=cosine,
                retarded_late_tail_amplitude=tail,analytic_positive_lower_bound=bound,
                operator_adapted_compact_source_tail_error=error))
    assert max_formula_error<1e-11 and max_quad_error<1e-11
    assert max_compact_repair_error<3e-8
    return dict(round=849,all_checks_passed=True,fresh_test_groups=1,
        diagnostic_is_original_matter_operator=False,
        ordinary_zero_moments=moments,all_finite_frequency_tails_nonzero=True,
        compact_operator_adapted_control_verified=True,
        maximum_tail_formula_error=max_formula_error,
        maximum_independent_quadrature_error=max_quad_error,
        maximum_compact_control_tail_error=max_compact_repair_error,rows=rows,
        original_matter_range_condition_numerically_verified=False,
        original_pure_metric_existence_disproved=False,
        arbitrary_compact_spacetime_test_obstruction_claimed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
