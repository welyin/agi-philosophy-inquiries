"""Working 857 probe: relabeling versus selecting a different material clock.
No numbered scientific group is published by this probe alone.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
HERE=Path(__file__).resolve().parent;TARGET=HERE/'clock_and_cell_dictionary_probe_results.json'

def run():
    y,z,c=F(-2),F(-3),F(3);d=c*c-y*z;eta=F(1)
    def q(la):return -(y+2*c*la+z*la*la)
    def density(la):
        assert q(la)>0
        return eta**2*d/q(la)
    base=density(F(0));rows=[]
    for la in (F(-1),F(-1,10),F(0),F(1,10)):
        active=density(la)
        # Same old clock expressed in the relabeled chart: T=T'-lambda Y1.
        dTold=(F(1),la-la)
        passive_q=-(y*dTold[0]**2+2*c*dTold[0]*dTold[1]+z*dTold[1]**2)
        passive=eta**2*d/passive_q;assert passive==base
        rows.append(dict(lambda_value=str(la),clock_q=str(q(la)),
            active_clock_density_squared=str(active),passive_relabel_density_squared=str(passive),
            active_difference=str(active-base)))
    detM=F(30);transported=(eta/detM)**2*(detM**2)*d/(-y)
    wrong=eta**2*(detM**2)*d/(-y)
    assert transported==base and wrong/base==900
    clock_scale=F(7,3)
    monotone=eta**2*clock_scale**2*d/(clock_scale**2*(-y));assert monotone==base
    # j^mu=(sqrt(D),0,0,0), g_00=-z/D; no floating square root needed.
    current_norm=d*(-z/d);assert current_norm==3
    slope=2*c*d/y**2;assert slope==F(9,2)
    return dict(kind='round_857_working_probe',formal_reports=856,
        newly_published_numbered_groups=0,all_checks_passed=True,
        base_density_squared=str(base),fixed_label_current_norm=str(current_norm),
        current_is_spacelike=True,rows=rows,clock_selection_first_derivative=str(slope),
        transported_label_density_squared=str(transported),untransported_label_error_factor=str(wrong/base),
        monotone_clock_relabel_density_squared=str(monotone),
        physical_counterterm_redundancy_or_scheme_equivalence_proved=False,
        original_on_shell_background_or_complete_continuum_mapping_proved=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
