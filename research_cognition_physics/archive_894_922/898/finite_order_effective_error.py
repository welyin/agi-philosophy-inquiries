"""898 formal entry: bounded finite candidate versus retained E coefficients."""
from pathlib import Path
import json
import finite_order_effective_error_probe as working
TARGET=Path(__file__).resolve().parent/'finite_order_effective_error_results.json'
def run():
    result=working.run()
    for key in ('working_round','status','formal_round'):result.pop(key)
    result['retained_parameter_degree']=result.pop('retained_time_degree')
    result['complex_parameter_radius']=result.pop('complex_time_radius')
    result.update(round=898,date='2026-10-06',cumulative_numbered_groups=3683,fresh_numbered_groups=1,
        argument_scope='Actual854 finite positive Q_eff versus the common retained-order E task/source polynomial, with explicit Cauchy remainder. No exact E(t), original-Q identification, physical hbar=1 accuracy, locality or full geometry/feedback completion is claimed.')
    return result
if __name__=='__main__':
    result=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},ensure_ascii=False,indent=2))
