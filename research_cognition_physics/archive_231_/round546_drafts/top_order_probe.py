"""Uncounted preliminary diagnostic; finite samples do not prove monotonicity."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
import math
import joint_singlet_common_mass_rg as model

def run():
    u=math.log(model.MU/173.34)
    q0=[0.,.2,.3,.4,.5,model.T/3]
    end=model.flow(q0,u,800)
    return dict(preliminary_round=546,completed_round=False,counted_scientific_checks=0,
        q0=q0,terminal_top_squared=end[0].tolist(),
        sampled_order=bool(all(end[0,i+1]>end[0,i] for i in range(len(q0)-1))),
        continuum_order_proved_by_this_probe=False)

if __name__=='__main__':
    result=run();target=Path(__file__).with_suffix('.json')
    if target.exists():assert result==json.loads(target.read_text('utf8'))
    else:
        with target.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))
