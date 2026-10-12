"""Unnumbered after-1148 phase-screen diagnostic; default readonly."""
from pathlib import Path
import json,math,argparse
import numpy as np
P=Path(__file__).resolve().parent
def compute():
    C=.2; delta=.01
    rows=[]
    for lam in (.05,.01,.001):
        T=math.pi/(C*lam*lam)
        # Leading coefficient mismatch persists; true quartic difference decays.
        theta_leading=-math.pi*delta/C
        theta_quartic=T*.03*lam**4
        theta=theta_leading+theta_quartic
        direct=np.linalg.norm(np.diag([1,1,1,np.exp(1j*theta)])-np.eye(4),2)
        analytic=2*abs(math.sin(theta/2))
        assert abs(direct-analytic)<1e-14
        assert abs(T*delta*lam**2-math.pi*delta/C)<1e-14
        rows.append(dict(lam=lam,T=T,leading_mismatch_phase=theta_leading,
                         quartic_phase=theta_quartic,code_norm=direct))
    assert rows[-1]["quartic_phase"]<rows[0]["quartic_phase"]/1000
    return dict(status="PASS",formal_round=False,formal_check_increment=0,
       new_axioms=0,diagnostic_only=True,actual_parent_field_constructed=False,rows=rows)
a=argparse.ArgumentParser();a.add_argument("--write-results",action="store_true")
args=a.parse_args();d=compute();p=P/'after1148_local_parent_results.json'
if args.write_results:
    assert not p.exists()
    p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8')
else:assert json.loads(p.read_text(encoding='utf-8'))==d
print(json.dumps(dict(status=d["status"],formal_round=False,rows=len(d["rows"]))))
