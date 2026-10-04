"""652 entry only: actual curved-target reference velocity vector fields."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_curved_quantum_source as model
import joint_gravity_material_coordinates as original


def vectors(x):
    phi=np.array([0.,x[0],0.,0.,x[1]])
    inverse=model.inverse(phi)
    return inverse[np.ix_([1,4],[1,4])]


def run():
    c=original.point_coefficients(original.fields(16))
    x=np.array([c['H'],c['S']]);phi=np.array([0.,x[0],0.,0.,x[1]])
    F=float(model.F(phi));X=vectors(x)
    exact=F**2/(6*model.M)*np.array([x[1],-x[0]])
    rows=[]
    for step in (.002,.001,.0005):
        deriv=np.stack([(vectors(x+step*d)-vectors(x-step*d))/(2*step) for d in np.eye(2)],axis=-1)
        bracket=deriv[:,1,:]@X[:,0]-deriv[:,0,:]@X[:,1]
        rows.append(dict(step=step,bracket=bracket.tolist(),max_error=float(np.max(abs(bracket-exact)))))
    assert rows[-1]['max_error']<rows[0]['max_error']/12
    assert np.linalg.norm(exact)>.1 and X[0,1]<0
    return dict(candidate_round=652,complete_round=False,original_point_radial_values=x.tolist(),
        F=F,gradient_vector_fields_columns=X.tolist(),gradient_Lie_bracket=exact.tolist(),
        commutator_h_velocity_s_coefficient_in_units_i_hbar_over_w=float(X[0,1]),rows=rows,
        scope='Formal differential expressions on compact smooth radial support h>0, F>0, with the original curved target. No selfadjoint joint reference or instrument is constructed here.')


if __name__=='__main__':
    r=run()
    with (HERE/'curved_reference_velocity_probe_results.json').open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
