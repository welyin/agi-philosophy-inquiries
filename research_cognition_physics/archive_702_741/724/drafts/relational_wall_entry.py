"""724 entry: smooth T wall agrees to first order, not as a full boundary."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_material_boundary_transport as old
import joint_squared_reference_readout as quantum
TARGET=HERE/'relational_wall_entry_results.json'


def run():
    c,lam,Y=old.base_data();H=c['H'];kappa=lam/H
    saved=json.loads((ARCHIVE/'joint_material_boundary_transport_results.json').read_text('utf8'))
    p=saved['original_wall']['rows'][-1]
    A=p['h_gradient_norm'];B=p['s_gradient_norm'];C=p['mixed_gradient']
    norm=B-2*kappa*H*C+kappa*kappa*H*H*A
    assert abs(norm-p['clock_relative_wall_normal_squared'])<1e-13 and norm>0
    j=old.chart_jet(0.);gi=j['gi'];g=j['g']
    # Full first jet of the inherited off-shell material-chart metric.
    dgi=np.array([-gi@d@gi for d in j['dg']])
    def normal(displacement,new):
        inverse=gi+np.einsum('a,aij->ij',displacement,dgi)
        df=np.eye(4)[1].copy()
        if new:df[0]-=kappa*displacement[0]
        return df/np.sqrt(df@inverse@df)
    cov=normal(np.zeros(4),False);step=2e-7
    matrices=[]
    for new in (False,True):
        dn=np.column_stack([(normal(step*e,new)-normal(-step*e,new))/(2*step) for e in np.eye(4)]).T
        covariant=dn-np.einsum('r,rab->ab',cov,j['Gamma'])
        matrices.append(covariant[np.ix_(old.TANGENT,old.TANGENT)])
    predicted=np.zeros((3,3));predicted[0,0]=-kappa/np.sqrt(gi[1,1])
    error=float(np.max(abs(matrices[1]-matrices[0]-predicted)))
    assert error<1e-6
    t=quantum.fixture(.017,'T');s=quantum.fixture(.017,'s');I=t['I']
    WT=t['W'];Ws=s['W'];VT=t['V'];Vs=s['V']
    Wst=np.exp(-4*.017)*(s['f']-.2*I)@(t['f']-.2*I)
    mixed=Wst-(Vs@VT+VT@Vs)/2
    V=Vs-kappa*VT;W=Ws-2*kappa*Wst+kappa*kappa*WT
    actual=W-V@V
    reconstructed=(Ws-Vs@Vs)-2*kappa*mixed+kappa*kappa*(WT-VT@VT)
    residual=float(np.linalg.norm(actual-reconstructed,2)/max(1,np.linalg.norm(actual,2)))
    omitted=float(np.linalg.norm(2*kappa*mixed,2))
    assert residual<1e-12 and omitted>1
    names=('research_note_651.md','research_note_652.md','research_note_723.md',
           'joint_material_boundary_transport.py','joint_material_boundary_transport_results.json',
           'joint_squared_reference_readout.py')
    return dict(entry_round=724,formal_round=False,
                classical_original_point=dict(h=H,old_lambda=lam,new_kappa=kappa,
                    same_positive_normal_squared=float(norm)),
                inherited_off_shell_boundary_jet=dict(extrinsic_curvature_shift=predicted.tolist(),
                    independent_normal_derivative_error=error,not_a_new_on_shell_solution=True),
                quantum_conditioned_neighbour_diagnostic=dict(mixed_identity_relative_error=residual,
                    omitted_mixed_operator_norm=omitted),
                dependencies={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in names},
                scope='Direct audit using old651/652. Linear T wall matches original h wall at one point through first derivatives, but changes second fundamental form. Its squared reference contains fixed mixed spatial/velocity terms. No new field equation, quantum wall, boundary source transport or new completed round.')


if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(r,ensure_ascii=False,indent=2))
