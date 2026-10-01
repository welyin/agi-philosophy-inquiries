"""Unfinished650: same-vacuum radial waves and the original boundary pairing."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_curved_quantum_source as original
import joint_geometric_boundary_matching as boundary


def run():
    mat,u,_=original.lattice.scalar.parameters()
    phi=np.zeros(5);phi[1]=np.sqrt(u[0]);phi[4]=np.sqrt(u[1])
    F=float(original.F(phi));K=original.metric(phi);idx=[1,4]
    kr=K[np.ix_(idx,idx)]
    Hess=2*np.diag(np.sqrt(u))@mat@np.diag(np.sqrt(u))/F**2
    e,V=np.linalg.eigh(kr);kinvhalf=(V*e**-.5)@V.T
    mass2,v=np.linalg.eigh(kinvhalf@Hess@kinvhalf)
    assert np.min(mass2)>0 and abs(original.node_potential(phi))<1e-25
    rows=[];k=.7;hE=np.diag([-1.,1.,1.]);hJ=hE/F
    for aindex in range(2):
        ar=kinvhalf@v[:,aindex];a=np.zeros(5);a[idx]=ar
        assert abs(a@K@a-1)<1e-14
        omega=np.sqrt(k*k+mass2[aindex])
        for t in (0.,.35,.9):
            amp=np.cos(omega*t)
            # At boundary x=0: delta1 phi=a*cos(omega t), normal derivative0;
            # delta2 phi=0, normal derivative=a*k*cos(omega t).
            dphi1=a*amp;vE2=a*k*amp;vJ2=np.sqrt(F)*vE2
            f=-phi/3;fdotv=float(f@vJ2)
            KJ2=-hJ*fdotv/(2*F)
            PJ,jJ=boundary.jordan_momenta(phi,hJ,KJ2,vJ2)
            rootJ=np.sqrt(abs(np.linalg.det(hJ)))
            PiJ=.5*rootJ*PJ;piJ=rootJ*jJ
            dhJ1=-float(f@dphi1)*hJ/F
            full=-float(np.einsum('ab,ab',PiJ,dhJ1)+piJ@dphi1)
            expected=k*amp**2
            # Deliberately inconsistent comparison freezes Jordan K while the
            # Einstein geometry is fixed: it drops the nonminimal normal term.
            _,jwrong=boundary.jordan_momenta(phi,hJ,np.zeros((3,3)),vJ2)
            wrong=-float(rootJ*jwrong@dphi1)
            assert abs(full-expected)<2e-14 and np.max(abs(PiJ))<2e-14
            rows.append(dict(mode=aindex,mass_squared=float(mass2[aindex]),time=t,
                complete_boundary_symplectic_flux=full,einstein_pairing=expected,
                frozen_Jordan_extrinsic_curvature_flux=wrong,omitted_term=full-wrong,
                first_order_Einstein_area_variation=0.))
    return dict(candidate_round=650,complete_round=False,original_F=F,
        original_radial_mass_squared=mass2.tolist(),rows=rows,
        scope='Local linearized waves at original constant vacuum; nonlinear asymptotically flat completion and boundary-generator integrability still require audit.')


if __name__=='__main__':
    r=run()
    with (HERE/'original_boundary_flux_probe_results.json').open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
