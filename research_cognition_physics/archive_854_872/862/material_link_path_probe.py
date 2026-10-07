"""862 working: field-dependent material paths retain curvature response.
Uses the original 753 non-Abelian color connection. The paths are an explicit
local kinematic calibration, not a solved h/Y coordinate grid or quantum graph.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'material_link_path_probe_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime

def exp_anti(A):
    v,U=np.linalg.eigh(1j*A)
    return (U*np.exp(-1j*v))@U.conj().T

def run():
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        color=old.color(1.)
    Ay=-1j*color['A'][1];Az=-1j*color['A'][2]
    curvature=Az@Ay-Ay@Az
    L=.7;shape=.4
    def transport(t,n):
        out=np.eye(3,dtype=complex)
        u=np.arange(n+1)/n;z=t*shape*u*(1-u)
        for dz in np.diff(z):out=exp_anti(-(Ay*L/n+Az*dz))@out
        return out
    points,weights=np.polynomial.legendre.leggauss(64)
    prediction=np.zeros((3,3),complex)
    for u,w in zip((points+1)/2,weights/2):
        prediction-=w*exp_anti(-Ay*L*(1-u))@curvature@exp_anti(-Ay*L*u)*L*shape*u*(1-u)
    norm=float(np.linalg.norm(prediction));assert norm>1e-5
    rows=[]
    for n in (128,256,512):
        step=1e-4
        actual=(transport(step,n)-transport(-step,n))/(2*step)
        err=float(np.linalg.norm(actual-prediction))
        assert err<norm*2/n**2
        rows.append(dict(segments=n,curvature_insertion_error=err,relative_error=err/norm,
                         path_omission_derivative_error=norm))
    assert rows[-1]['curvature_insertion_error']<rows[0]['curvature_insertion_error']/10
    # This is an invariant matrix-norm check; it is not by itself a physical
    # Gauss-state endpoint-bilinear expectation value.
    B=.23*Ay-.31*Az;G=exp_anti(B)
    norm_after=float(np.linalg.norm(G@prediction@G.conj().T))
    assert abs(norm_after-norm)<1e-14
    return dict(kind='round_862_working_material_path_probe',formal_reports=861,
        new_numbered_scientific_groups=0,all_checks_passed=True,
        original_753_color_connection_used=True,connection_field_held_fixed=True,
        endpoints_held_fixed=True,noncommuting_curvature_norm=float(np.linalg.norm(curvature)),
        full_path_derivative_norm=norm,frozen_path_derivative_norm=0.,rows=rows,
        global_color_conjugation_norm_error=abs(norm_after-norm),
        path_shape_is_declared_kinematic_calibration=True,
        full_constrained_material_chart_or_graph_algebra_proved=False,
        quantum_state_or_spectrum_transport_proved=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
