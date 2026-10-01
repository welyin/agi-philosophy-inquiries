"""Unfinished651: causal character of the original material-reference walls."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_gravity_material_coordinates as old


def run():
    rows=[]
    for N in (16,24,32):
        f=old.fields(N);idx=(0,N//4,N//8)
        A=float(f['rh'][idx]);B=float(f['rs'][idx])
        C=float(f['psi'][idx]**-4*(f['dh'][idx]@f['ds'][idx])-f['v'][idx][1]*f['v'][idx][4])
        lam=C/A;norm=B-C*C/A
        roots=sorted([(C-np.sqrt(C*C-A*B))/A,(C+np.sqrt(C*C-A*B))/A])
        assert A<0 and B<0 and norm>0
        v=(f['v'][idx][1],f['v'][idx][4]);dh=f['dh'][idx];ds=f['ds'][idx]
        direct=float(f['psi'][idx]**-4*np.sum((ds-lam*dh)**2)-(v[1]-lam*v[0])**2)
        assert abs(direct-norm)<1e-11
        rows.append(dict(N=N,point=[0.,float(np.pi/2),float(np.pi/4)],
            h_gradient_norm=A,s_gradient_norm=B,mixed_gradient=C,
            constant_s_wall_is_spacelike=True,
            clock_relative_slope=lam,clock_relative_wall_normal_squared=norm,
            timelike_wall_allowed_slope_interval=roots,
            direct_tensor_check=direct))
    return dict(candidate_round=651,complete_round=False,rows=rows,
        scope='Original constrained573 initial source; numerical diagnostics of material wall signature only. Need analytic bounds and the moving-boundary variational map before completion.')


if __name__=='__main__':
    result=run()
    with (HERE/'material_wall_signature_probe_results.json').open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))
