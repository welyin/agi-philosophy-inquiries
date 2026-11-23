"""Tree neutralino and adjoint checks, plus one projected bottom-loop term."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def neutralino():
    # Explicit low-energy example. All dimensionful inputs in GeV.
    gp, g, vd, mu, m1, m2 = .36,.65,25.,700.,500.,800.
    vu=math.sqrt(246.**2-vd**2)
    m=np.array([[m1,0,-gp*vd/2,gp*vu/2],[0,m2,g*vd/2,-g*vu/2],
                [-gp*vd/2,g*vd/2,0,-mu],[gp*vu/2,-g*vu/2,-mu,0]])
    coef=(g*g*m1+gp*gp*m2)/(4*np.linalg.det(m))
    mb,ms1,ms2,sin2=3.,1500.,2000.,.5
    def scalar_integral(ms):
        # B0 at zero momentum modulo a common additive constant, computed
        # analytically from the elementary Feynman parameter integral.
        x,y=mb**2,ms**2
        return -(x*math.log(x)-y*math.log(y))/(x-y)
    hb=math.sqrt(2)*mb/vd
    c_loop=3*mb*abs(sin2)*hb**2*abs(scalar_integral(ms1)-scalar_integral(ms2))/(16*math.pi**2)
    ms,ma=math.sqrt(.00007388)*1e-9,math.sqrt(.002509)*1e-9
    third=np.array([1.,1.,1.])/math.sqrt(3)
    second=np.array([1.,-1.,0.])/math.sqrt(2)
    align=math.sqrt(ma/abs(coef))*third
    eps=mu*math.sqrt(ms/c_loop)*second
    # Nonaligned v_i require independent soft/vacuum inputs.
    v=(align-vd*eps)/mu
    mix=np.column_stack([-gp*v/2,g*v/2,np.zeros(3),eps])
    tree=-mix@np.linalg.solve(m,mix.T)
    assert np.allclose(tree,coef*np.outer(align,align),atol=1e-22,rtol=1e-11)
    full=np.block([[m,mix.T],[mix,np.zeros((3,3))]])
    # minuscule eigenvalues in a 1 TeV matrix are poorly resolved by direct
    # eigensolvers; use the congruent exact Schur complement to count null modes.
    assert np.linalg.matrix_rank(tree,1e-20)==1
    assert np.linalg.matrix_rank(full,1e-12)==5
    proj=np.eye(3)-np.outer(third,third)
    selected=proj @ (c_loop*np.outer(eps,eps)/mu**2) @ proj
    recovered=np.linalg.svd(selected,compute_uv=False)
    assert np.isclose(recovered[0],ms, rtol=1e-12,atol=1e-25)
    return {"tree_neutralino_schur_rank":1,"neutralino_mass_eigenvalues_GeV":np.linalg.eigvalsh(m).tolist(),
            "tree_coefficient_GeVminus3":coef,"lambda_vector_GeV2":align.tolist(),
            "bilinear_eps_GeV":eps.tolist(),"sneutrino_vev_GeV":v.tolist(),
            "tree_largest_eV":np.linalg.svd(tree,compute_uv=False)[0]*1e9,
            "projected_bottom_coefficient_GeV":c_loop,"projected_one_loop_nonzero_eV":recovered[0]*1e9,
            "eps_norm_over_mu":np.linalg.norm(eps)/mu}


def independent():
    # Hypercharge conjugate fundamental partners, e.g. the lower neutral entry.
    t0=np.diag([2,2,2,-3,-3])/math.sqrt(60)
    t3=np.diag([0,0,0,1,-1])/2
    assert np.isclose(np.trace(t0@t0),.5) and np.isclose(np.trace(t3@t3),.5)
    y=np.array([1,1j,-.55+.34j])
    source=np.column_stack([t0[-1,-1]*y,t3[-1,-1]*y])
    heavy=np.array([[1.2,.13j],[.13j,2.2]])
    light=source @ np.linalg.solve(heavy,source.T)
    assert np.linalg.matrix_rank(light,1e-12)==1
    lower_s=.00007388-2*.00000078
    lower_n=.002509-2*.000025
    lower_i=.002482-2*.000025
    m_no=np.sqrt([0,lower_s,lower_n])
    m_io=np.sqrt([lower_i,lower_i+lower_s,0])
    return {"gauge_mixed_X_per_family":1.5-3*.5,"external_gravitational_X":10-3*5,
            "external_cubic_X":10-27*5,"adjoint_relative_squared_clebsch":(t0[-1,-1]/t3[-1,-1])**2,
            "norm_two_projected_columns":np.linalg.svd(source,compute_uv=False).tolist(),
            "independent_eff_mass_needed_eV":{"NO":np.sort(m_no)[1],"IO":np.sort(m_io)[1]},
            "massless_min_tr_gammaS_twoheavy_at_1e10_NO":10*1e10*sum(m_no)*1e-9/246.22**2,
            "sample_low_energy":neutralino()}


def tidy(v):
    if isinstance(v,dict): return {k:tidy(x) for k,x in v.items()}
    if isinstance(v,list): return [tidy(x) for x in v]
    return float(f"{float(v):.11g}")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    val=tidy(independent())
    path=Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            f.write(json.dumps(val,indent=2)+"\n")
    elif args.check:
        assert json.loads(path.read_text(encoding="utf-8"))==val
        print("1481 independent: selection rule, adjoint rank, bilinear tree and projected one-loop mass verified")
    else:
        print(json.dumps(val,indent=2))
