"""Joint light-kernel, quark-spectrum bounds and a unitary microscopic witness."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def load(relative):
    return json.loads((HERE/relative).read_text(encoding="utf-8-sig"))


def spectrum(A):
    return np.linalg.svd(A,compute_uv=False)


def ckm(s12,s13,s23,phase):
    c12,c13,c23=(math.sqrt(1-x*x) for x in (s12,s13,s23))
    z=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/z],
                     [-s12*c23-c12*s23*s13*z,c12*c23-s12*s23*s13*z,s23*c13],
                     [s12*s23-c12*c23*s13*z,-c12*s23-s12*c23*s13*z,c23*c13]])


def unitary_completion(T):
    left,s,right=np.linalg.svd(T)
    assert max(s)<=1
    d=np.diag(s)
    off=np.diag(np.sqrt(1-s*s))
    middle=np.block([[d,off],[off,-d]])
    unit=np.eye(6,dtype=complex)
    v=unit.copy()
    unit[:3,:3]=left
    v[:3,:3]=right
    result=unit@middle@v
    assert np.allclose(result[:3,:3],T) and np.allclose(result.conj().T@result,np.eye(6))
    return result


def witness(eps,up,down,V):
    q=eps*eps
    D=np.diag([q,q,1.])
    inverse=np.diag(1/np.diag(D))
    yu=np.diag(up)
    yd=V@np.diag(down)
    Cu=inverse@yu@inverse
    Cd=inverse@yd/eps**3
    # Keep the full up-type coefficient symmetric.
    assert np.linalg.norm(Cu,ord=2)<1
    off=np.diag(np.sqrt(1-np.diag(Cu)**2))
    Uu=np.block([[Cu,off],[off,-Cu]])
    Ud=unitary_completion(Cd)
    assert np.allclose(Uu,Uu.T) and np.allclose(spectrum(Uu),1)
    assert np.allclose(D@Uu[:3,:3]@D,yu)
    rebuilt=eps**3*D@Ud[:3,:3]
    assert np.allclose(rebuilt,yd) and np.allclose(spectrum(rebuilt),down[::-1])
    return {"kappa_projected_up":float(np.linalg.cond(Cu)),"kappa_projected_down":float(np.linalg.cond(Cd)),
            "norm_projected_up":float(spectrum(Cu)[0]),"norm_projected_down":float(spectrum(Cd)[0]),
            "both_full_coefficient_condition_numbers":1.0}


def geometry_checks():
    rng=np.random.default_rng(1584)
    counts=0
    for eps in (.03,.13,.45,.8):
        q=eps**2
        D=np.diag([1.]*3+[q]*3)
        for _ in range(24):
            A,B,C,F=[rng.normal(size=(3,k))+1j*rng.normal(size=(3,k)) for k in (2,2,1,1)]
            mass=np.block([[A, q*C],[q*B,F]])  # invertible rescaling of first two columns
            _,_,vh=np.linalg.svd(mass.T,full_matrices=True)
            P=vh.conj().T[:,3:]
            assert np.linalg.norm(mass.T@P)<1e-12
            L=np.linalg.norm(np.linalg.pinv(A.T)@B.T,ord=2)
            d=spectrum(D@P)
            assert d[1]/d[2]<=math.sqrt(1+L*L)+1e-9 and d[2]+1e-10>=q
            counts+=1
    return counts


def run():
    cfg=load("inputs.json")
    low=load(cfg["low"])
    high=load(cfg["high"])["high_scale"]
    sin_low=[a[0] for a in low["quark_sines_and_phase"]]
    angles=[a[0] for a in high["angles_theta12_theta13_theta23_delta"]]
    sin_high=[math.sin(a) for a in angles[:3]]+[angles[3]]
    cases={}
    for tag,raw,ps in [("MZ",low["yukawa"],sin_low),("GUT_SM",high,sin_high)]:
        up,down=raw["up"],raw["down"]
        n=cfg["widths"]
        def edge(pair,i):
            return pair[0]+i*n*pair[1]
        table=[]
        for tol in cfg["matching_relative_mass_tolerances"]:
            ug=(1-tol)**2/(1+tol)*edge(up[2],-1)*edge(up[1],-1)/edge(up[0],1)
            dg=(1-tol)/(1+tol)*edge(down[1],-1)/edge(down[0],1)
            L=cfg["heavy_alignment_cap"]
            norm=cfg["projected_kernel_norm_cap"]
            table.append({"relative_transfer":tol,"yt_yc_over_yu_floor":ug,"ys_over_yd_floor":dg,
                          "kappa_up_floor_L1_norm3":ug/((1+L*L)*norm),
                          "kappa_down_floor_L1":dg/math.sqrt(1+L*L),
                          "L_floor_if_all_kernel_norms_and_condition_numbers_leq3":
                              max(math.sqrt(max(0,ug/norm**2-1)),math.sqrt(max(0,(dg/norm)**2-1)))})
        cases[tag]={"bounds":table,"full_witness":witness(cfg["witness_epsilon"],
                     np.array([v[0] for v in up]),np.array([v[0] for v in down]),ckm(*ps))}
    return {"generic_complex_6x3_checks":geometry_checks(),"scales":cases}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    # Round floating linear algebra to 10 significant digits for cross-build replay.
    result=json.loads(json.dumps(result),parse_float=lambda x:float(format(float(x),".10g")))
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        assert result==json.loads(target.read_text())
    print(json.dumps(result))
