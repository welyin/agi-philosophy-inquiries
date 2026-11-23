"""Holomorphic flux multiplication and physical singular-value/CKM tests."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf8"))


def torus_matrix(degree, shape, cutoff):
    # A Fourier coefficient with total momentum k yields the entire k-th output.
    pairs=list(itertools.combinations_with_replacement(range(degree),2))
    matrix=np.zeros((2*degree,len(pairs)))
    for k in range(2*degree):
        for col,(i,j) in enumerate(pairs):
            if (k-i-j)%degree:
                continue
            n=i+degree*np.arange(-cutoff,cutoff+1)
            coeff=np.exp(-math.pi*shape*(2*n-k)**2/(2*degree)).sum()
            matrix[k,col]=coeff*(1 if i==j else math.sqrt(2))
    return matrix


def ckm(raw):
    a,b,c,phase=raw
    s12,s13,s23=map(math.sin,(a,b,c))
    c12,c13,c23=map(math.cos,(a,b,c))
    w=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/w],
                     [-s12*c23-c12*s23*s13*w,c12*c23-s12*s23*s13*w,s23*c13],
                     [s12*s23-c12*c23*s13*w,-c12*s23-s12*c23*s13*w,c23*c13]])


def vectorize(y):
    return np.array([y[i,j]*(1 if i==j else math.sqrt(2))
                     for i,j in itertools.combinations_with_replacement(range(3),2)])


def compute():
    cfg=read("inputs.json")
    data=cfg["one_TeV_DR"]
    ut=np.array(data["up"])
    dt=np.array(read(data["down_from"])["down_yukawa_dr_1TeV"])
    u,d=ut[:,0],dt[:,0]
    ratio=np.linalg.norm(u)/np.linalg.norm(d)
    V=ckm(data["angles_order_12_13_23_delta"])
    np.testing.assert_allclose(V.conj().T@V,np.eye(3),atol=3e-15)
    # Takagi factorization realizes both singular spectra in one left-handed basis.
    y1=vectorize(np.diag(u))
    y2=vectorize(V @ np.diag(d) @ V.T)
    np.testing.assert_allclose(np.linalg.svd(V@np.diag(d)@V.T,compute_uv=False)[::-1],d)
    cmax=float(u@abs(V)**2@d/(np.linalg.norm(u)*np.linalg.norm(d)))
    angular_floor=math.sqrt(1-cmax**2)
    n=cfg["half_widths"]
    low_ratio=np.linalg.norm(u-n*ut[:,1])/np.linalg.norm(d+n*dt[:,1])
    high_ratio=np.linalg.norm(u+n*ut[:,1])/np.linalg.norm(d-n*dt[:,1])
    lifetime=read(cfg["lifetime_bound"])["x_ceiling_over_two_widths_GeV_inv"]
    cases={}
    for shape in cfg["flat_rectangular_tau_imaginary"]:
        M=torus_matrix(3,shape,9)
        np.testing.assert_allclose(M,torus_matrix(3,shape,15),atol=1e-15)
        sv=np.linalg.svd(M,compute_uv=False)
        cond=sv[0]/sv[-1]
        hu,hd=[np.linalg.solve(M.T,y) for y in (y1,y2)]
        ru,rd=np.linalg.norm(hu),np.linalg.norm(hd)
        tan=ru/rd
        gain=math.hypot(ru,rd)
        beta=math.atan(tan)
        np.testing.assert_allclose(gain*math.sin(beta)*(M.T@(hu/ru)),y1,atol=1e-14)
        np.testing.assert_allclose(gain*math.cos(beta)*(M.T@(hd/rd)),y2,atol=1e-14)
        overlap=abs(np.vdot(hu,hd))/(ru*rd)
        sin_h=math.sqrt(max(0,1-overlap**2))
        assert sin_h+1e-12>=angular_floor/cond
        np.testing.assert_allclose(sv[:3],sv[0])
        np.testing.assert_allclose(sv[3:],sv[-1])
        cases[str(shape)]={"singular_values_common_scale_removed":sv.tolist(),
                           "equal_normalized_vertex_for_fit":gain/(3*shape)**.25,
                           "condition_number":float(cond),
                           "tan_beta_necessary_central":[ratio/cond,ratio*cond],
                           "tan_beta_necessary_two_width_box":[low_ratio/cond,high_ratio*cond],
                           "charged_mass_lower_GeV_if_same_TypeII_bc_bound":low_ratio/(cond*lifetime),
                           "min_sine_higgs_profile_angle_central":angular_floor/cond,
                           "one_fit_tan_beta":tan,"fit_profile_angle_sine":sin_h}
    square=torus_matrix(3,1,9)
    inv_gram=np.linalg.inv(square.T@square)
    blocks=[(0,4),(3,2),(5,1)]
    expected=np.zeros_like(inv_gram)
    base=inv_gram[np.ix_(blocks[0],blocks[0])]
    for i,j in blocks:
        expected[np.ix_([i,j],[i,j])]=base
    np.testing.assert_allclose(inv_gram,expected,atol=1e-13)
    pure_floor=min(base[0,0],base[1,1]-math.sqrt(2)*abs(base[0,1]))
    smallest=np.linalg.svd(square,compute_uv=False)[-1]
    def gain_floor(top,others):
        return (top*math.sqrt(pure_floor)-np.linalg.norm(others)/smallest)/3**.25
    low=ut[:,0]-n*ut[:,1]
    high=ut[:,0]+n*ut[:,1]
    g2,g2err=read("../1523/inputs.json")["sm_dr_at_one_TeV"]["g2"]
    conservative=gain_floor(low[2],high[:2])
    assert conservative>g2+n*g2err
    def solve_common_gain(shape):
        matrix=torus_matrix(3,shape,9)
        ru,rd=(np.linalg.norm(np.linalg.solve(matrix.T,y)) for y in (y1,y2))
        return math.hypot(ru,rd)/(3*shape)**.25,ru/rd
    lo,hi=1.,1.2
    endpoint=[solve_common_gain(t)[0] for t in (lo,hi)]
    assert endpoint[0]>g2>endpoint[1]
    for _ in range(48):
        mid=(lo+hi)/2
        if solve_common_gain(mid)[0]>g2:
            lo=mid
        else:
            hi=mid
    fitted=(lo+hi)/2
    matched,beta_ratio=solve_common_gain(fitted)
    np.testing.assert_allclose(matched,g2,atol=2e-13)
    physical={"inverse_gram_block":base.tolist(),"pure_input_inverse_norm_squared_floor":pure_floor,
              "unit_vertex_lower_central":gain_floor(u[2],u[:2]),
              "unit_vertex_lower_two_width_box":conservative,"weak_vertex_upper_two_width_box":g2+n*g2err,
              "minimum_group_or_matching_multiplier_over_weak_vertex":conservative/(g2+n*g2err),
              "fixed_frame_gains_at_shapes_1_and_1point2":endpoint,
              "one_rectangular_solution_for_weak_vertex":{"imag_tau":fitted,"gain":matched,"tan_beta":beta_ratio}}
    ranks={str(k):int(np.linalg.matrix_rank(torus_matrix(k,1,12))) for k in range(1,7)}
    assert ranks=={"1":1,"2":3,"3":6,"4":8,"5":10,"6":12}
    return {"ranks_flux_1_through_6":ranks,"up_to_down_norm_ratio":ratio,
            "two_width_Frobenius_norm_ratio": [low_ratio,high_ratio],
            "absolute_vertex_test":physical,
            "central_min_sine_between_Yukawa_directions":angular_floor,
            "absolute_J_at_example":abs(np.imag(V[0,0]*V[1,1]*V[0,1].conj()*V[1,0].conj())),
            "shape_cases":cases}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    else:
        np.testing.assert_allclose(a,b,rtol=3e-9,atol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    value=compute()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as out:
            json.dump(value,out,indent=2)
            out.write("\n")
    else:
        compare(value,read("results.json"))
    print("Verified flux multiplication, quark spectra, and normalized Higgs directions.")
