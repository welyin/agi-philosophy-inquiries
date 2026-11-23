"""Check symmetric Higgs-frame flavour transport and CPT kaon bounds."""
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def response(weights, mass, alpha, mu):
    loop = alpha/(4*math.pi)
    log = 2*math.log(mass/mu)
    lr1, lr2, s1, s2 = weights
    lr = -1.5*loop*lr1+(1-loop)*lr2
    ss = (1+loop*(4.5-3*log))*s1+loop*(.125-log/12)*s2
    return lr,ss


def action_frame(base, shape, gain):
    mat = gain*(3*shape)**.25*base.torus_matrix(3,shape,13)
    units = []
    for row in mat:
        y = np.zeros((3,3),complex)
        for val,(i,j) in zip(row,itertools.combinations_with_replacement(range(3),2)):
            y[i,j] = y[j,i] = val/(math.sqrt(2) if i != j else 1)
        units.append(y)
    return mat,np.asarray(units)


def verify_modes(base, shape, gain, data, responses):
    matrix, y = action_frame(base,shape,gain)
    v = base.ckm(data["angles_order_12_13_23_delta"])
    down = np.array(read("../1527/"+data["down_from"])["down_yukawa_dr_1TeV"])[:,0]
    light = np.linalg.solve(matrix.T,base.vectorize(v@np.diag(down)@v.T))
    light /= np.linalg.norm(light)
    massframe = np.array([v.conj().T@r@v.conj() for r in y])
    np.testing.assert_allclose(massframe,massframe.transpose(0,2,1),atol=2e-15)
    l = np.einsum("k,kij->ij",light,massframe)
    np.testing.assert_allclose(l,np.diag(np.diag(l)),atol=2e-15)
    pairs = [(0,1),(0,2),(1,2)]
    rows = np.array([massframe[:,i,j] for i,j in pairs])
    gram = rows@rows.conj().T
    minimum = np.linalg.svd(matrix,compute_uv=False)[-1]**2/2
    assert np.linalg.eigvalsh(gram)[0] >= minimum-2e-12
    # The remaining five complex fields become ten real neutral modes.
    _,_,vh = np.linalg.svd(light.conj()[None,:],full_matrices=True)
    heavy = vh.conj().T[:,1:]
    z = rows@heavy
    real_ports = np.concatenate((z,1j*z),axis=1)/math.sqrt(2)
    rng = np.random.default_rng(1529)
    rotation,_ = np.linalg.qr(rng.normal(size=(10,10)))
    kaon = (real_ports@rotation)[0]
    masses = np.linspace(530,1960,10)
    alpha,mu,gap,weights = responses
    terms = []
    for k,m in zip(kaon,masses):
        lr,ss = response(weights,m,alpha,mu)
        terms.append((lr*abs(k)**2+ss*(k*k).real)/m**2)
    assert min(terms)>0
    np.testing.assert_allclose(sum(abs(kaon)**2),gram[0,0].real,atol=1e-13)
    np.testing.assert_allclose(sum(kaon**2),0,atol=2e-15)
    assert sum(terms) >= gap*minimum/max(masses)**2
    # Dimensional count: the inverse image of all diagonal matrices has dimension three.
    _,_,vh_d = np.linalg.svd(rows,full_matrices=True)
    silent = vh_d.conj().T[:,3:]
    assert np.linalg.matrix_rank(rows,tol=1e-12)==3
    np.testing.assert_allclose(rows@silent,0,atol=2e-15)
    np.testing.assert_allclose(silent@silent.conj().T@light,light,atol=2e-15)
    return dict(K_Bd_Bs_source_norms_squared=gram.diagonal().real.tolist(),
                offdiagonal_gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
                coupled_real_modes=10, example_scalar_M12_abs_GeV=sum(terms),
                diagonal_complex_subspace=silent.shape[1])


def compute():
    cfg=read("inputs.json")
    previous=read(cfg["prior"])["absolute_vertex_test"]["one_rectangular_solution_for_weak_vertex"]
    spec=importlib.util.spec_from_file_location("flux_round_1527",HERE/"../1527/check.py")
    base=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    shape,gain=previous["imag_tau"],previous["gain"]
    g=read(cfg["gauge"])["sm_dr_at_one_TeV"]
    n=cfg["marginal_widths"]
    g_strong=g["g3"]
    alphas=[(g_strong[0]+sign*n*g_strong[1])**2/(4*math.pi) for sign in (-1,1)]
    mu=cfg["mu_GeV"]
    mass=cfg["mass_window_GeV"]
    weights=cfg["table1_K_LR1_LR2_SLL1_SLL2_GeV3"]
    frac=cfg["hadronic_relative_box"]
    gap=min(lr-abs(ss) for factors in itertools.product([1-frac,1+frac],repeat=4)
            for m,alpha in itertools.product(mass,alphas)
            for lr,ss in [response(np.array(factors)*weights,m,alpha,mu)])
    assert gap>0
    mat,_=action_frame(base,shape,gain)
    smin=np.linalg.svd(mat,compute_uv=False)[-1]
    # Positive theta coefficients decrease individually. Weyl with a Frobenius envelope.
    left=base.torus_matrix(3,1.,13)
    right=base.torus_matrix(3,1.2,13)
    assert np.all(left>=right-1e-15)
    uniform_promote=(g["g2"][0]-n*g["g2"][1])*3**.25
    uniform_s=uniform_promote*(np.linalg.svd(left,compute_uv=False)[-1]
                               -np.linalg.norm(left-right)-1e-12)
    assert uniform_s>0
    ktev=cfg["ktev"]
    freq,dfreq=ktev["delta_m_1e6_per_second"]
    tau,dt=ktev["tau_s_ps"]
    tl,*dlong=cfg["kloe_tau_L_ns"]
    delta_m=(freq+n*dfreq)*1e6
    ts=(tau-n*dt)*1e-12
    tl=(tl-n*sum(dlong))*1e-9
    measured=.5*cfg["hbar_GeV_s"]*math.sqrt(delta_m**2+1/(ts*tl))
    result={}
    for label,s in [("matched_rectangle",smin),("all_shapes_1_to_1p2",uniform_s)]:
        row_floor=.5*(cfg["frame_transfer_floor"]*s)**2
        lower=gap*row_floor/max(mass)**2
        result[label]=dict(geometric_smallest_singular_value=s,sd_total_weight_floor=row_floor,
            scalar_M12_abs_lower_GeV=lower, relative_to_data_cap=lower/measured,
            other_M12_min_GeV=lower-measured, relative_remainder_max=measured/lower,
            formal_response_scale_GeV=math.sqrt(gap*row_floor/measured))
    # If right and left source roles are distinct, an upper-triangular system hides all mixing.
    one_sided=np.zeros((6,3,3),complex)
    for k,(i,j) in enumerate(itertools.combinations_with_replacement(range(3),2)):
        one_sided[k,i,j]=1
    for i,j in [(0,1),(0,2),(1,2)]:
        z=one_sided[:,i,j]
        zz=np.r_[z,1j*z]/math.sqrt(2)
        assert np.vdot(z,z).real==1 and abs(sum(zz**2))<1e-14
        assert abs(np.vdot(one_sided[:,j,i],z))<1e-14
    [[s2,*ds2],[lr4,*dlr4]]=cfg["RBC_SUSY_R2_R4_3GeV"]
    ratio=(lr4-n*sum(dlr4))/(abs(s2)+n*sum(ds2))
    direct=verify_modes(base,shape,gain,read("../1527/inputs.json")["one_TeV_DR"],
        (g_strong[0]**2/(4*math.pi),mu,gap,weights))
    return dict(NLO_response_gap_GeV3=gap, two_width_total_M12_upper_GeV=measured,
                low_scale_SUSY_Q4_to_Q2_abs_ratio_floor=ratio, regions=result,
                mass_basis_and_nondegenerate_example=direct)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-9,atol=1e-20)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    options=parser.add_mutually_exclusive_group(required=True)
    options.add_argument("--check",action="store_true")
    options.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    out=compute()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        compare(out,read("results.json"))
        print("1529 verified: frame, chirality and NLO gap, finite kaon mass/width bound.")
