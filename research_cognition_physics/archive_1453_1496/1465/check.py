"""Two ordinary vector doublets: construct the common columns and test the physical currents."""
import argparse
import json
from math import hypot, pi, sqrt
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def hermitian_power(h, power):
    values, vectors=np.linalg.eigh(h)
    assert min(values)>0
    return (vectors*values**power)@vectors.conj().T


def complete_basis(last):
    complement=np.linalg.svd(last.conj().T,full_matrices=True)[2][len(last.T):].conj().T
    return np.column_stack((complement,last))


def make_currents(p, data):
    x,y,r,z=(float(data[k]["mean"]) for k in ("x","y","r","z"))
    m_u,m_c,m_t=p["light_masses_up_GeV"]
    down=np.array(p["down_masses_GeV"])
    a=np.sqrt(p["u_right_components_squared"])
    phase=np.exp(1j*p["light_left_phase"])

    def trial(t):
        v_d,v_s=(x+t)/2,(y+r*t)/2
        b_d,b_s=(x-t)/(2*a[0]),(y-r*t)/(2*a[1])
        b=np.array([[b_d,0.],[0.,b_s],[0.,0.],[sqrt(1-b_d**2),0.],[0.,sqrt(1-b_s**2)]])
        heavy=(m_u*a-down[:2]*np.array([v_d*b_d,v_s*b_s]))/(down[3:]*np.diag(b[3:]))
        first=np.array([v_d,v_s,z/phase,*heavy],dtype=complex)
        return first,b

    lo,hi=x,x+.01
    assert np.linalg.norm(trial(lo)[0])<1<np.linalg.norm(trial(hi)[0])
    for _ in range(80):
        mid=(lo+hi)/2
        if np.linalg.norm(trial(mid)[0])<1:
            lo=mid
        else:
            hi=mid
    first,b=trial((lo+hi)/2)
    target=down[:,None]*b
    gram=target.T@target
    complement=np.eye(2)-np.outer(a,a)
    root_inv=hermitian_power(complement,-.5)
    masses2,O=np.linalg.eigh(root_inv@(gram-m_u**2*np.outer(a,a))@root_inv)
    A=np.vstack([a,np.zeros((2,2)),O.T@hermitian_power(complement,.5)])
    up=np.array([m_u,m_c,m_t,*np.sqrt(masses2)])

    # Set the unused two light rows with a physical charged-current orientation.
    s13=abs(first[2])/np.linalg.norm(first[:3])
    c13=sqrt(1-s13**2)
    s12=abs(first[1])/np.linalg.norm(first[:2])
    c12=sqrt(1-s12**2)
    s23=p["light_left_s23"]
    c23=sqrt(1-s23**2)
    trial_rows=np.array([[-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13,0.,0.],
                         [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13,0.,0.]])
    normals=np.column_stack((target,first.conj()))
    q=np.linalg.qr(normals)[0]
    midrows=trial_rows@(np.eye(5)-q@q.conj().T)
    midrows=hermitian_power(midrows@midrows.conj().T,-.5)@midrows
    light=np.vstack([first,midrows])
    heavy_basis=np.linalg.svd(light,full_matrices=True)[2][3:]
    heavy=(up[3:,None]*A[3:])@np.linalg.solve(heavy_basis@target,heavy_basis)
    V=np.vstack([light,heavy])
    R=A@b.T
    Mu=np.diag(up)@complete_basis(A)
    Md=V@np.diag(down)@complete_basis(b)
    np.testing.assert_allclose(V@V.conj().T,np.eye(5),atol=2e-10)
    np.testing.assert_allclose(A.T@A,np.eye(2),atol=2e-12)
    np.testing.assert_allclose(b.T@b,np.eye(2),atol=2e-12)
    mismatch=float(np.max(abs(Mu[:,-2:]-Md[:,-2:])))
    assert mismatch<p["leak_check_GeV"]
    np.testing.assert_allclose(np.linalg.svd(Mu,compute_uv=False),up[::-1],atol=1e-8)
    np.testing.assert_allclose(np.linalg.svd(Md,compute_uv=False),down[::-1],atol=1e-8)
    np.testing.assert_allclose([abs(V[0,0]+R[0,0]),abs(V[0,1]+R[0,1]),
        abs((V[0,1]-R[0,1])/(V[0,0]-R[0,0])),abs(V[0,2])],[x,y,r,z],atol=1e-13)
    return V,R,A,b,up,mismatch


def run():
    p=read("inputs.json")
    V,R,A,b,up,err=make_currents(p,read(p["three_rates"]))
    fu,fd=A@A.T,b@b.T
    rd,rs=R[0,:2]
    used=rd**2/fd[0,0]+rs**2/fd[1,1]
    assert np.isclose(used,fu[0,0]) and fd[0,1]==0
    pv=read(p["weak_charges"])
    proton=pv["proton"]
    cesium=pv["caesium"]
    qp,qpe=proton["pv_electron_scattering_only"]
    cs,cse,csa=cesium["value_and_errors_experiment_atom"]
    dp,dc=2*fu[0,0]-fd[0,0],188*fu[0,0]-211*fd[0,0]
    pulls=[(dp-qp+proton["SM"][0])/hypot(qpe,proton["SM"][1]),
           (dc-cs+cesium["SM"][0])/hypot(cse,csa,cesium["SM"][1])]
    assert max(abs(v) for v in pulls)<1
    single=read(p["single_result"])
    needed=single["charged_windows"][1]["abs_RdRs_lower"]-single["joint_max_abs_RdRs"]
    assert 0<needed<abs(rd*rs)
    k=read(p["K_inputs"])
    lam=k["lambda_reference"]
    lc=V[1,1].conjugate()*V[1,0]
    lt=V[2,1].conjugate()*V[2,0]
    top,et=k["Xtop"]
    charm,ec=k["Pcharm"]
    t=lt*top/lam**5+lc.real*charm/lam
    kp,kpe=k["kplus_e11"]
    kl,kle=k["kL_e11"]
    kp0=(1+k["delta_EM"])*kp
    expK=read(p["kaon_experiment"])["B_rare"]
    b0=1e11*expK["central"]
    lower,upper=b0-2.e11*expK["minus"],b0+2.e11*expK["plus"]
    rmin=sqrt(lower/((1+k["delta_EM"])*(kp-2*kpe)))
    rmax=sqrt(upper/((1+k["delta_EM"])*(kp+2*kpe)))
    coefficient_budget=abs(lt)/lam**5*(2*et)+abs(lc.real)/lam*(2*ec)
    distance=min(abs(t)-rmin,rmax-abs(t))
    neutral_distance=sqrt(k["KOTO_90CL_e11"]/(kl+2*kle))-abs(t.imag)
    extra=min(distance,neutral_distance)-coefficient_budget
    assert extra>0
    return {"shared_two_mass_columns_max_residual_GeV":err,"heavy_up_GeV":up[3:].tolist(),
            "ordinary_left_u_row_moduli":abs(V[0]).tolist(),
            "up_right_light_budget":fu[0,0],"down_right_dd_ss_ds":[fd[0,0],fd[1,1],fd[0,1]],
            "right_ud_us":[rd,rs],"right_charged_product_magnitude":abs(rd*rs),
            "occupied_two_orthogonal_directions":used,"delta_Qp_QCs":[dp,dc],"PV_signed_widths":pulls,
            "minimum_second_direction_invariant_one_width":needed,
            "source_eigenvalues_common_bare_gram_GeV2":np.linalg.eigvalsh(b.T@np.diag(np.square(p["down_masses_GeV"]))@b).tolist(),
            "lambda_t_re_im":[lt.real,lt.imag],"lambda_c_re_im":[lc.real,lc.imag],
            "left_ordinary_top_charm_amplitude_re_im":[t.real,t.imag],
            "nominal_BR_plus_and_L_times_1e11":[kp0*abs(t)**2,kl*t.imag**2],
            "two_width_top_charm_amplitude_error_cap":coefficient_budget,
            "remaining_amplitude_for_every_additional_contribution":extra,
            "sufficient_extra_equivalent_rate_times_1e11":kp0*extra**2}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda x:float(f"{float(x):.11g}"))
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2)
            fp.write("\n")
    else:
        assert result==read("results.json")
        print("1465: common two mass columns, Gram budget and rare matching verified")
