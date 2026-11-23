"""Finite upper-block CP bound; exact asymptotic lower-block completion."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def low_band(values,n):
    return values[0]-n*math.hypot(*values[1:])


def upper_angles(eta,rd,rs,moduli,channel):
    leaks={"d":eta*rd**2/(1-rd**2),
           "s":eta/min(1-(1+eta)*rd**2,1/rs**2-1),
           "b":eta/(1-(1+eta)*rs**2)}
    products={"KK":("us","cs","ub","cb"),"db":("ud","cd","ub","cb")}
    values=[min(1.,leaks[k[1]]/moduli[k]) for k in products[channel]]
    return sum(math.asin(x) for x in values)


def invert(phi,channel,moduli,cfg):
    rd,rs=(cfg["template_broad_mass_ratios"][k] for k in ("d_over_s_max","s_over_b_max"))
    a,b=0.,.02
    assert upper_angles(b,rd,rs,moduli,channel)>phi
    for _ in range(65):
        mid=(a+b)/2
        if upper_angles(mid,rd,rs,moduli,channel)<phi:
            a=mid
        else:
            b=mid
    v=cfg["vev_GeV"]
    mb=cfg["vev_mass_eigenvalue_bounds_GeV"][-1]
    cap=cfg["p_norm_cap"]
    # N is the unbroken pole, and the extra phase shift is allowed at the weak port.
    mass_ceiling=math.hypot(mb,cap*v/math.sqrt(2*a))
    frac=a*(1-mb**2/100**2)
    return dict(minimum_phase_rad=phi,lower_eta=a,
                M_over_p_in_heavy_limit_GeV=v/math.sqrt(2*a),
                max_bare_mass_at_perturbative_p_GeV=mass_ceiling,
                lower_summed_Z_doublet_deficit_if_M_ge_100_GeV=frac/(1+frac))


def tree_data(cfg,n):
    raw=read(cfg["tree_moduli"])["direct_semileptonic_moduli"]
    charm=read(cfg["charm"])["charm_readings"][str(n)]
    mods={k[1:].lower():.99*low_band(values,n) for k,values in raw.items()}
    mods["cd"]=.99*charm["Vcd"][0]
    mods["cs"]=.99*charm["Vcd"][0]/charm["cd_over_cs"][1]
    mods["ud"]=.90
    bel=cfg["Belle_CP_even_KK"]
    amin=low_band(bel["A"],n)-bel["extra_absolute_A"]
    low=(low_band(bel["R"],n))*(1-bel["relative_R_transfer"])
    hi=(2*bel["R"][0]-low_band(bel["R"],n))*(1+bel["relative_R_transfer"])
    s=amin*low/(2*(1+math.sqrt(hi)))
    phase=math.asin(s)-bel["extra_weak_phase_radians"]
    assert phase>0
    out=dict(moduli_lower=mods,GLW_inputs={"A_min":amin,"R_min":low,"R_max":hi},
             raw_KK=invert(phase,"KK",mods,cfg),
             KK_quartet_imaginary_abs_lower=math.sin(phase)*math.prod(mods[i] for i in ("ub","us","cs","cb")))
    if n==2:
        gamma=cfg["gamma_charged_B"]
        bound=math.radians(gamma["range95_deg"][0]-gamma["extra_phase_translation_deg"])
        out["registered_gamma_db"]=invert(bound,"db",mods,cfg)
    return out


def empirical_matrix(inp):
    s12,s13,s23,phase=np.array(inp["quark_sines_and_phase"])[:,0]
    c12,c13,c23=np.sqrt(1-np.array([s12,s13,s23])**2)
    p=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/p],
                     [-s12*c23-c12*s23*s13*p,c12*c23-s12*s23*s13*p,s23*c13],
                     [s12*s23-c12*c23*s13*p,-c12*s23-s12*c23*s13*p,c23*c13]])


def lower_completion(inp,v,heavy):
    V=empirical_matrix(inp)
    eig=(np.array(inp["yukawa"]["down"])[:,0]*v/math.sqrt(2))**2
    H=(V*eig)@V.conj().T
    vals,vec=np.linalg.eigh(1j*H.imag)
    z=vec[:,0]*np.sqrt(-2*vals[0])
    K=H+np.outer(z,z.conj())
    assert np.max(abs(K.imag))<1e-14
    m=np.linalg.cholesky(K.real)
    s=np.linalg.solve(m,z)
    sn=float(np.linalg.norm(s))
    assert sn<1
    np.testing.assert_allclose(m@(np.eye(3)-np.outer(s,s.conj()))@m.T,H,atol=2e-14)
    J=float((V[0,1]*V[1,2]*V[0,2].conj()*V[1,1].conj()).imag)
    # Genuine finite check at a well-conditioned ratio; analytical bounds at the physical pole.
    def finite(L):
        M=np.block([[m,np.zeros((3,1))],
                    [L*s.conj()[None,:],np.array([[L*math.sqrt(1-sn**2)]])]])
        U,ss,_=np.linalg.svd(M)
        order=np.argsort(ss)
        r=U[:3,order[:3]]
        return float((r[0,1]*r[1,2]*r[0,2].conj()*r[1,1].conj()).imag)
    nz=float(np.vdot(z,z).real)
    bound=nz*eig[-1]/(heavy**2-eig[-1])
    return dict(target_light_mass_squared_GeV2=eig.tolist(),target_KK_quartet=J,
                lower_soft_row_over_total_mass=sn,bare_real_mass_over_total=math.sqrt(1-sn**2),
                real_down_yukawa_spectral_norm=math.sqrt(2)*float(np.linalg.norm(m,2))/v,
                J_at_L_300_GeV=finite(300.),J_at_L_3000_GeV=finite(3000.),
                schur_norm_difference_bound_for_lambda_le_target_bottom_GeV2=bound,
                total_weak_deficit_upper_for_3_light_at_true_pole=
                    3*nz*heavy**2/(heavy**2-eig[-1])**2)


def direct_matrix_test():
    rng=np.random.default_rng(1810)
    q,_=np.linalg.qr(rng.normal(size=(3,3)))
    m=q@np.diag([.004,.09,3.7])
    b=np.array([.3-.4j,.6+.4j,.15+.2j])
    M=14.
    matrix=np.block([[m,b[:,None]],[np.zeros((1,3)),np.array([[M]])]])
    U,ss,_=np.linalg.svd(matrix)
    Q=U[:,np.argsort(ss)]
    lam=np.sort(ss**2)
    J=float((Q[0,1]*Q[1,2]*Q[0,2].conj()*Q[1,1].conj()).imag)
    for j in range(3):
        val,x=lam[j],Q[:3,j]
        eff=m@m.T-val*np.outer(b,b.conj())/(M*M-val)
        np.testing.assert_allclose(eff@x,val*x,atol=2e-13)
    eta=np.vdot(b,b).real/(M*M-lam[2])
    for j in (1,2):
        coords=q.T@Q[:3,j]
        coords*=np.exp(-1j*np.angle(coords[j]))
        assert np.linalg.norm(coords.imag)<1.02*eta
    deficit=float(sum(abs(Q[:3,-1])**2))
    lb=float(np.vdot(b,b).real/(M*M+np.vdot(b,b).real))
    assert deficit>=lb and abs(J)>1e-9
    return dict(quartet=J,weak_deficit=deficit,deficit_lower=lb)


def run():
    cfg=read("inputs.json")
    datum=read(cfg["gauge_poles"])
    c=datum["common_MZ"]["two_doublets_one_triplet"]
    M=c["M_T_GeV"][0]
    cap=cfg["p_norm_cap"]
    old=datum["common_MZ"]["central"]
    # Same linear one-loop gauge solution as the central spectrum, with one low T.
    new_T=1000.
    ld=4*(math.log(c["matching_scale_GeV"]/new_T)-old["two_D_one_T_log_min_geometric_mean_D_over_T"])/7
    md=c["matching_scale_GeV"]*math.exp(-ld)
    assert 0<ld<old["log_scale"] and new_T<c["matching_scale_GeV"]
    eq=math.log(c["matching_scale_GeV"]/new_T)-7*ld/4
    assert abs(eq-old["two_D_one_T_log_min_geometric_mean_D_over_T"])<1e-12
    tables={str(n):tree_data(cfg,n) for n in (1,2)}
    v=cfg["vev_GeV"]
    eta=(cap*v)**2/(2*(M*M-25.))
    for x in (tables["2"]["raw_KK"],tables["2"]["registered_gamma_db"]):
        x["required_p_norm_at_registered_central_M"]=math.sqrt(2*x["lower_eta"]*(M*M-25.))/v
    pole_bands=read("../1809/results.json")
    both_floor=min(pole_bands[k]["original_window"]["B_T_lower_GeV"][0]
                   for k in ("common_MZ","independent_lattice"))
    return dict(bound_tables=tables,central_heavy=dict(M_GeV=M,eta_upper=eta,
        degenerate_doublets_in_1TeV_T_gauge_counterexample_GeV=md,
        KK_max_phase_for_perturbative_p=upper_angles(eta,.1,.05,tables["2"]["moduli_lower"],"KK")),
        two_complete_pairs_all_mass_bands=dict(minimum_T_GeV=both_floor,
            required_combined_p_operator_norm=math.sqrt(2*tables["2"]["raw_KK"]["lower_eta"]*(both_floor**2-25))/v),
        empirical_lower_block_competitor=lower_completion(read(cfg["yukawa_moduli"]),v,M),
        finite_upper_block_test=direct_matrix_test())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=run()
    raw=json.dumps(out,indent=2,allow_nan=False)+"\n"
    dest=HERE/"results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==out
    else:
        dest.write_bytes(raw.encode())
    print("1810",hashlib.sha256(raw.encode()).hexdigest()[:12])

