"""Use the archived actual weak matrices, with external clocks only as holdouts."""
import argparse
import cmath
import itertools
import json
import math
import runpy
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def distance(x,y):
    return max(0.,x[0]-y[1],y[0]-x[1])


def solve():
    p=read("inputs.json")
    obs=read(p["weak_states"])
    source=read(p["source_model"])
    single=read("../1465/"+source["three_rates"])
    # Reuse the archived mass construction, independently build the absolute mixing Hamiltonian.
    build=runpy.run_path(str((HERE/p["mass_builder"]).resolve()))["build"]
    helpers=runpy.run_path(str((HERE/"../1465/check.py").resolve()))
    weak,lattice,rare=[read(s) for s in p["absolute_calibrations"]]
    gf=weak["const_GeV"]["GF"]
    mw=weak["const_GeV"]["MW"]
    mb=weak["const_GeV"]["M_Bs"]
    hbar=lattice["hbar_GeV_ps"]
    s0,es=weak["short_distance"]["S0_published"]
    eta,eeta=weak["short_distance"]["etaB"]
    f,ef=lattice["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    x,z=p["eta_Z_RGI_range"]
    eta0=p["alpha_s_ref"]**(6/23)
    assert x<eta0<z
    cbox=gf**2*mw**2/(4*math.pi**2)
    cz=gf/math.sqrt(2)
    wall=p["extra_decay_phase_cap_rad"]
    dmean,*dw=p["omega_ps_inverse"]
    rows={}
    for label,base in obs.items():
        F=complex(*base["Fsb"])
        V,_,_,_=build(F,source,single,helpers["hermitian_power"],helpers["complete_basis"])
        lt=V[2,1].conjugate()*V[2,2]
        lc=V[1,1].conjugate()*V[1,2]
        om=2*mb/(3*hbar)
        def mixing(se,ez):
            return cbox*se*lt**2+cz*ez*F**2
        zero=cbox*eta*s0*lt**2
        total=mixing(eta*s0,eta0)
        tree=total-zero
        ph=cmath.phase(total/lc**2)
        freq=om*f*f*abs(total)
        theta=cmath.phase(total/zero)
        p_mu=complex(*base["P_mu"])
        mu_angle=cmath.phase(p_mu/lt)-theta/2
        rate=[1e9*rare["Hmu_heavy"]*abs(p_mu)**2*(1-(1-w)*math.sin(mu_angle)**2)
              for w in (rare["calculated_light_response"],rare["run3_central_light_response"])]
        trials={}
        for width in (1,2):
            extremes=[mixing(se,ej) for se,ej in itertools.product(
                ((eta-width*eeta)*(s0-width*es),(eta+width*eeta)*(s0+width*es)),(x,z))]
            # Positive real cross term throughout these five archived angles implies radial monotonicity.
            assert (zero.conjugate()*tree).real>0
            angles=[cmath.phase(v/lc**2) for v in extremes]
            phi_interval=[min(angles),max(angles)]
            omega_interval=[om*(f-width*ef)**2*min(map(abs,extremes)),
                            om*(f+width*ef)**2*max(map(abs,extremes))]
            def observed(v):
                return [v[0]-width*sum(v[1:]),v[0]+width*sum(v[1:])]
            data={key:observed(v) for key,v in p["phase_inputs_rad"].items()}
            # Only these two experiments measure the same J/psi phi final state.
            a,b=data["LHCb_KK_2015_18"],data["CMS_KK_13TeV_final"]
            kk=[max(a[0],b[0]),min(a[1],b[1])]
            assert kk[0]<kk[1]
            need=distance(phi_interval,kk)
            cobs=observed(p["omega_ps_inverse"])
            trials[str(width)]={"phi_range":phi_interval,"omega_range":omega_interval,
                "clock_min_gap_ps_inverse":distance(omega_interval,cobs),
                "min_phase_shift_for_individual":{key:distance(phi_interval,bounds) for key,bounds in data.items()},
                "same_KK_phase_box":kk,"common_KK_min_decay_shift_rad":need,
                "mixing_fraction_required_after_decay_cap":math.sin(max(0.,need-wall)),
                "central_decay_shift_for_KK":max(kk[0]-ph,0.)+min(kk[1]-ph,0.)}
        needed_f=math.sqrt(dmean/(om*abs(total)))
        rows[label]={"lambda_c":[lc.real,lc.imag],"lambda_t":[lt.real,lt.imag],
                     "tree_to_top_modulus":abs(tree/zero),"top_only_frequency":om*f*f*abs(zero),
                     "nominal_frequency":freq,"bag_sqrt_f_required_GeV":needed_f,
                     "top_only_phi":cmath.phase(zero/lc**2),"new_mixing_phase":theta,
                     "nominal_phi":ph,"updated_mu_e9":rate,"boxes":trials}
        if label=="free_-0.36":
            point=mixing(eta*s0,x)
            fw=math.sqrt(dmean/(om*abs(point)))
            phi=cmath.phase(point/lc**2)
            effective=phi+.0195
            assert abs(fw-f)<ef and abs(om*fw**2*abs(point)-dmean)<1.e-12
            assert all(v[0]-sum(v[1:])<=effective<=v[0]+sum(v[1:]) for v in p["phase_inputs_rad"].values())
            adjusted_angle=cmath.phase(p_mu/lt)-cmath.phase(point/zero)/2
            mu=[1e9*rare["Hmu_heavy"]*abs(p_mu)**2*(1-(1-w)*math.sin(adjusted_angle)**2)
                for w in (rare["calculated_light_response"],rare["run3_central_light_response"])]
            assert 3.18<mu[0]<4.54 and 2.88<mu[1]<3.77
            rows[label]["joint_conditional_witness"]={"eta_Z":x,"f_sqrtB_GeV":fw,"common_decay_offset_rad":.0195,
                "effective_phase_rad":effective,"clock_ps_inverse":om*fw**2*abs(point),"mu_e9":mu}
    assert rows["free_0.36"]["nominal_phi"]>rows["free_-0.36"]["nominal_phi"]
    return {"eta_RR_RGI_central_LL":eta0,"rows":rows}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write-new",action="store_true")
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=json.loads(json.dumps(solve()),parse_float=lambda x:float(f"{float(x):.12g}"))
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert result==read("results.json")
        print("1468 author OK")
