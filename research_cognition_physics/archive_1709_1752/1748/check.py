"""Use two actual neutral columns and a common real charged analyser for two heavy mass rows."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
OLD=runpy.run_path(str(HERE/CFG["matter"]))


def rational(x):
    # Exact certificates for the printed numerical mass/source input, not arbitrary UV values.
    return Fraction(str(float(x)))


def equations(bg):
    Y=bg["yW"][2,:,:3]
    L,y,_=np.linalg.svd(Y)
    L,y=L[:,::-1],y[::-1]
    N=L.conj().T@bg["N"]
    N=N/np.linalg.norm(N,axis=0)
    A,B=N[:,0],N[:,1]
    t=1/np.array([bg["cW"],bg["fW"]])**2-1
    gain=(bg["cW"]*bg["fW"])**(1/3)
    scale=sum(t)
    h=y**2
    left=1/h
    left/=sum(left)
    upper=h/sum(h)
    # w=t_ii/scale. The trace and inverse charged trace fix w_e,w_mu given w_tau=v.
    a=( (gain**-2-1)/scale -left[1])/(left[0]-left[1])
    b=(left[1]-left[2])/(left[0]-left[1])
    # Rational polynomial in v or odds r=z/(1-z), as appropriate.
    def p(*x):
        return Polynomial(np.array([rational(v) for v in x],dtype=object))
    w=[p(a,b),p(1-a,-1-b),p(0,1)]
    r=p(0,1)
    u=[p(abs(B[i])**2,abs(A[i])**2) for i in range(3)]
    denom=u[0]*u[1]*u[2]
    w_upper=(gain**2*np.prod(1+t)-1-scale)/scale
    E=[-rational(np.prod(t)/scale**2)*denom,p(0),p(0)]
    H=[-rational(w_upper)*denom,p(0),p(0)]
    for i,(j,k) in enumerate(((1,2),(2,0),(0,1))):
        cross=r*u[i]*rational(abs(A[j]*B[k]-A[k]*B[j])**2)
        for power,weight in enumerate((w[j]*w[k]).coef):
            E[power]+=cross*weight
            H[power]+=cross*weight*rational(scale*upper[i])
        for power,weight in enumerate(w[i].coef):
            H[power]-=denom*weight*rational(upper[i])
    c,b1,a1=E
    f,e,d=H
    resultant=(a1*f-c*d)**2-(a1*e-b1*d)*(b1*f-c*e)
    signs=[int(np.sign(q)) for q in resultant.coef if q!=0]
    assert all(isinstance(q,Fraction) for q in resultant.coef)
    assert sum(x!=y for x,y in zip(signs,signs[1:]))==2
    poly=Polynomial([float(x) for x in resultant.coef])
    roots=[]
    for lower,higher in ((Fraction(1,2),Fraction(1)),(Fraction(20),Fraction(23))):
        assert resultant(lower)*resultant(higher)<0
        lo,hi=float(lower),float(higher)
        for _ in range(62):
            mid=(lo+hi)/2
            if poly(lo)*poly(mid)<=0:
                hi=mid
            else:
                lo=mid
        odds=(lo+hi)/2
        E0,E1,E2=(float(q(odds)) for q in E)
        H0,H1,H2=(float(q(odds)) for q in H)
        v=(H2*E0-E2*H0)/(E2*H1-H2*E1)
        z=odds/(1+odds)
        x=scale*np.array([float(q(v)) for q in w])
        assert min(x)>0
        gains=np.sqrt(x/(z*abs(A)**2+(1-z)*abs(B)**2))
        roots.append(dict(odds=odds,source_power_A=z,w_tau=v,gains=gains))
    return L,y,N,roots,dict(descartes_signs=signs,root_brackets_odds=[[.5,1],[20,23]],
                            c=bg["cW"],fb=bg["fW"],common_charged_gain=gain)


def terms(bg, root=0, signs=(1,-1,1)):
    """Return the two-row dimensionless X and physical positive F in the 1744 weak basis."""
    L,y,N,roots,cert=equations(bg)
    sol=roots[root]
    z=sol["source_power_A"]
    X=(np.sqrt([z,1-z])[:,None]*N.conj().T)@np.diag(sol["gains"]*signs)@L.conj().T
    eig,U=np.linalg.eigh(X.conj().T@X)
    F=(U/np.sqrt(1+eig))@U.conj().T
    np.testing.assert_allclose(np.linalg.eigvalsh(F),[bg["cW"],bg["fW"],1.],atol=8e-11)
    eig2,V=np.linalg.eigh(X@X.conj().T)
    normalized=(V/np.sqrt(1+eig2))@V.conj().T@np.hstack([X,np.eye(2)])
    light=np.vstack([np.eye(3),-X])@F
    np.testing.assert_allclose(normalized@normalized.conj().T,np.eye(2),atol=2e-13)
    np.testing.assert_allclose(normalized@light,0,atol=2e-13)
    Y=bg["yW"][2,:,:3]
    _,ss,Vh=np.linalg.svd(F@Y)
    charged_gains=ss[::-1]/y
    np.testing.assert_allclose(charged_gains,cert["common_charged_gain"],atol=3e-10)
    Ue=Vh.conj().T[:,::-1]
    heavy=(Ue.conj().T@Y.conj().T@(np.eye(3)-F@F)@Y@Ue)/np.outer(ss[::-1],ss[::-1])
    return X,F,dict(**sol,signs=signs,charged_gains=charged_gains,
                    normalized_heavy_gram_magnitudes=np.abs(heavy),
                    F_real=F.real,F_imag=F.imag,physical_mass_rows_over_M_real=normalized.real,
                    physical_mass_rows_over_M_imag=normalized.imag),cert


def daughter(bg,F):
    a=1.
    target=bg["model"]["JOINT"]["LEPTONS"][2,0]
    for _ in range(3):
        b=dict(bg)
        b["yW"]=bg["yW"].copy()
        b["yW"][2]*=a
        row=OLD["daughter"](b,F)
        a*=target/row["MZ_leptons"][2]
    assert abs(target/row["MZ_leptons"][2]-1)<1e-11
    return a,row


def calculate(step):
    bg=OLD["mother"](step)
    examples=[]
    freq=bg["nu"]["FREQ"]
    for j in range(2):
        for signs in CFG["family_signs"]:
            _,F,port,cert=terms(bg,j,signs)
            a,row=daughter(bg,F)
            o=row["neutrinos"]
            def pull(x,ref):
                mid,low,high=freq[ref]
                return (x-mid)/(low if x<mid else high)
            clocks=[pull(o["dm31_eV2"],"juno_NO_dm31"),pull(o["dm32_eV2"],"accelerator_NO_dm32")]
            # Even when letting the solar normalization roam over its own two-width range.
            sol,sl,sh=freq["juno_NO_dm21"]
            ratio=o["dm31_eV2"]/sol
            lower=[sol-2*sl]
            upper=[sol+2*sh]
            for k,key in [(ratio,"juno_NO_dm31"),(ratio-1,"accelerator_NO_dm32")]:
                m,l,h=freq[key]
                lower.append((m-2*l)/k)
                upper.append((m+2*h)/k)
            free=max(lower)<=min(upper)
            examples.append(dict(root=j+1,port=port,charged_Mw_scale=a,
                                 angle_readings=o["angle_readings"],angle_in_boxes=o["angles_in_box"],
                                 solar_center_dm31=o["dm31_eV2"],solar_center_dm32=o["dm32_eV2"],
                                 atmospheric_marginal_offsets=clocks,
                                 clocks_both_with_solar_at_center=max(abs(np.array(clocks)))<=2,
                                 clocks_both_when_solar_varies_two_widths=free,
                                 admissible_solar_interval_if_nonempty=[max(lower),min(upper)] if free else None,
                                 solar_with_mother_normalization=o["dm21_with_mother_solar_gain"],
                                 delta_deg=o["delta_deg"],SK_T2K_68_90=o["SK_T2K_FC68_90"],
                                 mbb_meV=o["mbb_meV"],MZ_relative_e_mu_after_tau=row["MZ_relative_e_mu_after_tau_calibration"],
                                 quark_four_mass_radius=row["absolute_four_ratio_radius"],epsilonK_pull=row["epsilon_K_pull"],
                                 golden_S=row["S_Bd"],active_source_residual=row["matching_residual_neutral_columns"]))
    winners=[(a["root"],a["port"]["signs"]) for a in examples if
             a["clocks_both_when_solar_varies_two_widths"] and all(a["angle_in_boxes"])]
    assert winners==[(1,[1,-1,1])]
    assert max(abs(np.array(examples[1]["atmospheric_marginal_offsets"])))<1
    return dict(certificate=cert,physical_pole_GeV=bg["Mweak"],outputs=examples)


def tidy(o):
    if isinstance(o,dict):return {k:tidy(v) for k,v in o.items()}
    if isinstance(o,np.ndarray):return tidy(o.tolist())
    if isinstance(o,(list,tuple)):return [tidy(v) for v in o]
    if isinstance(o,(bool,np.bool_)):return bool(o)
    if isinstance(o,float):return 0 if abs(o)<1e-14 else float(f"{o:.10g}")
    return o


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate(.125)
    other=calculate(.25)
    for x,y in zip(result["outputs"],other["outputs"]):
        assert max(abs(np.array(x["angle_readings"])-y["angle_readings"]))<1e-6
        assert abs(x["solar_center_dm31"]-y["solar_center_dm31"])<2e-9
    result["step_half_dm31_witness_shift"]=abs(result["outputs"][1]["solar_center_dm31"]-other["outputs"][1]["solar_center_dm31"])
    out=tidy(result)
    path=HERE/"results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1748 OK: exact rational sign count, two roots, eight signed rows, split-spectrum tests.")
    else:
        path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")

