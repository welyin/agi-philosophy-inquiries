"""R1696. An independent pion+missing-energy test of the inferred gauge-boson scale."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PI = math.pi
SM = np.array([41/10, -19/6, -7.])
MSSM = np.array([33/5, 1., -3.])
PROJECTOR = np.array([5., -3., -2.])


def read(name):
    return json.loads((HERE/name).read_bytes())


def matching(gauge, mz, color, susy):
    b = MSSM if susy else SM
    normalization = 1. if susy else 2/3
    beta_t = normalization*np.array([.4, 0., 1.])
    beta_l = normalization*np.array([.6, 1., 0.])
    ell_s = math.log(1000./mz) if susy else 0.
    inverse = 4*PI/gauge**2 - (np.array([0., 2., 3.])/(12*PI) if susy else 0.)
    zt = sum(math.log(m/mz) for m in color)
    rhs = inverse - (SM-b)*ell_s/(2*PI) + beta_t*zt/(2*PI)
    mat = np.column_stack([np.ones(3), (b+2*beta_t+2*beta_l)/(2*PI), -beta_l/(2*PI)])
    a, ell, zl = np.linalg.solve(mat, rhs)
    M = mz * math.exp(ell)
    np.testing.assert_allclose(ell, PROJECTOR@(2*PI*inverse-(SM-b)*ell_s)/(PROJECTOR@b))
    # Gradients w.r.t. the measured inverse couplings, for certification over the box.
    da, de, _ = np.linalg.inv(mat)
    return dict(M=M, ainv=a, weak_geomean=mz*math.exp(zl/2), da=da, de=de)


def ckm(t12,t13,t23,delta):
    c,d,e=np.cos([t12,t13,t23]);s,t,u=np.sin([t12,t13,t23])
    phase=np.exp(1j*delta)
    return np.array([[c*d, s*d, t*phase.conjugate()],
                     [-s*e-c*u*t*phase, c*e-s*u*t*phase, u*d],
                     [s*u-c*e*t*phase,-c*u-s*e*t*phase,e*d]])


def baryon_transport(gauge,mz,center,color,light_weak,W,vud):
    """Leading dimension-six gauge transport on the same piecewise spectrum."""
    soft=1000.
    common=center["weak_geomean"]
    weak=np.array([light_weak,common**2/light_weak])
    M=center["M"]
    inv=4*PI/gauge**2
    below=np.array([11/20,9/4,2.])
    susy=np.array([11/30,3/2,4/3])
    log_transfer=0.
    grid=sorted({mz,soft,*color,*weak,M})
    for l,r in zip(grid,grid[1:]):
        if l==soft:
            inv=inv-np.array([0.,2.,3.])/(12*PI)
        mid=math.sqrt(l*r)
        b=(MSSM if mid>soft else SM).copy()
        b+=np.array([.4,0.,1.])*sum(color<mid)
        b+=np.array([.6,1.,0.])*sum(weak<mid)
        gamma=susy if mid>soft else below
        dt=math.log(r/l)/(2*PI)
        after=inv-b*dt
        assert min(after)>0
        integrand=np.repeat(dt,3)/inv
        nonzero=b!=0
        integrand[nonzero]=np.log(inv[nonzero]/after[nonzero])/b[nonzero]
        log_transfer+=gamma@integrand
        inv=after
    np.testing.assert_allclose(inv,center["ainv"],rtol=2e-12)
    factor=1.247*math.exp(log_transfer)
    mx=M;g2=4*PI/center["ainv"]
    width=.9382721/(32*PI)*(1-(.1395704/.9382721)**2)**2*(factor*g2/mx**2*W*vud)**2
    mass_matrix=np.array([[0.,common],[-common,weak[1]-weak[0]]])
    np.testing.assert_allclose(np.sort(np.linalg.svd(mass_matrix,compute_uv=False)),np.sort(weak))
    return dict(weak_poles_GeV=weak.tolist(),weak_bare_diagonal_GeV=weak[1]-weak[0],
                total_transfer=factor,partial_life_pi_nu_year=6.582119569e-25/width/(86400*365.25))

def compute():
    inp = read("../1369/inputs.json")
    gauge = np.array([inp[k] for k in ("g1_GUT","g2","g3")])
    color = np.array(read("../1685/results.json")["conditional_masses_TeV_descending"])*1000
    weak = np.array(read("../1397/inputs.json")["high_scale"]["angles_theta12_theta13_theta23_delta"])
    vud_floor=math.prod(math.cos(c+2*err) for c,err in weak[:2])
    vud_ceiling=math.prod(math.cos(c-2*err) for c,err in weak[:2])
    central_vud=math.prod(math.cos(c) for c in weak[:2,0])
    life_years, hbar, seconds_year=3.5e32,6.582119569e-25,365.25*86400
    mp, mpi = .9382721, .1395704
    W_floor, W_central=.098,.159
    constant = math.sqrt(32*PI*hbar / (life_years*seconds_year*mp*(1-(mpi/mp)**2)**2))
    out = {}
    for susy in (False,True):
        def at(g):
            return matching(g,inp["mu_ref_GeV"],color,susy)
        points = [at(gauge[:,0]+2*gauge[:,1]*sign)
                  for sign in itertools.product((-1,1), repeat=3)]
        center=at(gauge[:,0])
        ainv=[q["ainv"] for q in points]
        for val in (min(ainv),max(ainv)):
            slopes=2*center["de"] + center["da"]/val
            assert np.all(np.sign(slopes)==np.sign(PROJECTOR))
        def limit(row,W,vud):
            return constant*row["M"]**2*row["ainv"]/(4*PI*W*vud)
        # This envelope uses the same (conservative) W floor throughout.
        budget=[min(limit(q,W_floor,vud_ceiling) for q in points),
                max(limit(q,W_floor,vud_floor) for q in points)]
        central=limit(center,W_central,central_vud)
        out["MSSM_extra_unsplit_at_1TeV_soft" if susy else "SM_plus_extra_Dirac"] = dict(
            central_M_GeV=center["M"], central_gU_squared=4*PI/center["ainv"],
            central_formal_weak_geometric_mean_GeV=center["weak_geomean"],
            permission_at_reference_lattice_and_gauge_centers=central,
            two_width_gauge_box_max_allowed_transport_at_W_floor=budget,
            absolute_min_MX_over_crossing_if_unit_transport=1/math.sqrt(budget[1]),
            central_lifetime_years_if_MX_equals_crossing_and_unit_transport=life_years*central**2)
        if susy:
            key="MSSM_extra_unsplit_at_1TeV_soft"
            curves=[baryon_transport(gauge[:,0],inp["mu_ref_GeV"],center,color,l,W_central,central_vud)
                    for l in (center["weak_geomean"],1000.,inp["mu_ref_GeV"])]
            factors=[q["total_transfer"] for q in curves]
            assert factors[0]<factors[1]<factors[2]
            out[key]["explicit_gauge_only_spectra"]=curves
            out[key]["full_allowed_split_width_increase_at_fixed_product"]=(factors[-1]/factors[0])**2-1
    V=ckm(*weak[:,0]); vud=abs(V[0,0])
    rng=np.random.default_rng(1696)
    worst=0.
    for _ in range(48):
        q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        cvec=V[0,0]*q[0,:]
        worst=max(worst,abs(np.vdot(cvec,cvec).real-vud**2))
    P=np.array([[0,0,1],[1,0,0],[0,1,0]])
    W=P@V.conj().T
    assert abs(np.linalg.det(W)-1)<2e-15 and np.max(abs((W@V)[0,:2]))<2e-15
    assert worst<2e-15
    return dict(running_cases=out, symmetric_up_type_Vud_window=[vud_floor,vud_ceiling],
                inclusive_flavor_sum_tests=48,mixed_chirality_lattice_W_floor_GeV2=W_floor,
                recent_SuperK_pi_plus_neutrino_years=life_years,
                allowed_norm_of_low_energy_scalar_coefficient_GeV_minus2=constant/W_floor,
                competing_unit_determinant_right_rotation_first_row_abs=np.round(abs((W@V)[0,:]),12).tolist())


def rounded(d):
    if isinstance(d,dict):
        return {k: rounded(v) for k,v in d.items()}
    if isinstance(d,list):
        return [rounded(x) for x in d]
    return float(f"{d:.11g}")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    check=p.parse_args().check
    result=rounded(compute())
    target=HERE/"independent_results.json"
    if check:
        assert result==read(target)
        print("R1696 independent OK")
    else:
        target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
