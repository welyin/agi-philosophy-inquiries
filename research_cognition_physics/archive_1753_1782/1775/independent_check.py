"""Finite-window index certification under a bounded symmetry-breaking mass operator."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LIGHT_MAX = 5.0
SW2 = 0.23153                # Illustrative common tree-level reference, not a matched fit
YUKAWA_CAP = math.sqrt(4*math.pi)


def trace_error_bound(q, low, high, mass_shift):
    if high-low <= 2*mass_shift:
        raise ValueError("The symmetry-preserving reference cluster is not certified.")
    return 2*q*(mass_shift/(high-mass_shift-low))**2


def sign_test(data, error):
    eps = [data["heavy_modes"].get(k, {}).get("left_doublet_probability", 0) for k in "dsb"]
    em = SW2/3
    g_left = [-0.5*(1-x)+em for x in eps]
    assert min(abs(x) for x in g_left) > em
    branches = []
    for signs in itertools.product([1,-1],repeat=3):
        # For downward charges the left sign is fixed by 0<=doublet projection<=1.
        tr = 2*sum(s*em-a for s,a in zip(signs,g_left))
        dist = abs(tr-round(tr))
        branches.append(dict(signs=signs,visible_trace=tr,distance_to_integer=dist,
                             compatible_integers=list(range(math.ceil(tr-error),math.floor(tr+error)+1))))
    alternative = min(b["distance_to_integer"] for b in branches[1:])
    assert [b["compatible_integers"] for b in branches] == [[3]] + [[]]*7
    return dict(weak_trace=branches[0]["visible_trace"],by_number_flipped=[
                    branches[0]["visible_trace"] - n*4*em for n in range(4)],
                closest_wrong_integer_distance=alternative,
                uniform_absolute_vertex_error_open_upper_limit=(alternative-error)/12)


def product_unitary(size, rg):
    x = rg.normal(size=(size,size)) + 1j*rg.normal(size=(size,size))
    return np.linalg.qr(x)[0]


def stress():
    # Rectangular mass blocks at v=0. Include a light vector pair and four signed indices.
    rg = np.random.default_rng(1775)
    rows=[]
    for left, right in [(5,2),(4,2),(2,4),(4,4)]:
        dimension = left + right + 2
        sl, sr = dimension-left, dimension-right
        gauge = np.zeros((dimension,dimension),complex)
        for block, dims, light in [(np.s_[:left,:right], (left,right), True),
                                   (np.s_[left:,right:], (sl,sr), False)]:
            a,b = dims
            singular=np.zeros((a,b))
            values=np.linspace(1.2,2.0,min(a,b))
            if light:
                values[0]=0.04
            np.fill_diagonal(singular,values)
            gauge[block] = product_unitary(a,rg) @ singular @ product_unitary(b,rg)
        pl=np.diag(np.arange(dimension)<left)
        pr=np.diag(np.arange(dimension)<right)
        assert np.linalg.norm(pl@gauge-gauge@pr)<1e-12
        q=abs(left-right)+1
        worst=0.
        for _ in range(12):
            breaking=rg.normal(size=gauge.shape)+1j*rg.normal(size=gauge.shape)
            breaking[:left,:right]=0
            breaking[left:,right:]=0
            op=rg.uniform(0.04,0.30)
            breaking *= op/np.linalg.norm(breaking,2)
            mass=gauge+breaking
            ward=pl@mass-mass@pr
            np.testing.assert_allclose(ward,(2*pl-np.eye(dimension))@breaking,atol=1e-12)
            np.testing.assert_allclose(np.linalg.norm(ward,2),op,rtol=1e-12)
            u,s,vh=np.linalg.svd(mass)
            v=vh.conj().T
            np.testing.assert_allclose(u.conj().T@ward@v,
                (u.conj().T@pl@u)*s-s[:,None]*(v.conj().T@pr@v),atol=1e-12)
            mid = np.trace(u[:,-q:].conj().T@pl@u[:,-q:]).real
            mid -= np.trace(v[:,-q:].conj().T@pr@v[:,-q:]).real
            bound=trace_error_bound(q,s[-q],s[-q-1],op)
            err=abs(mid-left+right)
            assert err <=bound+1e-12
            worst=max(worst,err/bound)
        rows.append(dict(index=left-right,view_rank=q,cases=12,largest_error_over_bound=worst))
    return rows


def actual_ratio(central, stat, syst):
    width=math.hypot(stat,syst)
    ends=[min(1.0,central+width),max(-1.0,central-width)]
    return [math.sqrt((1-a)/(1+a)) for a in ends]


def calculate():
    saved=json.loads((HERE/"../1768/independent_results.json").read_bytes())
    high=saved["common_mass_GeV"]
    v=saved["vev_GeV"]
    cap=v*YUKAWA_CAP/math.sqrt(2)
    error=trace_error_bound(3,LIGHT_MAX,high,cap)
    # Spread one symmetry-breaking direction over individual Yukawa entries, preserving its norm.
    n=144
    direction=np.ones(n)/math.sqrt(n)
    projector=np.eye(n)-np.outer(direction,direction)
    basis=np.linalg.eigh(projector)[1][:,1:].T
    bare=np.vstack((np.zeros(n),high*basis))
    update=np.zeros_like(bare)
    update[0]=high*direction
    np.testing.assert_allclose(np.linalg.svd(bare+update,compute_uv=False),high,rtol=1e-12)
    assert np.max(abs(update)) < cap < np.linalg.norm(update,2)
    # Equal poles and equal total radial Higgs slope for different indices.
    rotation=np.array([[0.,1.],[-1.,0.]])
    def poles(t):
        return np.linalg.svd(np.eye(2)+t*rotation,compute_uv=False)
    delta=1e-5
    slopes=(np.log(poles(1+delta))-np.log(poles(1-delta)))/(2*delta)
    np.testing.assert_allclose(poles(1),[math.sqrt(2)]*2)
    np.testing.assert_allclose(slopes,[0.5,0.5],rtol=1e-8)
    return dict(q_light=3,light_mass_ceiling_GeV=LIGHT_MAX,assumed_all_other_poles_above_GeV=high,
                higgs_vev_GeV=v,assumed_operator_yukawa_cap=YUKAWA_CAP,mass_perturbation_cap_GeV=cap,
                maximum_visible_axial_trace_error=error,
                maximum_certifying_mass_perturbation_for_half_integer_GeV=(high-LIGHT_MAX)/(1+math.sqrt(12)),
                common_tree_angle_sw2=SW2,
                by_calibration={k:{name:sign_test(row,error) for name,row in rows.items()}
                                for k,rows in saved["by_calibration"].items()},
                measured_polarized_abs_right_over_left_marginal_one_width={
                    "s_SLD":actual_ratio(.895,.066,.062), "b_SLD_lepton_tag":actual_ratio(.919,.030,.024),
                    "same_tree_reference":(SW2/3)/(.5-SW2/3)},
                hidden_row_countermodel=dict(added_massive_lines=n,added_net_doublets=1,mass_each_GeV=high,
                    max_individual_yukawa=math.sqrt(2)*high/(v*math.sqrt(n)),operator_yukawa=math.sqrt(2)*high/v),
                scalar_slope_counterexample=dict(equal_masses_in_units_of_M=[math.sqrt(2)]*2,
                    vectorlike_index=0,vectorlike_slopes=[0.5,0.5],chiral_plus_singlet_index=1,
                    chiral_plus_singlet_slopes=[1.,0.],common_sum=1.),
                numerical_mixed_spectral_checks=stress())


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [clean(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,float) else x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=clean(calculate())
    payload=json.dumps(results,indent=2)+"\n"
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==results
    else:
        path.write_text(payload,encoding="utf-8")
    print("1775 finite mass-window chiral bound",hashlib.sha256(payload.encode()).hexdigest()[:12])
