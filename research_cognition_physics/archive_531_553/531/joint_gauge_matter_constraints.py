"""Round 531: exact joint gauge/matter constraints and a historical RG benchmark.

No spacetime dimension, full spectral triple, quantum gravity, or empirical fit
is derived. Known anomaly tools are reused under the report's explicit inputs.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_gauge_matter_constraints_results.json'


def add(*ps):
    out = {}
    for p in ps:
        for monomial, value in p.items():
            out[monomial] = out.get(monomial, F(0)) + value
    return {key: value for key, value in out.items() if value}


def mul(a, b):
    out = {}
    for p, u in a.items():
        for q, v in b.items():
            key = tuple(x+y for x, y in zip(p, q))
            out[key] = out.get(key, F(0)) + u*v
    return {key: value for key, value in out.items() if value}


def scale(a, c):
    return {p: c*v for p, v in a.items() if c*v}


def cube(a):
    return mul(mul(a, a), a)


def symbolic_anomalies():
    # Independent symbolic variables N,q,h, with no fixed color count.
    N, q, h = ({(1,0,0): F(1)}, {(0,1,0): F(1)}, {(0,0,1): F(1)})
    Nq = mul(N, q)
    u, d, l, e, n = add(scale(q,-1),scale(h,-1)), add(scale(q,-1),h), scale(Nq,-1), add(Nq,h), add(Nq,scale(h,-1))
    cubic = add(scale(mul(N,cube(q)),2), mul(N,cube(u)), mul(N,cube(d)),scale(cube(l),2),cube(e))
    grav = add(scale(Nq,2),mul(N,u),mul(N,d),scale(l,2),e)
    residual = add(h,scale(Nq,-1))
    assert grav == residual and cubic == cube(residual)
    assert add(grav,n) == {} and add(cubic,cube(n)) == {}
    assert add(scale(q,2),u,d) == {} and add(Nq,l) == {}
    return {'without_singlet': 'A_grav=h-Nq; A_Y3=(h-Nq)^3',
            'with_Dirac_singlet': 'both identically zero for n=Nq-h',
            'bare_Majorana': 'n=0 gives h=Nq', 'symbolic_all_N': True}


def fields(N):
    # name: (multiplicity, Y, quark-sector weight label)
    return {'Q': (2*N,F(1,2*N),True), 'uc': (N,F(-N-1,2*N),True),
            'dc': (N,F(N-1,2*N),True), 'L': (2,F(-1,2),False),
            'ec': (1,F(1),False), 'nuc': (1,F(0),False)}


def indices(N, wq=F(1), wl=F(1)):
    return (2*wq, (N*wq+wl)/2, wq*(F(N,2)+F(1,N))+3*wl/2)


def index_from_diagonal_generators(N, wq, wl):
    # Explicit representation diagonal entries, Tr_fund(T^2)=1/2.
    c = [F(1,2), F(-1,2)] + [F(0)]*(N-2)
    color = [x for x in c for _ in range(2)] + [-x for x in c]*2 + [F(0)]*4
    weak = [x for _ in range(N) for x in (F(1,2),F(-1,2))] + [F(0)]*(2*N) + [F(1,2),F(-1,2),F(0),F(0)]
    hyper = [y for multiplicity,y,_ in fields(N).values() for _ in range(multiplicity)]
    weights = [wq]*(4*N) + [wl]*4
    assert len(color) == len(weak) == len(hyper) == 4*N+4
    return tuple(sum(w*x*x for w,x in zip(weights, xs)) for xs in (color, weak, hyper))


def center_kernel(N):
    # Integer charges 2N*Y; theta/(2pi)=k/(2N), fixed by ec charge 2N.
    out = []
    for a in range(N):
        for b in range(2):
            for k in range(2*N):
                t = F(k,2*N)
                phases = (F(a,N)+F(b,2)+t,
                          F(-a,N)+(-N-1)*t, F(-a,N)+(N-1)*t,
                          F(b,2)-N*t, 2*N*t, F(b,2)+N*t)
                if all(p.denominator == 1 for p in phases):
                    out.append((a,b,k))
    return set(out)


def historical_running():
    # Fixed 2012 review benchmark, NOT current observations. Three generations,
    # N=3, one Higgs doublet, one-loop big-desert running, no thresholds.
    g0 = np.array([0.3575,0.6519,1.220], dtype=float)  # Y, weak, color
    beta = np.array([-41/6,19/6,7.0])
    x0 = 1/g0**2
    slope = beta/(8*math.pi**2)
    normal = np.array([-1.0,3.0,-4/3])
    t = -float(normal@x0)/float(normal@slope)
    x = x0+t*slope
    # Absorb overall kappa into wq and wl; exact inverse cone map.
    wq, wl = x[2]/2, 2*x[1]-3*x[2]/2
    expected = np.array(indices(3,wq,wl), dtype=float)[[2,1,0]]
    assert np.all(x>0) and wq>0 and wl>0
    assert np.max(np.abs(x-expected)) < 1e-12
    t23 = (x0[2]-x0[1])/(slope[1]-slope[2])
    x23 = x0+t23*slope
    # Same scale gc=gw does not also give xY=(5/3)xw in this fixed benchmark.
    mismatch = float(x23[0]-(5/3)*x23[1])
    assert abs(mismatch)>0.5 and abs(x23[1]-x23[2])<1e-12
    # Independent linear solve in unknown (t, wq, wl).
    columns = np.array([[11/6,3/2],[3/2,1/2],[2,0]],dtype=float)
    sol = np.linalg.solve(np.column_stack((slope,-columns)), -x0)
    assert np.max(np.abs(sol-np.array([t,wq,wl])))<1e-10
    def rd(v): return round(float(v),10)
    return {'source': '1204.0328 section 8.2, fixed historical 2012 central values',
            'input_g_Y_weak_color': g0.tolist(), 'MZ_GeV': 91.1876,
            'beta_minus_sign_convention': beta.tolist(),
            'weighted_trace_log_scale_ratio': rd(t),
            'weighted_trace_scale_GeV': rd(91.1876*math.exp(t)),
            'positive_weight_ratio_lepton_over_quark': rd(wl/wq),
            'matched_inverse_couplings': [rd(z) for z in x],
            'color_weak_equal_scale_GeV': rd(91.1876*math.exp(t23)),
            'hypercharge_mismatch_at_color_weak_equality': rd(mismatch),
            'current_empirical_fit': False, 'complete_spectral_action_realization': False}


def run():
    checks = []
    symbolic = symbolic_anomalies()
    checks.append('symbolic_anomaly_identities_all_N')
    for N in range(3,32):
        fs = fields(N)
        assert sum(m*y for m,y,_ in fs.values()) == 0
        assert sum(m*y**3 for m,y,_ in fs.values()) == 0
        q,u,d,l,e,n = [v[1] for v in fs.values()]
        assert 2*q+u+d == N*q+l == 0
        assert q+F(1,2)+u == q-F(1,2)+d == l-F(1,2)+e == l+F(1,2)+n == 0
    checks.append('explicit_chiral_multiplicities_yukawa_and_local_anomalies')
    assert all((N+1)%2 == 0 for N in range(3,32,2))
    assert (4+1)%2 == 1 and 2*(4+1)%2 == 0
    checks.append('usual_SU2_global_anomaly_and_even_generation_control')
    for N in (3,4,5,7,11):
        for wq,wl in ((F(1),F(1)),(F(2,3),F(7,5))):
            assert index_from_diagonal_generators(N,wq,wl) == indices(N,wq,wl)
    checks.append('independent_representation_generator_traces')
    for N in range(3,32):
        for wq,wl in ((F(1),F(1)),(F(2,3),F(7,5)),(F(9),F(1,7))):
            xc,xw,xy = indices(N,wq,wl)
            assert xy == 3*xw-F(N*N-1,2*N)*xc
            assert xw > F(N,4)*xc
            assert (xc/2,2*xw-F(N,2)*xc) == (wq,wl)
        boundary = indices(N,F(1),F(0))
        assert boundary[1] == F(N,4)*boundary[0]
    checks.append('positive_weight_cone_necessity_sufficiency_and_boundary')
    roots = [N for N in range(3,101) if 4-N>0]
    assert roots == [3] and indices(3,F(1),F(1))[0:2] == (F(2),F(2))
    assert indices(5)[0] != indices(5)[1]
    checks.append('additional_equal_coupling_selects_N3_within_ansatz')
    kernels = {}
    for N in (3,5,7,9):
        kernel = center_kernel(N)
        cyclic = {(((N-1)//2*k)%N,k%2,k) for k in range(2*N)}
        assert kernel == cyclic
        kernels[str(N)] = len(kernel)
    assert len(center_kernel(4)) == 4  # odd-N formula must not be used here.
    checks.append('exact_center_kernel_with_even_N_negative_control')
    rows = []
    for N in (3,5,7):
        ic,iw,iy = indices(N)
        assert ic>0 and iw>0 and iy>0
        rows.append({'N':N,'charges':{k:str(v[1]) for k,v in fields(N).items()},
            'indices_color_weak_Y':[str(v) for v in (ic,iw,iy)],
            'gc_squared_over_gw_squared':str(iw/ic),
            'gY_squared_over_gw_squared':str(iw/iy),
            'sine_squared_weak_angle_at_matching':str(iw/(iw+iy))})
    assert rows[1]['gc_squared_over_gw_squared']=='3/2'
    checks.append('odd_color_competitors_survive_without_extra_matching_equality')
    running = historical_running()
    checks.append('historical_RG_positive_trace_intersection_and_unweighted_mismatch')
    dependencies = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
        'unified_physics_condition_ledger.md','unified_physics_condition_audit_results.json',
        'research_note_339.md','research_note_344.md','research_note_360.md','research_note_490.md')}
    return {'round':531,'tests_run':len(checks),'failures':0,'errors':0,
            'checks':checks,'symbolic':symbolic,'known_color_family_reused':rows,
            'center_kernel_sizes':kernels,'historical_running':running,
            'dependency_hashes':dependencies,
            'scope':{'common_four_dimensional_EFT_ansatz_supplied':True,
                     'joint_positive_trace_compatibility_proved':True,
                     'complete_anomaly_classification':False,
                     'spacetime_dimension_or_GR_generated':False,
                     'complete_SM_selected_from_cognition':False,
                     'full_NCG_spectral_triple_or_new_empirical_fit':False}}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(result,ensure_ascii=False,indent=2))
