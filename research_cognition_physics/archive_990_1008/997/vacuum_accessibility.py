"""997: common effective vacuum density and finite causal round-trip domain.
All geometry is stipulated mean-field FLRW. No actual cosmological fit or instrument.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import argparse,hashlib,json,math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'vacuum_accessibility_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def exact(v):return dict(exact=str(v),value=float(v))
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-9,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)

def inverse_sqrt_bounds(d,scale=10**12):
    """Exact rational enclosure; d>0, floor(scale/sqrt(d)) by integer arithmetic."""
    k=isqrt(scale*scale*d.denominator//d.numerator)
    lo=F(k,scale);hi=F(k+1,scale)
    assert lo*lo*d<=1 and hi*hi*d>1
    return lo,hi

def run():
    old=read(STAGE/'992/cooling_availability_results.json')
    ip=old['inputs'];M=F(str(ip['material_mean_mass']));R=F(str(ip['initial_R']))
    V0=F(str(ip['V0']));L=M # New illustrative vacuum ENERGY in the initial cell.
    q=M/L;r=R/L;af=F(4)
    # H_Lambda*Delta eta = integral_{1/af}^1 dx/sqrt(1+q*x^3+r*x^4).
    n=512;dx=(1-1/af)/n;points=[1/af+j*dx for j in range(n+1)]
    boxes=[inverse_sqrt_bounds(1+q*x**3+r*x**4) for x in points]
    lower=dx*sum((lo for lo,hi in boxes[1:]),F(0))
    upper=dx*sum((hi for lo,hi in boxes[:-1]),F(0))
    assert 0<lower<upper<1
    x,w=np.polynomial.legendre.leggauss(96)
    xx=.25+.75*(x+1)/2
    value=float(.75/2*np.dot(w,1/np.sqrt(1+float(q)*xx**3+float(r)*xx**4)))
    assert float(lower)<value<float(upper)
    # Conditional future tail, using only positive matter/radiation and constant L.
    tl,tu=inverse_sqrt_bounds(1+q/af**3+r/af**4)
    future_lo=lower+tl/af;future_hi=upper+1/af
    assert future_hi<1 # Also bounded by 1 directly from H>=H_Lambda.
    close=F(1,5);oneway=F(11,20);separate=F(21,10)
    assert 2*close<lower
    assert oneway<lower and 2*oneway>1
    assert separate>2 # Disjoint future cones if both agents move and assumptions persist.
    # Same mean energy/source and material entropy preparation as 992.
    source_rows=[]
    for a in (F(1),F(2),F(4)):
        E=M+R/a+L*a**3
        pressure=R/(3*V0*a**4)-L/V0
        dE=-R/a**2+3*L*a**2
        assert dE+pressure*3*V0*a*a==0
        # (a_ddot/a)/H_Lambda^2, positive from a=1 in this chosen example.
        accel=1-q/(2*a**3)-r/a**4
        assert accel>0
        source_rows.append(dict(a=int(a),cell_energy=float(E),pressure=float(pressure),
            acceleration_over_HLambda_squared=float(accel),
            physical_cell_volume_ratio=int(a**3)))
    # Reuse the ORIGINAL material populations; a vacuum c-number does not create local athermality.
    energies=np.array(ip['material_energies'])+ip['mass_offset']
    probs=np.array(next(row['actual_populations'] for row in old['rows'] if row['a']==1.))
    entropy=-float(np.dot(probs,np.log(probs)))
    def availability(shift,temp):
        e=energies+shift;emin=float(min(e))
        return float(np.dot(probs,e))-temp*entropy-emin+temp*math.log(float(np.sum(np.exp(-(e-emin)/temp))))
    checks=[]
    for a in (1,2,4):
        T=ip['temp0']/a
        base=availability(0,T);shifted=availability(float(L)*a**3,T)
        inherited=next(row['availability'] for row in old['rows'] if row['a']==a)
        assert abs(base-inherited)<1e-12 and abs(shifted-base)<5e-12
        checks.append(dict(a=a,availability=base,vacuum_shift_residual=abs(shifted-base)))
    # Physical vacuum split is redundant only if the SAME density is compensated.
    delta=F(7,3);vb=F(1,4)*L/V0;vv=F(3,4)*L/V0
    assert (vb-delta)+(vv+delta)==L/V0
    # A material rest-energy shift has a^-3 density, unlike a vacuum-density shift.
    dm=F(2);dL=-dm
    assert dm+dL==0 and dm+dL*F(2)**3==-14
    inputs=[Path(__file__).resolve(),HERE/'drafts/selection.md',HERE/'drafts/STATUS.md',
        STAGE/'992/cooling_availability_results.json',STAGE/'research_note_992.md',
        STAGE/'993/common_candidate_v1.md',STAGE/'957/drafts/unified_operation_hypotheses_v0_2.md',
        STAGE/'996/mechanism_map_v0_8.md',STAGE/'research_note_983.md',
        STAGE.parent/'archive_301_341/research_note_304.md',
        STAGE.parent/'archive_301_341/research_note_319.md',
        STAGE.parent/'archive_531_553/research_note_553.md',
        STAGE.parent/'archive_585_628/research_note_601.md',
        STAGE.parent/'archive_585_628/research_note_613.md']
    return dict(round=997,all_scientific_checks_passed=True,
        kind='mean_source_and_causal_task_adoption_not_cosmological_prediction',
        inputs=dict(material_mass=exact(M),radiation_R=exact(R),V0=exact(V0),
            vacuum_cell_L=exact(L),q=exact(q),r=exact(r),scale_endpoint=int(af)),
        finite_conformal_window=dict(lower=exact(lower),upper=exact(upper),quadrature=value,
            interval_subdivisions=n,sqrt_integer_scale=10**12),
        conditional_future_window=dict(lower=exact(future_lo),upper=exact(future_hi),
            general_upper=1,depends_on_future_continuation=True),
        ideal_causal_tasks=dict(round_trip_distance=exact(close),
            one_way_only_comoving_distance=exact(oneway),
            disjoint_future_cones_distance=exact(separate),
            actual_signal_device_verified=False),
        source_rows=source_rows,availability_same_scale=checks,
        physical_vacuum_value_derived=False,early_boundary_generated=False,
        actual_cosmological_fit=False,full_quantum_geometry_verified=False,
        full_goal_completed=False,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','source_rows')},ensure_ascii=False,indent=2))
