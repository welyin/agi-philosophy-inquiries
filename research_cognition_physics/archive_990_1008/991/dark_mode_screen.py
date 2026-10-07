"""Bounded mechanism screen; no relic-abundance or cosmological simulation."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
TARGET = HERE / 'dark_mode_screen_results.json'

def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items(): out[k] = out.get(k, F(0)) + v
    return {k: v for k, v in out.items() if v}

def scale(p, c): return {k: v*c for k, v in p.items() if v*c}

def mul(a, b):
    out = {}
    for (i, j), c in a.items():
        for (k, l), d in b.items():
            key = i+k, j+l
            out[key] = out.get(key, F(0)) + c*d
    return {k: v for k, v in out.items() if v}

def derivative(p, axis):
    out = {}
    for k, v in p.items():
        if k[axis]:
            target = list(k); target[axis] -= 1
            out[tuple(target)] = v*k[axis]
    return out

def value(p, h, s): return sum(c*h**i*s**j for (i,j),c in p.items())

def pack(x): return dict(exact=str(x), value=float(x))

def run():
    v, lh, ld, portal = F(1), F(1,2), F(1,100), F(1,1000)
    bare = F(499,2000)
    mh2, ms2 = 2*lh*v*v, bare+portal*v*v/2
    assert mh2 > 0 and ms2 == F(1,4)
    omega, volume = F(1,2), F(1000)
    # Independent construction from the gauge-invariant x = H^dagger H.
    x = {(0,0):v*v/2, (1,0):v, (2,0):F(1,2)}
    y = add(x, {(0,0):-v*v/2})
    vh = scale(mul(y,y), lh)
    vd = {(0,2):bare/2, (0,4):ld/4}
    vp = scale(mul(x, {(0,2):F(1)}), portal/2)
    total = add(vh,vd,vp)
    expected = {(2,0):lh*v*v, (3,0):lh*v, (4,0):lh/4,
                (0,2):ms2/2, (0,4):ld/4, (1,2):portal*v/2,
                (2,2):portal/4}
    assert total == expected
    assert all(j % 2 == 0 for i,j in total)
    assert value(derivative(derivative(total,0),1), F(0),F(0)) == 0
    assert value(derivative(derivative(total,0),0), F(0),F(0)) == mh2
    assert value(derivative(derivative(total,1),1), F(0),F(0)) == ms2
    # One exact source/exchange identity for nonzero values of both fields.
    h, s, dh, ds = F(1,7),F(2,5),F(-1,3),F(3,8)
    jh = value(derivative(vp,0),h,s)
    js = value(derivative(vp,1),h,s)
    chain_direct = portal/2*(v+h)*s*s*dh + portal/2*(v+h)**2*s*ds
    exchange = jh*dh + js*ds - chain_direct
    assert exchange == 0
    # This finite matrix computes degree <= 4 moments on n=0,1,4 exactly;
    # it is not a dynamical Fock-space cutoff or a full interacting vacuum.
    dim = 8
    a = np.diag(np.sqrt(np.arange(1,dim,dtype=float)), 1)
    q = (a+a.T)/math.sqrt(float(2*omega*volume))
    q2, q4 = q@q, q@q@q@q
    parity = np.diag([(-1)**n for n in range(dim)])
    H = float(omega)*(np.diag(np.arange(dim))+.5*np.eye(dim)) + float(volume*ld/4)*q4
    J = float(portal*v/2)*q2
    base_energy = H[0,0]
    base_j = J[0,0]
    rows = []
    for n in (0,1,4):
        sn2 = F(2*n+1)/(2*omega*volume)
        sn4 = F(3*(2*n*n+2*n+1))/(4*omega**2*volume**2)
        delta_rho = omega*n/volume + F(3)*ld*n*(n+1)/(8*omega**2*volume**2)
        delta_j = portal*v*n/(2*omega*volume)
        residual = max(abs(q[n,n]),abs(q2[n,n]-float(sn2)),
            abs(q4[n,n]-float(sn4)),
            abs((H[n,n]-base_energy)/float(volume)-float(delta_rho)),
            abs(J[n,n]-base_j-float(delta_j)))
        assert residual < 1e-13
        assert n+2 < dim
        if n:
            assert delta_rho > 0 and 0 < delta_j/delta_rho <= portal*v/(2*ms2)
        rows.append(dict(n=n,mean_s=float(q[n,n]),s2=pack(sn2),s4=pack(sn4),
            delta_energy_density=pack(delta_rho),delta_higgs_source=pack(delta_j),
            source_to_energy_ratio=pack(delta_j/delta_rho) if n else None,
            max_moment_residual=residual))
    assert np.linalg.norm(parity@H-H@parity) == 0
    assert np.linalg.norm(parity@J-J@parity) == 0
    sources = [STAGE/'research_note_990.md', STAGE/'research_note_930.md',
        STAGE/'research_note_946.md', STAGE/'research_note_971.md',
        STAGE/'981/drafts/common_parent_contract_v1.md',
        HERE/'drafts/selection.md', HERE/'drafts/STATUS.md']
    sources += [STAGE.parent/'archive_301_341'/f'research_note_{n}.md' for n in (325,326,332,333)]
    sources += [STAGE.parent/'archive_342_369'/f'research_note_{n}.md' for n in (351,358,366)]
    return dict(round=991,all_scientific_checks_passed=True,
        research_kind='conditional_dark_sector_mechanism_screen',
        parameters={k:pack(val) for k,val in dict(v=v,lambda_h=lh,lambda_dark=ld,
            portal=portal,bare_mass_squared=bare,volume=volume).items()},
        masses_squared=dict(higgs=pack(mh2),dark=pack(ms2)),
        polynomial={f'h^{i} s^{j}':pack(c) for (i,j),c in sorted(total.items())},
        mixed_quadratic_coefficient=0,source_exchange_residual=pack(exchange),
        positivity_certificate='sum of nonnegative terms in x>=0 and s^2>=0',
        parity_commutator_norm=0.0,rows=rows,
        analytic_source_energy_ratio_upper=pack(portal*v/(2*ms2)),
        geometry_source_is_full_energy_not_higgs_source=True,
        calibration_state='free oscillator n=0; not interacting ground or universal vacuum subtraction',
        scalar_is_new_explicit_physical_input=True,
        scalar_model_is_original=False,source_moments_are_not_classical_geometry_certificate=True,
        relic_density_computed=False,current_experimental_viability_claimed=False,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float): assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-13),(a,b)
    else: assert a==b,(a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else: compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:out[k] for k in ('round','all_scientific_checks_passed','masses_squared',
        'analytic_source_energy_ratio_upper','rows','full_goal_completed')},ensure_ascii=False,indent=2))
