"""1289 independent finite tensor/overlap checks; no author imports.
Default is read-only comparison. --write exclusively creates the result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir() and (p/'scripts').is_dir())
OUT = HERE/'independent_results.json'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def exact_rank(rows):
    pivots = {}
    for row in rows:
        a = list(map(F, row))
        for j, old in sorted(pivots.items()):
            if a[j]:
                c = a[j]
                a = [x-c*y for x,y in zip(a,old)]
        j = next((j for j,x in enumerate(a) if x), None)
        if j is not None:
            c = a[j]
            pivots[j] = [x/c for x in a]
    return len(pivots)


def anti_generators(n):
    # Integer entries, reverse-pair order, and nonorthogonal Cartan basis.
    out = []
    for i in reversed(range(n-1)):
        for j in reversed(range(i+1,n)):
            a = np.zeros((n,n), complex)
            a[i,j], a[j,i] = 1, -1
            out.append(a)
            b = np.zeros((n,n), complex)
            b[i,j] = b[j,i] = 1j
            out.append(b)
    for i in reversed(range(n-1)):
        a = np.zeros((n,n), complex)
        a[i,i], a[n-1,n-1] = 1j, -1j
        out.append(a)
    return out


def integer_rows(columns):
    m = np.column_stack(columns)
    r = np.vstack((m.real,m.imag))
    assert np.max(np.abs(r-np.rint(r))) == 0
    return np.rint(r).astype(int).tolist()


def invariant_dimensions(n):
    xs = anti_generators(n)
    es = []
    for i in range(n):
        for j in range(n):
            e = np.zeros((n,n), complex); e[i,j] = 1
            es.append(e)
    result = {}
    for label, act in [('bilinear', lambda x,b:x.T@b+b@x),
                       ('sesquilinear', lambda x,b:x.conj().T@b+b@x)]:
        columns = [np.concatenate([act(x,b).ravel() for x in xs]) for b in es]
        columns += [1j*c for c in columns]
        rank = exact_rank(integer_rows(columns))
        null_real = 2*n*n-rank
        assert null_real % 2 == 0
        result[label+'_complex_dimension'] = null_real//2
    hs = []
    for h in anti_generators(n-1):
        a = np.zeros((n,n),complex); a[:n-1,:n-1] = h; hs.append(a)
    if hs:
        cols = [np.concatenate([(x@h-h@x).ravel() for h in hs]) for x in xs]
        centralizer = n*n-1-exact_rank(integer_rows(cols))
    else:
        centralizer = n*n-1
    result['stabilizer_adjoint_real_dimension'] = centralizer
    assert result['bilinear_complex_dimension'] == (1 if n == 2 else 0)
    assert result['sesquilinear_complex_dimension'] == 1
    assert centralizer == (3 if n == 2 else 1)
    return result


def hermitian_basis(n):
    # Gram-Schmidt from the above integer basis, rather than Gell-Mann formulas.
    out = []
    for x in anti_generators(n):
        t = -1j*x
        for b in out:
            t = t-2*np.trace(b@t).real*b
        t /= math.sqrt(2*np.trace(t@t).real)
        out.append(t)
    gram = np.array([[2*np.trace(a@b).real for b in out] for a in out])
    assert np.max(np.abs(gram-np.eye(n*n-1))) < 1e-13
    return out


def finite_overlap(n,g,v):
    ts = hermitian_basis(n)
    p = np.zeros(n,complex); p[-1] = v/math.sqrt(2)
    dp = np.array([-1j*g*t@p for t in ts])
    mass = 2*np.real(dp.conj()@dp.T)
    ev, vec = np.linalg.eigh(mass)
    massive = ev > 1e-10
    assert sum(massive) == 2*n-1
    j = np.array([float(np.real(1j*(np.vdot(p,d)-np.vdot(d,p)))) for d in dp])
    currents = [j]
    if n == 2:
        eps = np.array([[0,1],[-1,0]],complex)
        hol = np.array([p@eps@d for d in dp])
        currents += [2*hol.real,2*hol.imag]
    c = np.array(currents)
    assert np.linalg.matrix_rank(c,tol=1e-11) == (3 if n==2 else 1)
    assert np.linalg.matrix_rank(j.reshape(1,-1),tol=1e-11) == 1
    light = g*g*v*v/4
    heavy = g*g*v*v*(n-1)/(2*n)
    want = sorted([0]*((n-1)**2-1)+[light]*(2*n-2)+[heavy])
    assert np.max(np.abs(ev-np.array(want))) < 1e-12
    # Build centralizer projector by an independent SVD in this orthonormal basis.
    hs = hermitian_basis(n-1) if n > 2 else []
    if hs:
        zs = []
        for h in hs:
            hh = np.zeros((n,n),complex); hh[:n-1,:n-1]=h
            z = np.array([(t@hh-hh@t).ravel() for t in ts]).T
            zs.append(np.vstack((z.real,z.imag)))
        z = np.vstack(zs)
        _,s,vh = np.linalg.svd(z,full_matrices=False)
        null = vh[s < 1e-11].T
        central = null@null.T
    else:
        central = np.eye(len(ts))
    row_projector = c.T@np.linalg.pinv(c@c.T)@c
    assert np.max(np.abs(central-row_projector)) < 1e-12
    poles=[]
    for m in sorted(set([round(light,12),round(heavy,12)])):
        pick=np.abs(ev-m)<1e-10
        proj=vec[:,pick]@vec[:,pick].T
        residue=c@proj@c.T
        poles.append(dict(mass_squared=m,bare_multiplicity=int(sum(pick)),
                          residue_rank=int(np.linalg.matrix_rank(residue,tol=1e-11)),
                          residue_trace=float(np.trace(residue))))
    if n > 2:
        assert poles[0]['residue_rank'] == 0 and poles[1]['residue_rank'] == 1
    else:
        assert poles[0]['residue_rank'] == 3
    assert np.linalg.norm(c@vec[:,~massive]) < 1e-12
    return dict(N=n,g=g,v=v,bare_massive_count=int(sum(massive)),
        full_menu_rank=int(np.linalg.matrix_rank(c,tol=1e-11)),neutral_single_menu_rank=1,
        normalized_current_singular_values=np.linalg.svd(c/(g*v*v),compute_uv=False).tolist(),
        poles=poles,mass_spectrum_max_error=float(np.max(np.abs(ev-np.array(want)))),
        row_space_vs_stabilizer_max_error=float(np.max(np.abs(central-row_projector))))


def epsilon_checks():
    eps=np.array([[0,1],[-1,0]],complex)
    # Rational quaternion SU(2), exact after multiplication by its denominator.
    u3=np.array([[1,2+2j],[-2+2j,1]],complex)
    assert np.array_equal(u3.conj().T@u3,9*np.eye(2))
    assert np.array_equal(u3.T@eps@u3,9*eps)
    phase=1j*np.eye(2)
    assert np.array_equal(phase.T@eps@phase,-eps)
    p=np.array([2+1j,-1+3j]); d=np.array([1-2j,4+1j])
    h=p@eps@d
    assert p@eps@p == 0 and h != 0 and d@eps@p == -h
    return dict(local_SU2_invariance_exact=True,global_phase_i_action_on_epsilon=-1,
        holomorphic_global_charge=2,commuting_field_self_contraction=0,
        derivative_contraction=[int(h.real),int(h.imag)],
        derivative_on_opposite_slot_has_negative_sign=True,
        all_N_gt_2_center_obstruction='exp(4*pi*i/N) != 1, equivalently N does not divide 2')


def calculate():
    dims=[dict(N=n,**invariant_dimensions(n)) for n in (2,3,4)]
    cases=[finite_overlap(n,g,v) for n in (2,3,4) for g,v in ((1.,1.),(1.5,4/3))]
    gap=F(1,2)-F(1,4)
    assert gap == F(1,4)
    paths=[Path(__file__),HERE/'proof.md',HERE/'check.py',HERE/'results.json',
           HERE.parent/'_shared/round1289_admission.md']
    return dict(round=1289,status='PASS_WITH_SCOPE',independent_groups=1,
        implementation='Fraction elimination of integer intertwiner equations; independent Gram-Schmidt and spectral projectors; no author imports',
        exact_invariant_dimensions=dims,epsilon_and_derivative=epsilon_checks(),
        independent_mass_and_overlap_cases=cases,
        finite_matrix_conditioning=dict(normalization='currents divided by g*v^2',
            SU2_smallest_singular_value='1/2',operator_error_strict_upper='1/4',
            surviving_singular_value_strict_lower=str(gap),
            scope='matrix perturbation implication only, not a licensed detector error'),
        scope=dict(all_N_argument='center obstruction and residual SU(N-1) representation decomposition',
            finite_jets_smooth_same_background=True,
            expansion_background='A=0, phi=v*e_N/sqrt(2), all ordinary derivative jets zero',
            group_invariant_delta_epsilon_allowed=True,extra_colored_backgrounds_excluded=True,
            transverse_linear_overlap_only=True,tree_poles_only=True,
            global_charged_channels_require_full_compensation_account=True,
            neutral_pair_correlators_not_forbidden=True,
            arbitrary_partner_record_channel_proved=False,
            bare_modes_are_not_licensed_unknown_physical_codes=True,
            complete_detector_stress_tensor_constructed=False,
            nonperturbative_particle_spectrum_selected=False,
            weak_SU2_or_full_gauge_group_selected=False,new_axioms=0),
        source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in paths})


def same(a,b,path='root'):
    if isinstance(b,dict):
        assert a.keys()==b.keys(),path
        for k in b: same(a[k],b[k],path+'.'+k)
    elif isinstance(b,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): same(x,y,path+f'[{i}]')
    elif isinstance(b,float):
        assert math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-12),path
    else: assert a==b,path


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=calculate()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        mode='EXCLUSIVE_WRITE'
    else:
        same(json.loads(OUT.read_text('utf-8')),result); mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status=result['status'],mode=mode,
        exact_invariant_dimensions=result['exact_invariant_dimensions']),ensure_ascii=False))