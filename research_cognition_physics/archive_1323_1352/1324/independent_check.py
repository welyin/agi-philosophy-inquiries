"""Independent round 1324 checks. No author implementation is imported.

Small windows use all occupation states; the large window is propagated by
the minimal polynomial of explicit star Hamiltonians in its cyclic subspace.
Default --check is read-only and requires an existing complete result file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESULT = HERE / 'independent_results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def op(x):
    return float(np.linalg.norm(x, 2))


class Checks:
    def __init__(self):
        self.items = []

    def near(self, name, value, target=0., tol=3e-10):
        residual = float(np.linalg.norm(np.asarray(value) - np.asarray(target)))
        if not np.isfinite(residual) or residual > tol:
            raise AssertionError((name, residual, tol))
        self.items.append(dict(name=name, residual=residual, tolerance=tol))

    def bound(self, name, value, upper, tol=3e-10):
        value, upper = float(value), float(upper)
        if not (np.isfinite(value) and np.isfinite(upper) and value <= upper + tol):
            raise AssertionError((name, value, upper))
        self.items.append(dict(name=name, value=value, upper=upper, tolerance=tol))

    def truth(self, name, condition):
        if not condition:
            raise AssertionError(name)
        self.items.append(dict(name=name, passed=True))


def probabilities(nu, cutoff):
    logs = np.array([k * math.log(nu) - math.lgamma(k + 1)
                     for k in range(cutoff + 1)])
    values = np.exp(logs - logs.max())
    return values / values.sum()


def moments(nu, cutoff):
    p = probabilities(nu, cutoff)
    s = float(np.sqrt(p[:-1] * p[1:]).sum())
    c = float(np.dot(p, 1. / (np.arange(cutoff + 1) + 1)))
    mean = float(np.dot(p, np.arange(cutoff + 1)))
    a = float(np.sqrt(p[:-1] * p[1:] * np.arange(1, cutoff + 1)).sum())
    return dict(p=p, s=s, c=c, mean=mean, a=a)


def unitary(h, t):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j * t * vals)) @ vecs.conj().T


def effective(M, b):
    n = len(M)
    h = np.zeros((3*n, 3*n), complex)
    for alpha in range(2):
        h[(alpha+1)*n:(alpha+2)*n, :n] = b[alpha] * M
    return h + h.conj().T


def full_parent(cutoff, M, normalized):
    """Every single-fermion state in K <= cutoff+1, including dark states."""
    n = len(M)
    states = [(loc, j, k, total-k)
              for loc in range(3) for j in range(n)
              for total in range(cutoff+2-int(loc != 0))
              for k in range(total+1)]
    lookup = {state: row for row, state in enumerate(states)}
    conv = np.zeros((len(states), len(states)), complex)
    for col, (loc, j, k0, k1) in enumerate(states):
        if loc or not k0+k1:
            continue
        total = k0+k1
        divisor = math.sqrt(total) if normalized else 1.
        for alpha, count in enumerate((k0, k1)):
            if count:
                for i in range(n):
                    target = (alpha+1, i, k0-int(alpha == 0), k1-int(alpha == 1))
                    conv[lookup[target], col] = M[i, j] * math.sqrt(count) / divisor
    number = np.array([k0+k1 for _, _, k0, k1 in states], float)
    left = np.array([loc != 0 for loc, _, _, _ in states], float)
    return states, conv+conv.conj().T, number, left


def embed(states, n, p, phi, u):
    E = np.zeros((len(states), 3*n), complex)
    cutoff = len(p)-1
    for row, (loc, j, k0, k1) in enumerate(states):
        k = k0+k1
        if k <= cutoff:
            amplitude = (math.sqrt(p[k] * math.comb(k, k0)) *
                         np.exp(1j * phi * k) * u[0]**k0 * u[1]**k1)
            E[row, loc*n+j] = amplitude
    return E


def state_matrix(states, n, vector):
    """Uncompress the finite support solely for partial tracing reference."""
    refs = sorted(set((a, b) for _, _, a, b in states))
    lookup = {r: j for j, r in enumerate(refs)}
    out = np.zeros((3*n, len(refs)), complex)
    for row, (loc, j, a, b) in enumerate(states):
        out[loc*n+j, lookup[a, b]] = vector[row]
    return out


def pure_distance(a, b):
    return math.sqrt(max(0., 1. - abs(np.vdot(a, b))**2))


def small_checks(c):
    M = np.array([[.7+.1j, .17-.08j], [-.13j, .48-.04j]], complex)
    n, m = len(M), op(M)
    rotated = np.array([math.sqrt(.37)*np.exp(.27j),
                        math.sqrt(.63)*np.exp(-.51j)])
    directions = [np.array([1., 0.], complex), rotated]
    rng = np.random.default_rng(1324)
    rows = []
    for cutoff, nu in [(1, .4), (3, 2.), (5, 4.)]:
        mt = moments(nu, cutoff)
        p, s, cc = mt['p'], mt['s'], mt['c']
        c.near(f'L{cutoff}.finite_mean', mt['mean'], nu*(1-p[-1]))
        c.near(f'L{cutoff}.finite_c', cc, (1-p[0])/nu+p[-1]/(cutoff+1))
        c.near(f'L{cutoff}.finite_a', mt['a'], math.sqrt(nu)*(1-p[-1]))
        for direction, u in enumerate(directions):
            for norm in [True, False]:
                name = f'L{cutoff}.u{direction}.norm{int(norm)}'
                states, V, nb, pl = full_parent(cutoff, M, norm)
                d, phi = len(states), .31
                E = embed(states, n, p, phi, u)
                amp = s if norm else mt['a']
                h = effective(M, amp*np.exp(1j*phi)*u)
                c.truth(name+'.full_dimension', d == n*(cutoff+2)*(3*cutoff+5)//2)
                c.near(name+'.input_isometry', E.conj().T@E, np.eye(3*n))
                c.near(name+'.Hermitian', V, V.conj().T)
                c.near(name+'.K_conservation', ((nb+pl)[:, None]-(nb+pl)[None, :])*V)
                c.near(name+'.parent_norm', op(V), m*(1 if norm else math.sqrt(cutoff+1)))
                c.near(name+'.exact_calibration', E.conj().T@V@E, h)
                R = V@E-E@h
                right = (1-p[0]-s*s) if norm else (mt['mean']-amp*amp)
                parallel = (1-s*s) if norm else (1+mt['mean']-amp*amp)
                transverse = cc if norm else 1.
                G = parallel*np.outer(u, u.conj())+transverse*(np.eye(2)-np.outer(u,u.conj()))
                gram = np.zeros((3*n,3*n), complex)
                gram[:n,:n] = right*M.conj().T@M
                gram[n:,n:] = np.kron(G,M@M.conj().T)
                eta = m*math.sqrt(max(parallel, transverse))
                c.near(name+'.full_input_error_Gram', R.conj().T@R, gram)
                c.near(name+'.exact_error', op(R), eta)
                top = np.array([loc == 0 and a+b == cutoff+1 for loc,j,a,b in states])
                top_response = op((V@E)[top])
                c.truth(name+'.highest_right_is_reachable', top_response > 1e-8)
                zero_count = int(np.sum(abs(np.linalg.eigvalsh(V)) < 1e-9))
                c.truth(name+'.all_dark_dimensions', zero_count == n+n*cutoff*(cutoff+1)//2)
                # All inputs, all partners follow from this complete input-column bound.
                times = [.19/m, .71/m]
                last = None
                for it, t in enumerate(times):
                    actual = unitary(V,t)@E
                    ideal = E@unitary(h,t)
                    delta = op(actual-ideal)
                    c.bound(name+f'.t{it}.joint_Duhamel', delta, t*eta)
                    c.near(name+f'.t{it}.whole_isometry', actual.conj().T@actual,np.eye(3*n))
                    partner = rng.normal(size=(3*n,4))+1j*rng.normal(size=(3*n,4))
                    partner /= np.linalg.norm(partner)
                    dist = pure_distance(actual@partner, ideal@partner)
                    choi = pure_distance(actual/math.sqrt(3*n), ideal/math.sqrt(3*n))
                    c.bound(name+f'.t{it}.arbitrary_partner_sample', dist, min(1.,t*eta))
                    c.bound(name+f'.t{it}.full_Choi_partner', choi, min(1.,t*eta))
                    dn = actual.conj().T@(nb[:,None]*actual)-E.conj().T@(nb[:,None]*E)
                    dl = actual.conj().T@(pl[:,None]*actual)-E.conj().T@(pl[:,None]*E)
                    dl_eff = ideal.conj().T@(pl[:,None]*ideal)-E.conj().T@(pl[:,None]*E)
                    c.near(name+f'.t{it}.all_input_resource_identity',dn+dl)
                    c.bound(name+f'.t{it}.source_prediction_bound',op(dn+dl_eff),min(1.,t*eta))
                    H0 = .75*(1-pl)+1.15*pl+.4*nb
                    c.near(name+f'.t{it}.free_energy',actual.conj().T@(H0[:,None]*actual),
                           E.conj().T@(H0[:,None]*E))
                    c.near(name+f'.t{it}.interaction_energy',actual.conj().T@V@actual,h)
                    c.near(name+f'.t{it}.product_keeps_reference',
                           ideal.conj().T@(nb[:,None]*ideal),mt['mean']*np.eye(3*n))
                    last = dict(time=t,joint_operator_difference=delta,partner_distance=dist,
                                Choi_distance=choi,Duhamel_bound=min(1.,t*eta))
                psi=rng.normal(size=3*n)+1j*rng.normal(size=3*n)
                psi/=np.linalg.norm(psi)
                x=state_matrix(states,n,E@psi)
                vx=state_matrix(states,n,V@E@psi)
                derivative=-1j*(vx@x.conj().T-x@vx.conj().T)
                rho=np.outer(psi,psi.conj())
                c.near(name+'.reduced_initial_generator',derivative,-1j*(h@rho-rho@h))
                rows.append(dict(cutoff=cutoff,nu=nu,direction=direction,normalized=norm,
                                 full_dimension=d,exact_residual=eta,highest_right_response=top_response,
                                 dark_dimensions=zero_count,latest_evolution=last))
    return rows


def reachable(cutoff, g, normalized):
    """Construct closed occupation graph first, then independent finite stars.

    These are ALL states reachable from R/L0/L1 times the axial reference.
    Unreached states in the full invariant parent are not postselected away.
    """
    nodes={(0,0,0)}
    groups=[]
    for K in range(1,cutoff+2):
        pair=[(0,K,0),(1,K-1,0)]
        groups.append((K,pair))
        nodes.update(pair)
        triple=[(0,K-1,1),(2,K-1,0)]
        if K>1:
            triple.append((1,K-2,1))
        groups.append((K,triple))
        nodes.update(triple)
    nodes=sorted(nodes)
    index={x:i for i,x in enumerate(nodes)}
    blocks=[]
    for K,group in groups:
        block=np.zeros((len(group),len(group)),complex)
        r=group[0]
        for j,target in enumerate(group[1:],1):
            alpha=target[0]-1
            count=r[alpha+1]
            value=g*math.sqrt(count)/(math.sqrt(K) if normalized else 1.)
            block[j,0]=value
            block[0,j]=value
        blocks.append((np.array([index[x] for x in group]),block))
    return nodes,index,blocks


def block_apply(blocks, columns, t=None):
    out=np.zeros_like(columns) if t is None else columns.copy()
    for indices,h in blocks:
        x=columns[indices]
        if t is None:
            out[indices]=h@x
        else:
            h2=h@h
            lam=math.sqrt(float(np.trace(h2).real)/2)
            U=np.eye(len(indices))+(math.cos(t*lam)-1)/lam**2*h2-1j*math.sin(t*lam)/lam*h
            out[indices]=U@x
    return out


def large_checks(c):
    L,nu0,nu1,m0=400,64.,256.,.1
    phi=.37
    T=math.pi/(2*m0)
    cal=moments(nu0,L)
    gs={True:m0/cal['s'],False:m0/cal['a']}
    output=[]
    for normalized in [True,False]:
        g=gs[normalized]
        nodes,index,blocks=reachable(L,g,normalized)
        d=len(nodes)
        c.truth(f'large.norm{int(normalized)}.reachable_dimension',d==5*L+5)
        # A direct graph closure audit using number occupation transitions.
        closure_tests=[]
        for loc,k0,k1 in nodes:
            if loc==0:
                for alpha,count in enumerate((k0,k1)):
                    if count:
                        closure_tests.append((alpha+1,k0-int(alpha==0),k1-int(alpha==1)) in index)
            else:
                closure_tests.append((0,k0+int(loc==1),k1+int(loc==2)) in index)
        c.truth(f'large.norm{int(normalized)}.all_occupation_transitions_closed',all(closure_tests))
        nb=np.array([a+b for loc,a,b in nodes],float)
        pl=np.array([loc>0 for loc,a,b in nodes],float)
        H0=.75*(1-pl)+1.15*pl+.4*nb
        for nu in [nu0,nu1]:
            name=f'large.norm{int(normalized)}.nu{nu:g}'
            mt=moments(nu,L)
            p,s,cc=mt['p'],mt['s'],mt['c']
            coeff=np.sqrt(p)*np.exp(1j*phi*np.arange(L+1))
            E=np.zeros((d,3),complex)
            for loc in range(3):
                for k in range(L+1):
                    E[index[loc,k,0],loc]=coeff[k]
            amp=s if normalized else mt['a']
            h=effective(np.array([[g]]),np.array([amp*np.exp(1j*phi),0]))
            VE=block_apply(blocks,E)
            residual=VE-E@h
            eta=g*math.sqrt(max(1-s*s,cc) if normalized else 1+nu*p[-1]*(1-p[-1]))
            c.near(name+'.finite_normalization',p.sum(),1.)
            c.near(name+'.finite_mean',mt['mean'],nu*(1-p[-1]))
            c.near(name+'.finite_c',cc,(1-p[0])/nu+p[-1]/(L+1))
            c.near(name+'.finite_a',mt['a'],math.sqrt(nu)*(1-p[-1]))
            c.near(name+'.input_embedding',E.conj().T@E,np.eye(3))
            c.near(name+'.exact_mean_calibration',E.conj().T@VE,h)
            c.near(name+'.full_input_exact_residual',op(residual),eta)
            if nu==nu0:
                c.near(name+'.fixed_short_time_calibration',abs(h[1,0]),m0)
            actual=block_apply(blocks,E,T)
            ideal=E@unitary(h,T)
            duhamel=min(1.,T*eta)
            c.near(name+'.actual_isometry',actual.conj().T@actual,np.eye(3))
            joint_difference=op(actual-ideal)
            c.bound(name+'.whole_joint_Duhamel',joint_difference,T*eta)
            c.bound(name+'.full_Choi_partner',pure_distance(actual/math.sqrt(3),ideal/math.sqrt(3)),duhamel)
            actual_left=actual.conj().T@(pl[:,None]*actual)
            ideal_left=ideal.conj().T@(pl[:,None]*ideal)
            p_actual=float(actual_left[0,0].real)
            p_effective=float(ideal_left[0,0].real)
            weighted=(1-p[0])*math.sin(g*T)**2 if normalized else float(np.dot(p,np.sin(g*T*np.sqrt(np.arange(L+1)))**2))
            c.near(name+'.right_probability_finite_sum',p_actual,weighted)
            c.near(name+'.effective_probability',p_effective,math.sin(g*amp*T)**2)
            c.bound(name+'.observable_probability_error',abs(p_actual-p_effective),duhamel)
            dn=actual.conj().T@(nb[:,None]*actual)-E.conj().T@(nb[:,None]*E)
            dl=actual_left-E.conj().T@(pl[:,None]*E)
            dl_eff=ideal_left-E.conj().T@(pl[:,None]*E)
            c.near(name+'.all_input_reference_balance',dn+dl)
            c.bound(name+'.backreaction_from_effective_bounded_readout',op(dn+dl_eff),duhamel)
            c.near(name+'.free_energy',actual.conj().T@(H0[:,None]*actual),E.conj().T@(H0[:,None]*E),tol=2e-9)
            Vactual=block_apply(blocks,actual)
            c.near(name+'.interaction_energy',actual.conj().T@Vactual,h)
            c.near(name+'.total_energy',actual.conj().T@(H0[:,None]*actual+Vactual),
                   E.conj().T@(H0[:,None]*E+VE),tol=2e-9)
            c.near(name+'.frozen_product_reference',ideal.conj().T@(nb[:,None]*ideal),mt['mean']*np.eye(3),tol=2e-9)
            c.near(name+'.actual_reference_depletion',dn[0,0],-p_actual)
            transverse_up=float(1-actual_left[2,2].real)
            transverse_formula=float(np.dot(p,np.sin(g*T*(1 if normalized else np.sqrt(np.arange(L+1)+1)))**2/(np.arange(L+1)+1)))
            c.near(name+'.transverse_reverse_probability',transverse_up,transverse_formula)
            c.near(name+'.transverse_reference_gain',dn[2,2],transverse_up)
            # Actual unknown flavor and partner columns, not just integer labels.
            flavor_rows=[]
            for n in [3,4,5]:
                AJ=np.kron(actual,np.eye(n)); IJ=np.kron(ideal,np.eye(n))
                c.near(name+f'.n{n}.all_unknown_flavor_isometry',AJ.conj().T@AJ,np.eye(3*n))
                partner_distance=pure_distance(AJ/math.sqrt(3*n),IJ/math.sqrt(3*n))
                c.bound(name+f'.n{n}.full_flavor_Choi',partner_distance,duhamel)
                flavor_rows.append(dict(n=n,complete_input_dimension=3*n,partner_distance=partner_distance))
            output.append(dict(normalized=normalized,nu=nu,cutoff=L,g=g,
                               mean_reference=mt['mean'],p0=float(p[0]),p_top=float(p[-1]),
                               s=s,c=cc,a=mt['a'],time=T,effective_strength=g*amp,
                               residual=eta,Duhamel_bound=duhamel,joint_operator_difference=joint_difference,
                               actual_probability=p_actual,effective_probability=p_effective,
                               actual_reference_change=float(dn[0,0].real),
                               frozen_product_reference_change=0.,
                               transverse_reverse_probability=transverse_up,
                               flavor_partner_checks=flavor_rows))
    normal=next(r for r in output if r['normalized'] and r['nu']==nu1)
    bare=next(r for r in output if not r['normalized'] and r['nu']==nu1)
    numerical_allowance=1e-9
    low=max(0.,normal['effective_probability']-normal['Duhamel_bound']-numerical_allowance)
    high=min(1.,bare['effective_probability']+bare['Duhamel_bound']+numerical_allowance)
    certified_gap=low-high
    c.truth('heldout.certified_probability_gap',certified_gap>.70)
    c.truth('heldout.actual_probability_gap',normal['actual_probability']-bare['actual_probability']>.98)
    # The SAME >.70 conclusion also has a coarse rational certificate.
    # log(x) >= 2(x-1)/(x+1) gives both Poisson Chernoff exponents >30.
    # e > 8/3, and (8/3)^30 > 2*10^12, imply each conditional top weight <10^-12.
    tail=Fraction(1,10**12)
    c.truth('heldout.rational_tail_certificate',
            min(Fraction(7056,29),Fraction(1296,41))>30 and Fraction(8,3)**30>2*10**12)
    lo0=Fraction(127,128)-tail
    lo1=Fraction(511,512)-tail
    c.truth('heldout.rational_Jensen_bounds',
            Fraction(127,128)**2 < Fraction(64,65) and
            Fraction(511,512)**2 < Fraction(256,257))
    pi_hi=Fraction(22,7)
    deviation=pi_hi/2*max(1-lo1,1/lo0-1)
    normal_error=pi_hi/(2*lo0)*(Fraction(1,16)+16*tail)
    bare_error=pi_hi/(16*(1-tail))*(1+128*tail)
    bare_probability=(pi_hi*tail/(1-tail))**2
    rational_gap=1-deviation**2-normal_error-bare_error-bare_probability
    c.truth('heldout.rational_separation',rational_gap>Fraction(7,10))
    # Exact phase twirl kills the first moment, not the finite right-transfer task.
    phases=2*math.pi*np.arange(L+1)/(L+1)
    c.near('phase_twirl.mean_first_moment',np.exp(1j*phases).mean())
    c.truth('phase_twirl.nonzero_actual_task',normal['actual_probability']>.99)
    failure_bound=4*math.exp(-2*1000*.05**2)
    c.truth('registered_repetitions.failure_bound',failure_bound<.027)
    return dict(calibration_nu=nu0,heldout_nu=nu1,cutoff=L,calibration_strength=m0,
                time=T,reachable_occupation_dimension=5*L+5,
                data=output,certified_probability_gap=certified_gap,
                numerical_interval_allowance=numerical_allowance,
                rational_probability_gap_lower_bound=dict(numerator=str(rational_gap.numerator),
                                                         denominator=str(rational_gap.denominator),
                                                         decimal=float(rational_gap)),
                actual_probability_gap=normal['actual_probability']-bare['actual_probability'],
                conditional_repetition_failure_bound=failure_bound)


def compute():
    c=Checks()
    small=small_checks(c)
    large=large_checks(c)
    finite_residuals=[x['residual'] for x in c.items if 'residual' in x]
    return dict(schema='round1324_independent_v1',status='PASS_INDEPENDENT_CHECKS_ONLY',
                scientific_acceptance=False,author_code_imported=False,
                checks=len(c.items),max_absolute_residual=max(finite_residuals),
                small_full_occupation_cases=small,large_finite_reference=large,
                evidence=c.items,
                scope=dict(small_general_complex_flavors=2,
                           small_rotated_reference=True,
                           large_full_input_cyclic_subspace=True,
                           large_unreachable_complement_preserved_in_parent_definition=True,
                           large_flavors=[3,4,5],
                           initial_reference_independent_of_material_and_partner=True,
                           arbitrary_general_M_by_analytic_theorem=True,
                           physical_Higgs_source_not_certified=True,
                           new_axioms=0),
                source_bindings={str(p.relative_to(STAGE)).replace('\\','/'):sha(p)
                                 for p in [HERE/'independent_check.py',HERE/'proof.md',
                                           STAGE/'_admission'/'round1324_admission.md',
                                           STAGE/'1323'/'proof.md']})


def compare(actual,saved,path='result'):
    if isinstance(actual,dict):
        if not isinstance(saved,dict) or actual.keys()!=saved.keys():
            raise AssertionError(path+': keys/type')
        for key in actual:
            compare(actual[key],saved[key],path+'.'+key)
    elif isinstance(actual,list):
        if not isinstance(saved,list) or len(actual)!=len(saved):
            raise AssertionError(path+': sequence/type')
        for i,(a,b) in enumerate(zip(actual,saved)):
            compare(a,b,path+f'[{i}]')
    elif isinstance(actual,float):
        if isinstance(saved,bool) or not isinstance(saved,(float,int)) or not np.isfinite(saved):
            raise AssertionError(path+': numeric/type')
        if not math.isclose(actual,saved,rel_tol=5e-10,abs_tol=5e-12):
            raise AssertionError((path,actual,saved))
    elif type(actual) is not type(saved) or actual!=saved:
        raise AssertionError((path,actual,saved))


def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--check',action='store_true')
    mode.add_argument('--write-new',action='store_true')
    mode.add_argument('--preview',action='store_true',help='Compute without writing unaccepted evidence.')
    args=parser.parse_args()
    if args.write_new and RESULT.exists():
        raise FileExistsError(RESULT)
    if not args.write_new and not args.preview and not RESULT.is_file():
        raise FileNotFoundError(RESULT)
    result=compute()
    if args.write_new:
        with RESULT.open('x',encoding='utf-8',newline='\n') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2,allow_nan=False)
            handle.write('\n')
    elif not args.preview:
        compare(result,json.loads(RESULT.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',mode='preview' if args.preview else 'write-new' if args.write_new else 'check',
                          checks=result['checks'],max_residual=result['max_absolute_residual'],
                          certified_gap=result['large_finite_reference']['certified_probability_gap']),indent=2))


if __name__=='__main__':
    main()
