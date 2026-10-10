"""1060: check the global representation map, not a numeric bordism proof."""
from pathlib import Path
import argparse
import ast
import itertools
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OLD = ROOT/'archive_585_628/614/joint_spinor_subgroup_mass.py'
STATES = [tuple(c) for n in (0, 2, 4) for c in itertools.combinations(range(5), n)]
EPS = np.array([[0, 1], [-1, 0]], complex)
# Reuse only the old pure matrix functions, without its relocated physics imports.
tree = ast.parse(OLD.read_text(encoding='utf-8-sig'))
pure = [n for n in tree.body if isinstance(n, ast.FunctionDef)
        and n.name in {'op', 'exterior', 'dictionary'}]
assert len(pure) == 3
exec(compile(ast.Module(body=pure, type_ignores=[]), str(OLD), 'exec'))


def su(rng, n):
    z = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q = np.linalg.qr(z)[0]
    q[:,0] /= np.linalg.det(q)
    return q


def vector(k, C, W, theta):
    u = np.zeros((5,5), complex)
    u[:3,:3] = np.exp(-1j*(2+12*k)*theta)*C
    u[3:,3:] = np.exp(3j*theta)*W
    return u


def spinor(k, C, W, theta):
    return np.exp(18j*k*theta)*exterior(vector(k,C,W,theta))


def old_species_rep(k, C, W, theta):
    pieces = [np.exp(1j*(1+6*k)*theta)*np.kron(C,W),
              np.exp(1j*(-4-6*k)*theta)*C.conj(),
              np.exp(1j*(2-6*k)*theta)*C.conj(),
              np.exp(1j*(-3-18*k)*theta)*W,
              np.array([[np.exp(1j*(6+18*k)*theta)]]),
              np.array([[np.exp(18j*k*theta)]])]
    out = np.zeros((16,16), complex)
    start = 0
    for p in pieces:
        end = start+len(p)
        out[start:end,start:end] = p
        start = end
    return out


def run():
    J = dictionary()
    n_c = np.array([sum(i<3 for i in s) for s in STATES])
    n_w = np.array([sum(i>=3 for i in s) for s in STATES])
    q0 = -2*n_c+3*n_w
    q1 = 18-12*n_c
    x0 = 5-2*(n_c+n_w)
    b = 3-2*n_c
    assert np.array_equal(q1, 6*b)
    assert np.array_equal(3*x0+2*q0, 5*b)
    expected0 = np.array([1]*6+[-4]*3+[2]*3+[-3]*2+[6,0])
    expected1 = np.array([6]*6+[-6]*3+[-6]*3+[-18]*2+[18,18])
    assert op(J.conj().T@np.diag(q0)@J-np.diag(expected0)) == 0
    assert op(J.conj().T@np.diag(q1)@J-np.diag(expected1)) == 0
    assert np.all(x0 % 4 == 1)
    c = np.diag(np.exp(.5j*np.pi*x0))
    assert op(c-1j*np.eye(16)) < 1e-14
    assert op(c@c+np.eye(16)) < 1e-14
    # All-k cocharacter and central descent are integer-polynomial identities.
    assert 3*(-2)+2*3 == 0 and 3*(-12) == -36
    central0 = 2*np.array([1]*6+[-1]*6+[0]*4)
    weak = np.array([1]*6+[0]*6+[1]*2+[0]*2)
    assert np.all((central0+3*weak+expected0) % 6 == 0)
    assert np.all(expected1 % 6 == 0)

    rng = np.random.default_rng(1060)
    errs = dict(intertwiner=0., group_law=0., periodicity=0.,
                center_descent=0., higgs=0., unitarity=0.)
    ks = [-11,-2,-1,0,1,2,13]
    records = []
    for k in ks:
        for j in range(6):
            theta = j*np.pi/3
            C = np.exp(2j*theta)*np.eye(3)
            W = np.exp(-3j*theta)*np.eye(2)
            errs['center_descent'] = max(errs['center_descent'],
                                         op(spinor(k,C,W,theta)-np.eye(16)))
        for _ in range(6):
            C, D = su(rng,3), su(rng,3)
            W, V = su(rng,2), su(rng,2)
            theta, phi = rng.uniform(-1.,1.,size=2)
            R = spinor(k,C,W,theta)
            errs['intertwiner'] = max(errs['intertwiner'],
                op(R@J-J@old_species_rep(k,C,W,theta)))
            errs['group_law'] = max(errs['group_law'],
                op(R@spinor(k,D,V,phi)-spinor(k,C@D,W@V,theta+phi)))
            errs['periodicity'] = max(errs['periodicity'],
                op(R-spinor(k,C,W,theta+2*np.pi)))
            errs['unitarity'] = max(errs['unitarity'],op(R.conj().T@R-np.eye(16)))
            H = rng.normal(size=2)+1j*rng.normal(size=2)
            errs['higgs'] = max(errs['higgs'],float(np.linalg.norm(
                vector(k,C,W,theta)@np.r_[np.zeros(3),H]
                -np.r_[np.zeros(3),np.exp(3j*theta)*W@H])))
        records.append({'k':k, 'charges_Q_uc_dc_L_ec_N':
                        [1+6*k,-4-6*k,2-6*k,-3-18*k,6+18*k,18*k]})
    assert max(errs.values()) < 3e-11, errs
    return {'round':1060, 'passed':True,
            'all_k_exact_weight_coefficients_verified':True,
            'all_k_spin_cocharacter_even_sum_verified':True,
            'all_k_Z6_character_identity_verified':True,
            'same_round614_intertwiner_used':True,
            'nonAbelian_group_samples':42, 'central_samples':42,
            'diagnostic_max_operator_errors':errs, 'charge_examples':records,
            'parent_anomaly_theorem':'Wang-Wen-Witten 1810.00844v4 section 5.1.3',
            'all_background_anomaly_claim_based_on':'analytic structure-group pullback and adopted parent theorem',
            'bordism_computed_numerically':False,
            'physical_Spin10_gauge_group_added':False,
            'new_cognitive_axioms':0, 'new_empirical_groups':0,
            'whole_roadmap_completed':False}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=1e-12,abs_tol=2e-12)
    return a==b


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--save',action='store_true')
    result=run()
    if parser.parse_args().save:
        with (HERE/'results.json').open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert same(result,json.loads((HERE/'results.json').read_text(encoding='utf-8')))
    print(json.dumps(result,ensure_ascii=False))
