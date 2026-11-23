"""Discrete QCD measure, full determinant, and the coefficient price of leading powers."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def hierarchy(y,p,epsilon):
    c=np.array(y)/epsilon**np.array(p)
    return float(max(c)/min(c)),c


def optimize(y,p):
    # log(max|c|/min|c|) is piecewise affine in log(epsilon).
    # At infinity it diverges for any unequal exponent assignment.
    ys=np.log(y)
    knots=[0.]
    for i,j in itertools.combinations(range(3),2):
        if p[i]!=p[j]:
            crossing=(ys[i]-ys[j])/(p[i]-p[j])
            if crossing<0:
                knots.append(crossing)
    return min(hierarchy(y,p,math.exp(k))[0] for k in knots)


def solve():
    phi=-np.pi/6
    good=[]
    for powers in itertools.product(range(4),repeat=3):
        n=4+sum(powers)
        # The CP-real sign can be chosen once to put a real determinant on the positive axis.
        det_phase=np.exp(1j*n*phi)
        instanton=np.exp(1j*n*np.pi/2)
        if abs(det_phase.imag)<1e-12 and abs(instanton-1)<1e-12:
            good.append(list(powers))
    assert len(good)==3 and all(sorted(p)==[2,3,3] for p in good)
    inp=np.array(read('../1397/inputs.json')['high_scale']['up'])
    y=inp[:,0]
    by_perm={','.join(str(k) for k in p):round(optimize(y,p),6) for p in good}
    assert abs(min(by_perm.values())-y[1]/y[0])<6e-7
    low_y,high_y=y-2*inp[:,1],y+2*inp[:,1]
    min_two=low_y[1]/high_y[0]
    assert min_two>461
    lifted=np.array([5,3,0])
    eps=.095
    _,couplings=hierarchy(y,lifted,eps)
    assert np.all(low_y/eps**lifted>1/3) and np.all(high_y/eps**lifted<3)
    assert sum(lifted)+4==12 and sum(lifted//4)==1
    left=np.array([0,0,1])
    right=(lifted-left)%4
    assert np.all((lifted-left-right)%4==0)
    # Test a full complex texture in the CP-real convention. R and D
    # supply a nonaligned weak phase but drop out of the determinant.
    def matrix_pair(powers,angle):
        phase=lambda n: np.exp(1j*n*angle)
        down=np.array([[0,.011*phase(0),0],[-.019*phase(2),.023*phase(0),
                       -.17*phase(1)+.113*phase(5)],[0,.27*phase(1),.72*phase(2)]])
        up=np.diag(y*np.exp(1j*np.array(powers)*angle))
        return up,down
    for p in [*good,lifted,[0,0,2]]:
        N=4+sum(p)
        for angle in [phi,.132,phi+np.pi/2]:
            up,down=matrix_pair(p,angle)
            determinant=np.linalg.det(up)*np.linalg.det(down)
            np.testing.assert_allclose(determinant/abs(determinant),np.exp(1j*N*angle),atol=1e-13)
    frozen=np.array([0,0,2])
    assert (sum(frozen)+4)%4==2
    heavy=[]
    for count in (0,1,2):
        heavy=[p for p in itertools.combinations_with_replacement(range(4),count)
               if (6+sum(p))%12==0]
        if heavy:
            break
    assert heavy==[(3,3)]
    budget=read('../1637/results.json')['conditional_angle_budget_rad']
    slope=complex(np.exp(-4j*phi))
    tan=math.tan(budget)
    allowed=[-tan/(abs(slope.imag)+tan*slope.real),tan/(abs(slope.imag)-tan*slope.real)]
    for r in allowed:
        assert abs(abs(np.angle(1+r*slope))-budget)<1e-20
    return {'leading_allowed_up_powers_u_c_t':good,'optimized_leading_coefficient_ranges_by_permutation':by_perm,
            'minimum_coefficient_hierarchy_central':round(float(min(by_perm.values())),6),
            'minimum_possible_hierarchy_two_width_box':round(float(min_two),6),
            'optimal_epsilon_interval_central':[round(float(y[0]/y[2]),10),round(float(y[1]/y[2]),10)],
            'witness_epsilon':eps,'nonleading_up_powers':lifted.tolist(),'right_charge_residues':right.tolist(),
            'absolute_coefficients_u_c_t':couplings.round(9).tolist(),
            'two_width_coefficient_ranges':np.stack((low_y/eps**lifted,high_y/eps**lifted),axis=1).round(9).tolist(),
            'conditional_angle_budget_rad':budget,
            'allowed_c1_over_c5_interval_in_lifted_example':[float(f'{v*eps**4:.10g}') for v in allowed],
            'minimal_added_single_fundamental_Dirac_exponents_for_SU5_charges':[3,3]}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1681 independent PASS: instanton measure, determinant, optimal hierarchy, one lift and colored thresholds')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as output:
            output.write(json.dumps(result,indent=2)+'\n')
