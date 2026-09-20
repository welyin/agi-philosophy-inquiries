"""Optional multistart discovery for the coherent relation correction.

Numerical maxima and slopes are diagnostics only. Global certification is
performed separately from these searches, on the saved rational candidate.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import numpy as np

from four_angle_global_certificate import witness_data
from relational_history_decoder import PARAMETER, error_jet
from joint_history_certificate import stereographic_bloch
from soft_decoder_certificate import DENOMINATOR


def discover(parameter=PARAMETER, count=20):
    from scipy.optimize import minimize
    _, states, _ = witness_data()
    starts = [np.array(states,dtype=float).flatten()]
    for seed in range(count):
        point=np.random.default_rng(seed).normal(size=(2,3))
        point/=np.linalg.norm(point,axis=1)[:,None];starts.append(point.flatten())
    constraints=[{'type':'ineq','fun':lambda x:1-x[:3] @ x[:3],
                  'jac':lambda x:np.r_[-2*x[:3],np.zeros(3)]},
                 {'type':'ineq','fun':lambda x:1-x[3:] @ x[3:],
                  'jac':lambda x:np.r_[np.zeros(3),-2*x[3:]]}]
    def objective(point):
        jet=error_jet(point,parameter)
        return -jet['value'],-jet['gradient']
    clusters=[]
    for start in starts:
        result=minimize(objective,start,jac=True,method='SLSQP',constraints=constraints,
                        options={'ftol':1e-13,'maxiter':160})
        point=result.x.copy()
        point[:3]/=max(1.,np.linalg.norm(point[:3]));point[3:]/=max(1.,np.linalg.norm(point[3:]))
        if not any(np.linalg.norm(point-np.array(item['point']))<.01 for item in clusters):
            clusters.append({'error_diagnostic':error_jet(point,parameter)['value'],
                             'point':point.tolist(),'optimizer_message':str(result.message)})
    return sorted(clusters,key=lambda item:item['error_diagnostic'],reverse=True)


def save_witness(clusters):
    point=np.array(clusters[0]['point'])
    if point[1]<0: point[[1,4]]*=-1
    numbers=[int(round(x*10**6)) for x in np.r_[point[:2]/(1+point[2]),point[3:5]/(1+point[5])]]
    a=stereographic_bloch(*(F(n,10**6) for n in numbers[:2]))
    b=stereographic_bloch(*(F(n,10**6) for n in numbers[2:]))
    jet=error_jet(np.array(a+b,dtype=float))
    data={'parameter_exact':str(PARAMETER),'stereographic_numerators':numbers,
          'stereographic_denominator':10**6,'coefficient_denominator':DENOMINATOR,
          'coefficient_numerators':[int(n) for n in np.rint(jet['coefficients']*DENOMINATOR)],
          'error_diagnostic':jet['value'],'candidate_is_not_itself_a_proof':True}
    Path(__file__).with_name('relational_history_witness.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


def main():
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'.research_runtime'))
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-candidate',action='store_true')
    parser.add_argument('--scan',action='store_true');args=parser.parse_args()
    candidates={str(PARAMETER):discover()}
    if args.scan:
        for t in (F(0),F(1,2000),F(3,2000),F(1,400),F(1,200),F(1,100)):
            candidates[str(t)]=discover(t,10)
    report={'search_is_not_a_global_proof':True,'clusters_by_half_angle_tangent':candidates}
    if args.write_candidate:
        save_witness(candidates[str(PARAMETER)])
        Path(__file__).with_name('relational_history_search_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
