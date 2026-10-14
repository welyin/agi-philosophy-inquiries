"""Post-hoc stored-matrix guard; no assumed accuracy of eigh.

Replays precisely the RK4 state recurrence (without recomputing its residual
certificate), then requires bitwise agreement with the saved endpoint report.
--write is exclusive creation. Default read-only reproduction.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
UR=2.**-53
def up(x):return np.nextafter(x,np.inf)
def gamma(n):return float(up(n*UR/(1-n*UR)))
def frac_up(x):return float(up(float(x)))
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def float_leaves(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        return sum(float_leaves(a[k],b[k]) for k in a)
    assert a==b,(a,b)
    return int(isinstance(a,float))
def exact_matrix_l1(a,b):
    return sum(abs(F(float(x.real))-F(float(y.real)))+
               abs(F(float(x.imag))-F(float(y.imag))) for x,y in zip(a.ravel(),b.ravel()))
def run():
    m=load('author_endpoint_replay',HERE/'preview.py')
    cert=json.loads((HERE/'numerical_bound_checks.json').read_text(encoding='utf-8'))
    raw=m.evolve(m.COUPLING,cert['steps'])['u']
    gram=max(float(np.linalg.norm(a.conj().T@a-np.eye(4),'fro')) for a in raw.reshape(-1,4,4))
    gram_upper=float(up(gram+gamma(256)*16*float(np.max(abs(raw)))**2))
    shift=float(up(gram_upper/(1+np.sqrt(1-gram_upper))))
    assert gram_upper==cert['probe_gram_upper']
    assert shift==cert['right_polar_operator_correction']
    stored=np.empty_like(raw)
    for ii in np.ndindex(raw.shape[:-2]):
        a=raw[ii];ev,vec=np.linalg.eigh(a.conj().T@a)
        stored[ii]=a@(vec@np.diag(1/np.sqrt(ev))@vec.conj().T)
    menu=m.readouts(stored)
    matched=float_leaves(menu,cert['normalized_original_menu'])
    exact_l1=max(exact_matrix_l1(a,b) for a,b in zip(stored.reshape(-1,4,4),raw.reshape(-1,4,4)))
    matrix_l1=frac_up(exact_l1)
    frobenius=max(float(np.linalg.norm(a-b,'fro')) for a,b in zip(stored.reshape(-1,4,4),raw.reshape(-1,4,4)))
    eps_matrix=frac_up(F(matrix_l1)+F(shift))
    polar_probability=frac_up((2+F(eps_matrix))*F(eps_matrix))
    # Exact real witness target: 16 entries 1/4, or 8 entries 1/(2 sqrt(2)).
    massive=m.witness(False);zero=m.witness(True)
    exact_massive=np.zeros((4,3,4));exact_massive[:,1:,0]=.25;exact_massive[:,1:,3]=.25
    mw=sum(abs(F(float(x.real))-F(float(y)))+abs(F(float(x.imag)))
           for x,y in zip(massive,exact_massive.ravel()))
    rootlo=F('1.414213562373095');roothi=F('1.4142135623730954')
    assert rootlo**2<2<roothi**2
    alo=1/(2*roothi);ahi=1/(2*rootlo)
    zw=F(0)
    for p in range(4):
        for mass in range(3):
            for d in range(4):
                z=zero[12*p+4*mass+d]; zr=F(float(z.real));zi=F(float(z.imag))
                if mass==0 and d in (0,3):zw+=max(abs(zr-alo),abs(zr-ahi))+abs(zi)
                else:zw+=abs(zr)+abs(zi)
    witness_l1=frac_up(max(mw,zw))
    witness_probability=frac_up((1+F(eps_matrix))**2*(2+F(witness_l1))*F(witness_l1))
    # Same analytic Fourier coefficients and conservative fixed-size complex
    # product/sum budget as 1175 endpoint audit; no spectral-norm output certified.
    effect_probability=frac_up(F(64)*F(gamma(4096)))
    allowance=frac_up(F(polar_probability)+F(witness_probability)+F(effect_probability))
    assert allowance<cert['terminal_probability_arithmetic_padding']
    result={'status':'POSTHOC_POLAR_ENDPOINT_GUARD_PASS',
      'replay':'Same initial state, derivative and RK4 update as numerical_bound.run; extra residual evaluations do not mutate state.',
      'eigh_accuracy_assumed':False,'full_residual_recomputed':False,
      'normalized_menu_bitwise_matched_float_leaves':matched,
      'raw_Gram_and_polar_bound_bitwise_match':True,
      'stored_matrix_max_abs_entry':float(np.max(abs(stored))),
      'raw_matrix_max_abs_entry':float(np.max(abs(raw))),
      'max_stored_minus_raw_Frobenius_diagnostic':frobenius,
      'exact_stored_minus_raw_entry_l1_rational':str(exact_l1),
      'stored_minus_raw_operator_upper':matrix_l1,
      'raw_to_exact_polar_upper':shift,'stored_to_exact_polar_upper':eps_matrix,
      'polar_probability_allowance':polar_probability,
      'exact_witness_input_l1_upper':witness_l1,
      'witness_probability_allowance':witness_probability,
      'analytic_Fourier_probability_arithmetic_allowance':effect_probability,
      'total_probability_allowance':allowance,
      'published_endpoint_guard':cert['terminal_probability_arithmetic_padding'],
      'scope':'Exact rational post-hoc comparison of stored matrices, plus declared gamma_n basic arithmetic. Not formal IEEE/hardware verification; reported SVD effect distances are diagnostic.',
      'source_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
          for name in ('preview.py','numerical_bound.py','numerical_bound_checks.json')}}
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    r=run();path=HERE/'polar_endpoint_guard.json'
    if args.write:
        with path.open('x',encoding='utf-8',newline='\n') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS')
    else:
        assert r==json.loads(path.read_text(encoding='utf-8'));print('READ_ONLY_MATCH')
    print(json.dumps({k:r[k] for k in ('max_stored_minus_raw_Frobenius_diagnostic','total_probability_allowance','published_endpoint_guard')}))
