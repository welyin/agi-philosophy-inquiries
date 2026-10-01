"""654 entry: explicit neutral mass insertion, not an inherited Yukawa prescription.

Original complex coefficient is reused. The Euclidean physical-fermion mass
action and m=a*mu/2 normalization are declared choices to audit, not proven
equivalent to the full original matter Hamiltonian.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_fermion_gauss_completion as original
import joint_subgroup_measure_source as measure


def coefficient(s):
    phi=np.array([0.,0.,0.,0.,s])
    _,d=original.mass_matrices(phi)
    return complex(d[30,31])


def raw_partition(s,a,n):
    m=a*coefficient(s)/2
    shift=np.zeros((n,n))
    for j in range(n):shift[j,(j+1)%n]=-1 if j==n-1 else 1
    kinetic=np.kron((np.eye(n)-shift.T)/2,np.eye(2))
    eps=np.array([[0.,1.],[-1.,0.]])
    pair=np.kron(np.eye(n),m*eps)
    anti=np.block([[pair,-kinetic.T],[kinetic,-pair.conj()]])
    pf=measure.pfaffian(anti)*(-1)**n
    assert abs(pf.imag)<1e-12 and pf.real>0
    return float(pf.real)


def run():
    rows=[]
    for s in (.2,.6,1.2,2.8):
        mu=abs(coefficient(s))
        for a,n in ((1.,2),(.5,4),(.25,8)):
            beta=a*n;energy=2/a*np.arcsinh(a*mu/2)
            predicted=(2+2*np.cosh(beta*energy))/4**n
            actual=raw_partition(s,a,n)
            relative=abs(actual-predicted)/predicted
            assert relative<1e-11
            step=1e-5
            source=(np.log(raw_partition(s+step,a,n))-np.log(raw_partition(s-step,a,n)))/(2*step)
            f=original.original.F(np.array([0.,0.,0.,0.,s]))
            dmu=abs(original.Y['s'])*original.original.M/f**1.5
            de=dmu/np.sqrt(1+(a*mu/2)**2)
            expected_source=beta*np.tanh(beta*energy/2)*de
            original_source=beta*np.tanh(beta*mu/2)*dmu
            assert abs(source-expected_source)<1e-7
            rows.append(dict(s=s,time_step=a,time_sites=n,beta=beta,original_mass=mu,
                effective_energy=float(energy),mass_difference=float(energy-mu),
                Pfaffian_relative_error=float(relative),lattice_source=float(source),
                spectral_source=float(expected_source),original_H_source=float(original_source)))
    return dict(date='2026-10-02',status='entry probe;654 incomplete',
        candidate='A=[[m eps,-d.T],[d,-conj(m eps)]], d=(I-S_AP.T)/2, m=a*original_mu/2',
        rows=rows,dependency_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        for name in ('joint_fermion_gauss_completion.py','joint_subgroup_measure_source.py')},
        scope='Constant neutral Majorana background, trivial spatial sector, explicit candidate insertion; no full interacting/Gauss quantum equivalence.')


if __name__=='__main__':
    result=run()
    with (HERE/'neutral_mass_transfer_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=result['status'],cases=len(result['rows']),
        maximum_Pfaffian_error=max(r['Pfaffian_relative_error'] for r in result['rows']))))
