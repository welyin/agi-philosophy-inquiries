"""886: original unbroken electroweak direction and local-counterterm test.
Not a no-go for reference-inverse counterterms defined only on regular patches.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'885'))
import joint_static_thermal_weight as weight
old=weight.old
TARGET=HERE/'unbroken_holonomy_locality_results.json'

def em_charges(m,matter):
    qy,_=old.charges(m);p=qy[:32].copy()
    p[:12]+=np.tile(np.repeat([3.,-3.],2),3)
    p[24:28]+=np.repeat([3.,-3.],2)
    return np.r_[p,-p]

def joint_spectrum(m,q):
    M=np.exp(old.XI)*m.MASS
    charges=[];mass2=[]
    for value in np.unique(q):
        ids=np.flatnonzero(q==value)
        eig=np.linalg.eigvalsh((M@M)[np.ix_(ids,ids)])
        charges.extend([float(value)]*len(ids));mass2.extend(eig.tolist())
    assert min(mass2)>0
    return np.array(charges),np.array(mass2)

def theta(t,a,dual):
    j=np.arange(-20,21)
    if dual:
        return np.sqrt(np.pi/t)*float(np.sum(np.exp(-np.pi**2*j*j/t)*np.cos(2*np.pi*j*a)))
    return float(np.sum(np.exp(-t*(j+a)**2)))

def heat_difference(t,alpha,q,m2):
    if t<=1:
        j=np.arange(1,21)
        theta0=np.sqrt(np.pi/t)*(1+2*np.sum(np.exp(-np.pi**2*j*j/t)))
        delta=2*np.sqrt(np.pi/t)*np.sum(np.exp(-np.pi**2*j[:,None]**2/t)*
            (-2*np.sin(np.pi*j[:,None]*alpha*q[None,:])**2),axis=0)
    else:
        j=np.arange(-20,21)
        theta0=float(np.sum(np.exp(-t*j*j)))
        delta=np.sum(np.exp(-t*(j[:,None]+alpha*q[None,:])**2),axis=0)-theta0
    return float(theta0**2*np.dot(np.exp(-t*m2),delta))

def continuum_vacuum(alpha,q,m2,order):
    x,w=np.polynomial.legendre.leggauss(order)
    bounds=[np.log(.03),0.,np.log(80/min(m2))]
    total=0.
    for a,b in zip(bounds,bounds[1:]):
        u=(a+b)/2+(b-a)*x/2
        integrand=[np.exp(-v/2)*heat_difference(np.exp(v),alpha,q,m2) for v in u]
        total+=float(np.dot(w,integrand))*(b-a)/2
    return -old.BETA/(8*np.sqrt(np.pi))*total

def continuum_excitation(alpha,q,m2,N):
    k,r,c=weight.radial_counts(N)
    total=0.
    for x in k:
        ea=np.sqrt((x+alpha*q)[:,None]**2+m2[:,None]+r)
        eb=np.sqrt(x*x+m2[:,None]+r)
        total+=.5*float(np.sum((np.log1p(np.exp(-old.BETA*ea))-
                           np.log1p(np.exp(-old.BETA*eb)))*c))
    return total

def run():
    m=old.load();_,matter=old.charges(m);q=em_charges(m,matter)
    M=np.exp(old.XI)*m.MASS
    errors={}
    errors['mass_charge_commutator']=float(np.max(abs(M@np.diag(q)-np.diag(q)@M)))
    errors['Higgs_stabilizer']=float(np.linalg.norm(np.diag([6.,0.])@m.PHI[:2]))
    errors['group_representation']=0.
    for t in (.07,.39):
        R=matter.representation(np.eye(3),np.diag(np.exp(1j*t*np.array([3.,-3.]))),np.exp(1j*t))
        errors['group_representation']=max(errors['group_representation'],
            float(np.max(abs(R-np.diag(np.exp(1j*t*q[:32]))))))
    errors['matrix_square']=0.
    for k in ([0,0,0],[2,-1,3]):
        h=old.H(m,q,17,k,.04)
        sq=np.diag((k[0]+.04*q)**2+k[1]**2+k[2]**2)+M@M
        errors['matrix_square']=max(errors['matrix_square'],float(np.max(abs(h@h-sq))))
    assert max(errors.values())<3e-12
    coeff=weight.constants(q)['patch_N_squared_coefficient']
    rows=[]
    for N in (17,33,65,129):
        k,rad,count=weight.radial_counts(N)
        base=weight.x_squares(m,q,N,0.)
        d=weight.log_trace_difference(weight.x_squares(m,q,N,.04),base,rad,count)
        rows.append(dict(N=N,alpha=.04,**d,scaled_coefficient=d['logZ_difference']/(.04**2*N*N)))
    qq,m2=joint_spectrum(m,q)
    poisson=[]
    for t in (.3,1.,3.,8.):
        for a in (0.,.08,.24):
            err=abs(theta(t,a,False)-theta(t,a,True))
            poisson.append(dict(proper_time=t,shift=a,absolute_error=err))
            assert err<3e-14
    v96=continuum_vacuum(.04,qq,m2,96);v192=continuum_vacuum(.04,qq,m2,192)
    e33=continuum_excitation(.04,qq,m2,33);e49=continuum_excitation(.04,qq,m2,49)
    assert abs(v96-v192)<3e-10 and abs(e33-e49)<3e-10
    particle,counts=np.unique(q[:32],return_counts=True)
    return dict(round=886,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3671,
        argument_scope='For883 hopping, the original unbroken electroweak flat family has a divergent holonomy-dependent vacuum determinant. No finite-jet gauge-covariant counterterm regular on this family can subtract its relative divergence. The finite covariant heat-kernel relative target exists. This does not disprove the regular-reference-only773-796 branch or the full dynamic Gauss bridge.',
        original_species=32,Nambu_components=64,original_mass_and_Higgs_preserved=True,
        original_unbroken_particle_charges={str(int(v)):int(c) for v,c in zip(particle,counts)},
        full_charge_square_trace=float(np.sum(q*q)),original_mass_gap=float(np.sqrt(min(m2))),
        matrix_checks=errors,patch_leading_coefficient=coeff,rows=rows,Poisson_checks=poisson,
        continuum_relative=dict(alpha=.04,vacuum_heat_integral=v192,
            excitation_difference=e49,total_logZ_difference=v192+e49,
            vacuum_quadrature96_vs192=abs(v96-v192),
            excitation_cube33_vs49=abs(e33-e49),
            integration_proper_time_min=.03,integration_proper_time_max=float(80/min(m2)),
            integral_finiteness_is_analytic=True,numerical_value_is_calibration_only=True),
        regular_local_counterterm_obstruction_proved=True,
        all_reference_patch_counterterms_excluded=False,
        full_dynamic_Gauss_measure_derived=False,
        all_regulators_excluded=False,full_interacting_Q_E_bridge_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
