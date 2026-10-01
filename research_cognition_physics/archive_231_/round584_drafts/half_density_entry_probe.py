"""584 entry only: direct five-coordinate operators and a local record difference."""
from pathlib import Path
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'half_density_entry_results.json'


def density(phi):return np.sqrt(original.M)*original.F(phi)**-3


def qterm(phi):return 2.5-np.dot(phi,phi)/(3*original.M)


def divG(phi):
    f=original.F(phi);r2=np.dot(phi,phi)
    return phi*(-1/3+r2/(18*original.M)-f/original.M)


def derivatives(fn,point,step):
    eye=np.eye(5)*step;f=fn(point)
    grad=np.array([(fn(point+e)-fn(point-e))/(2*step) for e in eye])
    hessian=np.empty((5,5),dtype=complex)
    for a in range(5):
        hessian[a,a]=(fn(point+eye[a])-2*f+fn(point-eye[a]))/step**2
        for b in range(a):
            hessian[a,b]=hessian[b,a]=(fn(point+eye[a]+eye[b])-fn(point+eye[a]-eye[b])
                -fn(point-eye[a]+eye[b])+fn(point-eye[a]-eye[b]))/(4*step**2)
    return grad,hessian


def flat_minus_div(fn,point,step):
    gradient,hessian=derivatives(fn,point,step)
    return -np.sum(original.inverse(point)*hessian)-divG(point)@gradient


def operator_check():
    chi=lambda x:np.exp(-.21*np.dot(x,x)+1j*(.17*x[1]+.09*x[2]**2))*(1+.12*x[0]+.07*x[2]*x[4])
    psi=lambda x:chi(x)/np.sqrt(density(x))
    points=[np.zeros(5),np.array([.14,.53,-.17,.1,.42]),np.array([.7,-.8,.4,.9,1.1])]
    rows=[]
    for step in (.004,.002,.001):
        worst=0.;qerror=0.;unshifted=0.
        for point in points:
            f=original.F(point);G=original.inverse(point)
            dpsi,hpsi=derivatives(psi,point,step)
            drift=divG(point)+G@(point/f)
            transformed=-np.sqrt(density(point))*(np.sum(G*hpsi)+drift@dpsi)
            bare=flat_minus_div(chi,point,step)
            expected=bare+qterm(point)*chi(point)
            worst=max(worst,float(abs(transformed-expected)))
            unshifted=max(unshifted,float(abs(transformed-bare)))
            # Divergence of the actual G grad(log mu), independently differenced.
            vector=lambda x:original.inverse(x)@(x/original.F(x))
            divergence=sum((vector(point+step*e)[j]-vector(point-step*e)[j])/(2*step)
                for j,e in enumerate(np.eye(5)))
            gl=point/f
            derived=.5*divergence+.25*gl@G@gl
            qerror=max(qerror,abs(float(derived-qterm(point))))
        rows.append(dict(step=step,full_operator_error=worst,q_from_density_error=qerror,
            nonzero_error_without_q=unshifted))
    assert rows[-1]['full_operator_error']<rows[0]['full_operator_error']/12
    assert rows[-1]['q_from_density_error']<rows[0]['q_from_density_error']/12
    assert rows[-1]['nonzero_error_without_q']>1
    return dict(rows=rows,local_test_functions_only=True)


def record_check():
    point=np.array([.14,.53,-.17,.1,.42]);hbar=.7;weight=.8
    wave=lambda x:np.exp(-.21*np.dot(x,x))
    b=lambda x:np.sin(x[4]);Q=lambda x:hbar*hbar*qterm(x)/(2*weight)
    expected=hbar*hbar*original.F(point)**2*point[4]*np.cos(point[4])/(3*original.M**2*weight**2)
    rows=[]
    for step in (.004,.002,.001):
        H=lambda fn:hbar*hbar/(2*weight)*flat_minus_div(fn,point,step)+original.node_potential(point)*fn(point)
        commutator=(Q(point)*H(lambda x:b(x)*wave(x))-Q(point)*b(point)*H(wave)
            -H(lambda x:b(x)*Q(x)*wave(x))+b(point)*H(lambda x:Q(x)*wave(x)))
        actual=-commutator/(hbar*hbar*wave(point))
        error=float(abs(actual-expected))
        rows.append(dict(step=step,actual_double_commutator=float(actual.real),
            imaginary_error=float(abs(actual.imag)),error=error))
    assert rows[-1]['error']<rows[0]['error']/8 and expected>0
    return dict(expected_local_acceleration=float(expected),rows=rows,
        observable='b=sin(s); (D_correct^2-D_bare^2)b',
        not_a_full_wavepacket_or_time_evolution_simulation=True)


def run():
    return dict(status='584 entry; not completed round',checks_passed=2,
        operators=operator_check(),record_difference=record_check(),
        scope='pointwise full five-dimensional operator diagnostics on interior patches; finite graph domain, Gauss states, actual finite-time record theorem and geometric source still require formal completion')


if __name__=='__main__':
    result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if '--write' in sys.argv:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
