"""581 entry: actual metric/field substitution, finite remainder and length pullback."""
from pathlib import Path
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'field_redefinition_results.json'


def potential_derivatives(phi):
    # Full five-component original potential V/F^2, analytic derivatives.
    L,u,_=original.lattice.scalar.parameters()
    d=np.stack((np.sum(phi[...,:4]**2,axis=-1)-u[0],phi[...,4]**2-u[1]),axis=-1)
    z=d@L;V=np.sum(d*z,axis=-1)/4
    gV=phi*np.concatenate((np.repeat(z[...,0,None],4,axis=-1),z[...,1,None]),axis=-1)
    hV=np.zeros(phi.shape[:-1]+(5,5))
    hV[...,:4,:4]=z[...,0,None,None]*np.eye(4)+2*L[0,0]*phi[...,:4,None]*phi[...,None,:4]
    hV[...,4,4]=z[...,1]+2*L[1,1]*phi[...,4]**2
    hV[...,:4,4]=2*L[0,1]*phi[...,:4]*phi[...,4,None]
    hV[...,4,:4]=hV[...,:4,4]
    f=original.F(phi)
    gradient=gV/f[...,None]**2+2*V[...,None]*phi/(3*f[...,None]**3)
    hessian=(hV/f[...,None,None]**2
        +2*(gV[..., :,None]*phi[...,None,:]+phi[..., :,None]*gV[...,None,:])/(3*f[...,None,None]**3)
        +2*V[...,None,None]*np.eye(5)/(3*f[...,None,None]**3)
        +2*V[...,None,None]*phi[..., :,None]*phi[...,None,:]/(3*f[...,None,None]**4))
    inverse=original.inverse(phi)
    laplacian=np.einsum('...ab,...ab->...',inverse,hessian)-f*np.sum(phi*gradient,axis=-1)/(3*original.M)
    return gradient,laplacian


def derivative(values):
    freq=np.fft.fftfreq(len(values),1/len(values))
    if values.ndim>1:freq=freq[:,None]
    return np.fft.ifft(1j*freq*np.fft.fft(values,axis=0),axis=0).real


def run():
    x=np.linspace(0,2*np.pi,2048,endpoint=False)
    q=.73+.17*np.sin(x);qp=.17*np.cos(x);qpp=-.17*np.sin(x)
    direction=np.array([0.,np.sqrt(.6),0.,0.,np.sqrt(.4)])
    phi=np.sqrt(6*original.M)*np.tanh(q[:,None]/np.sqrt(6))*direction
    tangent=np.sqrt(original.M)/np.cosh(q[:,None]/np.sqrt(6))**2*direction
    gradient,trA=potential_derivatives(phi)
    inverse=original.inverse(phi)
    gradvector=np.einsum('...ab,...b->...a',inverse,gradient)
    S=qp*qp;U=original.node_potential(phi)
    C=(S+4*U)/24-trA/6-S/9
    D1=(S+U)/12+C;Dt=U/12+C
    trD=D1+3*Dt;T1=-2*D1+trD;Tt=-2*Dt+trD
    H1=S/4-U/2;Ht=-S/4-U/2
    Ephi_gradient=np.sum(gradient*(gradvector-tangent*qpp[:,None]),axis=1)
    expected=2*np.pi*np.mean(-H1*T1-3*Ht*Tt+Ephi_gradient/6)
    expected_length=2*np.pi*np.mean(T1)/2
    integral=lambda a:float(2*np.pi*np.mean(a))

    def action(t):
        changed_phi=phi+t*gradvector/6
        assert np.min(original.F(changed_phi))>1
        assert np.min(1-t*T1)>.9 and np.min(1-t*Tt)>.9
        a=-.5*np.log1p(-t*T1);b=-.5*np.log1p(-t*Tt)
        ap=derivative(a);bp=derivative(b);bpp=derivative(bp)
        R=-6*np.exp(-2*a)*(bpp+2*bp*bp-ap*bp)
        pprime=derivative(changed_phi)
        Sc=np.exp(-2*a)*np.einsum('...a,...ab,...b->...',pprime,original.metric(changed_phi),pprime)
        density=np.exp(a+3*b)*(-R/2+Sc/2+original.node_potential(changed_phi))
        return integral(density),integral(np.exp(a))

    initial,initial_length=action(0.)
    rows=[]
    for t in (.04,.02,.01,.005):
        new,length=action(t)
        remainder=new-initial-t*expected
        length_remainder=length-initial_length-t*expected_length
        rows.append(dict(parameter=t,action=new,linear_prediction=initial+t*float(expected),
            action_remainder=float(remainder),action_remainder_over_parameter_squared=float(remainder/t**2),
            pulled_back_length=length,length_linear_coefficient=float(expected_length),
            length_remainder=float(length_remainder)))
    ratios=[rows[i]['action_remainder']/rows[i+1]['action_remainder'] for i in range(3)]
    assert all(3.8<r<4.2 for r in ratios)
    assert abs(expected_length)>1e-3
    assert abs(rows[-1]['length_remainder'])<abs(rows[0]['length_remainder'])/50
    # Cross-check the analytic full potential gradient against original U directly.
    a=phi[333];g,_=potential_derivatives(a);step=1e-5
    numeric=np.array([(original.node_potential(a+step*np.eye(5)[j])-original.node_potential(a-step*np.eye(5)[j]))/(2*step) for j in range(5)])
    assert np.max(abs(g-numeric))<1e-8
    return dict(status='581 entry only; actual transformation test',checks_passed=2,
        initial_action=initial,first_action_variation=float(expected),rows=rows,
        successive_action_remainder_ratios=ratios,potential_gradient_error=float(np.max(abs(g-numeric))),
        untransformed_length=float(initial_length),
        scope='same original U and K on a periodic off-shell scalar profile; anisotropic positive Euclidean metric after substitution; verifies leading-action O(lambda^2) residual and nontrivial observable pullback, not full quantum equivalence or global invertibility')


if __name__=='__main__':
    result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if '--write' in sys.argv:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
