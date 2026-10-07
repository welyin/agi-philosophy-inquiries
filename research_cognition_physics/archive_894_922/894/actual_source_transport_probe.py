"""Working894: original shared-past pulse versus added horizontal path generator.
No formal round count yet. Prescribed physical-background path only; this does
not alone compare fully integrated autonomous gauge-coordinate experiments.
"""
from pathlib import Path
import json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'893'))
import curved_connection_heat_bridge as prior
old=prior.old
TARGET=HERE/'actual_source_transport_probe_results.json'

def pulse(t):
    z=2*t-1
    if abs(z)>=1:return 0.,0.
    f=np.exp(1-1/(1-z*z))
    return f,-4*z*f/(1-z*z)**2

def propagate(h0,G,eps,steps,transport):
    U=np.eye(4,dtype=complex);dt=1/steps
    for j in range(steps):
        f,fp=pulse((j+.5)*dt);h=h0+eps*f*G
        if transport:
            E=np.sqrt(np.trace(h@h).real/4)
            Ep=np.trace(h@G).real/(4*E)
            S=h/E;Sp=G/E-h*Ep/E**2
            A=(Sp@S-S@Sp)/4
            h=h+1j*eps*fp*A
        ev,V=np.linalg.eigh(h)
        U=((V*np.exp(-1j*dt*ev))@V.conj().T)@U
    return U

def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.prior.prior.prior.last.em_charges(m,matter)
    ids=[26,27,28,29];a0=-1/12
    h0=old.H(m,q,17,[0,0,0],a0,kind='continuum')[np.ix_(ids,ids)]
    G=(m.GAMMA[0]@np.diag(q))[np.ix_(ids,ids)]
    ev,V=np.linalg.eigh(h0);P=V[:,:2]@V[:,:2].conj().T;Q=np.eye(4)-P
    E=float(ev[-1])
    x,w=np.polynomial.legendre.leggauss(160);t=(x+1)/2
    integral=sum(ww/2*pulse(tt)[0]*np.exp(2j*E*tt) for tt,ww in zip(t,w))
    coefficient=float(np.trace(P@G@Q@G).real*abs(integral)**2)
    assert coefficient>0 and E<np.pi/2
    rows=[]
    for eps in (.004,.002,.001):
        for steps in (512,1024):
            row=dict(amplitude=eps,steps=steps)
            for transport in (False,True):
                U=propagate(h0,G,eps,steps,transport)
                pairs=float(np.trace(Q@U@P@U.conj().T).real)
                assert np.max(abs(U.conj().T@U-np.eye(4)))<1e-11
                row['horizontal_pairs' if transport else 'original_pairs']=pairs
            row['original_pairs_over_amplitude_squared']=row['original_pairs']/eps**2
            rows.append(row)
    assert abs(rows[-1]['original_pairs_over_amplitude_squared']-coefficient)<.01*coefficient
    assert max(abs(r['horizontal_pairs']) for r in rows)<1e-10
    return dict(status='working894_not_formal',date='2026-10-06',
        formal_round_still=893,cumulative_numbered_groups_still=3678,
        conditional_contract='Same prescribed smooth closed physical EM history, same past vacuum and final original particle readout.',
        pulse_duration=1.,base_alpha=a0,original_electron_gap=2*E,
        positive_original_source_pair_coefficient=coefficient,
        pulse_Fourier_integral_real=float(integral.real),pulse_Fourier_integral_imag=float(integral.imag),
        actual_evolution_rows=rows,
        scope='Original path h(a(t)) differs from h+i dot(a)[P_a,P]. The latter is the horizontal transitionless path kernel. Scalar vacuum phases cannot alter this pair count. Relation to fully autonomous integrated dynamics still requires proof.')
if __name__=='__main__':
    r=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
