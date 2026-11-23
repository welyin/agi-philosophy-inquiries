"""Replay the three-family formation test with NumPy. Run with --check or --write-new."""
import argparse
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).with_name("results.json")
I = np.eye(3)
Jplus = np.sqrt(2) * np.diag(np.ones(2), 1)
Jz = np.diag([1., 0., -1.])
X = np.block([[Jz, Jplus.T], [Jplus, -Jz]])
START = np.vstack([I, 0*I])
OTHER = np.vstack([0*I, I])

def portal(a, phase):
    return np.vstack([np.cos(a)*I, np.exp(1j*phase)*np.sin(a)*I])

def pulse(theta):
    a, v = np.linalg.eigh(X)
    return (v * np.exp(-1j*theta*a)) @ v.conj().T

def department(unitary, port, t, heavy):
    left = unitary @ START
    complement = unitary @ OTHER
    read = port.conj().T @ left
    effect = read.conj().T @ read
    mass = np.column_stack([t * port, heavy * complement])
    ev, vec = np.linalg.eigh(mass @ mass.conj().T)
    vals = np.linalg.svd(read, compute_uv=False)
    total = t*t + heavy*heavy
    disc = np.sqrt(np.maximum(0.,(heavy*heavy-t*t)**2 + 4*heavy*heavy*t*t*(1-vals**2)))
    high = (total + disc)/2
    low = (heavy*t*vals)**2 / high
    assert np.allclose(ev, np.sort(np.r_[low,high]), atol=3e-13,rtol=3e-13)
    assert vals[1] + 1e-12 >= vals[0]/3
    if ev[2]>1e-10:
        assert np.sqrt(max(ev[1],0)/ev[2])+1e-8 >= (np.sqrt(2)-1)/np.sqrt(3)
    return effect, ev.tolist(), vec[:,:3]

def measure(history, projector):
    # Record a distinct branch for each actual binary outcome.
    return {word+str(bit): p @ state for word,state in history.items()
            for bit,p in ((1,projector),(0,np.eye(6)-projector))}

def compare_frames():
    unitary=pulse(.63)
    u=portal(.47,.23)
    eu,mu,lu=department(unitary,u,.2,7.)
    rows=[]
    for label,shift in [('shared',0.),('independent_down_phase',np.pi)]:
        d=portal(.88,1.31+shift)
        ed,md,ld=department(unitary,d,.3,7.)
        hist=measure({'':unitary@START},u@u.conj().T)
        hist=measure({w:unitary.conj().T@v for w,v in hist.items()},START@START.T)
        hist=measure({w:unitary@v for w,v in hist.items()},d@d.conj().T)
        gram=sum(v.conj().T@v for v in hist.values())
        assert np.allclose(gram,I)
        probs={w:float(np.linalg.norm(v,'fro')**2/3) for w,v in hist.items()}
        p111=float(np.trace(eu@eu@ed).real/3)
        assert abs(probs['111']-p111)<1e-12
        v=lu.conj().T@ld
        comm=eu@ed-ed@eu
        rows.append(dict(label=label,up_squared_masses=mu,down_squared_masses=md,
                         outcome_probabilities=probs,p111=p111,
                         light_abs_current=np.abs(v).tolist(),
                         imaginary_commutator_trace=float(np.trace(comm@comm@comm).imag)))
    assert np.allclose(rows[0]['down_squared_masses'],rows[1]['down_squared_masses'])
    assert abs(rows[0]['p111']-rows[1]['p111'])>1e-3
    assert np.max(np.abs(np.array(rows[0]['light_abs_current'])-rows[1]['light_abs_current']))>1e-3
    return rows

def run():
    cases=0
    for theta in [0.,.13,.63,np.pi/3,np.pi,2*np.pi/3]:
        for a in [0.,.47,.88,np.pi/2]:
            for t,M in [(1.,1.001),(.3,7.),(2.,.5)]:
                department(pulse(theta),portal(a,.7),t,M)
                cases+=1
    exact=pulse(np.pi/3)
    s=np.linalg.svd(START.T@exact@START,compute_uv=False)
    assert np.allclose(s,[1,1/3,1/3])
    return dict(round=1335,scope="finite single exchange; three flavor roles supplied",status="PASS",
                parent_spectrum_cases=cases,sharp_leading_ratio=1/3,
                sharp_full_infimum=float((np.sqrt(2)-1)/np.sqrt(3)),comparison=compare_frames())

def check(actual, saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for k in actual: check(actual[k],saved[k])
    elif isinstance(actual,list):
        assert len(actual)==len(saved)
        for a,b in zip(actual,saved):check(a,b)
    elif isinstance(actual,float):
        assert np.isclose(actual,saved,atol=1e-11,rtol=1e-10)
    else: assert actual==saved

if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group(required=True)
    for option in ("--check","--write-new"):
        group.add_argument(option,action="store_true")
    args=p.parse_args(); data=run()
    if args.write_new:
        with OUT.open("x",encoding="utf8") as f: json.dump(data,f,indent=2)
    else:check(data,json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(status=data['status'],cases=data['parent_spectrum_cases'],
                          p111=[v['p111'] for v in data['comparison']])))