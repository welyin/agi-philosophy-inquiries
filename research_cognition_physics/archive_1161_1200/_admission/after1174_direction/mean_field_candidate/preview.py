"""Unnumbered feasibility preview. Full finite source-zone Hartree comparator.

Default compares saved JSON; --write creates once. No parent K-body simulation.
No claim that RK step doubling is a certified numerical error bound.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
GRID, LENGTH, TIME, STEPS = 65, 16.0, 2.0, 2048
EPS, NU, OMEGA, MS = .15, .35, .5, 2.0
MP = np.array([.3, .5])
DX = LENGTH / GRID
VOL = GRID ** 3
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.]).astype(complex)
BETA = np.kron(SZ, I2)
ALPHA = np.array([np.kron(SX, s) for s in [SX, SY, SZ]])
MODES = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]])
PMOM = 2 * np.pi * MODES / LENGTH
LMOM = np.sin(DX * PMOM) / DX
KIN = np.einsum('pb,bij->pij', LMOM, ALPHA)
CHI = np.array([1,0,0,1], complex) / np.sqrt(2)
CLOCK = (I2 - SY) / 2

def source_classes():
    s, weights = [], []
    for trip in itertools.combinations_with_replacement(range((GRID+1)//2), 3):
        counts = [trip.count(t) for t in set(trip)]
        perm = math.factorial(3)
        for c in counts:
            perm //= math.factorial(c)
        multiplicity = perm * 2 ** sum(t != 0 for t in trip)
        s.append(np.linalg.norm(np.sin(2*np.pi*np.array(trip)/GRID)) / DX)
        weights.append(multiplicity)
    assert sum(weights) == VOL and len(s) == 6545
    return np.array(s), np.array(weights, float)/VOL

SVAL, WEIGHTS = source_classes()

def exact_const(nclock, nspace, masses=MP):
    h = nspace*KIN[:,None,:,:] + nclock*masses[None,:,None,None]*BETA
    rad = np.sqrt((nspace**2*np.sum(LMOM**2, axis=1))[:,None]
                  +(nclock*masses[None,:])**2)
    # Arguments are already integrated coefficients, so this is exp(-i h).
    return np.cos(rad)[...,None,None]*np.eye(4) - 1j*np.sinc(rad/np.pi)[...,None,None]*h

def derivative(state, nu, vacuum):
    r, g, accum, u = state
    lapse, velo = 1+EPS*g[2], 1+nu*g[2]
    em = 0.0 if vacuum else 2*MS*np.dot(WEIGHTS, r[:,2])
    et = 0.0 if vacuum else 2*np.dot(WEIGHTS, SVAL*r[:,0])
    field = EPS*em + nu*et
    dr = np.empty_like(r)
    dr[:,0] = -2*lapse*MS*r[:,1]
    dr[:,1] = 2*(lapse*MS*r[:,0]-velo*SVAL*r[:,2])
    dr[:,2] = 2*velo*SVAL*r[:,1]
    if vacuum:
        dr.fill(0)
    dg = 2*np.cross(np.array([-OMEGA,0,field]), g)
    h = lapse*MP[None,:,None,None]*BETA + velo*KIN[:,None,:,:]
    du = -1j*(h@u)
    return dr, dg, np.array([lapse,velo]), du

def plus(x, y, factor):
    return tuple(a+factor*b for a,b in zip(x,y))

def evolve(nu, vacuum=False, steps=STEPS):
    r = np.zeros((len(SVAL),3)); r[:,2] = 0 if vacuum else 1
    state = (r, np.array([1.,0.,0.]), np.zeros(2),
             np.broadcast_to(np.eye(4,dtype=complex),(4,2,4,4)).copy())
    dt = TIME/steps
    history = []
    energy0 = 0 if vacuum else 2*MS
    energy0 -= OMEGA
    err_energy = 0.0
    for k in range(steps+1):
        r,g,accum,u=state
        em=0 if vacuum else 2*MS*np.dot(WEIGHTS,r[:,2])
        et=0 if vacuum else 2*np.dot(WEIGHTS,SVAL*r[:,0])
        n,v=1+EPS*g[2],1+nu*g[2]
        ene=n*em+v*et-OMEGA*g[0]
        err_energy=max(err_energy,abs(ene-energy0))
        if k%(steps//16)==0:
            history.append({'t':k*dt,'z':float(g[2]),'N':float(n),'v':float(v),
                            'a':float(n/v),'Em':float(em),'Et':float(et)})
        if k==steps:
            break
        k1=derivative(state,nu,vacuum)
        k2=derivative(plus(state,k1,dt/2),nu,vacuum)
        k3=derivative(plus(state,k2,dt/2),nu,vacuum)
        k4=derivative(plus(state,k3,dt),nu,vacuum)
        state=tuple(a+dt*(b+2*c+2*d+e)/6 for a,b,c,d,e in zip(state,k1,k2,k3,k4))
    r,g,accum,u=state
    errors={'source_bloch_norm':float(np.max(np.abs(np.sum(r*r,axis=1)-(0 if vacuum else 1)))),
            'G_bloch_norm':float(abs(np.dot(g,g)-1)),
            'energy_per_site':float(err_energy),
            'probe_unitarity':float(max(np.linalg.norm(v.conj().T@v-np.eye(4),2) for v in u.reshape(-1,4,4)))}
    return {'u':u,'tau':accum[0],'xi':accum[1],'history':history,'errors':errors}

def site_effect(axis):
    # Exact finite-site sampling effect, then Fourier compression to input band.
    j=np.arange(GRID)
    f=(1+np.sin(2*np.pi*j/GRID))/2
    mat=np.empty((4,4),complex)
    for p in range(4):
        for q in range(4):
            delta=MODES[q]-MODES[p]
            others=np.delete(delta,axis)
            mat[p,q]=0 if np.any(others) else np.mean(f*np.exp(2j*np.pi*delta[axis]*j/GRID))
    return mat

def effect_matrix(u, axis, clock=True):
    f=site_effect(axis)
    c=CLOCK if clock else I2
    out=np.zeros((32,32),complex)
    for p,q,m,n in itertools.product(range(4),range(4),range(2),range(2)):
        out[8*p+4*m:8*p+4*m+4,8*q+4*n:8*q+4*n+4]=f[p,q]*c[m,n]*(u[p,m].conj().T@u[q,n])
    return (out+out.conj().T)/2

WITNESS=np.tile(np.tile(CHI,2)/np.sqrt(2),4)/2

def readouts(u, reference):
    out={}
    for axis,name in enumerate('xyz'):
        for clock in [False,True]:
            e=effect_matrix(u,axis,clock); e0=effect_matrix(reference,axis,clock)
            key=name+('_joint_clock_minusY' if clock else '_position')
            pa=float(np.vdot(WITNESS,e@WITNESS).real)
            pb=float(np.vdot(WITNESS,e0@WITNESS).real)
            out[key]={'active':pa,'proper_time_flat':pb,'signed_gap':pa-pb,
                      'complete_input_effect_distance':float(np.linalg.norm(e-e0,2))}
    return out

def calculate():
    active=evolve(NU)
    fine=evolve(NU,steps=2*STEPS)
    lapse=evolve(EPS)
    vacuum=evolve(NU,vacuum=True)
    deltau=float(abs(active['tau']-fine['tau']))
    dU=float(max(np.linalg.norm(a-b,2) for a,b in zip(active['u'].reshape(-1,4,4),fine['u'].reshape(-1,4,4))))
    activeflat=exact_const(active['tau'],active['tau'])
    lapseflat=exact_const(lapse['tau'],lapse['tau'])
    massless=exact_const(0,active['xi'],np.zeros(2))
    masslessflat=exact_const(0,active['tau'],np.zeros(2))
    bS=2*EPS*MS+6*NU/DX
    bP=2*EPS*sum(MP)+12*NU/DX
    w=2*VOL*bS; Lbound=VOL*bP
    etaK=1e-4
    logK=2*(math.log10(4*Lbound*TIME/etaK)+6*w*TIME/math.log(10))
    lat=TIME*(1+NU)*DX**2*max(np.sqrt(np.sum(PMOM**6,axis=1)))/6
    lat_control=TIME*(1+EPS)*DX**2*max(np.sqrt(np.sum(PMOM**6,axis=1)))/6
    proper_flat_lat=float(active['tau']*DX**2*max(np.sqrt(np.sum(PMOM**6,axis=1)))/6)
    summaries={}
    for name,run in [('active',active),('pure_lapse',lapse),('vacuum',vacuum)]:
        summaries[name]={'tau':float(run['tau']),'xi':float(run['xi']),
                         'history':run['history'],'diagnostic_residuals':run['errors']}
    return {
        'kind':'UNNUMBERED_FEASIBILITY_PREVIEW',
        'parameters':{'grid':GRID,'V':VOL,'L':LENGTH,'ell':DX,'T':TIME,'steps':STEPS,
                      'epsilon':EPS,'nu':NU,'Omega':OMEGA,'source_mass':MS,'probe_masses':MP.tolist(),
                      'source_classes':len(SVAL),'probe_modes':MODES.tolist()},
        'scope':{'full_parent_Kbody_was_simulated':False,'Hartree_numeric_error_certified':False,
                 'actual_instruments_constructed':False,'GR_dynamics_claimed':False,
                 'source_full_CAR_not_Dirac_sea':True,'all_probe_frequency_branches_kept':True,
                 'terminal_sampling_menu_only':True},
        'background':summaries,
        'active_mobile_clock_and_site':readouts(active['u'],activeflat),
        'massless_site_and_control':readouts(massless,masslessflat),
        'pure_lapse_proper_time_control':{'full_prepared_unitary_error':float(max(np.linalg.norm(a-b,2) for a,b in zip(lapse['u'].reshape(-1,4,4),lapseflat.reshape(-1,4,4)))),
                                        'readouts':readouts(lapse['u'],lapseflat)},
        'vacuum_flat_control':{'full_prepared_unitary_error':float(max(np.linalg.norm(a-b,2) for a,b in zip(vacuum['u'].reshape(-1,4,4),exact_const(TIME,TIME).reshape(-1,4,4))))},
        'analytic_initial_tidal':{'z_second':-1.2,'a_second':.24,'FLRW_Kretschmann':12*.24**2},
        'error_parts':{'meanfield_target_eta':etaK,'w_bound':w,'probe_L_bound':Lbound,
                       'sufficient_log10_K':logK,'lattice_active_isometry_bound':float(lat),
                       'lattice_control_isometry_bound':float(lat_control),
                       'proper_time_flat_lattice_bound':proper_flat_lat,
                       'RK4_step_doubling_unitary_difference_DIAGNOSTIC_ONLY':dU,
                       'RK4_step_doubling_tau_difference_DIAGNOSTIC_ONLY':deltau,
                       'instrument_error':'Not yet certified; explicit future finite permission required'},
    }

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for k,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{k}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-10),(path,a,b)
    else: assert a==b,(path,a,b)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--write',action='store_true'); args=ap.parse_args()
    result=calculate(); out=HERE/'results_preview.json'
    if args.write:
        with out.open('x',encoding='utf-8') as f: json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
        status='WRITE_ONCE'
    else:
        compare(result,json.loads(out.read_text(encoding='utf-8'))); status='READ_ONLY_MATCH'
    print(json.dumps({'status':status,'joint_x':result['active_mobile_clock_and_site']['x_joint_clock_minusY'],
                      'massless_x':result['massless_site_and_control']['x_position'],
                      'errors':result['error_parts']},ensure_ascii=False,indent=2))
