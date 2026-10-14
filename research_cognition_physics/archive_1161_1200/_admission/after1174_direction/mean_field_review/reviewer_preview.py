"""Independent Taylor-jet audit from candidate.md, without author code imports.

Diagnostic floating-point reproduction, NOT a rigorous integration certificate.
Default is read-only; --write creates this directory's new JSON exclusively.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.]).astype(complex)
BETA = np.kron(SZ, np.eye(2))
ALPHA = np.array([np.kron(SX, a) for a in (SX, SY, SZ)])
PCLK = (np.eye(2)-SY)/2
MODES = np.array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]], int)
MASSES = np.array([.3, .5])
CHI = np.array([1.,0.,0.,1.], complex)/np.sqrt(2)


def shell_data(n):
    """Enumerate all signed triples, then group by sorted abs integer indices."""
    counts = {}
    signed = range(-(n//2), n//2+1)
    for xyz in itertools.product(signed, repeat=3):
        key = tuple(sorted(abs(x) for x in xyz))
        counts[key] = counts.get(key, 0)+1
    keys = np.array(sorted(counts), int)
    ell = 16/n
    s = np.sqrt(np.sum(np.sin(2*np.pi*keys/n)**2, axis=1))/ell
    weights = np.array([counts[tuple(k)] for k in keys], float)/n**3
    assert int(round(weights.sum()*n**3)) == n**3
    return s, weights


def taylor_run(n=65, nu=.35, order=12, steps=64, full_source=False):
    """Integrate source+G+full 4x4 probe blocks by quadratic Taylor recurrences."""
    eps, mass, omega, T = .15, 2., .5, 2.
    dt = T/steps
    s, weights = shell_data(n)
    ell = 16/n
    pv = 2*np.pi*MODES/16
    kp = np.array([np.sin(ell*pv)/ell, pv])
    ap = np.einsum('cmb,bij->cmij', kp, ALPHA)
    # Coordinate order: discretization, mode, flavor, Dirac row/column.
    uh = np.broadcast_to(np.eye(4, dtype=complex), (2,4,2,4,4)).copy()
    g = np.array([1.,0.,0.])
    ints = np.zeros(2)
    r = np.zeros((len(s), 3)); r[:,2] = 1
    if full_source:
        grid = np.array(list(itertools.product(range(-(n//2),n//2+1), repeat=3)))
        p = 2*np.pi*grid/16
        aks = np.einsum('mb,bij->mij', np.sin(ell*p)/ell, ALPHA)
        rho = np.broadcast_to((np.eye(4)+BETA)/2, (n**3,4,4)).copy()
    max_energy_error = 0.
    max_g_norm_error = 0.
    initial_energy = 4.-omega
    for _step in range(steps):
        gc = np.zeros((order+1,3)); gc[0]=g
        nc = np.zeros(order+1); vc = np.zeros(order+1)
        bc = np.zeros(order+1)
        uc = np.zeros((order+1,)+uh.shape, complex); uc[0]=uh
        if full_source:
            rc = np.zeros((order+1,)+rho.shape, complex); rc[0]=rho
        else:
            rc = np.zeros((order+1,)+r.shape); rc[0]=r
        for k in range(order):
            nc[k] = (1. if k==0 else 0.)+eps*gc[k,2]
            vc[k] = (1. if k==0 else 0.)+nu*gc[k,2]
            if full_source:
                em = mass*np.einsum('ij,mji->', BETA, rc[k]).real/n**3
                et = np.einsum('mij,mji->', aks, rc[k]).real/n**3
            else:
                em = 2*mass*np.dot(weights,rc[k,:,2])
                et = 2*np.dot(weights*s,rc[k,:,0])
            bc[k] = eps*em+nu*et
            gd = np.zeros(3)
            for l in range(k+1):
                u = k-l
                ox = -omega if l==0 else 0.
                gd += 2*np.array([-bc[l]*gc[u,1],
                                  bc[l]*gc[u,0]-ox*gc[u,2],
                                  ox*gc[u,1]])
                hs = nc[l]*mass*BETA + vc[l]*aks if full_source else None
                if full_source:
                    rc[k+1] += -1j*(hs@rc[u]-rc[u]@hs)/(k+1)
                else:
                    hx, hz = vc[l]*s, nc[l]*mass
                    rr=rc[u]
                    rc[k+1,:,0] += -2*hz*rr[:,1]/(k+1)
                    rc[k+1,:,1] += 2*(hz*rr[:,0]-hx*rr[:,2])/(k+1)
                    rc[k+1,:,2] += 2*hx*rr[:,1]/(k+1)
                hp = (nc[l]*MASSES[None,None,:,None,None]*BETA
                      +vc[l]*ap[:,:,None,:,:])
                uc[k+1] += -1j*(hp@uc[u])/(k+1)
            gc[k+1] = gd/(k+1)
        power = dt**np.arange(order+1)
        g = np.einsum('k,ki->i',power,gc)
        uh = np.einsum('k,kcmdij->cmdij',power,uc)
        ints += np.array([np.dot(nc[:-1], dt**np.arange(1,order+1)/np.arange(1,order+1)),
                          np.dot(vc[:-1], dt**np.arange(1,order+1)/np.arange(1,order+1))])
        if full_source:
            rho=np.einsum('k,kmij->mij',power,rc)
            em=mass*np.einsum('ij,mji->',BETA,rho).real/n**3
            et=np.einsum('mij,mji->',aks,rho).real/n**3
        else:
            r=np.einsum('k,kmi->mi',power,rc)
            em=2*mass*np.dot(weights,r[:,2]); et=2*np.dot(weights*s,r[:,0])
        energy=(1+eps*g[2])*em+(1+nu*g[2])*et-omega*g[0]
        max_energy_error=max(max_energy_error,abs(energy-initial_energy))
        max_g_norm_error=max(max_g_norm_error,abs(np.dot(g,g)-1))
    unitary=np.max(np.abs(np.swapaxes(uh.conj(),-1,-2)@uh-np.eye(4)))
    return {'g':g,'tau':ints[0],'xi':ints[1],'a':(1+eps*g[2])/(1+nu*g[2]),
            'U':uh,'source_classes':len(s),'max_energy_error':max_energy_error,
            'max_g_norm_error':max_g_norm_error,'max_unitary_error':unitary,
            'source_length_error':float(np.max(np.abs(np.sum(r*r,axis=1)-1))) if not full_source else
                float(np.max(np.abs(rho@rho-rho)))}


def flat_blocks(tau, continuum=False, massless=False):
    ell=16/65
    p=2*np.pi*MODES/16
    pp=p if continuum else np.sin(ell*p)/ell
    ap=np.einsum('mb,bij->mij',pp,ALPHA)
    masses=np.zeros(2) if massless else MASSES
    h=masses[None,:,None,None]*BETA+ap[:,None]
    energy=np.sqrt(np.sum(pp*pp,axis=1)[:,None]+masses[None,:]**2)
    sinc=np.sinc(energy*tau/np.pi)*tau
    return (np.cos(energy*tau)[...,None,None]*np.eye(4)
            -1j*sinc[...,None,None]*h)


def effects(blocks, n=65):
    """Direct site quadrature, retaining full 8-dimensional internal states."""
    nsites=n**3
    # Fixed witness propagated separately at each momentum and flavor.
    amplitudes=np.einsum('mdij,j->mid', blocks, CHI)/(2*np.sqrt(2))
    sites=np.array(list(itertools.product(range(n),repeat=3)))
    phase=np.exp(2j*np.pi*(sites@MODES.T)/n)/np.sqrt(nsites)
    psi=np.einsum('jm,mid->jid',phase,amplitudes)
    density=np.sum(np.abs(psi)**2,axis=(1,2))
    clock=np.einsum('jid,de,jie->j',psi.conj(),PCLK,psi).real
    fs=(1+np.sin(2*np.pi*sites/n))/2
    return {'position':fs.T@density,'joint':fs.T@clock,
            'clock':float(clock.sum()),'total':float(density.sum())}


def compressed_effects(blocks):
    """32x32 operators from exact Fourier coefficients; no dense site sum."""
    out=[]
    # blocks indexed momentum, flavor, Dirac; code order momentum,Dirac,flavor.
    code=np.zeros((32,32),complex)
    for m in range(4):
        for f in range(2):
            ix=np.array([m*8+d*2+f for d in range(4)])
            code[np.ix_(ix,ix)]=blocks[m,f]
    for axis in range(3):
        fourier=np.eye(4,dtype=complex)/2
        # <p|sin(x_b)|q>=(delta[q-p+e_b]-delta[q-p-e_b])/(2i).
        for p in range(4):
            for q in range(4):
                dif=MODES[q]-MODES[p]
                e=np.eye(3,dtype=int)[axis]
                if np.array_equal(dif,-e): fourier[p,q]+=-.25j
                if np.array_equal(dif,e): fourier[p,q]+=.25j
        for clock in (np.eye(2),PCLK):
            effect=np.kron(fourier,np.kron(np.eye(4),clock))
            out.append(code.conj().T@effect@code)
    return out


def encode(v):
    if isinstance(v,np.ndarray):
        if np.iscomplexobj(v): return {'real':v.real.tolist(),'imag':v.imag.tolist()}
        return v.tolist()
    if isinstance(v,np.generic): return v.item()
    if isinstance(v,dict): return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [encode(x) for x in v]
    return v


def run():
    print('Independent Taylor order12/64 and order14/80',flush=True)
    coarse=taylor_run(order=12,steps=64)
    fine=taylor_run(order=14,steps=80)
    print('Independent full small-grid Dirac density and shell comparison',flush=True)
    small=taylor_run(n=7,order=10,steps=64)
    small_full=taylor_run(n=7,order=10,steps=64,full_source=True)
    print('Lapse-only control and actual site quadrature',flush=True)
    lapse=taylor_run(n=65,nu=.15,order=12,steps=64)
    witness={}
    eigens={}
    for ix,name in enumerate(('lattice','continuum')):
        act=fine['U'][ix]
        flat=flat_blocks(fine['tau'],continuum=bool(ix))
        am=flat_blocks(fine['xi'],continuum=bool(ix),massless=True)
        fm=flat_blocks(fine['tau'],continuum=bool(ix),massless=True)
        records=[effects(u) for u in (act,flat,am,fm)]
        witness[name]={'active':records[0], 'flat_same_tau':records[1],
                       'massless_active':records[2], 'massless_flat_same_tau':records[3],
                       'position_gap':records[0]['position']-records[1]['position'],
                       'joint_gap':records[0]['joint']-records[1]['joint'],
                       'massless_gap':records[2]['position']-records[3]['position']}
        ae=compressed_effects(act); fe=compressed_effects(flat)
        eigens[name]=[float(np.max(np.abs(np.linalg.eigvalsh(a-f)))) for a,f in zip(ae,fe)]
        witness0=np.tile(np.kron(CHI,np.ones(2)/np.sqrt(2)),4)/2
        compressed=np.array([np.vdot(witness0,e@witness0).real for e in ae]).reshape(3,2)
        witness[name]['direct_site_vs_compressed_max']=float(np.max(np.abs(compressed-
             np.stack((records[0]['position'],records[0]['joint']),axis=1))))
    p=2*np.pi/16
    exact_delta=p-np.sin((16/65)*p)/(16/65)
    cubic_delta=(16/65)**2*p**3/6
    def lean(d): return {k:v for k,v in d.items() if k!='U'}
    result={'method':'independent quadratic Taylor jets; full 4x4 Dirac blocks; diagnostic floats',
            'independence':'Read candidate equations and parent-reported approximate outputs before coding. Did not read or import author preview.py.',
            'active':lean(fine),'lapse':lean(lapse),
            'refinement':{'g_max':float(np.max(np.abs(fine['g']-coarse['g']))),
                          'tau':float(abs(fine['tau']-coarse['tau'])),
                          'xi':float(abs(fine['xi']-coarse['xi'])),
                          'U_max':float(np.max(np.abs(fine['U']-coarse['U'])))},
            'small_grid_full_dirac':{'n':7,'g_max':float(np.max(np.abs(small['g']-small_full['g']))),
                                     'tau_error':float(abs(small['tau']-small_full['tau'])),
                                     'U_max':float(np.max(np.abs(small['U']-small_full['U']))),
                                     'projector_error':small_full['source_length_error']},
            'lapse_flat_max':float(np.max(np.abs(lapse['U'][0]-flat_blocks(lapse['tau'])))),
            'witness':witness,'unknown32_effect_difference_opnorm':eigens,
            'lattice_budget':{'exact_symbol_difference':exact_delta,'cubic_upper':cubic_delta,
                              'active_vector_estimate':cubic_delta*fine['xi'],
                              'flat_vector_estimate':cubic_delta*fine['tau'],
                              'sum_estimate':cubic_delta*(fine['xi']+fine['tau']),
                              'warning':'Integral values are diagnostic here; rigorous integral enclosures required.'},
            'geometry':{'addot_initial':.24,'kretschmann_initial':12*.24**2},
            'certification':'Independent numerical consistency only; no rigorous ODE interval or instrument-realization certificate.'}
    assert result['refinement']['U_max']<2e-10
    assert result['small_grid_full_dirac']['U_max']<2e-10
    assert result['lapse_flat_max']<2e-10
    assert max(witness[k]['direct_site_vs_compressed_max'] for k in witness)<2e-11
    return encode(result)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        return max((compare(a[k],b[k]) for k in a),default=0.)
    if isinstance(a,list):
        assert len(a)==len(b)
        return max((compare(x,y) for x,y in zip(a,b)),default=0.)
    if isinstance(a,(int,float)) and not isinstance(a,bool):
        err=abs(a-b)
        assert err<=2e-11*max(1,abs(a),abs(b)),(a,b)
        return err
    assert a==b
    return 0.


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); data=run(); path=ROOT/'reviewer_preview.json'
    if args.write:
        with path.open('x',encoding='utf-8') as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2);stream.write('\n')
        print('WROTE_NEW_REVIEWER_PREVIEW',flush=True)
    else:
        with path.open(encoding='utf-8') as stream: old=json.load(stream)
        print('READ_ONLY_MATCH',compare(data,old),flush=True)
    print(json.dumps({'g':data['active']['g'],'tau':data['active']['tau'],'xi':data['active']['xi'],
                      'joint_gap':data['witness']['lattice']['joint_gap'],
                      'massless_gap':data['witness']['lattice']['massless_gap'],
                      'refinement':data['refinement']},ensure_ascii=False),flush=True)
