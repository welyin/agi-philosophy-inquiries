"""Explicit same-parent mass-label extension of the frozen first preview.
No source, lattice, time or witness parameter changes. --write is create-once.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, itertools, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('candidate_preview',HERE/'preview.py')
p=importlib.util.module_from_spec(spec); spec.loader.exec_module(p)

def compute():
    original=json.loads((HERE/'results_preview.json').read_text(encoding='utf-8'))
    run=p.evolve(p.NU)
    u0=p.exact_const(0,run['xi'],np.zeros(1))
    u=np.concatenate([u0,run['u']],axis=1)
    uflat=p.exact_const(run['tau'],run['tau'],np.array([0,.3,.5]))
    big=np.zeros((48,48),complex)
    for mode,m in itertools.product(range(4),range(3)):
        i=12*mode+4*m; big[i:i+4,i:i+4]=u[mode,m]
    gram=float(np.linalg.norm(big.conj().T@big-np.eye(48),2))
    wmass=np.zeros((4,3,4),complex)
    wmass[:,1:,:]=p.WITNESS.reshape(4,2,4)
    wzero=np.zeros((4,3,4),complex);wzero[:,0,:]=p.CHI/2
    clock=np.zeros((3,3),complex);clock[1:,1:]=p.CLOCK
    max_reuse=0.0;actual={}
    for axis,name in enumerate('xyz'):
        for label,c,witness,oldname in [('massive_joint',clock,wmass,'active_mobile_clock_and_site'),
                                        ('massive_position',np.eye(3),wmass,'active_mobile_clock_and_site'),
                                        ('massless_position',np.eye(3),wzero,'massless_site_and_control')]:
            mats=[]
            for ut in [u,uflat]:
                e=np.zeros((48,48),complex);f=p.site_effect(axis)
                for j,k,m,n in itertools.product(range(4),range(4),range(3),range(3)):
                    e[12*j+4*m:12*j+4*m+4,12*k+4*n:12*k+4*n+4]=f[j,k]*c[m,n]*(ut[j,m].conj().T@ut[k,n])
                mats.append(e)
            psi=witness.ravel()
            vals=[float(np.vdot(psi,e@psi).real) for e in mats]
            key=name+('_joint_clock_minusY' if label=='massive_joint' else '_position')
            old=original[oldname][key]
            err=max(abs(vals[0]-old['active']),abs(vals[1]-old['proper_time_flat']))
            max_reuse=max(max_reuse,err)
            actual[name+'_'+label]={'active':vals[0],'flat':vals[1],'gap':vals[0]-vals[1],'old_menu_reuse_error':err}
    # Cross-mass coherence is retained in an explicitly non-block-diagonal state.
    coherent=np.arange(1,49,dtype=float)+1j*np.arange(48,0,-1,dtype=float)
    coherent/=np.linalg.norm(coherent)
    normerr=float(abs(np.linalg.norm(big@coherent)-1))
    bP=2*p.EPS*(.3+.5)+18*p.NU/p.DX
    Lbound=p.VOL*bP
    bS=2*p.EPS*p.MS+6*p.NU/p.DX; w=2*p.VOL*bS
    logK=2*(math.log10(4*Lbound*p.TIME/1e-4)+6*w*p.TIME/math.log(10))
    assert gram<1e-10 and max_reuse<1e-12 and normerr<1e-10
    return {'kind':'UNNUMBERED_EXPLICIT_DIRECT_SUM_EXTENSION',
            'previous_preview_unchanged':True,'mass_labels':[0,.3,.5],
            'full_CAR_modes_per_site':12,'unknown_prepared_dimension':48,
            'full_48_column_gram_error':gram,'cross_mass_coherent_norm_error':normerr,
            'maximum_old_menu_reuse_error':max_reuse,'menus':actual,
            'resources':{'bP_full_Fock_bound':bP,'L_full_Fock_bound':Lbound,'w_bound':w,
                         'meanfield_target_eta':1e-4,'sufficient_log10_K':logK},
            'sources':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['preview.py','results_preview.json']},
            'certification_scope':{'numerical_error_certified':False,'parent_was_simulated':False,
                                   'cross_mass_coherences_discarded':False}}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=compute();out=HERE/'results_unified_probe.json'
    if args.write:
        with out.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False)
        status='WRITE_ONCE'
    else:
        p.compare(data,json.loads(out.read_text(encoding='utf-8')));status='READ_ONLY_MATCH'
    print(json.dumps({'status':status,'Gram':data['full_48_column_gram_error'],
                     'reuse':data['maximum_old_menu_reuse_error'],'resources':data['resources']},indent=2))
