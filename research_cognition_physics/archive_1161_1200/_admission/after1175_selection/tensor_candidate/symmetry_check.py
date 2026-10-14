"""Read-only by default: original 48-code parity and passive menu identity.

Checks algebraic intertwiners, not a second integration or a numerical error bound.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import preview as p

HERE=Path(__file__).resolve().parent

def compute():
    data=json.loads((HERE/'results_preview.json').read_text(encoding='utf-8'))
    ud=p.BETA@np.kron(p.I2,-1j*(p.SX-p.SY)/np.sqrt(2))
    permutation=[0,2,1,3]
    trans=np.zeros((48,48),complex)
    for mode in range(4):
        for mass in range(3):
            src=12*mode+4*mass; dst=12*permutation[mode]+4*mass
            trans[dst:dst+4,src:src+4]=ud
    errs={}
    errs['UD_unitary']=float(np.linalg.norm(ud.conj().T@ud-np.eye(4),2))
    errs['beta']=float(np.linalg.norm(ud@p.BETA@ud.conj().T-p.BETA,2))
    errs['alpha_permutation']=float(max(np.linalg.norm(ud@p.ALPHA[b]@ud.conj().T-p.ALPHA[[1,0,2][b]],2) for b in range(3)))
    errs['full48_unitary']=float(np.linalg.norm(trans.conj().T@trans-np.eye(48),2))
    herrors=[]
    for line in data['background']['active']['history']:
        n=line['N']; v=np.array(line['v'])
        hs=[]
        for vel in (v,v[[1,0,2]]):
            hh=np.zeros((48,48),complex)
            for mode in range(4):
                for mass in range(3):
                    block=n*p.MP[mass]*p.BETA+sum(vel[b]*p.LMOM[mode,b]*p.ALPHA[b] for b in range(3))
                    start=12*mode+4*mass; hh[start:start+4,start:start+4]=block
            hs.append(hh)
        herrors.append(np.linalg.norm(hs[1]@trans-trans@hs[0],2))
    errs['history_generator_intertwiner']=float(max(herrors))
    eye=np.broadcast_to(np.eye(4,dtype=complex),(4,3,4,4))
    terminal_errors=[]
    for b in range(3):
        for c in (p.CLOCK,p.MASSLESS,np.eye(3)):
            eb=p.effect_matrix(eye,b,c)
            ea=p.effect_matrix(eye,[1,0,2][b],c)
            terminal_errors.append(np.linalg.norm(trans.conj().T@eb@trans-ea,2))
    errs['all_terminal_menu_transport']=float(max(terminal_errors))
    assert max(errs.values())<2e-13,errs
    return {'status':'ALGEBRAIC_PARITY_AND_COMMON_RECORD_PASS',
        'domain':'full original four-positive-axis modes x Dirac4 x mass3, tensor arbitrary R',
        'not_claimed':['new propagation integration','all instruments constructed','numerical certificate'],
        'errors':errs,'reflection':'swap x,y; z fixed; UD=beta kron(I,u), u=-i(sigma_x-sigma_y)/sqrt2',
        'consequence':'If input and all terminal records are transported together, probabilities agree. Source-only change with fixed original task is a different experiment.',
        'bindings':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['preview.py','results_preview.json']}}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--write',action='store_true'); args=ap.parse_args()
    out=compute(); dest=HERE/'results_symmetry.json'
    if args.write:
        with dest.open('x',encoding='utf-8') as f: json.dump(out,f,indent=2,ensure_ascii=False)
        mode='WRITE_ONCE'
    else:
        p.compare(out,json.loads(dest.read_text(encoding='utf-8'))); mode='READ_ONLY_MATCH'
    print(json.dumps({'mode':mode,**out},indent=2,ensure_ascii=False))
