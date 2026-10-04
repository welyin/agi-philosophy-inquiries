"""728 entry: inherited spectral filtering with the actual original CAR record."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_matter_ground_source as previous
TARGET=HERE/'macroscopic_source_entry_results.json'


def run():
    B,G=previous.star();d=previous.data(B,G);e=d['e'];v=d['v'];P=d['P']
    g=v.conj().T@G@v;de=e[:,None]-e[None,:];I=np.eye(len(e))
    # Nonselective component of the genuine sterile occupation instrument.
    # R=(-1)^n in Nambu form; this decomposition is only for moment evaluation.
    R=I.copy();R[30,30]=-1;R[128+30,128+30]=-1
    P1=R@P@R
    def moment(Q,A):
        mean=float(np.trace(Q@A).real/2)
        variance=float(np.trace(Q@A@(I-Q)@A).real/2)
        assert variance>-1e-9
        return mean,max(variance,0.)
    def mixture(A):
        m0,n0=moment(P,A);m1,n1=moment(P1,A)
        return (m0+m1)/2,(n0+n1)/2+(m0-m1)**2/4
    source_mean,source_var=moment(P,G)
    rows=[]
    for T in (1.,4.,16.,64.):
        f=np.exp(1j*de*T/2)*np.sinc(de*T/(2*np.pi))
        GT=v@(g*f)@v.conj().T
        mean,variance=moment(P,GT);pm,pv=mixture(GT)
        bound=source_var/(d['gap']**2*T*T)  # pair gaps >=2 one-particle gap
        assert abs(mean-source_mean)<1e-9 and variance<=min(source_var,bound)+1e-9
        rows.append(dict(T=T,ground_mean=mean,ground_variance=variance,
                         pair_gap_variance_bound=bound,
                         postrecord_mean=pm,postrecord_variance=pv))
    Ginf=v@(g*(abs(de)<1e-9))@v.conj().T
    Ginf2=v@(g*(abs(de)<1e-8))@v.conj().T
    unchanged=float(np.max(abs(Ginf-Ginf2)))
    im,iv=moment(P,Ginf);pm,pv=mixture(Ginf)
    em,ev=mixture(B)
    assert abs(im-source_mean)<1e-9 and iv<1e-9
    assert pv>.001 and abs(pm-source_mean)>.001 and ev>.001
    assert unchanged<1e-10
    deps=('research_note_591.md','research_note_602.md','research_note_630.md',
          'research_note_633.md','research_note_719.md','research_note_727.md',
          'joint_matter_ground_source.py')
    return dict(entry_round=728,new_formal_round=False,tests_run=2,failures=0,errors=0,
        results=dict(original_three_direction_star=True,original_ground_source_mean=source_mean,
                     original_ground_source_variance=source_var,finite_windows=rows,
                     infinite_window_ground_variance=iv,
                     infinite_window_postrecord_mean=pm,
                     infinite_window_postrecord_variance=pv,
                     postrecord_energy_increase=em-d['E'],postrecord_energy_variance=ev,
                     degeneracy_threshold_crosscheck=unchanged,
                     genuine_sterile_occupation_instrument_nonselective_moments=True),
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        scope='Entry reuse of591/602 spectral averages on727 original finite CAR star. One-sided time averages reduce stationary-ground off-band source noise. The genuine original local occupation record changes the state and leaves conserved source and energy fluctuations. An averaged operator is not an implemented instrument; no whole-Gauss interacting, dynamic-background or spatial continuum claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results'],ensure_ascii=False,indent=2))
