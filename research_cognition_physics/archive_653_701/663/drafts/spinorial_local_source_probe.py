"""663 entry: local lapse composition on original604 continuum symbol.

Does not reinstate the excluded naive lattice regulator, quantize gravity,
identify the full overlap state, or represent a physical many-body wavefunction.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_chiral_source_matching as original


def run():
    g=original.kinetic_matrices();b=original.mass();size=24
    axis=np.arange(size)*2*np.pi/size;x,y=np.meshgrid(axis,axis,indexing='ij')
    n=1+.2*np.sin(x);m=1+.3*np.sin(y)
    nx=.2*np.cos(x);my=.3*np.cos(y)
    zero=np.zeros_like(x)
    frequencies=np.fft.fftfreq(size,d=1/size)
    def derivative(v,axis):
        shape=[1,1,1];shape[axis]=size
        return np.fft.ifft(1j*frequencies.reshape(shape)*np.fft.fft(v,axis=axis),axis=axis)
    def act(mat,v):return np.einsum('ij,xyj->xyi',mat,v,optimize=True)
    def length(v):return float(np.sqrt(np.mean(np.sum(abs(v)**2,axis=-1))))
    rng=np.random.default_rng(663)
    amps=rng.normal(size=(3,64))+1j*rng.normal(size=(3,64))
    psi=np.exp(1j*x)[...,None]*amps[0]+np.exp(-1j*y)[...,None]*amps[1]+np.exp(1j*(x-y))[...,None]*amps[2]
    psi/=length(psi)
    def kinetic(v):return -1j*(act(g[0],derivative(v,0))+act(g[1],derivative(v,1)))
    def local(v,f,fx,fy,mass):
        return f[...,None]*(kinetic(v)+act(mass,v))-.5j*(fx[...,None]*act(g[0],v)+fy[...,None]*act(g[1],v))
    vx=-m*nx;vy=n*my
    div=.2*m*np.sin(x)-.3*n*np.sin(y)
    scalar=1j*(vx[...,None]*derivative(psi,0)+vy[...,None]*derivative(psi,1)+.5*div[...,None]*psi)
    comm=g[0]@g[1]-g[1]@g[0]
    spin=.25j*(nx*my)[...,None]*act(comm,psi)
    rows=[]
    for scale in (0.,1.,1.7):
        bm=scale*b
        hn=lambda v:local(v,n,nx,zero,bm)
        hm=lambda v:local(v,m,zero,my,bm)
        actual=(hn(hm(psi))-hm(hn(psi)))/1j
        error=length(actual-scalar-spin)
        assert error<2e-12
        rows.append(dict(mass_scale=scale,complete_local_bracket_error=error,
            scalar_shift_only_error=length(actual-scalar)))
    spin_norm=length(spin);assert spin_norm>.01
    q=.25j*comm
    noncentral=float(np.linalg.norm(q@g[0]-g[0]@q,2));assert noncentral>.5
    return dict(date='2026-10-02',status='663 actual probe, not completed round',
        inherited_model='604 original32CAR/64Nambu masses and continuum principal symbol',
        periodic_grid=[size,size],third_spatial_derivative_zero=True,
        positive_lapses='N=1+0.2sin(x), M=1+0.3sin(y)',
        complete_formula='[h[N],h[M]]/i = i(v.derivative+div(v)/2)+i[Gamma.gradN,Gamma.gradM]/4',
        vector_field='v=N grad(M)-M grad(N)',spin_term_norm=spin_norm,
        spin_generator_commutator_with_Clifford_norm=noncentral,
        rows=rows,Nambu_coefficient_not_many_body_state=True,
        no_naive_lattice_regulator_restoration=True,no_ADM_or_Kosmann_equivalence_yet=True,
        dependencies={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ('joint_chiral_source_matching.py','research_note_375.md','research_note_604.md',
             'research_note_612.md','research_note_662.md')})


if __name__=='__main__':
    result=run()
    with (HERE/'spinorial_local_source_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))
