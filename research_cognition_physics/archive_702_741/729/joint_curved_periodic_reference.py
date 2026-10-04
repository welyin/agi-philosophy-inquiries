"""729: original periodic links and conformal half-density fermion reference."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_matter_ground_source as old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_curved_periodic_reference_results.json'
QUARK=[0,1,2,3,12,13,18,19]
LEPTON=list(range(24,32))


def background(N):
    q=old.matter.original.shared_source(N)
    f=old.old.wall.old.fields(N)
    assert np.max(abs(q['phi']-f['q']['phi']))<1e-13
    W=old.matter.original.lattice.su2(q['a'],q['eps'])
    z=np.exp(1j*q['eps']*q['a0'])
    return q,f['psi'],W,z


def matrices(N,ids,gamma=0.,curved=True,source=None):
    q,psi,W,z=background(N);count=N**3;m=len(ids)
    c=1/psi.ravel() if curved else np.ones(count)
    if source is None:source=np.ones(count)
    c=c*np.exp(-gamma*np.asarray(source)/2)
    h=np.zeros((count*m,count*m),complex);d=np.zeros_like(h);K=np.zeros_like(h);G=np.zeros_like(h)
    ix=np.ix_(ids,ids)
    alphas=[a[:32,:32][ix] for a in old.chiral.kinetic_matrices()]
    for index in np.ndindex((N,N,N)):
        v=np.ravel_multi_index(index,(N,N,N));sl=slice(v*m,(v+1)*m)
        hh,dd=old.matter.mass_matrices(q['phi'][index]);h[sl,sl]=hh[ix];d[sl,sl]=dd[ix]
        for a in range(3):
            nxt=list(index);nxt[a]=(nxt[a]+1)%N;nxt=tuple(nxt)
            w=np.ravel_multi_index(nxt,(N,N,N));sr=slice(w*m,(w+1)*m)
            R=old.matter.representation(np.eye(3),W[index+(a,)],z[index+(a,)])[ix]
            edge=-1j*c[v]*c[w]*alphas[a]@R/(2*q['eps'])
            K[sl,sr]+=edge;K[sr,sl]+=edge.conj().T
            derivative=-(source[v]+source[w])*edge/2
            G[sl,sr]+=derivative;G[sr,sl]+=derivative.conj().T
    return h+K,d,K,G,psi


def certificate(B):
    e,v=np.linalg.eigh(B)
    X=(v/e)@v.conj().T
    # Every value in B is the constructed binary64 coefficient, not an exact
    # solution of the continuum elliptic equation for psi.
    residual=float(np.linalg.norm(np.eye(len(B))-B@X,'fro'))
    bnorm=float(np.linalg.norm(B,'fro'));xnorm=float(np.linalg.norm(X,'fro'))
    # Conservative standard scalar-operation roundoff envelope for complex
    # dot products (no under/overflow); not interval certification of B's input.
    eps=np.finfo(float).eps;n=len(B);u=32*n*eps
    outward=1+32*n*eps
    bupper=np.nextafter(bnorm*outward,np.inf)
    xupper=np.nextafter(xnorm*outward,np.inf)
    envelope=u/(1-u)*bupper*xupper+32*n*eps*(np.sqrt(n)+bupper*xupper)
    r=np.nextafter(residual*outward+envelope,np.inf)
    bound=np.nextafter((1-r)/xupper,0.)
    assert r<.01 and bound>0
    return e,v,dict(dimension=n,observed_gap=float(min(abs(e))),
                   inverse_residual_Frobenius=residual,roundoff_envelope=envelope,
                   inverse_Frobenius_norm=xnorm,inverse_norm_upper=xupper,
                   represented_matrix_gap_lower=bound,
                   coefficient_input_and_continuum_errors_not_included=True)


def probe(N):
    results={};total=0.;source_total=0.
    for name,ids,weight in [('one_color_quark',QUARK,3.),('complete_lepton',LEPTON,.5)]:
        h,d,K,G,psi=matrices(N,ids)
        B=h if name=='one_color_quark' else old.bdg(h,d)
        source=G if name=='one_color_quark' else old.bdg(G,np.zeros_like(G))
        e,v,cert=certificate(B);negative=e<0
        P=v[:,negative]@v[:,negative].conj().T
        energy=weight*float(e[negative].sum());mean=weight*float(np.trace(P@source).real)
        total+=energy;source_total+=mean
        cert.update(negative_levels=int(negative.sum()),weighted_energy=energy,
                    weighted_uniform_source=mean,color_or_Nambu_weight=weight)
        results[name]=cert
    return dict(N=N,physical_modes=32*N**3,psi_range=[float(psi.min()),float(psi.max())],
                complete_generation_ground_energy=total,uniform_scale_source=source_total,
                sectors=results)


def geometric_identity():
    rng=np.random.default_rng(729)
    A=old.chiral.kinetic_matrices();errors=[]
    for _ in range(8):
        psi=float(rng.uniform(.7,1.4));grad=rng.normal(size=3)
        chi=rng.normal(size=64)+1j*rng.normal(size=64)
        dchi=rng.normal(size=(3,64))+1j*rng.normal(size=(3,64))
        gradf=2*grad/psi
        omega=[sum((A[i]@A[j]-A[j]@A[i])*gradf[j]/4 for j in range(3)) for i in range(3)]
        contraction=sum(A[i]@omega[i] for i in range(3))
        errors.append(float(np.max(abs(contraction-sum(A[i]*gradf[i] for i in range(3))))))
        covariant=-1j*psi*sum(A[i]@(psi**-3*dchi[i]-3*psi**-4*grad[i]*chi
                                    +omega[i]@(psi**-3*chi)) for i in range(3))
        half=-1j/psi**2*sum(A[i]@dchi[i] for i in range(3)) \
              +1j/psi**3*sum(A[i]*grad[i] for i in range(3))@chi
        errors.append(float(np.max(abs(covariant-half))))
    # All three color blocks are retained, and are identical in this background.
    q,_,W,z=background(4);idx=(0,1,2)
    h,d=old.matter.mass_matrices(q['phi'][idx])
    matrices_to_check=[h,d]+[a[:32,:32] for a in A]+[
        old.matter.representation(np.eye(3),W[idx+(i,)],z[idx+(i,)]) for i in range(3)]
    sectors=[[4*c,4*c+1,4*c+2,4*c+3,12+2*c,13+2*c,18+2*c,19+2*c]
             for c in range(3)]+[LEPTON]
    for B in matrices_to_check:
        block=np.zeros_like(B)
        for ids in sectors:block[np.ix_(ids,ids)]=B[np.ix_(ids,ids)]
        errors.append(float(np.max(abs(B-block))))
        for ids in sectors[1:3]:errors.append(float(np.max(abs(B[np.ix_(QUARK,QUARK)]-B[np.ix_(ids,ids)]))))
    assert max(errors)<2e-12
    return dict(max_spin_connection_half_density_and_sector_error=max(errors),
        random_local_jets=8,all_original_physical_modes_retained=True,
        given_metric='g_ij=psi^4 delta_ij; sqrt(g)=psi^6; frozen lapse=1, shift=0',
        derived_canonical_kinetic='psi^-1 D_flat,A psi^-1',
        not_a_new_derivation_of_dimension_or_metric=True)


def ground_values(N,gamma=0.,source=None,curved=True):
    energy=0.;mean=0.;errors=[]
    for name,ids,weight in [('quark',QUARK,3.),('lepton',LEPTON,.5)]:
        h,d,K,G,_=matrices(N,ids,gamma=gamma,source=source,curved=curved)
        B=h if name=='quark' else old.bdg(h,d)
        BG=G if name=='quark' else old.bdg(G,np.zeros_like(G))
        e,v=np.linalg.eigh(B);negative=e<0
        energy+=weight*float(e[negative].sum())
        mean+=weight*float(np.sum(v[:,negative].conj()*(BG@v[:,negative])).real)
        r=np.ones(N**3) if source is None else source
        rr=np.repeat(r,len(ids))
        anti=-(rr[:,None]*K+K*rr[None,:])/2
        errors.extend([float(np.max(abs(G-anti))),float(np.max(abs(B-B.conj().T)))])
    return energy,mean,max(errors)


def source_identity():
    N=4;grid=np.indices((N,N,N))*2*np.pi/N
    local=(.7+.2*np.cos(grid[0])+.1*np.sin(grid[1]+grid[2])).ravel()
    step=2e-5;out={}
    for name,r in [('uniform',np.ones(N**3)),('local',local)]:
        energy,mean,anti_error=ground_values(N,source=r)
        plus=ground_values(N,step,r)[0];minus=ground_values(N,-step,r)[0]
        fd=(plus-minus)/(2*step);err=abs(fd-mean)
        h,d,K,G,_=matrices(N,LEPTON,source=r)
        Kplus=matrices(N,LEPTON,gamma=step,source=r)[2]
        Kminus=matrices(N,LEPTON,gamma=-step,source=r)[2]
        matrix_error=float(np.max(abs((Kplus-Kminus)/(2*step)-G)))
        assert anti_error<1e-12 and matrix_error<1e-8 and err<1e-5
        out[name]=dict(ground_energy=energy,source_mean=mean,energy_finite_difference=fd,
                      Hellmann_Feynman_error=err,edge_derivative_error=matrix_error,
                      anticommutator_and_Hermitian_error=anti_error)
    flat_energy,flat_source,_=ground_values(N,curved=False)
    out['comparison_only_flat_spatial_metric']=dict(ground_energy=flat_energy,source_mean=flat_source,
        original_internal_links_and_mass_retained=True,
        energy_change_from_original_geometry=flat_energy-out['uniform']['ground_energy'])
    assert abs(flat_energy-out['uniform']['ground_energy'])>1.
    return out


def run():
    results=dict(conformal_geometry_and_full_species=geometric_identity(),
                 original_periodic_reference=[probe(4),probe(6)],
                 same_ground_geometric_sources=source_identity())
    deps=('research_note_573.md','research_note_574.md','research_note_598.md',
          'research_note_604.md','research_note_663.md','research_note_667.md',
          'research_note_699.md','research_note_727.md','research_note_728.md',
          'joint_gravity_material_coordinates.py','joint_curved_quantum_source.py',
          'joint_matter_ground_source.py','joint_ground_gauss_lift.py',
          'round729_drafts/nonflat_reference_entry.py')
    return dict(round=729,tests_run=3,failures=0,errors=0,results=results,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Conformal half-density map derived for the inherited fixed 3D metric, connection, lapse=1 and shift=0. Original nonflat weak/circle links, all32 matter modes/site, positive numerical psi and full periodic graphs N=4,6 share an isolated reference and its geometric derivative. Inverse residual bounds apply to represented binary64 matrices under the stated arithmetic model, not certified exact continuum coefficients. No uniform continuum gap, one-generation chiral regulator, actual adiabatic dynamics, self-consistent Einstein source or generated GR claimed.604 doubling and699 candidate-specific obstruction remain.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results'],ensure_ascii=False,indent=2))
