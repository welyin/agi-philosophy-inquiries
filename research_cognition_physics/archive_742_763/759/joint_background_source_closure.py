"""759: original-matter source covariance versus a definite-background limit.

The no-go statement concerns a stated uniform finite-time moment contract,
not every effective theory. Code checks trace combinatorics, original CAR
coefficients and the exact noncommuting differential identity used in proof.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import numpy as np
import joint_record_mass_feedback as old
import joint_reference_constraint_strata as background
import joint_chiral_source_matching as chiral

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_background_source_closure_results.json'

def residue_count(n,r):
    return sum(comb(n,k) for k in range(r,n+1,3))

def trace_square(h,d,nq):
    q=Q(residue_count(nq-2,2),residue_count(nq,0))
    quark=float(q)*float(np.trace(h[:nq,:nq]@h[:nq,:nq]).real)
    lepton=.25*float(np.trace(h[nq:,nq:]@h[nq:,nq:]).real)
    pair=.25*float(np.vdot(d,d).real)
    return quark+lepton+pair,dict(quark=quark,lepton_Dirac=lepton,Majorana=pair,
                                 quark_pair_probability_exact=str(q))

def finite_trace_check():
    # Original one-node coefficients restricted solely to check the CAR trace
    # identity: three colors, one up spin, and all four neutral lepton modes.
    x=np.array([0.,.51,0.,0.,.39]);h,d=old.mass_x(x)
    ids=[0,4,8,12,14,16,24,25,30,31];ix=np.ix_(ids,ids);h,d=h[ix],d[ix]
    nq=6;states=[b for b in range(1<<len(ids)) if (b&((1<<nq)-1)).bit_count()%3==0]
    values=[];means=[]
    for b in states:
        image=old.car.quadratic({b:1.+0j},h,d)
        values.append(float(old.dot(image,image).real));means.append(image.get(b,0.).real)
        assert all((a&((1<<nq)-1)).bit_count()%3==0 for a in image)
    expected,parts=trace_square(h,d,nq);actual=float(np.mean(values))
    assert len(states)==residue_count(nq,0)*16==352
    assert abs(actual-expected)<2e-14 and max(abs(np.array(means)))==0.
    assert parts['Majorana']>0
    return dict(conditional_Fock_states=len(states),independent_enumerated_square=actual,
                analytic_square=expected,error=abs(actual-expected),parts=parts,
                restricted_coefficient_calibration_not_the_full_model=True)

def original_source_check():
    rng=np.random.default_rng(759)
    src,psi,*_=background.completed(12,1.)
    indices=((0,3,2),(1,3,2),(0,4,2))
    phi=np.array([src['phi'][i] for i in indices])
    x=phi/np.sqrt(old.geom.F(phi))[:,None]
    n=96;h=np.zeros((n,n),complex);d=np.zeros_like(h);hop=np.zeros_like(h)
    for v in range(3):
        sl=slice(32*v,32*v+32);h[sl,sl],d[sl,sl]=old.mass_x(x[v])
    eps=2*np.pi/12
    for v,w,axis in ((0,1,0),(1,2,1),(2,0,2)):
        R=old.matter.representation(
            old.matter.gauge.group_exp(rng.normal(size=8)*.2,3),
            old.matter.gauge.group_exp(rng.normal(size=3)*.2,2),np.exp(.11j))
        alpha=chiral.kinetic_matrices()[axis][:32,:32]
        # Original Weyl direction and conformal half-density endpoint weights.
        edge=-1j*alpha@R/(2*eps*float(psi[indices[v]]*psi[indices[w]]))
        a=slice(32*v,32*v+32);b=slice(32*w,32*w+32)
        hop[a,b]=edge;hop[b,a]=edge.conj().T
    quarks=np.array([i for i in range(n) if i%32<24])
    leptons=np.array([i for i in range(n) if i%32>=24]);order=np.r_[quarks,leptons]
    ix=np.ix_(order,order)
    square,parts=trace_square(h[ix],d[ix],len(quarks))
    p10=Q(residue_count(len(quarks)-2,2),residue_count(len(quarks),0))
    cross=float(p10)*np.trace(h[np.ix_(quarks,quarks)]@hop[np.ix_(quarks,quarks)])
    cross+=.25*np.trace(h[np.ix_(leptons,leptons)]@hop[np.ix_(leptons,leptons)])
    assert abs(cross)<1e-14
    lower=.5*abs(old.matter.Y['s'])**2*float(np.sum(x[:,4]**2))
    assert abs(lower-parts['Majorana'])<1e-14 and square>lower>0
    # Independent occupation columns: mass and hopping cannot return the same
    # occupation change. Test all32 species on this three-node coefficient fixture.
    max_cross=0.;samples=0
    while samples<18:
        bits=rng.integers(0,2,n)
        if int(bits[quarks].sum())%3:continue
        key=sum(int(b)<<i for i,b in enumerate(bits))
        m=old.car.quadratic({key:1.+0j},h,d)
        j=old.car.quadratic({key:1.+0j},hop,np.zeros_like(d))
        err=abs(old.dot(m,j));max_cross=max(max_cross,float(err));samples+=1
        assert all(sum((a>>int(i))&1 for i in quarks)%3==0 for a in set(m)|set(j))
    assert max_cross<1e-13
    # Vary all five mass coordinates, not only the Majorana coefficient.
    derivative_error=0.
    for r in (.8,1.,1.3):
        ds=1e-5
        def energy_square(a):
            hh=a*h+hop
            return trace_square(hh[ix],(a*d)[ix],len(quarks))[0]
        derivative=(energy_square(r+ds)-energy_square(r-ds))/(2*ds)
        derivative_error=max(derivative_error,abs(derivative-2*r*square))
    assert derivative_error<2e-8
    return dict(full_CAR_modes=n,original_Weyl_edges=3,mass_hopping_trace=float(cross.real),
                normalized_covariance_B_mass=square,strict_Majorana_lower=lower,parts=parts,
                independent_occupation_columns=samples,max_column_cross=max_cross,
                radial_spectral_second_moment_derivative_error=derivative_error,
                center_sector_dimension_exact=str(residue_count(len(quarks),0)*(1<<len(leptons))),
                complete_full_graph_statement_analytic=True,
                fixture_not_a_three_dimensional_periodic_dynamics_simulation=True)

# Exact rational matrix differential operators: keys are (derivative order,
# polynomial degree). Composition differentiates the right coefficient.
I=np.array([[Q(1),Q(0)],[Q(0),Q(1)]],dtype=object)
ZERO=I*0
def tidy(op):return {k:v for k,v in op.items() if any(v.ravel())}
def add(a,b,s=Q(1)):
    out={k:v.copy() for k,v in a.items()}
    for k,v in b.items():out[k]=out.get(k,ZERO)+s*v
    return tidy(out)
def scaled(a,s):return tidy({k:s*v for k,v in a.items()})
def compose(a,b):
    out={}
    for (i,p),A in a.items():
        for (j,q),B in b.items():
            for k in range(min(i,q)+1):
                derivative=1
                for r in range(k):derivative*=q-r
                key=(i+j-k,p+q-k)
                out[key]=out.get(key,ZERO)+comb(i,k)*derivative*(A@B)
    return tidy(out)
def comm(a,b):return add(compose(a,b),compose(b,a),Q(-1))
def sym(a,b):return scaled(add(compose(a,b),compose(b,a)),Q(1,2))

def differential_identity_check():
    M=np.array([[Q(1),Q(2)],[Q(2),Q(-1)]],dtype=object)
    N=np.array([[Q(0),Q(1)],[Q(1),Q(0)]],dtype=object)
    D={(1,0):I};B={(0,1):M,(0,0):N};Bp={(0,0):M}
    assert any(comm(B,Bp))
    rows=[]
    for hb in (Q(1,7),Q(2,5)):
        # Original H5 x5-only Laplacian and a scalar confining calibration.
        Hb={(2,0):-hb*hb*I/2,(2,2):-hb*hb*I/12,
            (1,1):-5*hb*hb*I/12,(0,4):I*Q(3,8),(0,2):I*Q(1,9)}
        H=add(Hb,B);A=sym(D,B)
        lhs=comm(H,A)
        rhs=add(add(sym(comm(Hb,D),B),sym(D,comm(Hb,B))),sym(Bp,B),Q(-1))
        residual=add(lhs,rhs,Q(-1));assert not residual
        assert not add(comm(B,D),Bp)
        assert not add(comm(H,D),add(comm(Hb,D),Bp,Q(-1)),Q(-1))
        # Replacing the symmetric force product by B'B would be wrong.
        wrong=add(rhs,sym(Bp,B));wrong=add(wrong,compose(Bp,B),Q(-1))
        assert add(lhs,wrong,Q(-1))
        rows.append(dict(hbar_exact=str(hb),identity_exact=True,wrong_unsymmetrized_force_rejected=True))
    return dict(rows=rows,noncommuting_mass_and_force=True,
                exact_commutator_moment_identity=True,
                differential_algebra_calibration_not_full_graph_evolution=True)

def run():
    a=finite_trace_check();b=original_source_check();c=differential_identity_check()
    deps=('research_note_350.md','research_note_525.md','research_note_574.md','research_note_643.md',
          'research_note_728.md','research_note_729.md','research_note_730.md','research_note_753.md',
          'research_note_757.md','research_note_758.md','joint_record_mass_feedback.py',
          'joint_fermion_gauss_completion.py','joint_reference_constraint_strata.py')
    return dict(round=759,tests_run=3,failures=0,errors=0,conditional_trace=a,
                original_mass_hopping_source=b,exact_moment_identity=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='A no-go for a precisely specified single deterministic material/gauge phase-orbit limit that retains a finite joint kinetic/source moment menu uniformly on a fixed physical-time interval, with continuous limiting matrix data at the initial instant. The original admissible maximally mixed stabilizer fiber gives a positive mass-source covariance contradicting the necessary moment identity. This is not a no-go for multibranch, noisy, memory-bearing or fully quantum backgrounds, finite-tolerance effective theories, or general relativity. The external spatial metric is not quantized here; no full dynamic replacement or continuum theorem is proved.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
