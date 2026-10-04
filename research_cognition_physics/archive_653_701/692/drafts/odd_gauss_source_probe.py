"""692 entry: original odd Gauss composite in the full source polynomial.

Retain every original auxiliary channel and the full unnormalized Pfaffian.
Fixed-background source identities are not an averaged positivity certificate.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
import joint_rational_physical_limit as sources
import joint_physical_source_reflection as reflection

TARGET=HERE/'odd_gauss_source_probe_results.json'
B=[((6,10,11),2.),((7,9,11),-2.),((8,9,10),2.)]


def matrices(e,phis):
    mat=list(base.fixed_matrices(e,phis))
    ev,vec=np.linalg.eigh(base.internal.spin.GAMMA[3])
    wp=vec[:,ev>.5];wm=base.internal.spin.G5@wp
    mat[0]=np.kron(np.kron(np.eye(4),(wp-wm)/np.sqrt(2)),np.eye(16))
    mat[1]=np.kron(np.kron(np.eye(4),(wp+wm)/np.sqrt(2)),np.eye(16))
    # This probe is strictly massless; no rotated mass convention is inferred.
    mat[-1]=np.zeros_like(mat[-1])
    return tuple(mat)


def source_rows(ry):
    j=base.mass.dictionary.dictionary()
    fields=[]
    for i in range(16):
        z=np.zeros(256,complex);z[64:80]=j.conj().T[i,:]
        fields.append(z)
    neutral=np.zeros(256,complex)
    neutral[64]=neutral[96]=1/np.sqrt(2)
    terms=[]
    for ia,ca in B:
        for ib,cb in B:
            left=[fields[k] for k in ia];right=[fields[k] for k in ib]
            odd=np.array([r.conj()@ry for r in reversed(left)]+right).T
            even=np.array([r.conj()@ry for r in reversed([neutral]+left)]+[neutral]+right).T
            terms.append((ca*cb,odd,even))
    return terms


def integral_coefficient(q,z):
    # Border Pfaffian convention differs from an ordered field product by (-1)^(k/2).
    return (-1)**(z.shape[1]//2)*sources.source_coeff(q,z)


def evaluate(links,e,phis,paired=False):
    u,v,d,h,gap=base.kernel(links)
    mat=matrices(e,phis)
    q=base.regular(u,v,d,mat,lam=0)
    _,ry,_=reflection.reflection_matrices(mat,np.eye(4)[[2,3,0,1]])
    terms=source_rows(ry)
    jm,jp,m,mb,_=mat
    rv=v.shape[1];r=128
    kl=jp.conj().T@d@v;w=jm.conj().T@v
    nw=np.block([[np.zeros((rv,rv)),-kl.T],[kl,np.zeros((r,r))]])
    ell=np.zeros((256,rv+r),complex)
    ell[:128,:rv]=w;ell[128:,rv:]=np.eye(128)
    sj=base.old.diag(np.column_stack((u,v)),np.column_stack((jm.conj(),jp.conj())))
    aux=base.pf(u.T@m@u)*base.pf(jm.conj().T@mb@jm.conj())/np.linalg.det(sj)
    physical=dict(N=nw,Phi=ell)
    direct=sum(c*integral_coefficient(q,z) for c,z,_ in terms)
    factored=aux*sum(c*integral_coefficient(physical,z) for c,z,_ in terms)
    denominator=max(abs(direct),abs(factored),1e-280)
    err=float(abs(direct-factored)/denominator)
    assert denominator>1e-22 and err<2e-9,(denominator,err)
    pair_value=None;pair_error=None
    if paired:
        pair_value=sum(c*integral_coefficient(q,z) for c,_,z in terms)
        pair_error=float(abs(pair_value-direct)/max(abs(direct),abs(pair_value),1e-280))
        assert pair_error<2e-9
    return dict(source=direct,factored=factored,factor_error=err,
                paired=pair_value,pair_error=pair_error,gap=gap,
                auxiliary=aux,weight=q['weight'],physical=q,terms=terms)


def run():
    links,e,phis=sources.fixture()
    # Reuse a known original nonflat full-group background, not an independent toy.
    first=evaluate(links,e,phis,paired=True)
    groups=[base.prior.group(69220+i,.19) for i in range(4)]
    lt,et,pt,_=base.transform(links,e,phis,groups)
    transformed=evaluate(lt,et,pt)
    rt=links[:,[2,3,0,1]].copy()
    rt[1]=np.array([a.conj().T for a in links[1]])
    reflected=evaluate(rt,e[[2,3,0,1]],phis[[2,3,0,1]])
    changed=evaluate(links,np.roll(e,1,axis=0),phis)
    gauge=float(abs(transformed['source']-first['source'])/abs(first['source']))
    conjugacy=float(abs(reflected['source']-first['source'].conjugate())/abs(first['source']))
    assert gauge<2e-9 and conjugacy<2e-9
    # E dependence belongs wholly to the inherited auxiliary factor. Dividing
    # here is a numerical check only at a resolved nonzero fixture, not a definition.
    e_identity=float(abs(changed['source']*first['auxiliary']-first['source']*changed['auxiliary'])/
        max(abs(changed['source']*first['auxiliary']),abs(first['source']*changed['auxiliary'])))
    assert e_identity<2e-9
    rows=[]
    for label,data in [('original',first),('local_G_transformed',transformed),
                       ('reflected',reflected),('changed_auxiliary',changed)]:
        rows.append(dict(label=label,source=base.old.cpair(data['source']),
            factored=base.old.cpair(data['factored']),factor_error=data['factor_error'],
            scalar_weight=base.old.cpair(data['weight']),Wilson_gap=data['gap']))
    deps=['research_note_675.md','research_note_676.md','research_note_679.md','research_note_680.md',
          'research_note_691.md','joint_gauss_boundary_functional.py','joint_rational_physical_limit.py',
          'joint_physical_source_reflection.py','joint_neutral_symmetry_reduction_results.json']
    return dict(date='2026-10-02',entry_round=692,latest_formal_round=691,
        not_formal_round=True,rows=rows,local_G_invariance_error=gauge,
        reflection_conjugacy_error=conjugacy,neutral_even_pair_error=first['pair_error'],
        exact_auxiliary_factor_identity_error=e_identity,
        source_definition_uses_no_inverse_or_weight_division=True,
        all_original_16_auxiliary_channels_retained=True,
        complete_S9_formula_inherited_from675_not_numerically_evaluated=True,
        complete_double_Haar_and_Hb_average_not_evaluated=True,
        no_RP_or_original_HF_identification=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    out=run()
    if TARGET.exists():assert out==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('entry_round','local_G_invariance_error',
        'reflection_conjugacy_error','neutral_even_pair_error','exact_auxiliary_factor_identity_error')}))
