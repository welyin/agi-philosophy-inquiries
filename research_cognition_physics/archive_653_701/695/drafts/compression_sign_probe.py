"""Explore exact complementary-subspace determinant reduction, not Haar sampling."""
import sys,json,importlib.util
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as b
import joint_rational_physical_limit as soft
import joint_critical_source_certificate as cert


def graph(links):
    u,v,d,h,gap=b.kernel(links)
    small=np.array([[1,0],[1,0],[0,1],[0,1]])/np.sqrt(2)
    frame=np.kron(np.kron(np.eye(4),small),np.eye(16))
    ev,vec=np.linalg.eigh(frame.conj().T@h@frame)
    if sum(ev<0)!=64:return None
    p=vec[:,ev<0]@vec[:,ev<0].conj().T
    top=np.array([32*i+j for i in range(4) for j in range(16)])
    bottom=top+16
    c=p[np.ix_(bottom,bottom)];w=-np.linalg.solve(c.T,p[np.ix_(top,bottom)].T).T
    k=b.internal.spin.G5@b.internal.spin.GAMMA[1]
    local=small.T@b.internal.B@k@small
    assert max(abs(local[0,0]),abs(local[1,1]),abs(local[0,1]-local[1,0]))<1e-14
    return dict(W=w,gap=gap,chart_min=float(np.linalg.eigvalsh(c)[0]),factor=local[0,1])


def pairing(e,factor):
    out=np.zeros((64,64),complex)
    for i,ei in enumerate(e):
        out[16*i:16*(i+1),16*i:16*(i+1)]=factor*sum(x*t for x,t in zip(ei,b.internal.T))
    return out


def coefficient(data,e):
    t=pairing(e,data['factor']);w=data['W'];x=w.T@t
    a=np.eye(64)-x@x.conj()
    sg,log=np.linalg.slogdet(a)
    return sg,float(log),a


def run():
    spec=importlib.util.spec_from_file_location('e694',ROOT/'round694_drafts/critical_source_entry.py')
    entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
    fixtures=[('full_group_fixture',soft.fixture()[0]),
              ('flux_anchor',entry.links_at(4*np.arctan(float(cert.LEFT))/np.pi))]
    for seed in (67351,69571,69573):
        links=np.array([[b.prior.rep(*b.prior.group(seed+8*mu+i,.9)) for i in range(4)] for mu in range(2)])
        fixtures.append((f'full_group_{seed}',links))
    results=[];negative=[]
    e0=np.tile(np.eye(10)[0],(4,1))
    for label,links in fixtures:
        data=graph(links)
        if data is None:results.append(dict(label=label,unbalanced=True));continue
        ref,logref,_=coefficient(data,e0);minrel=1.;worst=None
        for seed in range(69500,69564):
            rng=np.random.default_rng(seed);e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
            phase,log,a=coefficient(data,e);rel=phase/ref
            assert abs(rel.imag)<1e-7
            if rel.real<0:
                sv=float(np.linalg.svd(a,compute_uv=False)[-1])
                negative.append(dict(label=label,seed=seed,relative_phase=[float(rel.real),float(rel.imag)],
                    log_absolute_relative=log-logref,minimum_singular=sv))
                worst=seed;break
            minrel=min(minrel,float(np.exp(log-logref)))
        results.append(dict(label=label,chart_min=data['chart_min'],W_norm=float(np.linalg.norm(data['W'],2)),
            minimum_sample_relative=minrel,negative_seed=worst))
    return dict(rows=results,negative=negative,
        samples_are_only_counterexample_search_not_integral_or_positivity_certificate=True)


if __name__=='__main__':
    result=run();print(json.dumps(result))
    target=HERE/'compression_sign_probe_results.json'
    if target.exists():assert result==json.loads(target.read_text('utf8'))
    else:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,indent=2)+'\n')
