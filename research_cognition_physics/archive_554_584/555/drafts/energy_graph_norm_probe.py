"""Preliminary 555 graph-norm identity on the inherited two-scalar one-cell potential.

Exact Gaussian polynomial quadrature, not a full-lattice theorem or a thermal
preparation result. No completed research round is registered here.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'energy_graph_norm_probe_results.json'


def run():
    source=HERE.parent/'joint_singlet_common_mass_rg_results.json'
    row=json.loads(source.read_text('utf8'))['examples'][2]['state']
    L=np.array([[row['lambda_H'],row['p']],[row['p'],row['lambda_s']]])
    C=.25*np.array([row['x'],row['y']]);u=np.linalg.solve(L,C);mu=np.sqrt(u)
    ell=np.linalg.eigvalsh(L)[0];r=1.;dimension=2
    coeff=np.array([3*L[0,0]+L[0,1],3*L[1,1]+L[0,1]])
    alpha=max(0.,*coeff)
    A=2*alpha/(ell*r)
    B=alpha*(sum(u)+dimension*r/2)+sum(abs(C))
    hbar=v=1.;a=hbar*hbar/(2*v);variance=.2
    nodes,weights=np.polynomial.hermite.hermgauss(12)
    grid=np.array(np.meshgrid(nodes,nodes,indexing='ij')).reshape(2,-1).T
    w=np.outer(weights,weights).ravel()/np.pi
    q=mu+np.sqrt(2*variance)*grid
    x=q*q;delta=x-u
    W=np.einsum('ni,ij,nj->n',delta,L,delta)/4
    lapW=x@coeff-sum(C)
    assert np.max(lapW-A*W-B)<0
    def moments(momentum):
        loggrad=-(q-mu)/(2*variance)+1j*momentum/hbar
        T=-a*np.sum(loggrad**2-1/(2*variance),axis=1)
        E1=float(np.real(w@(T+W)))
        E2=float(w@abs(T+W)**2)
        T2=float(w@abs(T)**2);W2=float(w@(W*W))
        cross=float(2*a*(w@(W*np.sum(abs(loggrad)**2,axis=1)))-a*(w@lapW))
        residual=abs(E2-(T2+W2+cross))/max(1.,E2)
        assert residual<1e-12
        bound=E2+a*A*E1+a*B
        assert T2+W2<=bound+1e-10*max(1.,bound)
        return np.array([E1,E2,T2,W2,bound,residual])
    base=moments(np.zeros(2));rows=[]
    for M in (2.,4.,16.,64.):
        weight=1/(M*M)
        plus=moments(np.array([M,0.]));minus=moments(np.array([-M,0.]))
        mixed=(1-weight)*base+weight*(plus+minus)/2
        assert abs(mixed[0]-(base[0]+.5))<1e-12
        rows.append(dict(M=M,mixture_weight=weight,E1=float(mixed[0]),E2=float(mixed[1]),
            T2_plus_W2=float(mixed[2]+mixed[3]),graph_norm_bound=float(mixed[4]),
            largest_component_identity_relative_residual=float(max(base[5],plus[5],minus[5]))))
    assert all(rows[i+1]['E2']>rows[i]['E2'] for i in range(3))
    return dict(status='preliminary_555_not_a_completed_round',L=L.tolist(),C=C.tolist(),u=u.tolist(),
        laplacian_potential_majorant=dict(A=float(A),B=float(B)),examples=rows,
        exact_polynomial_Gaussian_quadrature=True,operator_domain_extension_not_yet_reviewed=True,
        finite_lattice_full_readout_bound_not_yet_verified=True,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())


if __name__=='__main__':
    result=run()
    if TARGET.exists():
        assert json.loads(TARGET.read_text('utf8'))==result
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
