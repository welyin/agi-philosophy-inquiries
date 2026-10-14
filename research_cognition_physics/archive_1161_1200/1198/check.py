"""1198: same-parent Hartree response audit, not a finite-K gap certificate.

Group 1: exact rational bounds conditional on the proved spectral identities.
Group 2: independent position-space versus momentum-block source response.
No image, installation, fitting, or frozen result overwrite.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def rational_bounds():
    eps, nu, g, om, m, ell = map(F, ('3/20','7/20','1/5','1/2','2','16/65'))
    zlo, w = F(49,50), F(1,4)
    sq3, xs = F(17321,10000), F(199,1000)
    assert sq3*sq3>3 and xs*xs>1-zlo*zlo
    n, v = 1+eps*zlo, 1+nu*zlo
    # For every three-dimensional torus with n_b>=3, <sin^2>=1/2.
    # At zlo, Emax<10 and the same fixed-point map satisfies F(zlo)>zlo.
    emax2=(m*n)**2+3*(v/ell)**2
    blo=2*(eps*m*m*n+nu*v*3/(2*ell**2))/10+eps*om
    root_margin=blo*blo*(1-zlo*zlo)-(om*zlo)**2
    assert emax2<100 and root_margin>0
    # Reuse 1196's volume-independent full-interval contraction inequality.
    abs_b_min=2*eps*(1-nu)*m/(1+nu)+eps*om
    bp_no_sqrt3=4*m*(eps-nu)**2/(3*(1-nu)**3)
    lipschitz2=om**4*bp_no_sqrt3**2/(3*(om**2+abs_b_min**2)**3)
    assert lipschitz2<F(6,25)**2
    bs = 2*eps**2*m/n+2*nu**2*sq3/(v*ell)
    bt = 2*g**2*sq3/(v*ell)
    cs, ct = xs**3/om, 1/(om*n)
    rs2, rt2, gap = om**2/(1-zlo*zlo), (om*n)**2, 2*m*n
    static = cs*bs+ct*bt
    feedback = (cs*bs/(1-w*w/(4*rs2))+ct*bt/(1-w*w/(4*rt2)))/(1-w/gap)
    assert static<F(753,1000) and feedback<F(835,1000)
    old=HERE.parent/'1196/certificate_results.json'
    root=json.loads(old.read_text(encoding='utf-8'))['fixed_point']['root_enclosure']
    assert F(root['lo'])>zlo and F(root['hi'])<1
    def pair(x): return {'exact':str(x),'decimal':float(x)}
    return {'root_certificate_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
            'z_lower':str(zlo), 'sqrt3_upper':str(sq3), 'xs_upper':str(xs),
            'uniform_periodic_boxes':{'minimum_sites_each_axis':3,
                'Emax_squared_at_zlo':str(emax2),'absolute_b_lower_at_zlo':str(blo),
                'root_bracket_squared_margin':str(root_margin),
                'inherited_global_lipschitz_squared_upper':str(lipschitz2),
                'all_same_spacing_filled_fixed_points_above_zlo':True},
            'Bs_upper':pair(bs),'Bt_upper':pair(bt),'source_pair_gap_lower':str(gap),
            'static_feedback_norm_upper':pair(static),
            'low_frequency_feedback_norm_upper':pair(feedback),
            'weighted_inverse_margin_lower':pair(1-feedback),
            'frequency_interval':['-1/4','1/4']}


def source_response_check():
    # A small three-dimensional torus audits the full matrix construction.
    # It does NOT replace the 65^3 reference or certify a new preparation.
    n=3; coords=np.array(list(itertools.product(range(n),repeat=3)))
    volume=len(coords); ell=16/65
    eps,nu,g,m,z=.15,.35,.2,2.,.989835554078098
    lapse,speed=1+eps*z,1+nu*z
    eye=np.eye(2); sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1,-1])
    beta=np.kron(sz,eye); alpha=[np.kron(sx,a) for a in (sx,sy,sz)]
    index={tuple(x):i for i,x in enumerate(coords)}
    # Scalar-position bond generators, built without Fourier vertex formulas.
    hopping=[]
    for b in range(3):
        mat=np.zeros((volume,volume),complex)
        for i,x in enumerate(coords):
            y=x.copy(); y[b]=(y[b]+1)%n; j=index[tuple(y)]
            mat[i,j]+=-.5j/ell; mat[j,i]+=.5j/ell
        hopping.append(np.kron(mat,alpha[b]))
    mass=np.kron(np.eye(volume),m*beta)
    h=lapse*mass+speed*sum(hopping)
    ev,u=np.linalg.eigh(h); fill=(ev<0).astype(float)
    assert np.min(np.abs(ev))>m*lapse-1e-12
    fourier=np.exp(2j*np.pi*(coords@coords.T)/n)/np.sqrt(volume)
    ft=np.kron(fourier,np.eye(4)); hf=ft.conj().T@h@ft
    kval=np.sin(2*np.pi*coords/n)/ell
    blocks=[m*lapse*beta+speed*sum(k[b]*alpha[b] for b in range(3)) for k in kval]
    hdiag=np.zeros_like(hf)
    for i,block in enumerate(blocks): hdiag[4*i:4*i+4,4*i:4*i+4]=block
    vertex_error=float(np.max(np.abs(hf-hdiag)))
    response_error=0.; hermitian_error=0.; minimum_qeig=1e9
    checks=[]
    for q in [(0,0,0),(1,0,0),(1,1,0),(1,2,1)]:
        phase=np.exp(2j*np.pi*(coords@np.array(q))/n)
        phase_op=np.kron(np.diag(phase),np.eye(4))
        # Each edge belongs half to each endpoint.
        half=[(phase_op@t+t@phase_op)/2 for t in hopping]
        ds=eps*phase_op@mass+nu*sum(half); dt=g*(half[0]-half[1])
        dreal=[ds,dt]; dft=[ft.conj().T@d@ft for d in dreal]
        dblocks=[]; destinations=[]
        for j,p in enumerate(coords):
            dest=index[tuple((p+np.array(q))%n)]; destinations.append(dest)
            k=(kval[j]+kval[dest])/2
            dblocks.append([eps*m*beta+nu*sum(k[b]*alpha[b] for b in range(3)),
                            g*(k[0]*alpha[0]-k[1]*alpha[1])])
        for a in range(2):
            expected=np.zeros_like(h)
            for j,dest in enumerate(destinations):
                expected[4*dest:4*dest+4,4*j:4*j+4]=dblocks[j][a]
            vertex_error=max(vertex_error,float(np.max(np.abs(dft[a]-expected))))
        eig_d=[u.conj().T@d@u for d in dreal]
        for w in [0.,.25,-.25]:
            # Full position-space Liouville linearization:
            # w*dC = [h,dC]+[D,C]. Terms within an occupied/empty band vanish.
            numerator=fill[None,:]-fill[:,None]
            mask=np.abs(numerator)>0
            coeff=np.zeros_like(numerator)
            coeff[mask]=numerator[mask]/(w-(ev[:,None]-ev[None,:])[mask])
            pi_real=np.array([[np.vdot(eig_d[a],coeff*eig_d[b])/volume
                               for b in range(2)] for a in range(2)])
            # Independent four-dimensional blocks and exact spectral projectors.
            pi_block=np.zeros((2,2),complex)
            for j,dest in enumerate(destinations):
                e0=np.sqrt((m*lapse)**2+speed**2*np.dot(kval[j],kval[j]))
                e1=np.sqrt((m*lapse)**2+speed**2*np.dot(kval[dest],kval[dest]))
                for s0 in [-1,1]:
                    p0=(np.eye(4)+s0*blocks[j]/e0)/2
                    for s1 in [-1,1]:
                        if s0==s1: continue
                        p1=(np.eye(4)+s1*blocks[dest]/e1)/2
                        factor=((s0<0)-(s1<0))/(w-(s1*e1-s0*e0))
                        for a in range(2):
                            for b in range(2):
                                pi_block[a,b]+=factor*np.trace(dblocks[j][a].conj().T@p1@dblocks[j][b]@p0)/volume
            response_error=max(response_error,float(np.max(np.abs(pi_real-pi_block))))
            hermitian_error=max(hermitian_error,float(np.max(np.abs(pi_real-pi_real.conj().T))))
            minimum_qeig=min(minimum_qeig,float(np.linalg.eigvalsh(-pi_real)[0]))
            checks.append({'q':list(q),'omega':w,'Pi_real':pi_real.real.tolist()})
    assert vertex_error<2e-13 and response_error<2e-12
    assert hermitian_error<2e-13 and minimum_qeig>-2e-13
    return {'diagnostic_torus':[n,n,n], 'full_single_particle_dimension':4*volume,
            'vertex_max_error':vertex_error,'response_max_error':response_error,
            'hermitian_max_error':hermitian_error,'min_Q_eigenvalue':minimum_qeig,
            'tests':checks,'scope':'Numerical algebra audit only, not a finite-K physical gap or detector certificate.'}


def run():
    return {'schema':'original_source_ir_response_1198_v1','status':'PASS','check_groups':2,
            'rational_certificate':rational_bounds(),'position_momentum_audit':source_response_check(),
            'scope':'Hartree retarded source-relation response of the specified filled reference only; no full finite-K gap, finite-window GR no-go, or dark-sector identity.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    result=run(); result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
    if args.write:
        with (HERE/'results.json').open('x',encoding='utf-8') as stream:stream.write(output)
    print(output,end='')
