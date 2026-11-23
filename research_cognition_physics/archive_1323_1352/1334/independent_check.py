"""1334 independent full-mass, charged-current and invariant diagnostics.
No author imports. Exact rational local-rank certificate uses stdlib Fraction.
Only --write-new writes; --check is read-only. No apparatus is constructed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import permutations
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
S = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
I = np.eye(3)
TOL = 3e-9


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(a):
    return float(np.linalg.norm(a, 2))


def make(m, r, q, g, heavy, h):
    d = r * (I + q * S)
    k = np.zeros((3*m, 3*(m+1)), complex)
    for j in range(m):
        k[3*j:3*j+3, 3*j:3*j+3] = heavy * I
        k[3*j:3*j+3, 3*j+3:3*j+6] = -heavy * d
    p = np.zeros((3, 3*(m+1)))
    p[:, -3:] = I
    powers = [np.linalg.matrix_power(d, j) for j in range(m+1)]
    stack = np.vstack(powers)
    z = stack.conj().T @ stack
    # Stacked QR avoids cancellation in the strongly hierarchical Gram spectrum.
    _, rr = np.linalg.qr(stack, mode='reduced')
    _, sqrtz, vh = np.linalg.svd(rr)
    zi_half = (vh.conj().T / sqrtz) @ vh
    zero = np.vstack(powers[::-1]) @ zi_half
    mass = np.vstack([h*g*p, k])
    left, sig, right_h = np.linalg.svd(mass)
    order = np.arange(len(sig)-1, -1, -1)
    left, sig, right = left[:, order], sig[order], right_h.conj().T[:, order]
    # For q,g real, numerical SVD is real up to storage dtype.
    if abs(np.imag(q)) + abs(np.imag(g)) == 0:
        left = left.real
        right = right.real
    return dict(d=d, k=k, p=p, z=z, sqrtz=sqrtz, zi_half=zi_half, zero=zero,
                mass=mass, left=left, sig=sig, right=right, g=g, h=h,
                m=m, r=r, q=q, heavy=heavy)


def moments(m, r):
    kk = np.arange(m+1, dtype=float)
    cc = kk*(kk-1)/2
    w = r**(2*kk)
    return np.array([sum(w), sum(w*kk), sum(w*kk**2),
                     sum(w*cc), sum(w*kk*cc), sum(w*cc**2)])


def moment_z(m, r, q):
    a, b, c, d, e, f = moments(m, r)
    x = abs(q)**2
    return np.array([[a+c*x+f*x*x, q.conjugate()*(b+e*x), q.conjugate()**2*d],
                     [q*(b+e*x), a+c*x, q.conjugate()*b], [q*q*d, q*b, a]], complex)


def current(u, d, phase):
    qq = np.diag(np.exp(1j*phase*np.arange(3)))
    return u['left'][:3, :].conj().T @ qq @ d['left'][:3, :]


def phase_coeff(a, b):
    t = a*b
    return np.array([sum(t*t)-2*t[0]*t[2],
                     2*(t[0]*t[1]+t[1]*t[2]), 4*t[0]*t[2]])


def qrange(coeff):
    a, b, c = coeff
    xs = [-1., 1.]
    if c != 0 and -1 <= -b/(2*c) <= 1:
        xs.append(-b/(2*c))
    values = [a+b*x+c*x*x for x in xs]
    return min(values), max(values)


def response(obj, z):
    _, ss, vh = np.linalg.svd(obj['k'], full_matrices=False)
    ends = obj['p'] @ vh.conj().T
    return abs(obj['g'])**2 * (ends/(ss*ss-z)) @ ends.conj().T


def determinant(a):
    a = [[F(x) for x in row] for row in a]
    ans, n = F(1), len(a)
    for col in range(n):
        pivot = next((j for j in range(col, n) if a[j][col]), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            a[pivot], a[col] = a[col], a[pivot]
            ans = -ans
        p = a[col][col]
        ans *= p
        for j in range(col+1, n):
            fac = a[j][col]/p
            for k in range(col+1, n):
                a[j][k] -= fac*a[col][k]
            a[j][col] = F(0)
    return ans


def det_derivative(a, da):
    result = F(0)
    for perm in permutations(range(3)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        for row in range(3):
            prod = F(sign)
            for i in range(3):
                prod *= da[i, perm[i]] if i == row else a[i, perm[i]]
            result += prod
    return result


def rational_rank():
    # Invariants of H_f^{-1}=Z_f/|g_f|^2, so their ranks are physical spectral ranks.
    si = np.array([[F(int(x)) for x in row] for row in S], dtype=object)
    eye = np.array([[F(int(x)) for x in row] for row in I], dtype=object)
    def data(t, x, amp):
        z = np.zeros((3,3), dtype=object)
        dt = np.zeros((3,3), dtype=object)
        dx = np.zeros((3,3), dtype=object)
        for k in range(3):
            c = k*(k-1)//2
            b = eye+k*x*si+c*x*x*(si@si)
            db = k*si+2*c*x*(si@si)
            z += t**k*(b.T@b)
            if k:
                dt += k*t**(k-1)*(b.T@b)
            dx += t**k*(db.T@b+b.T@db)
        tt = z/amp
        deriv = [dt/amp, dx/amp, -z/(amp*amp)]
        rows = []
        for func in range(3):
            rows.append([np.trace(e) if func == 0 else
                         2*np.trace(tt@e) if func == 1 else det_derivative(tt,e)
                         for e in deriv])
        return tt, rows
    # Coordinates (t=r²,x=|qu|,y=|qd|,A=|gu|²,B=|gd|²,c=cos delta).
    tu, ru = data(F(1), F(1), F(1))
    td, rd = data(F(1), F(2), F(1))
    jac = [[F(0) for _ in range(5)] for __ in range(5)]
    for row in range(3):
        jac[row][0], jac[row][1], jac[row][3] = ru[row]
    for row in range(2):
        jac[3+row][0], jac[3+row][2], jac[3+row][4] = rd[row]
    five = determinant(jac)
    # Mixed inverse-mass invariant tr(Tu Q Td Q†), derivative at cos delta=0.
    six_factor = 2*(tu[0,1]*td[0,1]+tu[1,2]*td[1,2])
    return dict(point=dict(m=2,r_squared=1,qu_absolute=1,qd_absolute=2,
                           gu_squared=1,gd_squared=1,cos_delta=0),
                five_by_five_jacobian=[[str(x) for x in row] for row in jac],
                mass_minor=str(five), mixed_phase_derivative=str(six_factor),
                six_minor=str(five*six_factor), rank_six=bool(five*six_factor))


def run():
    count, maxerr = 0, 0.
    def close(a,b,label):
        nonlocal count,maxerr
        err = float(np.max(abs(np.asarray(a)-np.asarray(b))))/(1+float(np.max(abs(np.asarray(b)))))
        count += 1
        maxerr = max(maxerr,err)
        if err > TOL:
            raise AssertionError((label,err))
    def positive(x,label):
        nonlocal count
        count += 1
        if x < -TOL:
            raise AssertionError((label,x))

    moment_rows=[]
    for m,r,q in [(1,.7,.3+.2j),(2,1.,1j),(2,1.3,.7),
                   (3,.8,1.2-.4j),(5,1.1,.5+.3j)]:
        o=make(m,r,q,.2,10.,1.)
        close(o['z'],moment_z(m,r,complex(q)),'six moments')
        close(o['zero'].conj().T@o['zero'],I,'zero normalization')
        close(o['k']@o['zero'],np.zeros((3*m,3)),'zero kernel')
        close(o['p']@o['zero'],o['zi_half'],'fixed Rm portal')
        close(np.linalg.det(o['mass']),(-1)**(3*m)*(o['h']*o['g'])**3*o['heavy']**(3*m),'full determinant')
        moment_rows.append(dict(m=m,r=r,q=[float(np.real(q)),float(np.imag(q))],moments=moments(m,r).tolist()))

    rank=rational_rank()
    positive(int(rank['rank_six'])-.5,'exact rank six')
    # Cross-check the nonzero mixed derivative through ordinary matrices.
    zu=moment_z(2,1.,1+0j); zd=moment_z(2,1.,2+0j)
    def mixed(c):
        q=np.diag(np.exp(1j*np.arccos(c)*np.arange(3)))
        return np.trace(zu@q@zd@q.conj().T).real
    close((mixed(.125)-mixed(-.125))/.25,float(F(rank['mixed_phase_derivative'])),'exact quadratic phase derivative')

    finite=[]
    for label,r,xu,xd,gu,gd,mu,md in [('independent',1.3,.7,1.1,.2,.3,10.,13.),
                                     ('author_scope',1.,1.,1.,.2,.3,10.,13.)]:
        u=make(2,r,xu,gu,mu,1.);d=make(2,r,xd,gd,md,1.)
        pp=np.zeros((9,9));pp[:3,:3]=I
        phase_rows=[]
        for delta in (0.,.3,.8,np.pi/2,2.7,np.pi):
            q=np.diag(np.exp(1j*delta*np.arange(3)))
            qall=np.kron(np.eye(3),q)
            dc=make(2,r,xd*np.exp(1j*delta),gd,md,1.)
            close(dc['mass'],qall@d['mass']@qall.conj().T,'full phase dictionary')
            close(dc['sig'],d['sig'],'phase-independent full masses')
            mapped=qall@d['left']
            close(dc['mass']@dc['mass'].conj().T@mapped,mapped*(d['sig']**2),'complete eigendictionary')
            v=current(u,d,delta)
            close(v@v.conj().T,u['left'].conj().T@pp@u['left'],'full up source projection')
            close(v.conj().T@v,d['left'].conj().T@pp@d['left'],'full down source projection')
            nu=u['left'][:3,:3].T@u['left'][:3,:3]
            nd=d['left'][:3,:3].T@d['left'][:3,:3]
            close(nu,v[:3,:3]@v[:3,:3].conj().T+v[:3,3:]@v[:3,3:].conj().T,'light-heavy up completion')
            close(nd,v[:3,:3].conj().T@v[:3,:3]+v[3:,:3].conj().T@v[3:,:3],'light-heavy down completion')
            for i in range(3):
                for j in range(3):
                    coeff=phase_coeff(u['left'][:3,i],d['left'][:3,j])
                    c=np.cos(delta)
                    pred=coeff@[1,c,c*c]
                    close(abs(v[i,j])**2,pred,'one common phase polynomial')
                    lo,hi=qrange(coeff)
                    positive(pred-lo,'global phase lower');positive(hi-pred,'global phase upper')
            cp=np.imag(v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate())
            phase_rows.append(dict(delta=delta,CP_quartet=float(cp),weak_light_norm_defect=norm(I-nu)))
        positive(max(abs(x['CP_quartet']) for x in phase_rows)-1e-5,'nonzero finite CP witness')
        for o in (u,d):
            gap=np.linalg.svd(o['k'],compute_uv=False)[-1]
            leading=np.sort(abs(o['g'])/o['sqrtz'])
            zeta=leading[-1]**2/gap**2
            positive(1-zeta,'three-light safe gap')
            positive(o['sig'][3]-gap,'heavy lower gap')
            c6=response(o,0.)
            gram=o['k']@o['k'].conj().T
            direct=abs(o['g'])**2*o['p']@o['k'].conj().T@np.linalg.solve(gram,np.linalg.solve(gram,o['k']))@o['p'].conj().T
            close(c6,direct,'same-parent heavy response')
            for j,z in enumerate(o['sig'][:3]**2):
                a=o['left'][:3,j]
                cz=response(o,z)
                schur=abs(o['g'])**2*(o['zi_half']@o['zi_half'])-z*(I+cz)
                close(schur@a,np.zeros(3),'exact physical Schur root')
                positive(np.linalg.eigvalsh(cz-c6)[0],'response lower')
                positive(np.linalg.eigvalsh(c6/(1-zeta)-cz)[0],'response upper')
            # Analytic source matrices, differentiated before matching/diagonalization.
            for name in ('h','g','heavy','r','qreal','qimag'):
                eps=1e-5
                pars=dict(m=o['m'],r=o['r'],q=o['q'],g=o['g'],heavy=o['heavy'],h=o['h'])
                plus=pars.copy();minus=pars.copy()
                key='q' if name.startswith('q') else name
                direction=1j if name=='qimag' else 1.
                plus[key]+=eps*direction;minus[key]-=eps*direction
                delta_m=(make(**plus)['mass']-make(**minus)['mass'])/(2*eps)
                analytic=np.zeros_like(o['mass'])
                if name=='h':analytic[:3]=o['g']*o['p']
                elif name=='g':analytic[:3]=o['h']*o['p']
                elif name=='heavy':analytic[3:]=o['k']/o['heavy']
                else:
                    dd=(I+o['q']*S) if name=='r' else o['r']*direction*S
                    for j in range(o['m']):
                        analytic[3+3*j:6+3*j,3*j+3:3*j+6]=-o['heavy']*dd
                close(delta_m,analytic,'full matrix source '+name)
                sqp=make(**plus)['mass'];sqm=make(**minus)['mass']
                numeric=(sqp.conj().T@sqp-sqm.conj().T@sqm)/(2*eps)
                close(numeric,analytic.conj().T@o['mass']+o['mass'].conj().T@analytic,'source square '+name)
        finite.append(dict(label=label,m=2,r=r,qu=xu,qd_absolute=xd,g_u=gu,g_d=gd,M_u=mu,M_d=md,h=1.,
                           up_masses=u['sig'].tolist(),down_masses=d['sig'].tolist(),phase_rows=phase_rows,
                           up_zeta=float((abs(gu)/u['sqrtz'][-1]/np.linalg.svd(u['k'],compute_uv=False)[-1])**2),
                           down_zeta=float((abs(gd)/d['sqrtz'][-1]/np.linalg.svd(d['k'],compute_uv=False)[-1])**2)))

    minmax=[]
    for heavy in (.3,1.,10.):
        for g in (.2,1.5):
            o=make(2,1.3,.7+.4j,g,heavy,1.)
            upper=abs(g)/o['sqrtz'][-1]
            positive(upper-o['sig'][2],'all-M third light bound')
            minmax.append(dict(M=heavy,g=g,third_mass=float(o['sig'][2]),upper=float(upper)))

    strong=[]
    for flavor,q,target in [('u',10.,171.7),('d',1.6,2.89)]:
        o=make(5,8.,q,1.,1.,174.)
        lower=target*o['sqrtz'][-1]/174.
        rigorous=target/174.*8.**5/(1+5*q+15*q*q)
        positive(lower-rigorous,'analytic conservative portal bound')
        if flavor=='u':positive(rigorous-20.,'strict high portal obstruction')
        strong.append(dict(flavor=flavor,m=5,r=8.,q=q,h_GeV=174.,target_GeV=target,
                           sqrt_z_min=float(o['sqrtz'][-1]),g_necessary_lower=float(lower),
                           conservative_analytic_lower=rigorous,loop_diagnostic_lower=float(lower**2/(16*np.pi**2))))

    limits=[]
    for r,q in [(1.3,10.),(1.3,30.),(1.3,100.),(2.,1.),(5.,1.),(10.,1.)]:
        o=make(2,r,q,1.,10.,1.)
        limits.append(dict(r=r,q=q,leading_highest=1/float(o['sqrtz'][-1]),
                           sqrtz=o['sqrtz'].tolist()))
    # m2 exact moment Schur pivots (q -> infinity at fixed r): r^4,r^2,1.
    for r in (.7,1.,1.3,2.):
        a,b,c,d,e,f=moments(2,r)
        pivot2=c-e*e/f
        pivot3=a-d*d/f-(b-e*d/f)**2/pivot2
        close(pivot2,r*r,'fixed-r second asymptotic pivot')
        close(pivot3,1.,'fixed-r third asymptotic pivot')

    return dict(status='PASS',checks=count,max_normalized_residual=maxerr,
                provenance={'independent_check.py':digest(Path(__file__)),
                            'round1334_admission.md':digest(HERE.parent/'_admission/round1334_admission.md')},
                moment_cases=moment_rows,exact_rational_local_rank=rank,finite_parent_cases=finite,
                all_M_minmax_cases=minmax,strong_portal_case=strong,nonuniform_limit_diagnostics=limits,
                scope='Finite tree-level matrices and conditional source identities; no apparatus, loop control theorem or realistic global fit.',
                fixed_m_two_rank_is_local=True,three_generations_input=True,strong_CP_not_solved=True,
                new_axioms=0,stage_complete=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),(path,'keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert np.isclose(a,b,rtol=1e-8,atol=3e-10),(path,a,b)
    else:assert a==b,(path,a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write-new',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();result=run();out=HERE/'independent_results.json'
    if args.write_new:
        with out.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(result,json.loads(out.read_text(encoding='utf-8')))
    print(json.dumps({k:result[k] for k in ('status','checks','max_normalized_residual')}))
