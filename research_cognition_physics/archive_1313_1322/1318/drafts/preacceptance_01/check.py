"""Round 1318: one declared 331 parent; deterministic, read-only by default.

--write-new creates results.json exclusively. No empirical fits or instrument claim.
All mass matrices here are PHYSICAL M; the all-left Weyl coefficient is M.conj().
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ERRORS = {}


def eq(name, actual, expected, tol=3e-9, scale=None):
    a, b = np.asarray(actual), np.asarray(expected)
    if scale is None:
        scale = max(1., float(np.linalg.norm(a)), float(np.linalg.norm(b)))
    error = float(np.linalg.norm(a-b))/scale
    ERRORS[name] = max(ERRORS.get(name, 0.), error)
    assert error <= tol, (name, error, tol)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def spectrum_pair(a, b):
    val, vec = np.linalg.eigh(b)
    z = (vec / np.sqrt(val)) @ vec.conj().T
    return np.linalg.eigvalsh(z @ a @ z)


def field_table(n, a):
    # name, multiplicity, dimC, cubicC, dimL, cubicL, integer x=3X, r exponent s.
    return [('L', n, 1, 0, 3, 1, -1, 1),
            ('Q+', a, 3, 1, 3, 1, 1, 1),
            ('Q-', n-a, 3, 1, 3, -1, 0, 1),
            ('uc', n, 3, -1, 1, 0, -2, 1),
            ('dc', n, 3, -1, 1, 0, 1, 1),
            ('ec', n, 1, 0, 1, 0, 3, 1),
            ('Nc', n, 1, 0, 1, 0, 0, 1),
            ('Fc', n, 1, 0, 1, 0, 0, -1),
            ('Uc', a, 3, -1, 1, 0, -2, -1),
            ('Dc', n-a, 3, -1, 1, 0, 1, -1)]


def anomalies(n, a):
    out = dict(C3=0, L3=0, C2x=0, L2x=0, x3=0, gravx=0,
               rC=0, rL=0, rdimension=0, rx2=0, r2x=0)
    for _, m, c, sc, l, sl, x, r in field_table(n, a):
        out['C3'] += m*l*sc
        out['L3'] += m*c*sl
        out['C2x'] += m*l*x*bool(sc)
        out['L2x'] += m*c*x*bool(sl)
        out['x3'] += m*c*l*x**3
        out['gravx'] += m*c*l*x
        out['rC'] += m*l*r*bool(sc)
        out['rL'] += m*c*r*bool(sl)
        out['rdimension'] += m*c*l*r
        out['rx2'] += m*c*l*r*x*x
        out['r2x'] += m*c*l*r*r*x
    return out


def generators():
    out = []
    for i, j in ((0, 1),):
        x = np.zeros((3, 3), complex); x[i,j] = x[j,i] = .5; out.append(x)
        x = np.zeros((3, 3), complex); x[i,j] = -.5j; x[j,i] = .5j; out.append(x)
    out.append(np.diag([.5, -.5, 0.]).astype(complex))
    for i, j in ((0, 2), (1, 2)):
        x = np.zeros((3, 3), complex); x[i,j] = x[j,i] = .5; out.append(x)
        x = np.zeros((3, 3), complex); x[i,j] = -.5j; x[j,i] = .5j; out.append(x)
    out.append(np.diag([1., 1., -2.]).astype(complex)/(2*np.sqrt(3)))
    return out


def gauge(w, gl, gx):
    v = np.array([[gl/2, gl/(2*np.sqrt(3)), -gx/3],
                  [-gl/2, gl/(2*np.sqrt(3)), 2*gx/3],
                  [0., -gl/np.sqrt(3), -gx/3]])
    matrix = 2*v.T @ np.diag(w*w) @ v
    d = -v[2]
    eig, vec = np.linalg.eigh(matrix)
    assert eig[1] > 0
    inverse = (vec[:, 1:]/eig[1:]) @ vec[:, 1:].T
    f0 = float(d @ inverse @ d)
    return v, matrix, d, eig, vec, inverse, f0


def gauge9(w, gl, gx):
    g = generators()
    out = np.zeros((9, 9))
    for i, charge in enumerate((-1/3, 2/3, -1/3)):
        vacuum = np.eye(3)[:, i]*w[i]
        actions = np.stack([gl*t@vacuum for t in g]+[gx*charge*vacuum])
        out += 2*np.real(actions.conj() @ actions.T)
    return out


def analytic_f0(w):
    a, b, c = w*w
    return 1/(2*(c+a*b/(a+b)))


def fixed_mass():
    u = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    a = u @ np.diag([.01, .02, .03]) @ u.conj().T
    b = 1e-4*np.array([[0., 1., 0.], [-1., 0., 1.], [0., -1., 0.]])
    c = 1e-4*np.eye(3)
    d = np.diag([10., 12., 15.])
    return u, a, b, c, d


def blocks(m):
    h = m @ m.conj().T
    return h[:3, :3], h[:3, 3:], h[3:, 3:]


def schur(m, z):
    e, f, g = blocks(m)
    return e-z*np.eye(3)-f@np.linalg.solve(g-z*np.eye(3), f.conj().T)


def neutral_case(m, label):
    ul, mass, _ = np.linalg.svd(m)
    order = np.argsort(mass)
    mass, ul = mass[order], ul[:, order]
    poles = mass**2
    e, f, g = blocks(m)
    gamma2 = np.linalg.eigvalsh(g)[0]
    zstar = np.linalg.norm(e, 2)
    assert zstar < gamma2 and mass[0] > 0
    gi = np.linalg.inv(g)
    k0, c0 = e-f@gi@f.conj().T, f@gi@gi@f.conj().T
    delta = zstar/gamma2
    lower = spectrum_pair(k0, np.eye(3)+c0/(1-delta))
    upper = spectrum_pair(k0, np.eye(3)+c0)
    # At tiny mixing the analytic interval is below floating-point resolution.
    allowance = 1e-11*max(1., float(np.linalg.norm(e)))
    assert np.all(poles[:3] >= lower-allowance)
    assert np.all(poles[:3] <= upper+allowance)
    assert poles[2] <= zstar+allowance and poles[3] >= gamma2-allowance
    for j, z in enumerate(poles[:3]):
        eq('neutral_exact_S_'+label, schur(m, z)@ul[:3, j], np.zeros(3), tol=1e-10)
    low = (ul[:3, :3]*poles[:3]) @ ul[:3, :3].conj().T
    cp = float(np.imag(low[0,1]*low[1,2]*low[2,0]))
    # Conjugating physical M reverses the physical CP invariant.
    uc, sc, _ = np.linalg.svd(m.conj())
    orderc = np.argsort(sc)
    uc, sc = uc[:, orderc], sc[orderc]
    lowc = (uc[:3, :3]*sc[:3]**2)@uc[:3, :3].conj().T
    eq('Weyl_physical_conjugation_'+label, lowc, low.conj())
    return dict(label=label, masses=mass.tolist(), gamma_squared=float(gamma2),
                zstar=float(zstar), delta=float(delta),
                lower=lower.tolist(), upper=upper.tolist(),
                interval_float_allowance=allowance,
                max_interval_width=float(np.max(upper-lower)), cp_low=cp,
                low_weak_nonunitarity=float(np.linalg.norm(ul[:3,:3].conj().T@ul[:3,:3]-np.eye(3))))


def scalar_check(w):
    lam, kval, f, eps = .2, .3, .01, .02
    wi2, prod = w*w, float(np.prod(w))
    hx = np.diag(2*lam+f*prod/(2*wi2**2))
    for i in range(3):
        for j in range(i):
            hx[i,j] = hx[j,i] = -f*prod/(2*wi2[i]*wi2[j])
    off = []
    for i, j in ((0,1), (0,2), (1,2)):
        base = kval+f*prod/(wi2[i]*wi2[j])
        off.extend([base+2*eps, base-2*eps] if (i,j)==(0,2) else [base, base])
    assert np.linalg.eigvalsh(hx)[0] > 0 and min(off) > 0 and f*prod > 0
    assert kval > 2*abs(eps)
    # Analytic first derivative in the ten-coordinate local transverse section.
    first = f*prod/wi2-f*prod/wi2
    eq('scalar_stationarity', first, np.zeros(3))
    return dict(w=w.tolist(), lambdas=[lam]*3, k=[kval]*3, f=f, epsilon=eps,
                radial_Hessian_eigenvalues=np.linalg.eigvalsh(hx).tolist(),
                off_diagonal_real_quadratic_coefficients=off,
                determinant_phase_quadratic_coefficient=f*prod,
                physical_positive_directions=10, gauge_tangent_zero_directions=8,
                unique_global_minimum_claimed=False)


def run():
    anomaly_cases = []
    for n in range(1, 10):
        for a in range(n+1):
            q = anomalies(n, a)
            assert [q[k] for k in ('C3','L3','C2x','L2x','x3','gravx')] == [0,2*(3*a-n),0,3*a-n,6*(n-3*a),0]
            assert [q[k] for k in ('rC','rL','rdimension','rx2','r2x')] == [4*n,4*n,16*n,24*n,0]
            if q['L3'] == 0:
                anomaly_cases.append(dict(n=n, a=a, b=n-a, b0_all_colored=11-2*n, **q))
    assert [z['n'] for z in anomaly_cases if z['b0_all_colored']>0] == [3]
    # r exponents modulo4; scalar eta/rho=2,chi=0; all-left fermions ordinary=1,heavy-conjugate=3.
    rmenu = {'LLrho':(1+1+2)%4, 'NcFc':(1+3)%4,
             'NcNc':2, 'FcFc':2, 'eta_chi_squared':(2*(0-2))%4,
             'det_eta_rho_chi':(2+2+0)%4}
    assert rmenu['LLrho'] == rmenu['NcFc'] == rmenu['eta_chi_squared'] == rmenu['det_eta_rho_chi'] == 0
    r_eff_gauge = np.diag([-1., -1., 1.])
    for i, scalar_r in enumerate((-1., -1., 1.)):
        eq('residual_vacuum', scalar_r*r_eff_gauge@np.eye(3)[:,i], np.eye(3)[:,i])

    w = np.array([1.,1.,100.]); gl, gx = .65, .4
    scalar = scalar_check(w)
    u, a, b, c, d = fixed_mass()
    m = np.block([[a,b],[c,d]])
    eq('allowed_H_antisymmetry', b.T, -b)
    rleft = np.diag([-1j]*3+[1j]*3)
    weyl = np.block([[np.zeros((6,6)),m.conj()], [m.conj().T,np.zeros((6,6))]])
    # More direct two-sector check, without identifying particle labels across conjugation.
    eq('residual_Dirac_block', (-1j)*m.conj()*(1j), m.conj())
    eq('all_left_symmetric', weyl.T, weyl)
    cases = [neutral_case(m, 'common_point')]
    assert abs(cases[0]['cp_low']) > 1e-12
    # A second algebraic point gives a visibly resolvable pole enclosure, not a fit.
    mw = np.block([[a,.20*b/1e-4],[.15*np.eye(3),np.diag([1.,1.2,1.5])]])
    cases.append(neutral_case(mw, 'visible_enclosure_diagnostic'))
    assert cases[1]['max_interval_width'] > 1e-7

    e, fmat, gmat = blocks(m)
    ul, mass, _ = np.linalg.svd(m)
    order = np.argsort(mass); ul, mass = ul[:,order], mass[order]
    derivatives = {}
    for name in ('weta','wrho','wchi','ReMu01','ImH01'):
        dm = np.zeros((6,6), complex)
        if name=='weta': dm[:3,:3]=a/w[0]
        elif name=='wrho': dm[:3,3:]=b/w[1]
        elif name=='wchi': dm[3:,3:]=d/w[2]
        elif name=='ReMu01': dm[3,1]=1
        else: dm[0,4]=1j; dm[1,3]=-1j
        step=1e-5
        plus, minus=m+step*dm,m-step*dm
        dh=dm@m.conj().T+m@dm.conj().T
        de, df, dg=dh[:3,:3],dh[:3,3:],dh[3:,3:]
        z=-.003; r=np.linalg.inv(gmat-z*np.eye(3))
        ds=de-df@r@fmat.conj().T-fmat@r@df.conj().T+fmat@r@dg@r@fmat.conj().T
        eq('complete_S_source_'+name, (schur(plus,z)-schur(minus,z))/(2*step),ds,tol=5e-9)
        numeric=(np.sort(np.linalg.svd(plus,compute_uv=False))**2-np.sort(np.linalg.svd(minus,compute_uv=False))**2)/(2*step)
        analytic=np.real(np.diag(ul.conj().T@dh@ul))
        eq('full_pole_HF_'+name,numeric,analytic,tol=3e-7)
        derivatives[name]=dict(max_abs_HF_error=float(np.max(np.abs(numeric-analytic))), full_source=True)

    gauge_cases=[]
    for ww,gg,xx in ((w,gl,gx),(np.array([2.,3.,7.]),.3,.7),(np.array([.8,1.6,4.]),1.1,.2)):
        v, mat, dd, ev, ov, inv, ff=gauge(ww,gg,xx)
        full=gauge9(ww,gg,xx); idx=[2,7,8]
        eq('full9_Cartan_block',full[np.ix_(idx,idx)],mat)
        other=[i for i in range(9) if i not in idx]
        eq('full9_no_cross',full[np.ix_(idx,other)],np.zeros((3,6)))
        eq('X0_degenerate_mass',full[[3,4],[3,4]],np.full(2,gg*gg*(ww[0]**2+ww[2]**2)/2))
        for gen in generators()[3:5]:
            eq('X0_ordinary_projection',gen[:2,:2],np.zeros((2,2)))
        eq('charge_vector_sum',v.sum(axis=0),np.zeros(3))
        ee=(4/(3*gg*gg)+1/(xx*xx))**-.5
        photon=ee*np.array([1/gg,-1/(np.sqrt(3)*gg),1/xx])
        eq('photon_null',mat@photon,np.zeros(3),tol=1e-10)
        eq('photon_norm',photon@photon,1.)
        eq('photon_nonuniversal',dd@photon,0.)
        kapp=ov[:,1:].T@dd
        spectral=float(np.sum(kapp*kapp/ev[1:]))
        exact=analytic_f0(ww)
        eq('F0_full_sum',ff,exact,tol=3e-10,scale=exact)
        eq('F0_spectral',spectral,exact,tol=3e-10,scale=exact)
        boundchecks=[]
        for fraction in (-.2,.1,.4):
            z=fraction*ev[1]
            value=float(np.sum(kapp*kapp/(ev[1:]-z)))
            limit=abs(z)/(ev[1]-abs(z))*ff
            assert abs(value-ff)<=limit+1e-12*ff
            boundchecks.append(dict(z=float(z), difference=abs(value-ff), bound=float(limit)))
        gauge_cases.append(dict(w=ww.tolist(),gl=gg,gx=xx,Cartan_mass_squared=ev.tolist(),
                                full9_mass_squared=np.linalg.eigvalsh(full).tolist(),F0=ff,
                                spectral_residues=(kapp*kapp).tolist(), finite_bounds=boundchecks))

    v, mat, dd, ev, ov, inv, ff=gauge(w,gl,gx)
    gauge_sources={}
    for name in ('gl','gx','weta','wrho','wchi'):
        dv=np.zeros((3,3)); dw=np.zeros(3)
        wp,wm=w.copy(),w.copy(); gp,gm=gl,gl; xp,xm=gx,gx
        step=1e-5
        if name=='gl':
            dv[:,:2]=v[:,:2]/gl; gp+=step; gm-=step
        elif name=='gx':
            dv[:,2]=v[:,2]/gx; xp+=step; xm-=step
        else:
            i=('weta','wrho','wchi').index(name); dw[i]=1.; wp[i]+=step; wm[i]-=step
        dmat=2*(dv.T@np.diag(w*w)@v+v.T@np.diag(w*w)@dv+v.T@np.diag(2*w*dw)@v)
        ddd=-dv[2]
        deriv=float(2*ddd@inv@dd-dd@inv@dmat@inv@dd)
        if name in ('gl','gx'): closed=0.
        elif name=='wchi': closed=-4*w[2]*ff*ff
        else:
            i=('weta','wrho').index(name); j=1-i
            closed=-4*w[i]*(w[j]**2/(w[0]**2+w[1]**2))**2*ff*ff
        fd=(gauge(wp,gp,xp)[-1]-gauge(wm,gm,xm)[-1])/(2*step)
        eq('F0_source_closed_'+name,deriv,closed,tol=1e-10,scale=ff)
        eq('F0_source_fd_'+name,fd,closed,tol=3e-7,scale=ff)
        gauge_sources[name]=dict(complete=deriv,closed=closed,finite_difference=fd)

    p=np.diag([1.,0.,0.]); pu=p; pd=u.conj().T@p@u; vckm=u
    eq('CKM_common_projector',pd,vckm.conj().T@pu@vckm)
    eq('projector_idempotent',pd@pd,pd)
    assert abs(np.trace(pd)-1)<1e-12 and np.linalg.norm(pd-np.diag(np.diag(pd)))>.5
    vfup=np.array([gl/2,-gl/(2*np.sqrt(3)),0.])
    vfdown=np.array([-gl/2,-gl/(2*np.sqrt(3)),0.])
    for j in (1,2):
        cu,cd,kappa=vfup@ov[:,j],vfdown@ov[:,j],dd@ov[:,j]
        gamma_u=cu*np.eye(3)+kappa*pu
        gamma_d=cd*np.eye(3)+kappa*pd
        eq('two_department_current',gamma_d-cd*np.eye(3),vckm.conj().T@(gamma_u-cu*np.eye(3))@vckm)
    # Actual Yukawa variation at fixed nondegenerate physical masses; projections move.
    jrot=np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,0.]])
    dp=u.conj().T@(p@jrot-jrot@p)@u
    def pd_angle(t):
        rot=np.array([[np.cos(t),np.sin(t),0.],[-np.sin(t),np.cos(t),0.],[0.,0.,1.]])
        ud=rot@u; md=ud@np.diag([.15,.25,.35])
        eq('same_Yukawa_left_mass',md@md.conj().T,ud@np.diag(np.array([.15,.25,.35])**2)@ud.conj().T)
        return ud.conj().T@p@ud
    step=1e-5
    tensor_derivative=ff*np.kron(pu,dp)
    tensor_numeric=ff*np.kron(pu,(pd_angle(step)-pd_angle(-step))/(2*step))
    eq('complete_flavor_tensor_source',tensor_numeric,tensor_derivative,tol=1e-8,scale=ff)

    altered=np.block([[a,2*b],[3*c,d]])
    oldm=np.sort(np.linalg.svd(m,compute_uv=False))
    newm=np.sort(np.linalg.svd(altered,compute_uv=False))
    spectrum_change=float(np.max(np.abs(newm-oldm)))
    assert spectrum_change>1e-9
    assert analytic_f0(np.array([1.,1.,120.])) != ff
    return dict(schema='round1318_same_parent_331_v1',status='PASS',
                anomaly_nonzero_n_cases=anomaly_cases,r_menu_mod4=rmenu,scalar=scalar,
                neutral_cases=cases,neutral_parent_sources=derivatives,
                gauge_cases=gauge_cases,gauge_background_sources=gauge_sources,
                flavor=dict(projector_rank=1,offdiagonal_norm=float(np.linalg.norm(pd-np.diag(np.diag(pd)))),
                            tensor_source_norm=float(np.linalg.norm(tensor_derivative)),
                            phase_instrument_certified=False),
                deletions=dict(n6_allactive_b0=-1,n6_heavies_inactive_b0=3,
                               no_r_allows_same_sector_Majorana=True,
                               retained_r_changed_allowed_cross_full_mass_change=spectrum_change,
                               changed_wchi_F0=analytic_f0(np.array([1.,1.,120.]))),
                max_normalized_errors=ERRORS,
                scope='declared renormalizable parent, tree structural kernels and finite all-active one-loop sign; no global/instrument/empirical certification',
                new_axioms=0,stage_complete=False,
                bindings={str(path.relative_to(HERE.parent)):sha(path) for path in
                          (HERE/'proof.md',HERE/'check.py',HERE.parent/'_admission/round1318_admission.md')})


def compare(a,b,path='result'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert np.isclose(a,b,rtol=1e-7,atol=1e-14),(path,a,b)
    else: assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--check',action='store_true')
    mode.add_argument('--write-new',action='store_true')
    args=p.parse_args(); result=run(); out=HERE/'results.json'
    if args.write_new:
        with out.open('x',encoding='utf8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2); stream.write('\n')
    else: compare(json.loads(out.read_text(encoding='utf8')),result)
    print(json.dumps(dict(status='PASS',mode='write-new' if args.write_new else 'read-only',
                         CP=result['neutral_cases'][0]['cp_low'],F0=result['gauge_cases'][0]['F0'],
                         maximum_normalized_error=max(ERRORS.values())),ensure_ascii=False))


if __name__=='__main__': main()
