"""573: material coordinates on the exact common-source Einstein initial branch.

Analytic barriers and a maximum principle prove local rank; collocation checks
are diagnostics. The last two time derivatives cancel from the determinant.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gauss_einstein_initial_data as geo

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gravity_material_coordinates_results.json'


def rational_bounds():
    """Conservative rational bounds; no sampled maximum is used."""
    hmin=Q('.95')*Q('.65');hmax=Q('1.05')*Q('.67');smax=Q('1.06')*Q('.55')
    amin=Q('.05')*Q('.65');bmin=Q('.06')*Q('.54');bmax=Q('.06')*Q('.55')
    amax=Q('.05')*Q('.67');kap=Q('.0056')/amin;eta=Q('.0056')/bmin
    u=Q('.38');aw=Q('.41');bw=Q('.29');e=Q('.17')
    r2w=2*e*u*(1+u/2)
    r3w=aw*e*(2*u-3*u**3/4+3*u*u/2+u**3/4)
    q0w=6*r3w
    e0w=Q('.13')+6*aw*e*(2*u-3*u**3/4+3*u*u/4+u**3/12)
    mx=Q('.0056')+q0w*Q('.03')+Q('.07')*Q('.11')
    my=Q('.0056')+Q('.025')*bmax+Q('.0056')+r2w*Q('.348')+q0w*Q('.08')+e0w*Q('.11')
    shear=Q('3.16')*(mx+my)
    p2=kap**2+(Q('.025')+eta)**2+4*((2*e*u*(1+u))**2+r3w**2)/hmin**2
    Amax=shear**2+2*p2
    Fmin=2-(hmax*hmax+smax*smax)/6
    ax=aw*(1+u);ay=bw*Q('1.2');az=Q('.12')
    D2=2*amax**2+bmax**2+hmax**2*((ax/2+Q('.09'))**2+(ay/2+Q('.24'))**2+(az/2)**2)
    radial2=(hmax*amax)**2+(hmax*amax+smax*bmax)**2
    Bmax=D2/Fmin+radial2/(6*Fmin**2)
    ex=e*(2*u+3*u*u/2)
    Ymax=Q('.421')*(ex**2+Q('.11')**2+Q('.09')**2)/2
    Ymax+=Q('.129')*(e0w**2+Q('.07')**2)/72
    Ymax+=(ax**2*ay**2+az**2*(ax**2+ay**2))/(2*Q('.420'))
    Ymax+=18*Q('.11')**2/Q('.128')
    upper=Q('1.5')
    reaction=2*upper**5-upper-2*upper**-7-4*upper**-3
    Rmin=3*Q('1.8')**3*Q('.97')*Q('.14')**2/(2*bmax**2)
    assert Fmin>Q('1.8') and Amax<2 and Bmax<1 and Ymax<2
    assert r2w<Q('.155') and q0w<Q('.40') and e0w<Q('.48')
    assert reaction>0 and Rmin>150>upper**8
    values=dict(F_lower=Fmin,A_upper=Amax,B_upper=Bmax,Y_upper=Ymax,
                momentum_Wiener_bound=mx+my,shear_bound=shear,
                supersolution_reaction_lower=reaction,R_lower=Rmin,psi_eighth_upper=upper**8)
    return {k:dict(exact=str(v),decimal=float(v)) for k,v in values.items()}


@lru_cache(maxsize=3)
def fields(N):
    q=geo.make_source(N);psi,tensor,stats=geo.solve_hamiltonian(q)
    phi,p,F=q['phi'],q['p'],q['F']
    v=psi[...,None]**-6*F[...,None]*(p-phi*np.sum(phi*p,axis=-1)[...,None]/12)
    x,y,z=np.moveaxis(q['grid'],-1,0)
    hs=np.sqrt(geo.old.PAR['h2']);_,u,_=geo.old.scalar.parameters();ss=np.sqrt(u[1])
    dh=np.zeros(q['grid'].shape);ds=np.zeros_like(dh)
    dh[...,0]=dh[...,1]=.05*hs*np.cos(x+y);ds[...,1]=-.06*ss*np.sin(y)
    vh,vs=v[...,1],v[...,4]
    rh=psi**-4*np.sum(dh*dh,axis=-1)-vh*vh
    rs=psi**-4*np.sum(ds*ds,axis=-1)-vs*vs
    grh=np.stack([geo.derivative(rh,i) for i in range(3)],axis=-1)
    grs=np.stack([geo.derivative(rs,i) for i in range(3)],axis=-1)
    psiz=geo.derivative(psi,2)
    return dict(q=q,psi=psi,tensor=tensor,stats=stats,v=v,dh=dh,ds=ds,
                rh=rh,rs=rs,grh=grh,grs=grs,psiz=psiz,hstar=hs,sstar=ss)


def point_coefficients(f):
    H=1.05*f['hstar'];S=f['sstar'];b=.06*S
    kap=.0056/(.05*f['hstar']);eta=.0056/b;p=.025-eta
    F=2-(H*H+S*S)/6;alpha=1-H*H/12;beta=1-S*S/12;d=H*S/12
    Ch=-F*d*p;Cs=F*beta*p;Chx=F*alpha*kap;Csx=-F*d*kap
    R=3*F**3*beta*p*p/(2*alpha*b*b)
    return dict(H=H,S=S,b=b,kap=kap,eta=eta,p=p,F=F,alpha=alpha,beta=beta,
                Ch=Ch,Cs=Cs,Chx=Chx,Csx=Csx,R=R)


def coefficient_bounds_check():
    hs=np.sqrt(geo.old.PAR['h2']);_,u,_=geo.old.scalar.parameters();ss=np.sqrt(u[1])
    gw2=2*geo.old.PAR['b'][1];gy2=72*geo.old.PAR['b'][2]
    assert .65<hs<.67 and .54<ss<.55 and .420<gw2<.421 and .128<gy2<.129
    bounds=rational_bounds()
    return dict(hstar=hs,sstar=ss,gw_squared=gw2,gY_squared=gy2,rational_bounds=bounds)


def source_structure_check():
    f=fields(16);q=f['q'];grid=q['grid'];x,y,z=np.moveaxis(grid,-1,0)
    cp=point_coefficients(f);ph=-cp['kap']*np.cos(x+y);ps=.025*np.cos(x)-cp['eta']*np.sin(y)
    err=max(float(np.max(abs(q['PX'][...,1].real-ph))),float(np.max(abs(q['Ps']-ps))))
    assert err<1e-14
    assert np.max(abs(q['mom'][...,2]))<1e-14
    assert np.max(abs(q['mom']-q['mom'][:,:,0:1,:]))<1e-14
    A=np.sum(f['tensor']**2,axis=(-1,-2))+q['pKp']
    assert np.max(abs(A-A[:,:,0:1]))<1e-13
    bcoef=.12**2*q['f']['h']**2/(4*q['F'])
    ycoef=.12**2*(q['f']['a'][...,0,0]**2+q['f']['a'][...,1,1]**2)/(2*(2*geo.old.PAR['b'][1]))
    bbase=q['B']-bcoef*np.cos(z)**2;ybase=q['Y']-ycoef*np.cos(z)**2
    polyerr=max(float(np.max(abs(bbase-bbase[:,:,0:1]))),float(np.max(abs(ybase-ybase[:,:,0:1]))))
    assert polyerr<1e-13
    assert np.max(A)<2 and np.max(q['B'])<1 and np.max(q['Y'])<2
    return dict(radial_formula_error=err,z_structure_error=polyerr,
                Mz_max=float(np.max(abs(q['mom'][...,2]))),b_min=float(np.min(bcoef)),y_min=float(np.min(ycoef)))


def monotonic_geometry_check():
    rows=[]
    for N in (16,24):
        f=fields(N);psi=f['psi'];q=f['q'];idx=(-np.arange(N))%N
        reflected=np.take(psi,idx,axis=2)
        sym=max(np.max(abs(psi-reflected)),np.max(abs(psi-np.roll(psi,N//2,axis=2))))
        interior=f['psiz'][:,:,1:N//4]
        assert sym<1e-12 and np.max(interior)<0 and np.max(psi)<1.5
        A=np.sum(f['tensor']**2,axis=(-1,-2))+q['pKp']
        weight=4*q['C']*psi**4+8*A*psi**-8+8*q['Y']*psi**-4
        assert np.min(weight)>0
        rows.append(dict(N=N,symmetry_error=float(sym),interior_max_dpsi_z=float(np.max(interior)),
                         psi_max=float(np.max(psi)),weighted_zero_order_min=float(np.min(weight))))
    return dict(rows=rows,strict_continuum_sign_from_maximum_principle_not_sampling=True)


def reference_minor_check():
    rows=[]
    for N in (16,24,32):
        f=fields(N);c=point_coefficients(f);i=(0,N//4,N//8)
        psi=float(f['psi'][i]);psiz=float(f['psiz'][i]);vh=c['Ch']*psi**-6
        ax,az=f['grh'][i][[0,2]];bx,bz=f['grs'][i][[0,2]]
        minor=-vh*(-c['b'])*(ax*bz-az*bx)
        bracket=c['Chx']*c['b']**2*(c['R']-psi**8)
        analytic=8*c['Ch']**2*(-c['b'])*psi**-31*psiz*bracket
        relative=abs(minor-analytic)/abs(analytic)
        assert analytic>0 and minor>0 and relative<.012
        assert abs(f['v'][i][1]-vh)<1e-14
        rows.append(dict(N=N,psi=psi,dpsi_z=psiz,normal_h=vh,clock_norm=-vh**2,
                         exact_formula_on_collocation=analytic,independent_spatial_minor=float(minor),
                         relative_derivative_discrepancy=float(relative),R=c['R']))
    assert rows[-1]['relative_derivative_discrepancy']<1e-5
    return dict(point=[0.,float(np.pi/2),float(np.pi/4)],rows=rows,
                last_two_time_derivatives_cancel_exactly=True)


def determinant_algebra_check():
    # Exact rational check with arbitrary accelerations/unused time derivatives.
    from itertools import permutations
    def det(m):
        total=Q(0)
        for perm in permutations(range(4)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
            term=Q(sign)
            for i,j in enumerate(perm):term*=m[i][j]
            total+=term
        return total
    vh,vs,b=Q(2,7),Q(-3,5),Q(1,4)
    ax,ay,az=Q(2,3),Q(-2,9),Q(5,7);bx,by,bz=Q(-1,3),Q(2,5),Q(7,8)
    expected=vh*b*(ax*bz-az*bx)
    for at,bt in ((Q(0),Q(0)),(Q(71),Q(-91)),(Q(2,3),Q(8,11))):
        m=[[vh,0,0,0],[vs,0,-b,0],[at,ax,ay,az],[bt,bx,by,bz]]
        assert det(m)==expected
    return dict(exact_determinant=str(expected),arbitrary_last_time_entries=3)


def chart_boundary_and_frame_check():
    f=fields(24);N=24
    wall=max(float(np.max(abs(f['psiz'][:,:,0]))),float(np.max(abs(f['psiz'][:,:,N//4]))))
    assert wall<1e-12
    i=(0,N//4,N//8);mirror=(0,N//4,(-N//8)%N)
    menu=np.array([f['q']['f']['h'][i],f['q']['f']['s'][i],f['rh'][i],f['rs'][i]])
    reflected=np.array([f['q']['f']['h'][mirror],f['q']['f']['s'][mirror],f['rh'][mirror],f['rs'][mirror]])
    assert np.max(abs(menu-reflected))<1e-12
    c=point_coefficients(f);F=c['F']
    transform=np.eye(4);transform[2:,:2]=np.outer(menu[2:],[-c['H']/3,-c['S']/3])
    transform[2:,2:]=F*np.eye(2)
    assert abs(np.linalg.det(transform)-F*F)<1e-13
    return dict(wall_dpsi_max=wall,mirror_menu_error=float(np.max(abs(menu-reflected))),
                Jordan_menu_determinant_factor=float(np.linalg.det(transform)),
                same_menu_not_a_global_chart=True)


def run():
    checks=('coefficient_bounds_check','source_structure_check','monotonic_geometry_check',
            'reference_minor_check','determinant_algebra_check','chart_boundary_and_frame_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_gauss_einstein_initial_data.py','joint_gauss_einstein_initial_data_results.json',
          'joint_gauss_continuum_sampling.py','research_round_572_checks.json')
    return dict(round=573,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='same changed source and prescribed Einstein-Yang-Mills-sigma model; analytic local material chart, no added source or field; not global chart, actual apparatus, dimension or GR generation')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
