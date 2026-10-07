"""919: weak Gauss transport of the original finite Fourier/Hermite fields.
N is the fixed evolution grid; M is solely the integration grid. No continuum
error certificate is inferred from the two M values or selected weak moments.
"""
from pathlib import Path
import sys,json,hashlib,time,argparse,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'918'))
import euler_density_contacts as density
c=density.c;stress=density.stress;mb=c.aj.mb
TARGET=HERE/'weak_gauss_transport_results.json'

class GridField:
    """Evaluate exactly the saved finite interpolant on a new quadrature grid."""
    def __init__(self,field,t,M,spinor=False):
        self.field=field;self.t=t;self.M=M;self.spinor=spinor;self.data={}
        k=np.rint(field.k).astype(int)
        assert np.max(abs(field.k-k))<1e-13 and np.max(abs(k))<(M/2)
        self.ids=tuple((k[:,i]%M) for i in range(3))
    def get(self,axes=()):
        key=tuple(sorted(axes))
        if key not in self.data:
            field=self.field;co=field.co if self.spinor else field.c
            ntime=key.count(0);powers=np.zeros(4)
            for n in range(ntime,4):
                powers[n]=math.factorial(n)/math.factorial(n-ntime)*(self.t/field.T)**(n-ntime)/field.T**ntime
            v=np.einsum('p,kpf->kf',powers,co)
            for axis in key:
                if axis:v=v*(1j*field.k[:,axis-1,None])
            grid=np.zeros((self.M,self.M,self.M,v.shape[-1]),complex)
            grid[self.ids]=v
            value=np.fft.ifftn(grid,axes=(0,1,2)).reshape((-1,v.shape[-1]))*self.M**3
            self.data[key]=value if self.spinor else value.real
        return self.data[key]
    def jets(self,lo,hi):
        g,phi,A=mb.unpack(self.get()[lo:hi])
        dg,dphi,dA=mb.unpack(np.stack([self.get((m,))[lo:hi] for m in range(4)],axis=1))
        ddg,ddphi,ddA=mb.unpack(np.stack([np.stack([self.get((m,n))[lo:hi] for n in range(4)],axis=1) for m in range(4)],axis=1))
        return dict(g=g,phi=phi,A=A,dg=dg,dphi=dphi,dA=dA,ddg=ddg,ddphi=ddphi,ddA=ddA)

def grid_points(t,M,ids):
    idx=np.array(np.unravel_index(np.asarray(ids),(M,M,M))).T
    return np.column_stack((np.full(len(ids),t),2*np.pi*idx/M))

def validate_adapter(bg,field,u):
    errors={}
    for M in (17,25):
        ids=np.array([0,1,M+3,M**2+M+1,M**3//2,M**3-1]);t=.173*bg.T
        x=grid_points(t,M,ids)
        for name,f in (('background',bg),('response',field),('receiver',u)):
            gf=GridField(f,t,M,name=='receiver');diff=[]
            axes=((),(0,),(2,)) if name=='receiver' else ((),(0,),(2,),(0,0),(0,2),(1,3))
            for a in axes:
                direct=f.evaluate(x,None if not a else a[0]) if name=='receiver' else c.aj.previous.derivative(f,x,a)
                diff.append(c.maximum(gf.get(a)[ids]-direct))
            errors[f'{name}_M{M}']=max(diff)
    assert max(errors.values())<2e-10,errors
    return errors

def row(bg,field,u,par,t,M,batch=256):
    bggrid=GridField(bg,t,M);vgrid=GridField(field,t,M);ugrid=GridField(u,t,M,True)
    names=('test_time','spatial_derivative','spatial_connection','scalar_exchange','background_contact')
    sums={k:np.zeros(3) for k in ('moment',)+names}
    absolute={k:np.zeros(3) for k in ('moment',)+names}
    maxima={'Gauss':0.,'internal_contact':0.,'scalar_calibration':0.}
    initial_gauss_identity=0.;covariant_pairing_error=0.
    for lo in range(0,M**3,batch):
        hi=min(lo+batch,M**3);p=bggrid.jets(lo,hi);v=vgrid.jets(lo,hi)
        base=c.euler(p,par);z=base['z'];J0=density.density(base,z);h=1e-24
        shifted=c.euler({k:p[k].astype(complex)+1j*h*v[k] for k in p},par)
        Q={k:a.imag/h for k,a in density.density(shifted,shifted['z']).items()}
        # Neutral receiver has no direct internal current: Q_YM,Q_phi are J'.
        C=sum(c.bracket(v['A'][:,mu,:],J0['YM'][:,mu,:]) for mu in range(4))
        C+=np.einsum('bA,aAB,bB->ba',J0['scalar'],c.REP,v['phi'])
        scalar=np.einsum('bA,aAB,bB->ba',Q['scalar'],c.REP,p['phi'])
        uu=ugrid.get()[lo:hi];du=np.stack([ugrid.get((mu,))[lo:hi] for mu in range(4)],axis=1)
        ell=.5*np.einsum('bij,bmji->bm',z['invgamma'],p['dg'][:,:,1:,1:])
        bil=lambda a,b:np.einsum('bi,ij,bj->b',a.conj(),stress.r.BETA,b)
        a=bil(uu,uu).real/z['vol']
        da=np.stack([2*bil(du[:,mu],uu).real/z['vol']-ell[:,mu]*a for mu in range(4)],axis=1)
        F=base['F'];A=p['A']
        dF=p['ddA']-p['ddA'].swapaxes(2,3)+c.bracket(p['dA'][:,:,:,None,:],A[:,None,None,:,:])+c.bracket(A[:,None,:,None,:],p['dA'][:,:,None,:,:])
        values={k:[] for k in ('moment',)+names}
        for i,j in ((1,2),(1,3),(2,3)):
            theta=a[:,None]*F[:,i,j,:]
            dtheta=da[:,:,None]*F[:,None,i,j,:]+a[:,None,None]*dF[:,:,i,j,:]
            values['moment'].append(np.sum(theta*Q['YM'][:,0,:],axis=1))
            values['test_time'].append(np.sum(dtheta[:,0,:]*Q['YM'][:,0,:],axis=1))
            values['spatial_derivative'].append(np.sum(dtheta[:,1:,:]*Q['YM'][:,1:,:],axis=(1,2)))
            connection=sum(c.bracket(A[:,mu,:],Q['YM'][:,mu,:]) for mu in range(1,4))
            values['spatial_connection'].append(-np.sum(theta*connection,axis=1))
            alternative=sum(np.sum(c.bracket(A[:,mu,:],theta)*Q['YM'][:,mu,:],axis=1) for mu in range(1,4))
            covariant_pairing_error=max(covariant_pairing_error,c.maximum(alternative-values['spatial_connection'][-1]))
            values['scalar_exchange'].append(-np.sum(theta*scalar,axis=1))
            values['background_contact'].append(-np.sum(theta*C,axis=1))
        for k in values:
            array=np.stack(values[k],axis=1);sums[k]+=np.sum(array,axis=0);absolute[k]+=np.sum(abs(array),axis=0)
        direct=z['mu'][:,None]*(shifted['YM'][:,0,:].imag/h)+(shifted['z']['mu'].imag/h)[:,None]*base['YM'][:,0,:]
        initial_gauss_identity=max(initial_gauss_identity,c.maximum(direct-Q['YM'][:,0,:]))
        maxima['Gauss']=max(maxima['Gauss'],c.maximum(Q['YM'][:,0,:]))
        maxima['internal_contact']=max(maxima['internal_contact'],c.maximum(C))
        maxima['scalar_calibration']=max(maxima['scalar_calibration'],c.maximum(a))
    factor=(2*np.pi)**3/M**3
    sums={k:factor*v for k,v in sums.items()};absolute={k:factor*v for k,v in absolute.items()}
    rhs=sum(sums[k] for k in names)
    out=dict(M=M,time=t,moment=sums.pop('moment').tolist(),rhs=rhs.tolist(),contributions={k:v.tolist() for k,v in sums.items()},
      absolute_integrand_quadrature_not_bound={k:v.tolist() for k,v in absolute.items()},sample_maxima=maxima,
      same_Gauss_density_error=initial_gauss_identity,ad_invariant_pairing_error=covariant_pairing_error)
    assert initial_gauss_identity<1e-16 and covariant_pairing_error<1e-15
    return out

def run():
    start=time.time();par=c.parameters();bg,field,initial=c.ex.previous.rebuilt.build_pair()
    u,_,_=stress.weighted.receiver_modes();adapter=validate_adapter(bg,field,u);rows=[];transport=[]
    for M in (17,25):
        current=[]
        for t in (-bg.T,0.,bg.T):
            rr=row(bg,field,u,par,t,M);rows.append(rr);current.append(rr)
            print(json.dumps(dict(M=M,t=t,moment=rr['moment'],rhs=rr['rhs'],elapsed=round(time.time()-start,2))),flush=True)
        first,mid,last=current;dt=2*bg.T
        integrate=lambda key:dt/6*(np.array(first[key])+4*np.array(mid[key])+np.array(last[key]))
        change=np.array(last['moment'])-np.array(first['moment']);driver=integrate('rhs')
        terms={k:(dt/6*(np.array(first['contributions'][k])+4*np.array(mid['contributions'][k])+np.array(last['contributions'][k]))).tolist() for k in mid['contributions']}
        defect=change-driver
        drop_initial=np.array(last['moment'])-driver
        drop_spatial=defect+np.array(terms['spatial_derivative'])+np.array(terms['spatial_connection'])
        transport.append(dict(M=M,actual_initial=first['moment'],actual_final=last['moment'],actual_change=change.tolist(),Simpson_rhs=driver.tolist(),
          Simpson_contributions=terms,weak_balance_defect=defect.tolist(),drop_initial_defect=drop_initial.tolist(),drop_spatial_terms_defect=drop_spatial.tolist(),
          quadrature_and_roundoff_error_certified=False))
    differences={k:(np.array(transport[1][k])-np.array(transport[0][k])).tolist() for k in ('actual_initial','actual_final','actual_change','Simpson_rhs','weak_balance_defect')}
    return dict(round=919,date='2026-10-06',physical_interpolant_N=17,time_half_window=bg.T,quadrature_sizes=[17,25],
      tests='theta_ij^a=(bar psi psi) F_ij^a for spatial pairs12,13,23; original receiver and connection',
      same_background_response_and_receiver=True,new_physical_preparation=False,finite_constraint_tests_not_complete_observable_basis=True,
      adapter_direct_Fourier_checks=adapter,rows=rows,transport=transport,quadrature_differences_not_error_bounds=differences,
      actual_finite_observable_error_certified=False,all_constraints_certified=False,full872_response_computed=False,full_goal_completed=False,
      elapsed_seconds=round(time.time()-start,3),source_hashes={str(q.relative_to(STAGE)):hashlib.sha256(q.read_bytes()).hexdigest() for q in
       (Path(__file__),STAGE/'918/euler_density_contacts.py',STAGE/'916/covariant_joint_residual.py',STAGE/'914/weighted_receiver_pairing.py',STAGE/'913/independent_record_rebuild.py')})
if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--write',action='store_true');a=pa.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'saved':a.write,'elapsed':out['elapsed_seconds'],'transport':out['transport']},ensure_ascii=False,indent=2))
