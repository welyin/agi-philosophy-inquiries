"""1190: exact Cartan source identities and full-Fock Gaussian certificate.
Finite matrices diagnose implementation only. Default mode checks saved results.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,math
import numpy as np

class D:
    def __init__(self,x,dx=0):self.x,self.dx=F(x),F(dx)
    def __add__(self,o):
        o=o if isinstance(o,D) else D(o);return D(self.x+o.x,self.dx+o.dx)
    __radd__=__add__
    def __neg__(self):return D(-self.x,-self.dx)
    def __sub__(self,o):return self+-o if isinstance(o,D) else self+D(-o)
    def __rsub__(self,o):return D(o)+-self
    def __mul__(self,o):
        o=o if isinstance(o,D) else D(o);return D(self.x*o.x,self.dx*o.x+self.x*o.dx)
    __rmul__=__mul__
    def __pow__(self,n):
        return D(self.x**n,n*self.x**(n-1)*self.dx) if n else D(1)
    def __truediv__(self,o):return self*(o**-1 if isinstance(o,D) else D(o)**-1)
    def __rtruediv__(self,o):return D(o)*self**-1

def calculate():
    # Same fixed q and its coordinate derivative; tetrad is a common scale.
    zeta,v,m,q,X=F(2),F(1),F(3,2),F(2,3),F(-1,2)
    rows=[]
    for ph,hs in ((F(1),F(1)),(F(9,4),F(3,2)),(F(4),F(2))):
        for c in (F(0),F(1,2),F(1)):
            def lag(phid,scale):
                f=1-c+c*phid**2/v**2
                return -F(1,2)*scale**4*zeta*phid**3*(X/scale**2+m*m*f*q*q)
            pd=D(ph,1);fd=1-c+c*pd**2/v**2
            Z=zeta*ph**3;Zp=3*zeta*ph**2
            W=Z*m*m*fd.x;Wp=m*m*(Zp*fd.x+Z*fd.dx)
            trace=-Z*X-2*W*q*q
            je=-(Zp*X+Wp*q*q)/2
            ja=je+trace/ph
            je_test=lag(pd,D(1)).dx
            ja_test=lag(pd,D(1,1/ph)).dx
            jh_density_test=lag(pd,D(1,-1/(2*ph))).dx
            # Independently differentiated Einstein-frame action at fixed h.
            Lh=-F(1,2)*hs**4*(zeta*v*pd**2*X/hs**2+zeta*v*v*pd*m*m*fd*q*q)
            assert je_test==je and ja_test==ja
            assert jh_density_test==je-trace/(2*ph)==Lh.dx
            jh=(v/ph)**2*(ja-3*trace/(2*ph))
            assert jh==Lh.dx/hs**4
            rows.append(dict(phi=str(ph),c=str(c),trace=str(trace),j_e=str(je),j_A=str(ja),j_h=str(jh)))
    # Algebraic auxiliary elimination, with spacelike/timelike components.
    eta=[-1,1,1,1];u=[F(2),F(1),F(3,2),F(1,2)]
    for ph in (F(1),F(3,2)):
        aux=[ph*eta[i]*u[i]/2 for i in range(4)]
        norm=sum(eta[i]*aux[i]**2 for i in range(4))
        firstorder=2*ph*sum(aux[i]*u[i] for i in range(4))-2*norm
        assert firstorder==ph**2*sum(eta[i]*u[i]**2 for i in range(4))/2
    # Complete oscillator quadratic compressions: do not square truncated q,p.
    N=16;a=np.diag(np.sqrt(np.arange(1,N)),1);ad=a.T
    num=np.diag(np.arange(N));I=np.eye(N);pair=a@a+ad@ad
    Q2=(2*num+I+pair)/2;P2=(2*num+I-pair)/2
    H=(P2+Q2)/2;src=-P2+Q2/2
    exact=-H/2+3*pair/4
    assert np.max(np.abs(src-exact))<1e-14
    # Zero-to-two survives exactly in the untruncated algebra.
    assert abs(src[2,0]-3*math.sqrt(2)/4)<1e-14
    # At phi=v=1 all c give the same H but different full source.
    family=[]
    for c in (F(0),F(1,2),F(1)):
        cp=float(c);J=-P2+(.5+cp)*Q2
        alpha=c-F(1,2)
        for n0 in range(5):assert abs(J[n0,n0]-float(alpha)*(n0+.5))<1e-13
        family.append(dict(c=str(c),gap_at_reference='1',d_log_gap_at_reference=str(alpha),positive_for_all_phi_gt_zero=True))
    # Full-Fock theorem: Bogoliubov transformation gives exact P2.
    ph=F(121,100);Momega=F(1331,1000)
    n=((Momega-1/Momega)/2)**2
    assert F(2,25)<n<F(1,10)
    assert F(11,10)**3<F(16,9)
    # P2>3/100 follows directly by positive squaring as well.
    assert n*n>F(3,50)**2*(1+n)**3
    Pexact=float(n)/(2*(1+float(n))**1.5)
    t=11*math.pi/20
    diagnostics=[]
    for size in (12,24,48):
        a=np.diag(np.sqrt(np.arange(1,size)),1);pair=a@a+(a.T)@(a.T)
        ni=np.diag(np.arange(size));eye=np.eye(size)
        q2=(2*ni+eye+pair)/2;p2=(2*ni+eye-pair)/2
        Hph=p2/(2*float(ph)**2)+float(ph)*q2/2
        eigen,U=np.linalg.eigh(Hph)
        psi=U@(np.exp(-1j*t*eigen)*U[0,:])
        p2out=float(abs(psi[2])**2)
        energy=float(np.vdot(psi,Hph@psi).real)
        initial=float(Hph[0,0])
        diagnostics.append(dict(dimension=size,P2=p2out,analytic_difference=abs(p2out-Pexact),norm_error=abs(float(np.vdot(psi,psi).real)-1),fixed_H_energy_error=abs(energy-initial)))
    assert diagnostics[-1]['analytic_difference']<1e-12
    return dict(status='PASS_WITH_SCOPE',check_groups=3,
        dynamic_auxiliary_and_double_sources=dict(exact_rows=rows,all_residuals='0',material_q_fixed=True),
        complete_quantum_source=dict(zero_to_two=3*math.sqrt(2)/4,identity_error=float(np.max(np.abs(src-exact))),energy_only_source_sufficient=False),
        full_Fock_probability=dict(phi=str(ph),time='11*pi/20',occupation_exact=str(n),P2=Pexact,strict_lower='3/100',isospectral_frozen_basis_P2='0',certificate='exact Bogoliubov identity plus rational inequalities',finite_matrix_diagnostics=diagnostics),
        allowed_material_family=family,
        scope=dict(full_unknown_input_R_preserved_by_declared_unitaries=True,finite_preparation_and_readout_implemented=False,dynamic_geometry_solved=False,clock_driving_resource_certified=False,new_axioms=0,GR_generated=False,dark_sector_identified=False))

def compare(a,b):
    if isinstance(a,float):return math.isclose(a,b,rel_tol=1e-9,abs_tol=5e-13)
    if isinstance(a,dict):return a.keys()==b.keys() and all(compare(a[k],b[k]) for k in a)
    if isinstance(a,list):return len(a)==len(b) and all(compare(x,y) for x,y in zip(a,b))
    return a==b

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--write',action='store_true');args=pa.parse_args()
    result=calculate();out=Path(__file__).with_name('results.json')
    if args.write:
        with out.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert compare(result,json.loads(out.read_text(encoding='utf-8')))
    print(json.dumps(dict(status=result['status'],check_groups=3,P2=result['full_Fock_probability']['P2'],strict_lower='3/100',actual_instrument=False)))
if __name__=='__main__':main()
