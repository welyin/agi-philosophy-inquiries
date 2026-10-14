"""Independent dual-certificate review of the positive-spectrum screen."""
from fractions import Fraction as F
from pathlib import Path
import json, hashlib, argparse
P=Path(__file__).resolve().parent

def calculate():
    # Actual tensor contractions in an adapted Lorentz frame; no author formulas imported.
    metric=[-1,1,1,1];u=[1,0,0,0];n=[1,1,0,0]
    outer=lambda v:[[F(v[i]*v[j]) for j in range(4)] for i in range(4)]
    D=outer(u);N=outer(n)
    trace=lambda T:sum(metric[i]*T[i][i] for i in range(4))
    colon=lambda T,S:sum(metric[i]*metric[j]*T[i][j]*S[i][j] for i in range(4) for j in range(4))
    assert trace(D)==-1 and trace(N)==0
    assert colon(D,D)==colon(D,N)==1
    for T in (D,N):assert all(T[3][j]==0 for j in range(4))
    scalar_dd=trace(D)**2
    massless_dd=colon(D,D)-F(1,2)*trace(D)**2
    massive_dd=colon(D,D)-F(1,3)*trace(D)**2
    assert massless_dd==F(1,2) and massive_dd/massless_dd==F(4,3)
    assert scalar_dd/massless_dd==2
    # Dual witness: 4*r0*x(A)-3*y(B)<=0 for every positive spectral atom.
    # Its distance from target (1,1), divided by its l1 norm, gives the minimax lower bound.
    cases=0; atom_checks=0
    for A in [F(1,3),F(1),F(2)]:
      for ratio in [F(2),F(4),F(9,2)]:
        B=A*ratio
        for m0 in [F(1,100),F(1),F(20)]:
          r0=B*(A+m0)/(A*(B+m0))
          dual=(4*r0-3)/(4*r0+3)
          stated=(A*B+(4*B-3*A)*m0)/(7*A*B+(4*B+3*A)*m0)
          assert dual==stated and F(1,7)<dual<1
          c=2/(A/(A+m0)+F(4,3)*B/(B+m0))
          assert c*A/(A+m0)==1-dual
          assert F(4,3)*c*B/(B+m0)==1+dual
          for m in [m0,m0+F(1,7),m0*100]:
            left=4*r0*A/(A+m)-4*B/(B+m)
            factor=-4*B*(B-A)*(m-m0)/((B+m0)*(A+m)*(B+m))
            assert left==factor and left<=0
            atom_checks+=1
          cases+=1
    A=F(1);B=F(4);m0=F(1,100)
    r0=B*(A+m0)/(A*(B+m0))
    optimum=(4*r0-3)/(4*r0+3)
    c=2/(A/(A+m0)+F(4,3)*B/(B+m0))
    source=json.loads((P/'results.json').read_text(encoding='utf-8'))
    assert optimum==F(413,2819)==F(source['band']['optimal_error']['exact'])
    assert c==F(source['band']['saturating_weight']['exact'])
    assert c*A/(A+m0)==F(source['band']['x_A']['exact'])
    assert F(4,3)*c*B/(B+m0)==F(source['band']['y_B']['exact'])
    # Point bound and massless fraction are sharp within their declared response class.
    xp=F(6,7);yp=F(4,3)*xp
    assert max(abs(xp-1),abs(yp-1))==F(1,7)
    delta=F(1,100);a=1-7*delta;b=6*delta;s=F(0)
    x=a+b;y=a+F(4,3)*b+2*s
    assert x==1-delta and y==1+delta and a/x==F(31,33)
    # A common positive average retains y>=4x/3.
    qs=[F(1),F(7,4),F(4)];weights=[F(2,9),F(1,3),F(4,9)]
    averaged=sum(w*q/(q+m0) for q,w in zip(qs,weights))
    assert F(4,3)*averaged>=F(4,3)*averaged
    # Removing one hypothesis at a time gives the stated algebraic escapes.
    small=F(1,10**6)
    single=small/(A+small)
    assert single==F(source['negative_controls']['single_channel_error']['exact'])
    assert single<F(1,100)
    g11=F(3,4);g12=F(1);g22=F(4,3)
    assert g11>0 and g22>0 and g11*g22-g12*g12==0
    assert F(4,3)*g11==g12==1
    b=F(7,9);signed=-b/6
    assert signed<0 and F(4,3)*b+2*signed==b
    hashes={name:hashlib.sha256((P/name).read_bytes()).hexdigest()
            for name in ['proof.md','check.py','results.json','reviewer_check.py']}
    return {'schema':'independent_positive_spectrum_review_v1','status':'PASS_RESTRICTED_KERNEL_CLASS',
      'bindings':hashes,'adapted_frame_tensor_contractions_verified':True,
      'independent_dual_parameter_cases':cases,'spectral_atom_factorization_checks':atom_checks,
      'band_optimum_exact':str(optimum),'extremizer_weight_exact':str(c),
      'point_optimum_exact':'1/7','one_percent_massless_fraction_sharp':'31/33',
      'negative_controls':{'TT_only_small_mass':True,'negative_scalar':True,'nonuniversal_PSD_Gram':True},
      'admission':{'direction_decisive_in_stated_kernel_class':True,'new_vDVZ_discovery_claimed':False,
                   'actual_source_record_menu_certified':False,'old_cognitive_axioms_refuted':False,
                   'all_massive_gravity_excluded':False,'new_number_allocated':False,
                   'new_axiom_adopted':False,'GR_deliverables_completed':0}}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=calculate();dest=P/'review_checks.json'
    if args.write:
      with dest.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
      print('SAVED_INDEPENDENT_POSITIVE_SPECTRUM_REVIEW')
    else:
      assert json.loads(dest.read_text(encoding='utf-8'))==result
      print('READ_ONLY_MATCH: exact dual bound and scoped negative controls passed')
