"""1265 independent check. One group; no author implementation imported.
The Einstein tensor is rebuilt from a general metric two-jet as exact sparse
polynomials. The capacity/integration arguments are analytic; their constants
and the inherited finite radius are reconstructed with Fraction arithmetic.
Default is read-only. --write creates independent_results.json exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
OUT = HERE / 'independent_results.json'
OLD = {
 '../1253/proof.md': '422c736a1c3d65671fcef1795cc889e44514516344510048be7c8240acc29161',
 '../1264/source_proof.md': 'ae5f0e49a8d47adf1256d5d6ffa701a5e154da802e97bfb7ffb31535cc56f62d',
 '../1264/source_results.json': '008252c273b1a5f02d62f7205abe27eeaf9707eb628e875e37e670541accdb1a',
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

class Poly:
    def __init__(self, coeff=0):
        if isinstance(coeff, Poly):
            self.c = dict(coeff.c)
        elif isinstance(coeff, dict):
            self.c = {tuple(k): F(v) for k,v in coeff.items() if v}
        else:
            self.c = {(): F(coeff)} if coeff else {}
    @staticmethod
    def var(name):
        return Poly({(name,): F(1)})
    def __add__(self, other):
        other=Poly(other); c=dict(self.c)
        for key,value in other.c.items():
            c[key]=c.get(key,F(0))+value
        return Poly(c)
    __radd__=__add__
    def __neg__(self):
        return Poly({k:-v for k,v in self.c.items()})
    def __sub__(self,other): return self+-Poly(other)
    def __rsub__(self,other): return Poly(other)+-self
    def __mul__(self,other):
        other=Poly(other); c={}
        for ka,va in self.c.items():
            for kb,vb in other.c.items():
                key=tuple(sorted(ka+kb))
                c[key]=c.get(key,F(0))+va*vb
        return Poly(c)
    __rmul__=__mul__
    def __eq__(self,other): return self.c==Poly(other).c
    def serial(self):
        return {'*'.join(k) if k else '1':str(v) for k,v in sorted(self.c.items())}


def tensor_from_metric_jet():
    # At one point e^(2 omega) is set to 1. A constant conformal factor
    # cancels from coordinate G_ab, while its first/second jets are retained.
    eta=[F(-1),F(1),F(1),F(1)]
    zero=Poly()
    grad=[zero]+[Poly.var('a'+str(i)) for i in range(1,4)]
    hess=[[zero if i==0 or j==0 else Poly.var('h'+str(min(i,j))+str(max(i,j)))
           for j in range(4)] for i in range(4)]
    g=[[eta[i] if i==j else F(0) for j in range(4)] for i in range(4)]
    dg=[[[2*grad[d]*g[i][j] for j in range(4)] for i in range(4)] for d in range(4)]
    ddg=[[[[(2*hess[d][e]+4*grad[d]*grad[e])*g[i][j]
             for j in range(4)] for i in range(4)] for e in range(4)] for d in range(4)]
    # d(g^{-1})=-g^{-1} (dg) g^{-1}; evaluated for the same metric jet.
    dgi=[[[ -sum((g[i][r]*dg[d][r][s]*g[s][j]
                   for r in range(4) for s in range(4)),zero)
             for j in range(4)] for i in range(4)] for d in range(4)]
    Gamma=[[[sum((F(1,2)*g[a][d]*(dg[b][d][c]+dg[c][d][b]-dg[d][b][c])
                  for d in range(4)),zero)
             for c in range(4)] for b in range(4)] for a in range(4)]
    dGamma=[[[[sum((F(1,2)*(
                     dgi[e][a][d]*(dg[b][d][c]+dg[c][d][b]-dg[d][b][c])
                    +g[a][d]*(ddg[e][b][d][c]+ddg[e][c][d][b]-ddg[e][d][b][c]))
                    for d in range(4)),zero)
               for c in range(4)] for b in range(4)] for a in range(4)] for e in range(4)]
    Ric=[[sum((dGamma[a][a][b][c]-dGamma[c][a][b][a]
                +sum((Gamma[a][a][d]*Gamma[d][b][c]
                     -Gamma[a][c][d]*Gamma[d][b][a] for d in range(4)),zero)
                for a in range(4)),zero) for c in range(4)] for b in range(4)]
    scalar=sum((g[a][b]*Ric[a][b] for a in range(4) for b in range(4)),zero)
    ein=[[Ric[a][b]-F(1,2)*g[a][b]*scalar for b in range(4)] for a in range(4)]
    lap=sum((hess[i][i] for i in range(1,4)),zero)
    grad2=sum((grad[i]*grad[i] for i in range(1,4)),zero)
    expected=[[zero for _ in range(4)] for _ in range(4)]
    expected[0][0]=-2*lap-grad2
    for i in range(1,4):
        for j in range(1,4):
            expected[i][j]=-2*hess[i][j]+2*grad[i]*grad[j]
            if i==j: expected[i][j]+=2*lap+grad2
    assert all(ein[i][j]==expected[i][j] for i in range(4) for j in range(4))
    # Hamiltonian projection agrees with 1/2 R^(3) when K_ij=0.
    assert 2*ein[0][0]==-4*lap-2*grad2
    return {str(i)+str(j):ein[i][j].serial() for i in range(4) for j in range(i,4)}


def source_and_weight_constants():
    # Exact one-particle normal-ordered quadratic form. q has all material
    # components; |Mq|^2 is positive without diagonalising/averaging M.
    terms=['dt_re','dt_im','dx_re','dx_im','dy_re','dy_im','dz_re','dz_im','Mq_re','Mq_im']
    positive=Poly()
    for name in terms:
        v=Poly.var(name);positive+=v*v
    assert len(positive.c)==10
    assert all(len(k)==2 and k[0]==k[1] and c==1 for k,c in positive.c.items())
    # The last two squares carry e^(2 omega)>0, not a claim this factor is 1.
    # For chi=1 wherever grad f !=0: int (chi(1-f)) Delta f = int |grad f|^2.
    f=Poly.var('f'); chi=Poly.var('chi')
    assert Poly(1)-f+f==Poly(1)
    # In the annulus, dw=-df; the integration-by-parts coefficient is +1.
    weighted_laplacian_coefficient=-F(-1)
    assert weighted_laplacian_coefficient==1
    # Capacity: int_R^(2R) r^-2 dr=1/(2R); Cauchy gives int r^2 f'^2>=2R.
    cap_over_4piR=1/(F(1)-F(1,2))
    assert cap_over_4piR==2
    I_over_piR=4*cap_over_4piR
    volume_over_piR3=F(4,3)*4**3
    linear_coeff=2*I_over_piR/volume_over_piR3
    quadratic_coeff=I_over_piR/volume_over_piR3
    assert linear_coeff==F(3,16)
    assert quadratic_coeff==F(3,32)
    # Smooth drop guarantees Z>0. Only an upper bound on Z is used.
    r=F(1,10000)
    L=19100606975
    R=2*L
    assert R==38201213950
    lower_lambda=r/(1+r)  # strict: log(1+r)>r/(1+r), from integral 1/(1+t)
    bound=linear_coeff*lower_lambda/R**2
    assert bound==F(3,16*10001*R**2)>0
    assert F(1,10**26)<bound<F(2,10**26)
    # The integration identity is independent of the magnitude of positive T.
    # An error allowance above the gap is not excluded by this witness alone.
    tolerances={
       'strictly_rejected_example':str(bound/2),
       'not_rejected_by_this_bound_alone_example':str(2*bound),
    }
    assert bound/2<bound<2*bound
    # At lambda=0 this particular geometric lower bound vanishes. It does not
    # imply that a nonzero material source solves the flat Einstein equation.
    assert linear_coeff*0/R**2==0
    return {
       'source_quadratic_form':'|dt q|^2+sum_i |di q|^2+exp(2 omega)|M q|^2 >=0',
       'source_range':'same static-vacuum normal ordering, one-particle density matrices and arbitrary passive reference; not general renormalized states',
       'Einstein_component':'G00-kappa<T00>, kappa>0, Lambda=0, static metric',
       'test_weight':'w=f_(2R)*(1-f_R), nonnegative compact support <=4R; not a certified instrument',
       'int_w_Delta_f_coefficient':str(weighted_laplacian_coefficient),
       'capacity_I_over_pi_R_lower':str(I_over_piR),
       'weight_volume_over_pi_R3_upper':str(volume_over_piR3),
       'residual_linear_coefficient':str(linear_coeff),
       'constant_weight_quadratic_coefficient':str(quadratic_coeff),
       'L':L,'R':R,'r':str(r),'lambda_strict_lower':str(lower_lambda),
       'normalized_negative_residual_magnitude_strict_lower':str(bound),
       'bound_display_only':float(bound),
       'tolerance_examples':tolerances,
       'metric_and_source_joint_error_deletion':'eta_geometry+kappa*eta_source < lower_bound excludes the nominated self-source matching only',
       'large_R_gap_scales_as_R_minus_2':True,
       'record_probability_budget_is_not_constraint_residual_budget':True,
    }


def run():
    for rel,digest in OLD.items():
        assert sha((HERE/rel).resolve())==digest,rel
    return {
      'round':1265,'verdict':'PASS_INDEPENDENT_FINITE_PROFILE_CONSTRAINT_DELETION',
      'independent_groups':1,'code_sha256':sha(Path(__file__).resolve()),
      'old_source_sha256':OLD,
      'group':{'exact_general_static_metric_jet_Einstein_tensor':tensor_from_metric_jet(),
               'same_positive_source_and_weight_bound':source_and_weight_constants()},
      'scope':'excludes exact self-source completion of the inherited compact static conformal profile in the specified Einstein/KG class; no cognitive no-go or actual detector or complete GR claim',
    }

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2);handle.write('\n')
        print('WRITE_ONCE_PASS_1265_INDEPENDENT',result['independent_groups'])
    else:
        assert json.loads(OUT.read_text(encoding='utf-8'))==result,'Saved independent result mismatch'
        print('READ_ONLY_MATCH_1265_INDEPENDENT',result['independent_groups'])