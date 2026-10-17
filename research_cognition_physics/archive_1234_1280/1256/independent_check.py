"""Independent round 1256 audit. Standard library; author code is not imported.

Finite rational fixtures check implementations of the analytic symplectic/source
identities. They do not prove the continuous-field theorem or a quantum model.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent


def dot(a, b):
    return sum((x*y for x,y in zip(a,b)), F(0))


def mv(a, b):
    return [dot(row,b) for row in a]


def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def inv(a):
    b = [list(row)+[F(i==j) for j in range(3)] for i,row in enumerate(a)]
    for k in range(3):
        pivot = next(i for i in range(k,3) if b[i][k])
        b[k], b[pivot] = b[pivot], b[k]
        scale = b[k][k]
        b[k] = [x/scale for x in b[k]]
        for i in range(3):
            if i!=k:
                scale = b[i][k]
                b[i] = [x-scale*y for x,y in zip(b[i],b[k])]
    out = [row[3:] for row in b]
    assert all(dot(a[i], [out[k][j] for k in range(3)]) == (i==j)
               for i in range(3) for j in range(3))
    return out


def pack(x):
    if isinstance(x,F):
        return str(x)
    if isinstance(x,list):
        return [pack(v) for v in x]
    if isinstance(x,dict):
        return {k:pack(v) for k,v in x.items()}
    return x


def canonical_source_transport():
    # DZ[A][i] is not a diagonal volume rescaling; the density is not unity.
    B = [[F(2),F(1,3),F(0)], [F(0),F(3,2),F(1,5)], [F(1,7),F(0),F(4,3)]]
    A = inv(B)
    jac, label_n = det(B), F(3,7)
    eta = label_n*jac
    assert jac > 0 and eta > 0
    I, p = F(5,4), F(-1,3)
    P = [F(2,9),F(-1,8),F(1,6)]
    phase_grad = [F(1,5),F(-2,7),F(3,11)]
    Pi_phase = eta*I
    weighted = [P[i]+I*phase_grad[i] for i in range(3)]
    Pi_Z = [-eta*sum(weighted[i]*A[i][a] for i in range(3)) for a in range(3)]
    dZ = [F(1,13),F(-2,17),F(3,19)]
    dphase = F(-2,23)
    dX = [-x for x in mv(A,dZ)]
    dphase_L = dphase+dot(phase_grad,dX)
    original_theta = eta*(dot(P,dX)+I*dphase_L)
    transformed_theta = dot(Pi_Z,dZ)+Pi_phase*dphase
    assert original_theta == transformed_theta
    D = [sum(Pi_Z[a]*B[a][i] for a in range(3))+Pi_phase*phase_grad[i] for i in range(3)]
    assert D == [-eta*v for v in P]
    recovered_P = [-v/eta for v in D]
    assert recovered_P == P and Pi_phase/eta == I
    hi = [[F(2),F(1,4),F(0)], [F(1,4),F(1),F(1,5)], [F(0),F(1,5),F(3,2)]]
    assert hi[0][0] > 0 and hi[0][0]*hi[1][1]-hi[0][1]**2 > 0 and det(hi)>0
    M = (1+(1+p)*I)/10
    sigma = (eta+(1+p)*Pi_phase)/10
    assert M>0 and sigma==eta*M
    e2 = dot(P,mv(hi,P))+M*M
    E2 = dot(D,mv(hi,D))+sigma*sigma
    assert E2 == eta*eta*e2
    # Compare densitized stress directly to the original moving material.
    # Multiplication by E removes the harmless positive square root.
    stress_times_E = [[D[i]*D[j] for j in range(3)] for i in range(3)]
    original_stress_times_E = [[eta*eta*P[i]*P[j] for j in range(3)] for i in range(3)]
    assert stress_times_E == original_stress_times_E
    assert D == [-eta*x for x in P]  # sqrt(h) T_ni, with the original shift convention
    field_source_times_E = sigma*Pi_phase/10
    assert field_source_times_E == eta*eta*M*I/10
    # Choose a=N/e>0, allowing rational fixtures without confusing e and E=eta e.
    a = F(7,5)
    beta = [F(1,9),F(-1,11),F(2,13)]
    v = [a*x-b for x,b in zip(mv(hi,P),beta)]
    hamilton_velocity = [a*x/eta+b for x,b in zip(mv(hi,D),beta)]
    assert hamilton_velocity == [-x for x in v]
    Zdot = mv(B,hamilton_velocity)
    assert all(Zdot[j]+dot(B[j],v)==0 for j in range(3))
    phase_dot = dot(hamilton_velocity,phase_grad)+(a/eta)*sigma*(1+p)/10
    assert phase_dot+dot(v,phase_grad) == a*M*(1+p)/10
    p_source = -(a/eta)*sigma*Pi_phase/10
    assert p_source == -eta*a*M*I/10
    # Nonzero second jets of Z and a nonconstant label density test the actual
    # Jacobian continuity identity; relative modes have not been set constant.
    Z_hess = [[[F((j+1)*(i+k+1)*(-1 if (j+i+k)%2 else 1),31+7*j+3*(i+k))
                for k in range(3)] for i in range(3)] for j in range(3)]
    density_label_grad = [F(1,13),F(-1,17),F(1,19)]
    v_grad = [[F((j+1)*(i+2)*(-1 if i==j else 1),37+5*j+7*i)
               for i in range(3)] for j in range(3)]  # [j][i]=partial_i v^j
    Bdot = [[-sum(v_grad[j][i]*B[k][j]+v[j]*Z_hess[k][i][j] for j in range(3))
             for i in range(3)] for k in range(3)]
    jac_dot = jac*sum(A[i][k]*Bdot[k][i] for i in range(3) for k in range(3))
    n_dot = dot(density_label_grad,Zdot)
    eta_dot = n_dot*jac+label_n*jac_dot
    eta_grad = []
    for i in range(3):
        n_i = sum(density_label_grad[k]*B[k][i] for k in range(3))
        jac_i = jac*sum(A[j][k]*Z_hess[k][i][j] for j in range(3) for k in range(3))
        eta_grad.append(n_i*jac+label_n*jac_i)
    div_v = sum(v_grad[i][i] for i in range(3))
    assert eta_dot+dot(v,eta_grad)+eta*div_v == 0
    I_grad = [F(2,13),F(-1,17),F(1,19)]
    Pi_grad = [I*x+eta*y for x,y in zip(eta_grad,I_grad)]
    Pi_dot = -dot(v,Pi_grad)-Pi_phase*div_v
    I_dot = (Pi_dot-I*eta_dot)/eta
    assert I_dot+dot(v,I_grad)==0
    # Recover the original mechanical force with arbitrary gradients of I, P,
    # metric, lapse, shift and p.  e*grad(N) is rational for the chosen fixtures.
    p_grad = [F(1,29),F(-2,31),F(1,37)]
    N_grad_over_e = [F(1,41),F(-1,43),F(2,47)]
    P_grad = [[F((j+1)*(i+1),53+3*j+5*i) for i in range(3)] for j in range(3)]
    beta_grad = [[F((-1 if j==i else 1)*(j+i+1),61+3*j+7*i)
                  for i in range(3)] for j in range(3)]
    hi_grad = [[[F((i+1)*(j+k+1),71+5*i+3*(j+k)) for k in range(3)]
                for j in range(3)] for i in range(3)]
    force_residual = []
    raised_P = mv(hi,P)
    a0 = a*M/10
    for i in range(3):
        Mi = ((1+p)*I_grad[i]+I*p_grad[i])/10
        metric_term = dot(P,mv(hi_grad[i],P))
        Ei_times_e = dot(raised_P,[P_grad[j][i] for j in range(3)])+metric_term/2+M*Mi
        a0i = (N_grad_over_e[i]*M+a*Mi-a*M*Ei_times_e/e2)/10
        bi = (1+p)*a0i+p_grad[i]*a0
        grad_raised_P = [dot(hi_grad[i][j],P)+dot(hi[j],[P_grad[k][i] for k in range(3)])
                         for j in range(3)]
        vi = [N_grad_over_e[i]*raised_P[j]+a*grad_raised_P[j]
              -a*raised_P[j]*Ei_times_e/e2-beta_grad[j][i] for j in range(3)]
        force_eulerian = -dot(P,vi)-a0i-I*bi
        force_original = (-e2*N_grad_over_e[i]-a*metric_term/2
                          +dot(P,[beta_grad[j][i] for j in range(3)])-a*M*I*p_grad[i]/10)
        force_residual.append(force_eulerian-force_original)
    assert force_residual == [0,0,0]
    return pack({
        'status':'PASS','DZ':B,'det_DZ':jac,'label_density':label_n,'eta':eta,
        'I':I,'p':p,'mechanical_P':P,'phase_gradient':phase_grad,
        'Pi_phase':Pi_phase,'Pi_Z':Pi_Z,'momentum_density_D':D,
        'one_form_value':original_theta,'one_form_residual':F(0),
        'mass':M,'mass_density_sigma':sigma,'energy_density_squared':E2,
        'energy_identity_residual':F(0),'spatial_stress_density_times_E':stress_times_E,
        'matter_H_p_derivative_times_E_over_N':field_source_times_E,
        'physical_coordinate_velocity':v,'actual_phase_advective_rate':a*M*(1+p)/10,
        'eta_continuity_residual':F(0),'I_advection_residual':F(0),'mechanical_force_residual':force_residual,
        'jet_scope':'Nonzero finite Z Hessians, density gradients and velocity gradients; not a PDE proof'
    })


def deletion_periodicity():
    eta,I,p = F(1,10),F(1),F(0)
    B = [[F(i==j) for j in range(3)] for i in range(3)]
    P = [F(0)]*3
    grad_phase = [F(1),F(0),F(0)]
    Pi_phase = eta*I
    Pi_Z = [-eta*I*x for x in grad_phase]
    correct_D = [Pi_Z[i]+Pi_phase*grad_phase[i] for i in range(3)]
    bad_D = list(Pi_Z)  # omits the compensated phase-gradient contribution
    bad_P = [-x/eta for x in bad_D]
    sigma = (eta+(1+p)*Pi_phase)/10
    correct_E2 = dot(correct_D,correct_D)+sigma*sigma
    bad_E2 = dot(bad_D,bad_D)+sigma*sigma
    assert correct_D == [0,0,0] and sigma==F(1,50)
    assert correct_E2==F(1,2500) and bad_E2==F(26,2500)
    threshold = F(2,25)
    # E_bad > E_correct + 2/25, both positive; no floating square root.
    strict_squared_margin = bad_E2-(sigma+threshold)**2
    assert strict_squared_margin == F(1,2500) and strict_squared_margin>0
    assert bad_P == [1,0,0]
    return pack({
        'status':'PASS','eta':eta,'Pi_phase':Pi_phase,'Pi_Z':Pi_Z,
        'correct_D':correct_D,'bad_D':bad_D,'bad_mechanical_P':bad_P,
        'bad_physical_momentum_density':[eta*x for x in bad_P],
        'correct_E':sigma,'bad_E_squared':bad_E2,
        'strict_energy_gap_lower':threshold,'positive_squared_margin':strict_squared_margin,
        'periodic_realization':{
            'coordinate_period':'L=10',
            'phase':'theta(x)=(5/pi) sin(pi*x1/5)',
            'phase_at_zero':'0','gradient_at_zero':['1','0','0'],
            'periodicity_certificate':'theta(x1+10)=theta(x1); elementary sine period 2*pi',
            'label_map':'Z=x on the same torus; identity winding, positive unit Jacobian',
            'full_original_state':'I=1 and P=0 everywhere; Pi_Z1=-(1/10) cos(pi*x1/5)',
            'scope':'Smooth periodic material phase; no extra spatial direction, energy, or measurement permission'
        },
        'scope':'Tests the same constitutive/source identity; it does not exclude different materials or geometries'
    })


def calculate():
    files = [HERE/'admission.md',HERE/'proof.md',HERE/'NEXT.md',HERE/'check.py',
             HERE/'results.json',HERE.parent/'research_note_1256.md',Path(__file__).resolve()]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise SystemExit('Await final author files before first binding: '+', '.join(missing))
    return {'round':1256,'status':'PASS','independent_groups':2,
            'canonical_source_transport':canonical_source_transport(),
            'deletion_periodicity':deletion_periodicity(),
            'input_sha256':{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--write',action='store_true',help='exclusive first save only')
    args=ap.parse_args()
    result=calculate()
    dest=HERE/'independent_results.json'
    if args.write:
        with dest.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        print('SAVED_1256_INDEPENDENT: two groups; PASS')
    else:
        assert result==json.loads(dest.read_text(encoding='utf-8')), 'Results or bound source versions differ'
        print('READ_ONLY_MATCH_1256_INDEPENDENT: two groups; PASS')


if __name__=='__main__':
    main()