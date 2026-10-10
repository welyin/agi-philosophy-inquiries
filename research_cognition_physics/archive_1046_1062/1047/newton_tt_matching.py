"""1047: one bound Hamiltonian, matched TT source, two finite-pulse outputs.

Only Python standard library and NumPy. Default is read-only reproduction.
--write exclusively creates the initial results. No Fock or atomic truncation
is used to assert a full interacting evolution: computed probabilities are
linear-response coefficients, and the nonzero-pulse remainder is not certified.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, pi, sqrt
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / 'newton_tt_matching_results.json'
HISTORY = (
    'archive_935_955/research_note_940.md',
    'archive_956_989/research_note_985.md',
    'archive_956_989/research_note_986.md',
    'archive_956_989/research_note_982.md',
    'archive_956_989/research_note_977.md',
    'archive_956_989/961/drafts/deferred_material_gravity_interface.md',
    'archive_1009_1043/research_note_1017.md',
    'archive_1009_1043/research_note_1018.md',
    'archive_1009_1043/research_note_1027.md',
    'archive_1009_1043/research_note_1028.md',
    'archive_1009_1043/research_note_1029.md',
    'archive_1009_1043/research_note_1033.md',
    'archive_1044_/research_note_1045.md',
)


def radial_exact(n):
    """Rational angular-averaged integrals BEFORE 4*pi*N_s*N_d.

    Units a=mu=1 imply alpha=1 only for integral normalization, not the
    physical nonrelativistic parameters. Final expressions restore mu*alpha^2.
    """
    b = F(1, 3)
    c, polynomial = (F(4, 3), (F(1),)) if n == 1 else (F(5, 6), (F(2), F(-1)))
    def rint(power):
        return sum(v * F(factorial(power+j)) / c**(power+j+1)
                   for j, v in enumerate(polynomial))
    kinetic = -(4*rint(2)-F(12,5)*b*rint(3)+F(4,15)*b*b*rint(4))
    binding = -F(4,15)*rint(3)
    quad = F(4,15)*rint(6)
    gap = F(1,2*n*n)-F(1,18)
    assert kinetic+binding == -gap*gap*quad/2
    return dict(kinetic=str(kinetic), binding=str(binding), full=str(kinetic+binding),
                quadrupole=str(quad), gap=str(gap),
                amplitude_full_over_kinetic=str((kinetic+binding)/kinetic))


def angular_grid(order=40):
    z, wz = np.polynomial.legendre.leggauss(order)
    phi = 2*pi*np.arange(2*order)/(2*order)
    x = np.sqrt(1-z[:,None]**2)*np.cos(phi)[None,:]
    y = np.sqrt(1-z[:,None]**2)*np.sin(phi)[None,:]
    weights = np.broadcast_to(wz[:,None]/(4*order), x.shape)
    assert abs(weights.sum()-1) < 1e-14
    return x*x, y*y, np.broadcast_to(z[:,None]**2,x.shape), weights


def spatial_integrals(n, order=56):
    """Independent Cartesian differentiated wavefunction plus angular/radial quadrature."""
    u, w = np.polynomial.laguerre.laggauss(order)
    c = 1/n+1/3
    r = u/c
    final_poly = np.ones_like(r) if n == 1 else 2-r
    Ng = 1/sqrt(pi) if n == 1 else 1/(4*sqrt(2*pi))
    Nd = 1/(81*sqrt(2*pi))
    pref = 4*pi*Ng*Nd
    nx2, ny2, nz2, wa = angular_grid()
    p = nx2-ny2
    b = 1/3
    # Cartesian D=(partial_x^2-partial_y^2) acting on P exp(-b r).
    # D/e^(-br) = 4 -4 b r(nx^2+ny^2)+(b^2 r^2+b r)(nx^2-ny^2)^2.
    differentiated = (4-4*b*r[:,None,None]*(nx2+ny2)
        +(b*b*r[:,None,None]**2+b*r[:,None,None])*p*p)
    derivative_angular = np.sum(differentiated*wa,axis=(1,2))
    avg_p2 = float(np.sum(wa*p*p))
    kinetic = -pref*np.sum(w*r*r*final_poly*derivative_angular)/c
    binding = -pref*avg_p2*np.sum(w*r**3*final_poly)/c
    quad = pref*avg_p2*np.sum(w*r**6*final_poly)/c
    # Direct finite metric variation, not source inserted by hand.
    radial_potential = pref*np.sum(w*r**3*final_poly)/c
    def potential_matrix(s):
        distance_factor = np.sqrt(np.exp(2*s)*nx2+np.exp(-2*s)*ny2+nz2)
        return -radial_potential*float(np.sum(wa*p/distance_factor))
    errors = []
    for h in (1e-3, 3e-4, 1e-4):
        # Sum px^2+py^2 matrix vanishes by x/y antisymmetry against the s state.
        kinetic_dH = -kinetic*np.sinh(2*h)/(2*h)
        dH = kinetic_dH+(potential_matrix(h)-potential_matrix(-h))/(2*h)
        errors.append(dict(step=h, derivative=float(dH),
                           error=float(abs(dH+kinetic+binding))))
    assert errors[-1]['error'] < 3e-9
    return dict(kinetic=float(kinetic), binding=float(binding), full=float(kinetic+binding),
                quadrupole=float(quad), metric_derivative_checks=errors)


def norms():
    u,w = np.polynomial.laguerre.laggauss(56)
    Ng=1/sqrt(pi); N2=1/(4*sqrt(2*pi)); Nd=1/(81*sqrt(2*pi))
    a1=4*pi*Ng*Ng*np.sum(w*(u/2)**2)/2
    a2=4*pi*N2*N2*np.sum(w*u*u*(2-u)**2)
    # Integral P^2 over angles = 16*pi*r^4/15.
    ad=16*pi/15*Nd*Nd*np.sum(w*(u/(2/3))**6)/(2/3)
    return [float(a1),float(a2),float(ad)]


def fourier(q):
    def F0(x):
        return 2*pi if abs(x)<1e-14 else -np.expm1(-2j*pi*x)/(1j*x)
    return .25*(F0(q-1)+F0(q+1)-F0(q))-.125*(F0(q-2)+F0(q+2))


def pulse_checks():
    z,w=np.polynomial.legendre.leggauss(180)
    u=pi*(1+z); w=pi*w
    chi=np.sin(u/2)**2; chip=.5*np.sin(u); chipp=.5*np.cos(u)
    shape=chi*np.cos(u)
    shapepp=chipp*np.cos(u)-2*chip*np.sin(u)-chi*np.cos(u)
    rows=[]
    for q,tau,ratio in ((1,-1/(8*sqrt(2)),F(16,25)),(5/32,1024/15625,F(1,64))):
        phase=np.exp(-1j*q*u)
        J=complex(np.sum(w*shape*phase))
        assert abs(J-fourier(q))<1e-12
        lhs=-q*q*J
        rhs=complex(np.sum(w*shapepp*phase))
        naive=complex(np.sum(w*(-chi*np.cos(u))*phase))
        assert abs(lhs-rhs)<2e-12
        coeff=abs(tau*J/(4/9))**2
        rows.append(dict(gap_over_first_gap=q, fourier_real=J.real, fourier_imag=J.imag,
                         coefficient_full=float(coeff), coefficient_kinetic=float(coeff*float(ratio)),
                         bare_over_full=str(ratio), tidal_identity_error=float(abs(lhs-rhs)),
                         dropped_envelope_terms_error=float(abs(naive-lhs))))
    assert abs(rows[0]['coefficient_full']-81*pi*pi/8192)<1e-13
    assert rows[1]['coefficient_full']>0
    assert rows[1]['dropped_envelope_terms_error']>1
    q=F(5,32)
    R=-(q*q+2)/(2*q*(q*q-1)*(q*q-4))
    assert abs(fourier(float(q))-(1-np.exp(-2j*pi*float(q)))*complex(0,-float(R)))<1e-12
    # Exact two-output response calibration obstruction; x=|lambda|^2 >=0.
    a,b=F(16,25),F(1,64)
    x=2/(a+b); error=(a-b)/(a+b)
    assert x==F(3200,1049) and error==F(999,1049)
    assert a*x-1==1-b*x==error
    return dict(duration_in_inverse_first_gap='2*pi',same_pulse_for_both_outputs=True,
                rows=rows,second_fourier_rational_factor=str(R),
                optimal_common_intensity_scale=str(x),
                sharp_relative_response_error=str(error),
                first_line_calibrated_second_response_ratio=str(F(25,1024)))


def compute():
    exact={str(n):radial_exact(n) for n in (1,2)}
    numerical={str(n):spatial_integrals(n) for n in (1,2)}
    expected={1:(-1/(10*sqrt(2)),-1/(40*sqrt(2)),81/(64*sqrt(2))),
              2:(128/15625,896/15625,-10616832/390625)}
    residuals=[]
    for n in (1,2):
        row=numerical[str(n)]
        for name,value in zip(('kinetic','binding','quadrupole'),expected[n]):
            residuals.append(abs(row[name]-value))
            assert abs(row[name]-value)<5e-11
    norm=norms()
    assert max(abs(v-1) for v in norm)<1e-12
    return dict(round=1047,new_science_groups=1,new_cognitive_axioms=0,
                full_roadmap_completed_here=False,probabilities_are_leading_coefficients=True,
                finite_nonzero_amplitude_remainder_certified=False,
                exact_radial_factors=exact,spatial_quadrature=numerical,
                state_norms=norm,spatial_max_error=float(max(residuals)),
                pulse=pulse_checks(),all_scientific_checks_passed=True,
                historical_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in HISTORY})


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert np.isclose(a,b,rtol=2e-9,atol=2e-11),(path,a,b)
    else: assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=compute()
    if args.write:
        with OUTPUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else: compare(result,json.loads(OUTPUT.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1047,all_scientific_checks_passed=True,
        mode='exclusive_write' if args.write else 'read_only_compare',
        spatial_max_error=result['spatial_max_error'],pulse=result['pulse']),ensure_ascii=False))


if __name__=='__main__': main()
