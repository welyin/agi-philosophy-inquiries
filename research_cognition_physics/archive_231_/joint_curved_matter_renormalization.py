"""553: joint matter/curvature one-loop interfaces, not quantum-gravity closure.

Euclidean +bW W^2; fixed classical metric; mass-independent MS; all listed
matter active, no hard fermion mass or graviton loops. Fractions certify the
displayed algebra; floating inherited endpoints are diagnostics, not ODE bounds.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_singlet_common_mass_rg as rg

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_curved_matter_renormalization_results.json'


def matter_B(q, S, z, lh, p, ls, gy, gw):
    return [[12*lh+6*q+2*S-(3*gy+9*gw)/2, 2*p], [8*p, 6*ls+4*z]]


def scalar_tensor(a, b, c, d, lh, p, ls):
    indices = (a, b, c, d)
    if all(i < 4 for i in indices):
        return 2*lh*((a == b)*(c == d)+(a == c)*(b == d)+(a == d)*(b == c))
    if all(i == 4 for i in indices):
        return 6*ls
    if indices.count(4) == 2:
        h = [i for i in indices if i != 4]
        return 2*p*(h[0] == h[1])
    return Q(0)


def counts(N, ng):
    ns, nw, nv = 5, 4*(N+1)*ng, N*N+3
    c = Q(ns+3*nw+12*nv, 120)
    a = Q(2*ns+11*nw+124*nv, 720)
    b2 = -Q(22, 3)+Q(2, 3)*Q((N+1)*ng, 2)+Q(1, 6)
    return ns, nw, nv, c, -a, b2


def poly_sub(a, b):
    return {key: value for key in a.keys() | b.keys()
            if (value := a.get(key, Q(0))-b.get(key, Q(0))) != 0}


def run():
    checks = []
    pars = list(map(Q, ('0.2','0.3','0.4','0.25','0.15','0.7','0.16','0.32')))
    q,S,z,lh,p,ls,gy,gw = pars
    B = matter_B(*pars)
    gamma = [3*q+S-(3*gy+9*gw)/4]*4+[2*z]
    tensor_columns = []
    for delta in ([Q(1)]*4+[Q(0)], [Q(0)]*4+[Q(1)]):
        beta = [sum(scalar_tensor(a,a,c,c,lh,p,ls)*delta[c] for c in range(5))
                +2*gamma[a]*delta[a] for a in range(5)]
        assert len(set(beta[:4])) == 1
        tensor_columns.append([beta[0], beta[4]])
    assert [[tensor_columns[j][i] for j in range(2)] for i in range(2)] == B
    for dh,ds in ((Q(0),Q(0)), (Q(1,7),Q(-1,11)), (Q(0),Q(1))):
        v = [B[0][0]*dh+B[0][1]*ds, B[1][0]*dh+B[1][1]*ds]
        inherited = rg.beta_numerator([*map(float,pars[:6]),float(dh),float(ds)],
                                     [float(gy),float(gw),.9])[6:]
        assert np.allclose(inherited, list(map(float,v)), atol=1e-14, rtol=0)
        if dh == ds == 0:
            assert v == [0,0]
    assert B[0][1] != 0 and B[1][0] != 0  # One conformal field alone does not suffice.
    checks.append('quartic_tensor_and_mass_spurion_map_including_portal_and_wavefunction')

    # Half trace of (-C + delta*R)^2, in Euclidean V0-M02*R/2+bR*R^2.
    Ch,Cs,dh,ds = Q(2,3),Q(3,5),Q(1,7),Q(-2,9)
    heat = [Q(0),Q(0),Q(0)]
    for C,delta,multiplicity in ((Ch,dh,4),(Cs,ds,1)):
        for j,val in enumerate((C*C/2,-C*delta,delta*delta/2)):
            heat[j] += multiplicity*val
    assert heat == [2*Ch*Ch+Cs*Cs/2, -(4*Ch*dh+Cs*ds), (4*dh*dh+ds*ds)/2]
    assert -2*heat[1] == 2*(4*Ch*dh+Cs*ds)
    checks.append('scalar_heat_kernel_vacuum_Einstein_and_R_squared_normalizations')

    # Independent published SM alpha basis, then add s and three Weyl singlets.
    xi = Q(1,6)
    a1,a2,a3 = 2*xi*xi-2*xi/3-Q(277,144),Q(571,90),-Q(293,720)
    sm = [a1+(a2+a3)/3, a2/2+2*a3, -a2/2-a3]
    assert sm == [0,Q(283,120),-Q(1991,720)]
    ns,nw,nv,c,e,b2 = counts(3,3)
    assert (ns,nw,nv) == (5,48,12)
    assert c == sm[1]+Q(1+3*3,120) == Q(293,120)
    assert e == sm[2]-Q(1,360)-3*Q(11,720) == -Q(1013,360)
    assert b2 == -Q(19,6)
    assert c != Q(2+3*48+12*12,120)  # h,s radial count misses Goldstones.
    assert c != Q(5+3*47+12*12,120)  # A zero-Yukawa nu_R still gravitates.
    checks.append('Euclidean_curvature_basis_and_complete_unweighted_particle_count')

    # Exact polynomial identity in N,ng (monomial key is their pair of powers).
    cpoly = {(2,0):Q(12,120),(1,1):Q(12,120),(0,1):Q(12,120),(0,0):Q(41,120)}
    spectral_slope = {(1,1):Q(12,120),(0,1):Q(12,120),(0,0):-Q(258,120)}
    gap_poly = poly_sub(cpoly,spectral_slope)
    assert gap_poly == {(2,0):Q(12,120),(0,0):Q(299,120)}
    # Weights enter the bare trace but cancel before any loop count is taken.
    weight_rows = []
    for N,ng,wq,wl in ((3,3,Q(1),Q(2)), (5,1,Q(2,3),Q(7,4)), (4,2,Q(1,5),Q(9))):
        Iw = (N*wq+wl)/2
        n = 16*ng*Iw
        # Ratio bW / g_w^(-2), cancelling the common f0/pi^2.
        ratio = (-n/Q(320))/(ng*Iw/3)
        assert ratio == -Q(3,20)
        *_,cn,en,bn = counts(N,ng)
        gap = cn-3*bn/10
        assert gap == Q(12*N*N+299,120) > 0
        weight_rows.append(dict(N=N,ng=ng,wq=str(wq),wl=str(wl),ratio=str(ratio),gap=str(gap)))
    checks.append('all_N_generation_cancellation_as_exact_polynomial_and_spectral_weight_elimination')

    # tau=ln(mu/mu0)/(16*pi^2); solve this decoupled one-loop pair exactly.
    flow_rows = []
    for N,ng in ((3,3),(3,20),(5,1),(4,2)):
        *_,cn,en,bn = counts(N,ng)
        x0 = Q(4); bw0 = -3*x0/20
        for tau in (Q(-1,100), Q(0), Q(1,100)):
            x = x0-2*bn*tau
            bw = bw0+cn*tau
            assert x > 0
            deviation = bw+3*x/20
            expected = Q(12*N*N+299,120)*tau
            assert deviation == expected
            assert (deviation == 0) == (tau == 0)
            flow_rows.append(dict(N=N,ng=ng,tau=str(tau),inverse_weak_squared=str(x),
                                  bW=str(bw),deviation=str(deviation)))
    checks.append('exact_flow_leaves_spectral_surface_but_single_scale_matching_survives')

    # Polynomial in X=h^2,Y=s^2: insert C=L(X,Y) into beta Delta.
    bH,bP,bS = rg.beta_numerator([*map(float,pars[:6]),1.,1.], [float(gy),float(gw),.9])[3:6]
    exact_b = [24*lh*lh-(3*gy+9*gw)*lh+Q(3,8)*(gy*gy+2*gy*gw+3*gw*gw)
               +4*(3*q+S)*lh-2*(3*q*q+S*S)+2*p*p,
               p*(12*lh+6*ls+8*p+6*q+2*S+4*z-(3*gy+9*gw)/2)-4*z*S,
               18*ls*ls+8*ls*z+8*p*p-8*z*z]
    assert np.allclose([bH,bP,bS], list(map(float,exact_b)), atol=1e-14, rtol=0)
    tree_coeff = [2*lh*lh+p*p/2-(B[0][0]*lh+2*p*p)/2+exact_b[0]/4,
                  4*lh*p+p*ls-(B[0][0]*p+2*p*ls+8*p*lh+B[1][1]*p)/2+exact_b[1]/2,
                  2*p*p+ls*ls/2-(8*p*p+B[1][1]*ls)/2+exact_b[2]/4]
    spectrum_coeff = [2*lh*lh+3*gw*gw/16+Q(3,32)*(gy+gw)**2-Q(3,2)*q*q-S*S/2,
                      4*p*p-2*z*S, 2*ls*ls-2*z*z]
    assert tree_coeff == spectrum_coeff
    spec = importlib.util.spec_from_file_location('round553_vacuum_probe',HERE/'round553_drafts/vacuum_identity_probe.py')
    probe = importlib.util.module_from_spec(spec); spec.loader.exec_module(probe)
    diagnostic = probe.run()
    assert diagnostic == json.loads((HERE/'round553_drafts/vacuum_identity_probe_results.json').read_text('utf8'))
    checks.append('same_matter_tree_vacuum_supertrace_and_explicit_scale_cancellation')

    dependencies = ('research_note_531.md','research_note_533.md','research_note_544.md',
                    'research_note_552.md','joint_singlet_common_mass_rg.py',
                    'joint_singlet_common_mass_rg_results.json',
                    'round553_drafts/vacuum_identity_probe.py','round553_drafts/vacuum_identity_probe_results.json')
    return dict(round=553,tests_run=len(checks),failures=0,errors=0,checks=checks,
        conformal_matrix=[[str(a) for a in row] for row in B],
        SM_plus_s_plus_three_RH_curvature=dict(real_scalars=ns,Weyl_fermions=nw,vectors=nv,
            beta_bW_numerator=str(c),beta_bE_numerator=str(e),b2=str(b2),
            spectral_deviation_numerator=str(c-3*b2/10)),
        general_gap_polynomial={str(k):str(v) for k,v in sorted(gap_poly.items())},
        positive_weight_proxies=weight_rows,exact_joint_flow=flow_rows,
        vacuum_identity_polynomial_coefficients=list(map(str,tree_coeff)),
        frozen_endpoint_diagnostics=diagnostic['examples'],
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in dependencies},
        scope=dict(classical_metric_matter_only_one_loop_MS=True,all_listed_fields_active=True,
            no_hard_fermion_masses=True,no_graviton_loops=True,
            inherited_beta_functions_not_new_discoveries=True,
            new_result_is_joint_all_scale_spectral_relation_obstruction=True,
            single_scale_matching_allowed=True,quantum_vacuum_not_solved=True,
            dimension_and_GR_not_derived=True,unified_goal_completed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args(); result = run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors',
                    'SM_plus_s_plus_three_RH_curvature')},ensure_ascii=False))
