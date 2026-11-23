"""Independent round 1328 parent-matrix, CG, sources and preparation audit.

No author imports. --write-new creates output exclusively; --check is read-only.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def opnorm(matrix):
    return float(np.linalg.norm(matrix, ord=2))


def evolution(H, time):
    values, vectors = np.linalg.eigh(H)
    return (vectors * np.exp(-1j * time * values)) @ vectors.conj().T


def parent(N, MD, MU, p=0.0):
    n = len(MD)
    states = []
    for role, photons in [('D', N), ('L0', N-1), ('L1', N-1), ('U', N-2)]:
        states += [(role, blue, flavor) for blue in range(photons+1) for flavor in range(n)]
    index = {state: i for i, state in enumerate(states)}
    T = np.zeros((len(states), len(states)), complex)
    for col, (role, blue, flavor) in enumerate(states):
        if role == 'D':
            for i in range(n):
                if blue < N:
                    T[index[('L0', blue, i)], col] += MD[i, flavor] * math.sqrt((N-blue)/N)
                if blue > 0:
                    T[index[('L1', blue-1, i)], col] += MD[i, flavor] * math.sqrt(blue/N)
        elif role == 'U':
            for i in range(n):
                T[index[('L0', blue+1, i)], col] += MU[i, flavor] * math.sqrt((blue+1)/N)
                T[index[('L1', blue, i)], col] -= MU[i, flavor] * math.sqrt((N-1-blue)/N)
    PD = np.diag([float(s[0] == 'D') for s in states])
    PU = np.diag([float(s[0] == 'U') for s in states])
    PL = np.diag([float(s[0].startswith('L')) for s in states])
    NB = N*PD + (N-2)*PU + (N-1)*PL
    V = T + T.conj().T
    H = p*(PD+PU-PL)+V
    return dict(states=states, index=index, T=T, V=V, H=H, PD=PD, PU=PU, PL=PL, NB=NB)


def dictionary(N, n, physical):
    # Dictionary order: reference blue, direction, role(right,left), flavor.
    D = np.zeros((len(physical['states']), 4*N*n), complex)
    Pplus = np.zeros((2*N, 2*N), complex)
    for blue in range(N+1):
        vector = np.zeros(2*N)
        if blue < N:
            vector[2*blue] = math.sqrt((N-blue)/N)
        if blue > 0:
            vector[2*(blue-1)+1] = math.sqrt(blue/N)
        Pplus += np.outer(vector, vector)
    for blue in range(N):
        for direction in (0, 1):
            for role in (0, 1):
                for i in range(n):
                    col = (((blue*2+direction)*2+role)*n+i)
                    if role == 1:
                        D[physical['index'][(f'L{direction}', blue, i)], col] = 1
                    elif direction == 0:
                        D[physical['index'][('D', blue, i)], col] = math.sqrt((N-blue)/N)
                        if blue:
                            D[physical['index'][('U', blue-1, i)], col] = math.sqrt(blue/N)
                    else:
                        D[physical['index'][('D', blue+1, i)], col] = math.sqrt((blue+1)/N)
                        if blue < N-1:
                            D[physical['index'][('U', blue, i)], col] = -math.sqrt((N-1-blue)/N)
    E = D[:, :4*n]
    return D, E, Pplus


def block(M, p=0.0):
    n = len(M)
    return np.block([[p*np.eye(n), M.conj().T], [M, -p*np.eye(n)]])


def directsum(A, B):
    return np.block([[A, np.zeros_like(A)], [np.zeros_like(B), B]])


def compute():
    checks = []
    cases = []
    spectra = []
    source_samples = []

    def close(name, value, expected, tolerance=4e-10):
        residual = float(np.linalg.norm(np.asarray(value)-np.asarray(expected)))
        checks.append(dict(name=name, residual=residual, tolerance=tolerance))
        if residual > tolerance:
            raise AssertionError((name, residual, tolerance))

    def below(name, value, bound, tolerance=4e-10):
        residual = max(0.0, float(value-bound))
        checks.append(dict(name=name, residual=residual, tolerance=tolerance, value=float(value), bound=float(bound)))
        if residual > tolerance:
            raise AssertionError((name, value, bound))

    MD = np.array([[0.37+0.04j, 0.05-0.03j], [0.08+0.11j, 0.24-0.02j]])
    MU = np.array([[0.19-0.01j, 0.09+0.07j], [-0.04+0.02j, 0.31+0.03j]])
    assert opnorm(MD@MU-MU@MD) > 0.01
    n = 2
    zero = np.zeros_like(MD)
    PR = np.diag([1.0]*n+[0.0]*n)
    PL = np.eye(2*n)-PR
    spin_D = np.diag([1.0, 0.0])
    spin_U = np.diag([0.0, 1.0])
    N0 = 3
    calibration_parent = parent(N0, MD, MU)
    _, E0, _ = dictionary(N0, n, calibration_parent)
    compression = E0.conj().T @ calibration_parent['H'] @ E0
    A_up = compression[n:2*n, :n]
    A_down = compression[3*n:4*n, 2*n:3*n]
    calibrated_D = A_up
    calibrated_U = (N0*A_down-A_up)/(N0-1)
    close('one_calibration.MD', calibrated_D, MD)
    close('one_calibration.MU', calibrated_U, MU)
    # These matrices are held fixed in every later test.
    MD, MU = calibrated_D, calibrated_U
    for N in (2, 3, 7, 17):
        for p in (0.0, 0.13):
            physical = parent(N, MD, MU, p)
            D, E, Pplus = dictionary(N, n, physical)
            Pminus = np.eye(2*N)-Pplus
            hD, hU = block(MD, p), block(MU, p)
            h = directsum(hD, hU)
            factor = np.kron(Pplus, hD) + np.kron(Pminus, hU)
            tag = f'N{N}.p{p}'
            assert D.shape == (4*N*n, 4*N*n)
            close(tag+'.whole_dictionary', D.conj().T@D, np.eye(4*N*n))
            close(tag+'.whole_parent_factor', physical['H']@D, D@factor)
            close(tag+'.projector', Pplus@Pplus, Pplus)
            close(tag+'.rank_plus', np.trace(Pplus), N+1)
            close(tag+'.rank_minus', np.trace(Pminus), N-1)
            close(tag+'.embedding', E.conj().T@E, np.eye(4*n))
            close(tag+'.exact_generator_error', opnorm(physical['H']@E-E@h), opnorm(hD-hU)/math.sqrt(N))
            compressed_expected = directsum(hD, hU+(hD-hU)/N)
            close(tag+'.actual_compression_not_target', E.conj().T@physical['H']@E, compressed_expected)
            close(tag+'.D_count_source', physical['PD']@D, D@np.kron(Pplus, PR))
            close(tag+'.U_count_source', physical['PU']@D, D@np.kron(Pminus, PR))
            close(tag+'.L_count_source', physical['PL']@D, D@np.kron(np.eye(2*N), PL))
            close(tag+'.reference_identity', physical['NB']+physical['PL']+2*physical['PU'], N*np.eye(4*N*n))
            close(tag+'.momentum_source', (physical['PD']+physical['PU']-physical['PL'])@E,
                  E@np.kron(np.eye(2), PR-PL))
            sources = []
            for department in ('D', 'U'):
                for i in range(n):
                    for j in range(n):
                        for component, z in [('real', 1.0), ('imag', 1j)]:
                            derivative = np.zeros_like(MD)
                            derivative[i,j] = z
                            smD = derivative if department == 'D' else zero
                            smU = derivative if department == 'U' else zero
                            source = parent(N, smD, smU)['V']
                            sD, sU = block(smD), block(smU)
                            s = directsum(sD, sU)
                            name = f'{tag}.source{department}{i}{j}.{component}'
                            close(name+'.full_dictionary', source@D, D@(np.kron(Pplus,sD)+np.kron(Pminus,sU)))
                            close(name+'.vector_error', opnorm(source@E-E@s), opnorm(sD-sU)/math.sqrt(N))
                            close(name+'.compression_error', opnorm(E.conj().T@source@E-s), opnorm(sD-sU)/N)
                            if department == 'D' and i == 0 and j == 1 and component == 'imag':
                                sources.append((source,s,sD,sU))
            # Existing physical preparation: unknown natural left, then complete equal-g swap.
            left_columns = list(range(n,2*n))+list(range(3*n,4*n))
            right_columns = list(range(n))+list(range(2*n,3*n))
            EL, ER = E[:,left_columns], E[:,right_columns]
            prep = parent(N, np.eye(n)*0.7, np.eye(n)*0.7, 0)
            prep_time = math.pi/(2*0.7)
            close(tag+'.actual_right_preparation', evolution(prep['H'],prep_time)@EL, -1j*ER)
            close(tag+'.natural_left_only_reference', physical['NB']@EL, (N-1)*EL)
            # Source and preparation are held fixed across the three prediction times.
            time_rows = []
            for time in (0.37, 2.1, 9.0):
                U = evolution(physical['H'],time)
                u = evolution(h,time)
                delta = opnorm(U@E-E@u)
                exact = opnorm(evolution(hD,time)-evolution(hU,time))/math.sqrt(N)
                close(tag+f'.time{time}.exact_process_error', delta, exact)
                below(tag+f'.time{time}.uniform_bound', delta, min(2,time*opnorm(MD-MU))/math.sqrt(N))
                for source,s,sD,sU in sources:
                    actual_source = (U@E).conj().T@source@(U@E)
                    effective_source = u.conj().T@s@u
                    observable_error = opnorm(actual_source-effective_source)
                    observable_bound = 2*opnorm(source)*delta+opnorm(sD-sU)/N
                    below(tag+f'.time{time}.evolved_source_bound', observable_error, observable_bound)
                time_rows.append(dict(time=time, whole_state_error=delta, exact_formula=exact))
            # Maximal entanglement with a reference checks every prepared left/right input column.
            H0 = 3.0*np.eye(4*N*n)+0.6*physical['NB'] + 0.6*physical['PL'] + 1.2*physical['PU']
            H0 += p*(physical['PD']+physical['PU']-physical['PL'])
            Hfull = H0+physical['V']
            close(tag+'.independent_energy_sources', Hfull, (3.0+0.6*N)*np.eye(4*N*n)+physical['H'])
            for input_name, initial in [('natural_left',EL),('swap_right',-1j*ER)]:
                initial = initial/math.sqrt(2*n)
                out = evolution(physical['H'],2.1)@initial
                work_off = -float(np.vdot(out,physical['V']@out).real)
                source_after_off = float(np.vdot(out,H0@out).real)
                energy_before_off = float(np.vdot(out,Hfull@out).real)
                close(tag+'.'+input_name+'.quench_work', source_after_off-energy_before_off, work_off)
                branches = [physical[key]@out for key in ('PD','PU','PL')]
                close(tag+'.'+input_name+'.record_norm', sum(np.vdot(b,b).real for b in branches),1)
                close(tag+'.'+input_name+'.record_postquench_energy',sum(np.vdot(b,H0@b).real for b in branches),source_after_off)
                close(tag+'.'+input_name+'.full_energy', np.vdot(out,Hfull@out),np.vdot(initial,Hfull@initial))
                delta_reference = np.vdot(out,physical['NB']@out)-np.vdot(initial,physical['NB']@initial)
                delta_left = np.vdot(out,physical['PL']@out)-np.vdot(initial,physical['PL']@initial)
                delta_up = np.vdot(out,physical['PU']@out)-np.vdot(initial,physical['PU']@initial)
                close(tag+'.'+input_name+'.actual_reference_balance',delta_reference+delta_left+2*delta_up,0)
                source_samples.append(dict(N=N,p=p,input=input_name,work_off=work_off,
                                           probabilities=[float(np.vdot(b,b).real) for b in branches]))
            # A D/U/L record itself has a dictionary error, separate from propagation error.
            actual_records = np.concatenate([physical[key]@E for key in ('PD','PU','PL')],axis=0)
            target_records = np.concatenate([E@np.kron(spin_D,PR),E@np.kron(spin_U,PR),E@np.kron(np.eye(2),PL)],axis=0)
            close(tag+'.record_dictionary_extra_error',opnorm(actual_records-target_records),math.sqrt(2/N))
            cases.append(dict(N=N,n=n,p=p,parent_dimension=4*N*n,all_time_bound=2/math.sqrt(N),times=time_rows))

    for N in (2,3,7):
        d,u=0.8,0.3
        physical=parent(N,np.array([[d]]),np.array([[u]]))
        values,vectors=np.linalg.eigh(physical['H'])
        expected=np.repeat([-d,-u,u,d],[N+1,N-1,N-1,N+1])
        close(f'spectral.N{N}.parent_multiplicity',values,expected)
        optimal=np.repeat([-d,-(d+u)/2,(d+u)/2,d],N)
        fitted=(vectors*optimal)@vectors.conj().T
        close(f'spectral.N{N}.attained_minimax',opnorm(physical['H']-fitted),(d-u)/2)
        fixed=np.repeat([-d,-u,u,d],N)
        fixed_fit=(vectors*fixed)@vectors.conj().T
        close(f'spectral.N{N}.fixed_spectrum_gap',opnorm(physical['H']-fixed_fit),d-u)
        _,E,_=dictionary(N,1,physical)
        initial=E[:,3:4]  # natural left, down
        for time in (0.4,1.7,4.0):
            out=evolution(physical['H'],time)@initial
            actual_D=float(np.vdot(out,physical['PD']@out).real)
            actual_U=float(np.vdot(out,physical['PU']@out).real)
            close(f'spectral.N{N}.t{time}.D_flow',actual_D,math.sin(d*time)**2/N)
            close(f'spectral.N{N}.t{time}.U_flow',actual_U,(N-1)*math.sin(u*time)**2/N)
            close(f'spectral.N{N}.t{time}.reference_flow',np.vdot(out,physical['NB']@out).real-(N-1),actual_D-actual_U)
            h=directsum(block(np.array([[d]])),block(np.array([[u]])))
            close(f'spectral.N{N}.t{time}.uniform_exact',opnorm(evolution(physical['H'],time)@E-E@evolution(h,time)),2*abs(math.sin((d-u)*time/2))/math.sqrt(N))
        spectra.append(dict(N=N,d=d,u=u,minimax_error=(d-u)/2,fixed_target_error=d-u))

    return dict(round=1328,status='PASS',stage_complete=False,new_axioms=0,
                independent_construction='Direct physical occupation transition matrix and inverse CG coefficients; no author imports.',
                check_count=len(checks),max_residual=max(c['residual'] for c in checks),
                calibration=dict(N0=N0,MD_real=MD.real.tolist(),MD_imag=MD.imag.tolist(),MU_real=MU.real.tolist(),MU_imag=MU.imag.tolist()),
                cases=cases,spectral_cases=spectra,source_record_samples=source_samples,checks=checks)


def compare(a,b,path='root'):
    if isinstance(b,dict):
        assert a.keys()==b.keys(),path
        for key in b: compare(a[key],b[key],path+'.'+key)
    elif isinstance(b,list):
        assert len(a)==len(b),path
        for i,(aa,bb) in enumerate(zip(a,b)): compare(aa,bb,path+f'[{i}]')
    elif isinstance(b,float):
        assert math.isfinite(a) and abs(a-b)<=1e-9*(1+abs(b)),path
    else:
        assert a==b,path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-new',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=compute()
    output=HERE/'independent_results.json'
    if args.write_new:
        with output.open('x',encoding='utf-8') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2)
            handle.write('\n')
    else:
        compare(json.loads(output.read_text(encoding='utf-8')),result)
    print(json.dumps({key:result[key] for key in ('round','status','check_count','max_residual','spectral_cases')},ensure_ascii=False))


if __name__=='__main__':
    main()
