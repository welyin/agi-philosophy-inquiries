"""595 entry: positive-energy battery fibre dilation diagnostics.

The finite spectral fixture tests the general construction, not eigenvalues
of the original full joint H0. Its application there is analytic. No new
matter species or local autonomous device is inferred from these matrices.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'positive_battery_entry_results.json'


def profile(e,R,L):
    x=(e-R)/L
    return np.where((x>0)&(x<1),np.sqrt(2/L)*np.sin(np.pi*x),0.)


def integrate_segments(points,n=32):
    z,w=np.polynomial.legendre.leggauss(n)
    p=np.unique(points)
    x=np.concatenate([(a+b)/2+(b-a)*z/2 for a,b in zip(p[:-1],p[1:])])
    weights=np.concatenate([(b-a)*w/2 for a,b in zip(p[:-1],p[1:])])
    return x,weights


def sine_profile_checks():
    rows=[]
    R,L=3.,8.
    for shift in (0.,.4,1.7,3.):
        x,w=integrate_segments([0,R-shift,R,R+L-shift,R+L,R+L+shift])
        f=profile(x,R,L);g=profile(x+shift,R,L)
        overlap=float(w@(f*g))
        ratio=shift/L
        expected=(1-ratio)*np.cos(np.pi*ratio)+np.sin(np.pi*ratio)/np.pi
        assert abs(overlap-expected)<2e-14
        error=float(np.sqrt(w@((f-g)**2)))
        assert error<=np.pi*shift/L+1e-13
        rows.append(dict(shift=shift,overlap=overlap,analytic_overlap=float(expected),
                         vector_translation_error=error,derivative_bound=np.pi*shift/L))
    x,w=integrate_segments([R,R+L])
    f=profile(x,R,L)
    mean=float(w@(x*f*f));var=float(w@((x-mean)**2*f*f))
    assert abs(mean-(R+L/2))<2e-14
    assert abs(var-L*L*(1/12-1/(2*np.pi*np.pi)))<3e-14
    derivative=np.sqrt(2/L)*np.pi/L*np.cos(np.pi*(x-R)/L)
    assert abs(float(w@(derivative*derivative))-(np.pi/L)**2)<1e-14
    return dict(rows=rows,mean=mean,variance=var,positive_energy_support=[R,R+L])


def energy_fibre_fixture():
    # Deliberately non-equally spaced diagnostic energies. These are not
    # numerical eigenvalues of the research model's joint H0.
    energies=np.array([.2,.9,2.1]);E=np.repeat(energies,2)
    A=np.array([[.70,.065,.015],[.065,.83,.035],[.015,.035,.94]])
    assert np.linalg.eigvalsh(A)[0]>np.pi/6
    assert np.linalg.eigvalsh(A)[-1]<np.pi/3
    sy=np.array([[0,-1j],[1j,0]])
    vals,vecs=np.linalg.eigh(np.kron(A,sy))
    U=(vecs*np.exp(-1j*vals))@vecs.conj().T
    assert np.linalg.norm(U.conj().T@U-np.eye(6))<5e-15
    psi=np.zeros((6,2),complex)
    psi[0,0]=np.sqrt(.8)
    psi[2,1]=np.sqrt(.15)
    psi[4,0]=1j*np.sqrt(.05)
    target=U@psi
    input_E=float(np.sum(E[:,None]*abs(psi)**2))
    input_E2=float(np.sum(E[:,None]**2*abs(psi)**2))
    ideal_E=float(np.sum(E[:,None]*abs(target)**2))
    R=3.;rows=[]
    for L in (4.,8.,16.,32.):
        shifts=E[:,None]-E[None,:]
        bounds=[0.,R+L+R]
        bounds.extend((R-shifts).ravel());bounds.extend((R+L-shifts).ravel())
        x,w=integrate_segments(bounds,40)
        f=profile(x[:,None,None]+shifts[None,:,:],R,L)
        out=np.einsum('mn,nr,xmn->xmr',U,psi,f)
        probability=np.sum(abs(out)**2,axis=2)
        norm=float(np.sum(w[:,None]*probability))
        out_E=float(np.sum(w[:,None]*probability*E[None,:]))
        out_battery=float(np.sum(w[:,None]*probability*x[:,None]))
        total_square=float(np.sum(w[:,None]*probability*(E[None,:]+x[:,None])**2))
        battery_mean=R+L/2
        battery_var=L*L*(1/12-1/(2*np.pi*np.pi))
        input_total_square=input_E2+2*input_E*battery_mean+battery_var+battery_mean**2
        assert abs(norm-1)<8e-15
        assert abs(out_E+out_battery-input_E-battery_mean)<3e-13
        assert abs(total_square-input_total_square)<8e-12
        direct=profile(x,R,L)[:,None,None]*target[None,:,:]
        vector_error=float(np.sqrt(np.sum(w[:,None,None]*abs(out-direct)**2)))
        assert vector_error<np.pi*R/L
        flat=out.reshape(len(x),-1)
        reduced=np.einsum('x,xi,xj->ij',w,flat,flat.conj())
        target_vec=target.ravel()
        delta=reduced-np.outer(target_vec,target_vec.conj())
        distance=float(np.sum(abs(np.linalg.eigvalsh(delta)))/2)
        assert distance<=vector_error+1e-13
        rows.append(dict(width=L,norm=norm,system_energy_gain=out_E-input_E,
                         ideal_system_energy_gain=ideal_E-input_E,
                         battery_energy_change=out_battery-battery_mean,
                         total_energy_error=out_E+out_battery-input_E-battery_mean,
                         total_energy_square_error=total_square-input_total_square,
                         vector_error=vector_error,reference_retaining_trace_distance=distance,
                         uniform_battery_error_bound=np.pi*R/L))
    assert rows[-1]['vector_error']<rows[0]['vector_error']/6
    assert rows[-1]['reference_retaining_trace_distance']<rows[0]['reference_retaining_trace_distance']/30
    return dict(diagnostic_energies=energies.tolist(),rows=rows,
                spectrum_not_claimed_to_be_original_H0=True,
                half_line_boundary_completed_fibrewise=True)


def run():
    names=('research_note_594.md','research_note_592.md','research_note_401.md',
           'round594_drafts/premeasurement_energy_entry.py')
    return dict(status='595 entry; not completed numbered round',checks_passed=True,
                sine_profile=sine_profile_checks(),energy_fibre=energy_fibre_fixture(),
                original_model_application_is_analytic=True,
                full_original_H0_spectrum_and_propagation_not_computed=True,
                half_infinite_clock_is_declared_extra_hardware=True,
                geometry_local_realization_and_resource_preparation_not_derived=True,
                dependency_hashes={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
