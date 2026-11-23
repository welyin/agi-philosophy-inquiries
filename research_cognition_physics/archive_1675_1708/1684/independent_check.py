"""Resolve which heavy-port information survives diagonal or gauge-only measurements."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def transfer(B,left,right,t):
    n=len(B)
    h=np.block([[np.zeros_like(B),B],[B.conj().T,np.zeros_like(B)]])
    # right- to left-handed heavy propagator in a Euclidean momentum frame
    return left.conj()@np.linalg.inv(h+1j*np.sqrt(t)*np.eye(2*n))[n:,:n]@right


def unitary(seed):
    rng=np.random.default_rng(seed)
    return np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]


def main():
    M=np.array([1.,2.,3.])
    lam=.2*np.sqrt(M)
    signs=np.array([[1,1,-1],[1,-1,1]])
    B=np.diag(M)
    t=np.array([0.,1.,4.])
    raw=[]
    diagonal=[]
    residues=[]
    gauge=[]
    tests=[]
    for s in signs:
        eta=lam*s
        z=(lam*eta*M)[None,:]/(t[:,None]+M**2)
        amplitudes=z.sum(axis=1)
        raw.append(amplitudes)
        np.testing.assert_allclose([transfer(B,lam,eta,k) for k in t],amplitudes,atol=1.e-14)
        currents=np.vstack([lam,eta])
        spectrum=np.einsum('ai,bi->iab',currents,currents.conj())
        assert np.linalg.eigvalsh(spectrum).min()>-1.e-15
        residues.append(spectrum)
        diagonal.append(np.array([((abs(currents)**2)/(v+M**2)).sum(axis=1) for v in t]))
        gauge.append([sum(1/(k+M**2)) for k in t])
        # Basis change in the heavy multiplicity space: preserve also the
        # physical port vectors. The nontrivial signs in eta cannot follow
        # from such a basis change with everything else held fixed.
        ul,ur=unitary(442),unitary(643)
        rotated=np.array([transfer(ul@B@ur.conj().T,ur@lam,ul@eta,k) for k in t])
        error=float(max(abs(rotated-amplitudes)))
        assert error<2.e-16
        tests.append(error)
    raw=np.array(raw)
    diagonal=np.array(diagonal)
    np.testing.assert_array_equal(diagonal[0],diagonal[1])
    np.testing.assert_array_equal(gauge[0],gauge[1])
    np.testing.assert_allclose(raw[:,0],[.04,.04],atol=1.e-15)
    reduced=raw/raw[:,:1]
    np.testing.assert_allclose(reduced[:,1],[.4,.6],atol=3.e-15)
    neutral=1/(t+1)
    D=(reduced/neutral)**2
    np.testing.assert_allclose(D[:,1],[.64,1.44],atol=3.e-15)
    # Three signed amplitude data at distinct t reconstruct all three mixed
    # residues when the poles have been determined independently.
    C=1/(t[:,None]+M**2)
    recovered=np.linalg.solve(C,raw.T).T
    np.testing.assert_allclose(recovered,lam**2*M*signs,atol=1.e-14)
    # Conservative worst case when each of four positive reduced RATES can
    # independently vary by eps relative to its true value.
    DA,DB=1.07383,.84014
    odds=DA/DB
    critical=lambda p:(odds**(1/p)-1)/(odds**(1/p)+1)
    erate,eamp=critical(4),critical(8)
    H=(1+erate)/(1-erate)
    np.testing.assert_allclose(DA/H**2,DB*H**2)
    return dict(masses=M.tolist(),common_left=lam.tolist(),two_sink_signs=signs.tolist(),
        equal_port_norms_and_gauge_trace=True,matrix_spectral_min_eigenvalue=0,
        canonical_rotation_errors=tests,spacelike_t=t.tolist(),mixed_amplitudes=raw.tolist(),
        normalized_charged_shapes=reduced.tolist(),shared_neutral_shape=neutral.tolist(),
        observable_double_ratio=D.tolist(),Cauchy_signed_inversion_condition=float(np.linalg.cond(C)),
        guaranteed_separation_independent_fractional_rate_error=erate,
        guaranteed_separation_independent_amplitude_contamination=eamp,
        worst_case_touching_ratio=DA/H**2)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    out=main()
    path=HERE/'independent_results.json'
    if args.check:
        for k,v in json.loads(path.read_bytes()).items():
            np.testing.assert_allclose(v,out[k],rtol=1.e-10,atol=1.e-12)
        print('PASS: basis-invariant mixed ports, equal diagonal data, Cauchy inversion and error margins')
    else:
        with path.open('x',encoding='utf8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
