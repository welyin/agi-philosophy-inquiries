"""Round 1353. Full source evolution and the three-pack permutation test."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def generators(group):
    result = []
    for i, j in itertools.combinations(range(3), 2):
        E = np.zeros((3, 3))
        E[i, j] = 1
        if group == "SU3":
            result.extend([(E + E.T) / 2, 1j * (E.T - E) / 2])
        else:
            result.append(1j * (E.T - E))
    if group == "SU3":
        result.extend([np.diag([1., -1., 0.]) / 2, np.diag([1., 1., -2.]) / np.sqrt(12)])
    return result


def source_blocks(theta, group, second_arm_phase=0.):
    T = generators(group)
    cr = sum(np.trace(t @ t).real for t in T) / 3
    overlap = np.exp(1j * second_arm_phase)
    coupling = np.concatenate([(np.kron(t, np.eye(3)) + overlap * np.kron(np.eye(3), t)) /
                               np.sqrt(cr) for t in T])
    H = np.zeros((9 + coupling.shape[0],) * 2, complex)
    H[9:, :9] = coupling
    H[:9, 9:] = coupling.conj().T
    w, vec = np.linalg.eigh(H)
    all_blocks = (vec * np.exp(-1j * theta * w)) @ vec[:9].conj().T
    assert np.allclose(all_blocks.conj().T @ all_blocks, np.eye(9))
    return all_blocks[:9], all_blocks[9:]


def antisymmetrize(vectors):
    # Orthonormal coordinates for their joint span give an exact small, pure-state
    # compression of the 3 x 81 tensor product. Never discard a partner leg.
    basis = np.linalg.qr(vectors)[0]
    small = basis.conj().T @ vectors
    assert np.allclose(basis @ small, vectors)
    word = np.einsum("a,b,c->abc", *small.T)
    result = np.zeros_like(word)
    for permutation in itertools.permutations(range(3)):
        inversions = sum(a > b for a, b in itertools.combinations(permutation, 2))
        result += (-1)**inversions * word.transpose(permutation) / 6
    return float(np.vdot(result, result).real)


def marginal_antisymmetry(rhos):
    a, b, c = rhos
    value = (np.trace(a)*np.trace(b)*np.trace(c) -
             np.trace(a@b)*np.trace(c) - np.trace(a@c)*np.trace(b) -
             np.trace(b@c)*np.trace(a) + 2*np.trace(a@b@c).real)
    return float(value.real/6)


def compute():
    swap = np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
    baseline = np.eye(9) / 3
    exchanged = swap @ baseline
    first_two = np.column_stack([baseline.ravel(), exchanged.ravel()])
    projector = first_two @ np.linalg.solve(first_two.conj().T @ first_two, first_two.conj().T)
    outputs = {}
    for group in ["SU3", "SO3", "mu3SO3"]:
        times = {}
        for name, theta in [("first",np.pi/np.sqrt(3)), ("holdout",np.pi/(2*np.sqrt(3)))]:
            base, exc = source_blocks(theta, group)
            success = base @ baseline
            accepted = float(np.linalg.norm(success)**2)
            rejected = float(np.linalg.norm(exc @ baseline)**2)
            assert np.isclose(accepted + rejected, 1.)
            vectors = np.column_stack([first_two, success.ravel()])
            joint = antisymmetrize(vectors)
            residual = float(np.linalg.norm((np.eye(81)-projector) @ vectors[:,2])**2)
            expected = 0. if group=="SU3" else 10/729*(1-np.cos(np.sqrt(3)*theta))**2
            assert abs(joint - expected) < 1e-13
            assert abs(joint - (4/27)*residual) < 1e-13
            mats = [baseline, exchanged, success]
            active_only = marginal_antisymmetry([q @ q.conj().T for q in mats])
            assert np.isclose(active_only,28*accepted/243)
            tag_effects = {}
            for overlap_sq in [0.,.5,1.]:
                tagged = np.column_stack([
                    np.kron([1.,0.], vectors[:,0]), np.kron([1.,0.], vectors[:,1]),
                    np.kron(np.sqrt([overlap_sq,1-overlap_sq]),vectors[:,2])])
                p_tag = antisymmetrize(tagged)
                assert abs(p_tag - (overlap_sq*joint + (1-overlap_sq)*4*accepted/27))<1e-13
                tag_effects[str(overlap_sq)] = p_tag
            wrong_base,_ = source_blocks(theta,group,second_arm_phase=np.pi/2)
            wrong_signal=antisymmetrize(np.column_stack([first_two,(wrong_base @ baseline).ravel()]))
            assert wrong_signal < 1e-24
            times[name] = dict(returned=accepted, source_excited=rejected, joint=joint,
                               returned_other=accepted-joint, residual_norm_squared=residual,
                               gram_determinant=float(np.linalg.det(vectors.conj().T@vectors).real),
                               reduced_active_legs_false_signal=active_only,
                               third_tag_overlap_squared_to_joint=tag_effects,
                               incoincident_arm_phase_joint=wrong_signal)
        outputs[group] = times
    nshots=math.ceil(math.log(200)/(2*.0025**2))
    return dict(groups=outputs, shots_for_pointwise_two_tail=nshots,
                smallest_predicted_signal=10/729,
                source_contract="Source overlaps and independent partner access required")


def compare(actual,saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for key in actual:
            compare(actual[key],saved[key])
    elif isinstance(actual,str):
        assert actual == saved
    else:
        assert np.isclose(actual,saved,atol=1e-12,rtol=1e-12),(actual,saved)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=Path(__file__).with_name("independent_results.json")
    result=compute()
    if args.check:
        compare(result,json.loads(out.read_text(encoding="utf8")))
    else:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    print("PASS",result["groups"]["SO3"]["first"]["joint"],
          result["groups"]["SO3"]["holdout"]["joint"])
