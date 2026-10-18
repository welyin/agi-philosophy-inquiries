"""1288: exact cyclic Gauss parents of the fixed 958/959 material contract.

Default: read-only recomputation against results.json.  --write creates only a
previously absent result.  Old scientific files and their historical hashes
are never changed.  General N/time/partner quantifiers are proved in proof.md.
"""
from pathlib import Path
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
ROOT = RESEARCH.parent
OLD = RESEARCH / "archive_956_989"
TARGET = HERE / "results.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def fro(a):
    """Frobenius residual; no uncomputed diamond norm is claimed."""
    return float(np.linalg.norm(a))


def require(error, tolerance=2e-11):
    assert error < tolerance, (error, tolerance)
    return error


def evolution(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * time * values)) @ vectors.conj().T


def load_old_material():
    path = OLD / "956/native_material_interface.py"
    spec = importlib.util.spec_from_file_location("round956_material_1288", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cyclic_link(n):
    shift = np.zeros((n, n))
    for e in range(n):
        shift[(e + 1) % n, e] = 1
    electric = np.array([e if e <= n // 2 else e - n for e in range(n)], float)
    return shift, electric


def dressing(charges, n):
    charges = np.asarray(charges, int)
    indices = np.arange(len(charges)) * n + charges % n
    j = np.zeros((len(charges) * n, len(charges)))
    j[indices, np.arange(len(charges))] = 1
    return j, indices


def lift_operator(m, charges, n, shift):
    d = len(charges)
    out = np.zeros((d * n, d * n), complex)
    for s in range(d):
        for t in range(d):
            if m[s, t] != 0:
                out[s*n:(s+1)*n, t*n:(t+1)*n] = (
                    m[s, t] * np.linalg.matrix_power(shift, int(charges[s] - charges[t])))
    return out


def parent(n, charges, tp, onsite, hopping, coupling):
    shift, electric = cyclic_link(n)
    dl = len(charges) * n
    j, indices = dressing(charges, n)
    local_t = np.kron(tp, shift)
    local_hop = -hopping * (local_t + local_t.T)
    e = np.tile(electric, len(charges))
    q = np.repeat(charges, n)
    residue = np.tile(np.arange(n), len(charges))
    gauss = np.exp(2j * np.pi * (residue - q) / n)
    source_u = e[:, None]**2 + e[None, :]**2
    source_k = e[:, None] * e[None, :]
    energy = onsite * source_u + coupling * source_k

    def apply_h(x, electric_part=True):
        tensor = np.asarray(x).reshape(dl, dl, -1)
        out = np.einsum("ia,ajk->ijk", local_hop, tensor)
        out += np.einsum("jb,ibk->ijk", local_hop, tensor)
        if electric_part:
            out += energy[:, :, None] * tensor
        return out.reshape(dl * dl, -1)

    joint_j = np.kron(j, j)
    return dict(n=n, dim=dl*dl, local_dim=dl, shift=shift, electric=e,
                j=j, joint_j=joint_j, gauss=gauss, local_hop=local_hop,
                source_u=source_u, source_k=source_k, apply_h=apply_h,
                local_indices=indices)


def old_record_data(material, saved):
    b = np.eye(16)
    code = np.stack([(b[5]-b[6]-b[9]+b[10])/2,
                     (2*b[3]+2*b[12]-b[5]-b[6]-b[9]-b[10])/np.sqrt(12)], axis=1)
    local = np.kron(material["W"], np.eye(4)) @ code
    time = saved["exact_interface"]["controlled_phase_time"]
    target = (np.exp(1j * material["J"] * time) * local[:, 0] + local[:, 1]) / np.sqrt(2)
    return local, target


def act_material(unitary, x):
    a = x.reshape(6, 4, 6, 4, -1).transpose(0, 2, 1, 3, 4)
    a = (unitary @ a.reshape(36, -1)).reshape(6, 6, 4, 4, -1)
    return a.transpose(0, 2, 1, 3, 4).reshape(576, -1)


def old_probabilities(h, time, local, target):
    unitary = evolution(h, time)
    pair_code = np.kron(local, local)
    plus = np.ones(2) / np.sqrt(2)
    probabilities = []
    for bit in (0, 1):
        initial = pair_code @ np.kron(np.eye(2)[:, bit], plus)
        final = act_material(unitary, initial).reshape(24, 24)
        probabilities.append(float(np.linalg.norm(final @ target.conj())**2))
    return probabilities


def check_instrument(charges, n, target, rng):
    """The old actual terminal effect and its complement, no new read port."""
    q = np.repeat(charges, 4)
    shift, _ = cyclic_link(n)
    j, indices = dressing(q, n)
    effect = np.outer(target, target.conj())
    lifted = lift_operator(effect, q, n, shift)
    complement = np.eye(24*n) - lifted
    local_gauss = np.exp(2j*np.pi*(np.tile(np.arange(n),24)-np.repeat(q,n))/n)
    completeness = require(fro(lifted.conj().T @ lifted + complement.conj().T @ complement - np.eye(24*n)))
    intertwining = require(fro(lifted @ j - j @ effect))
    invariance = require(fro((local_gauss[:, None]-local_gauss[None, :])*lifted))
    # A sample supplement to the all-R operator proof.  No random experiment
    # is used to claim its universal quantifier.
    psi = rng.normal(size=(24,24,3)) + 1j*rng.normal(size=(24,24,3))
    psi /= np.linalg.norm(psi)
    encoded = np.zeros((24*n,24*n,3), complex)
    encoded[indices[:, None], indices[None, :], :] = psi
    old_rho_r = psi.reshape(576,3).conj().T @ psi.reshape(576,3)
    recovered_rho_r = np.zeros((3,3), complex)
    errors = []
    for old_k, new_k in ((effect,lifted),(np.eye(24)-effect,complement)):
        after = np.einsum("jb,ibk->ijk", new_k, encoded)
        expected = np.einsum("jb,ibk->ijk", old_k, psi)
        pulled = after[indices[:,None],indices[None,:],:]
        errors.append(require(fro(pulled-expected)))
        require(abs(float(np.linalg.norm(after)**2)-float(np.linalg.norm(expected)**2)))
        flat = after.reshape(-1,3)
        recovered_rho_r += flat.conj().T @ flat
    marginal = require(fro(recovered_rho_r-old_rho_r))
    return dict(completeness_fro=completeness,branch_intertwining_fro=intertwining,
                gauss_fro=invariance,complete_QR_branch_sample_fro=errors,
                nonselective_passive_R_marginal_fro=marginal,
                conditional_R_or_system_R_unchanged_claimed=False)


def run():
    old = load_old_material()
    saved = read(OLD / "958/capacitive_material_write_results.json")
    onsite = saved["parameters"]["U"]
    hopping = saved["parameters"]["hopping"]
    coupling = saved["parameters"]["kappa"]
    assert (onsite,hopping,coupling) == (1.0,0.1,0.2)
    material = old.material(onsite, hopping)
    cs = [old.annihilate(i) for i in range(4)]
    ns = [c.T@c for c in cs]
    sector = np.eye(16)[:, [s for s in range(16) if s.bit_count()==2]]
    charge = sector.T @ (ns[0]+ns[1]-ns[2]-ns[3]) @ sector / 2
    charges = np.diag(charge).astype(int)
    tp = sector.T @ sum(cs[i].T@cs[i+2] for i in (0,1)) @ sector
    require(fro(charge@tp-tp@charge-tp))
    require(fro(material["D"]-charge@charge))
    im = np.eye(6)
    hop_source = -np.kron(tp+tp.T,im)-np.kron(im,tp+tp.T)
    u_source = np.kron(charge@charge,im)+np.kron(im,charge@charge)
    k_source = np.kron(charge,charge)
    hm = hopping*hop_source+onsite*u_source+coupling*k_source
    require(fro(hm-(np.kron(material["H"],im)+np.kron(im,material["H"])+coupling*k_source)))
    local, target = old_record_data(material, saved)
    rng = np.random.default_rng(1288)
    rows = []
    for n in (3,4,5,7):
        p = parent(n, charges, tp, onsite, hopping, coupling)
        j = p["joint_j"]
        hj = p["apply_h"](j)
        hp = j.T @ hj
        gauss = p["gauss"]
        physical_dim = int(np.count_nonzero(np.abs(gauss-1)<1e-12))**2
        assert physical_dim == 36
        source_errors = dict(
            U=require(fro(p["source_u"].reshape(-1,1)*j-j@u_source)),
            v=require(fro(p["apply_h"](j,False)/hopping-j@hop_source)),
            kappa=require(fro(p["source_k"].reshape(-1,1)*j-j@k_source)))
        generator = require(fro(hj-j@hm))
        gauss_error = require(fro((gauss[:,None]-gauss[None,:])*p["local_hop"]))
        isometry = require(fro(j.T@j-np.eye(36)))
        hermiticity = require(fro(hp-hp.T))
        table = []
        for old_row in saved["exact_interface"]["finite_time_rows"]:
            probabilities = old_probabilities(hp,old_row["time"],local,target)
            old_p = old_row["receiver_probabilities_for_logical_0_1"]
            err = max(abs(a-b) for a,b in zip(probabilities,old_p))
            require(err,1e-8)
            table.append(dict(time=old_row["time"],receiver_probabilities=probabilities,
                              difference_from_original=err))
        full_propagator = None
        if n == 3:
            h_full = p["apply_h"](np.eye(p["dim"]))
            require(fro(h_full-h_full.T))
            full_propagator = require(fro(evolution(h_full,0.7)@j-j@evolution(hm,0.7)),1e-10)
        rows.append(dict(N=n,kinematical_dimension=p["dim"],physical_dimension=physical_dim,
                         isometry_fro=isometry,generator_fro=generator,gauss_fro=gauss_error,
                         hermiticity_fro=hermiticity,source_fro=source_errors,
                         full_ambient_propagator_fro_at_t_07=full_propagator,
                         original_read_window=table,
                         terminal_instrument=check_instrument(charges,n,target,rng)))
    # N=2 still has an exact isometric physical code and a complete *-lift.
    p2 = parent(2,charges,tp,onsite,hopping,coupling)
    j2 = p2["joint_j"]
    n2_source_error = p2["source_k"].reshape(-1,1)*j2-j2@k_source
    n2_h_error = p2["apply_h"](j2)-j2@hm
    p20 = parent(2,charges,tp,onsite,hopping,0.)
    zero_k_error = require(fro(p20["apply_h"](j2)-j2@(hopping*hop_source+onsite*u_source)))
    n2 = dict(physical_dimension=int(np.count_nonzero(np.abs(p2["gauss"]-1)<1e-12))**2,
              isometry_fro=require(fro(j2.T@j2-np.eye(36))),
              onsite_source_fro=require(fro(p2["source_u"].reshape(-1,1)*j2-j2@u_source)),
              same_hopping_and_onsite_at_kappa_zero_fro=zero_k_error,
              actual_cross_source_operator_error=float(np.linalg.norm(n2_source_error,2)),
              actual_H_operator_error=float(np.linalg.norm(n2_h_error,2)),
              best_any_pure_link_source_collision_lower=str(Fraction(1)),
              best_any_pure_link_H_collision_lower_at_original_kappa=str(Fraction(1,5)),
              terminal_instrument=check_instrument(charges,2,target,rng),
              all_Z2_realizations_excluded=False)
    assert n2["physical_dimension"]==36
    assert abs(n2["actual_cross_source_operator_error"]-2)<1e-12
    assert abs(n2["actual_H_operator_error"]-.4)<1e-12
    # Integer diagnostics only; the proposition for all finite S is analytic.
    menu_rows = []
    for menu in ((-1,0,1),(0,3,7),(-3,0,3)):
        diam = max(menu)-min(menu)
        for n in range(1,diam+4):
            injective = len({q % n for q in menu})==len(menu)
            if n > diam:
                assert injective
            menu_rows.append(dict(charges=list(menu),N=n,injective=injective,
                                  diameter_sufficient=n>diam))
    sources = [Path(__file__),OLD/"956/native_material_interface.py",
               OLD/"958/capacitive_material_write_results.json",
               OLD/"959/gauss_material_parent.py",OLD/"research_note_958.md",
               OLD/"research_note_959.md"]
    return dict(round=1288,date="2026-10-18",all_scientific_checks_passed=True,
                fixed_parameters=dict(U=onsite,v=hopping,kappa=coupling,
                    original_read_time=saved["exact_interface"]["controlled_phase_time"],
                    inherited_unknown_input_error_upper=saved["exact_interface"]["all_time_unknown_input_error_upper"],
                    inherited_read_window_contrast_lower=saved["exact_interface"]["analytic_receiver_contrast_lower"]),
                cyclic_parent_checks=rows,N2_source_collision=n2,
                finite_charge_menu_checks=menu_rows,
                scope=dict(all_N_ge3_and_all_times_are_analytic=True,
                    whole_declared_36d_physical_sector=True,arbitrary_passive_partner_by_map_identity=True,
                    original_effective_preparation_and_read_permissions_only=True,
                    only_original_parameter_sources_certified=True,
                    all_cyclic_groups_excluded=False,U1_uniquely_selected=False,
                    full_stress_tensor_or_Lorentz_QFT=False,new_microscopic_instrument_built=False,
                    all_current_axioms_countermodel=False,stage_goal_completed=False),
                source_hashes={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sources})


def compare(actual, saved, location="root"):
    if isinstance(actual,dict):
        assert set(actual)==set(saved),location
        for key,value in actual.items():
            compare(value,saved[key],location+"."+key)
    elif isinstance(actual,list):
        assert len(actual)==len(saved),location
        for i,(a,b) in enumerate(zip(actual,saved)):
            compare(a,b,location+f"[{i}]")
    elif isinstance(actual,float):
        assert math.isclose(actual,saved,rel_tol=1e-11,abs_tol=1e-9),(location,actual,saved)
    else:
        assert actual==saved,(location,actual,saved)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    if args.write and TARGET.exists():
        raise FileExistsError("results.json exists; refusing to overwrite scientific output")
    result=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
        status="CREATED"
    else:
        compare(result,read(TARGET))
        status="MATCH"
    print(json.dumps(dict(status=status,round=1288,
        N=[row["N"] for row in result["cyclic_parent_checks"]],
        N2_cross_source_error=result["N2_source_collision"]["actual_cross_source_operator_error"],
        N2_kappa_zero_error=result["N2_source_collision"]["same_hopping_and_onsite_at_kappa_zero_fro"],
        all_scientific_checks_passed=True),ensure_ascii=False))
