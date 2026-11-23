"""Partial trace of two difference interactions and a test of literal quark/lepton sharing."""
import argparse
import json
from pathlib import Path
import runpy
import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
PARENT = runpy.run_path(str(HERE / CFG["parent"]))
I = np.eye(3)


def spectral(h):
    values, frame = np.linalg.eigh(h)
    values = values / values[-1]
    assert min(np.diff(values)) > 0
    return values, frame


def projector(line):
    return np.outer(line, line.conj())


def spectator_test(h):
    values, frame = spectral(h)
    norm_h = (frame * values) @ frame.conj().T
    swap = np.eye(9).reshape(3, 3, 9).transpose(1, 0, 2).reshape(9, 9)
    anti = (np.eye(9) - swap) / 2
    sym = np.eye(9) - anti
    h0 = np.kron(norm_h, np.eye(9))
    h1 = np.kron(I, np.kron(norm_h, I))
    h2 = np.kron(np.eye(9), norm_h)
    pair = ((h0-h1) @ (h0-h2)).reshape(3, 9, 3, 9)

    def read(rho):
        return np.einsum("ts,isjt->ij", rho, pair)

    polynomial = 3*norm_h@norm_h - 2*sum(values)*norm_h + sum(
        values[i]*values[j] for i,j in ((0,1),(0,2),(1,2))) * I
    np.testing.assert_allclose(read(anti), polynomial, atol=4e-15)
    # Positive definite spectator state on wedge^2 with unequal weights and coherences.
    idx = np.arange(9) + 1.
    v = np.exp(1j*idx) * np.sqrt(idx)
    rho = anti @ (np.diag(idx) + projector(v)) @ anti
    rho /= np.trace(rho)
    restricted = frame.conj().T @ read(rho) @ frame
    np.testing.assert_allclose(restricted, np.diag(np.diag(restricted)), atol=2e-15)
    sign = [1, -1, 1]
    assert list(np.sign(restricted.diagonal().real)) == sign
    a, b = np.diff(values)
    critical = 2*a*b/(a*a + b*b + a*b)
    for scale, expected in ((.99, -1), (1.01, 1)):
        leak = scale*critical
        channel = (1-leak)*anti/3 + leak*sym/6
        assert np.sign(np.vdot(frame[:,1], read(channel) @ frame[:,1]).real) == expected
    return dict(normalized_eigenvalues=values, varied_ancilla_response=restricted.diagonal().real,
                equal_mixed_symmetric_leakage_flip=critical,
                polynomial_max_error=np.max(abs(polynomial-read(anti))))


def solve():
    bg = PARENT["OLD"]["mother"](CFG["matching_step"])
    y = bg["yW"][2,:,:3]
    charged = y @ y.conj().T
    _, le = spectral(charged)
    down = bg["yW"][1,:,:3] @ bg["right_quark_to_GJ"]
    # The five-bar identification acts on conjugate right-handed down flavours.
    down_h = (down.conj().T @ down).T
    _, ld = spectral(down_h)
    pd, pe = projector(ld[:,1]), projector(le[:,1])
    np.testing.assert_allclose(pd, projector(bg["tagged_light_columns"][:,1]), atol=5e-12)
    sd, se = I-2*pd, I-2*pe
    overlap = abs(np.vdot(le[:,1],ld[:,1]))
    angle = np.arccos(overlap)
    distance = np.linalg.norm(sd-se,2)
    np.testing.assert_allclose(distance,2*np.sin(angle),atol=2e-13)
    _, middle_axes = np.linalg.eigh(pd+pe)
    midpoint = I-2*projector(middle_axes[:,-1])
    radius=2*np.sin(angle/2)
    np.testing.assert_allclose([np.linalg.norm(midpoint-q,2) for q in (sd,se)],radius,atol=2e-13)
    _, positive, _, _=PARENT["terms"](bg,CFG["source_root_zero_based"],(1,1,1))
    _, win, _, _=PARENT["terms"](bg,CFG["source_root_zero_based"],(1,-1,1))
    np.testing.assert_allclose(se@positive@se,win,atol=5e-14)
    windows=bg["nu"]["FREQ"]
    rows={}
    for key, gate in zip(CFG["candidates"],(I,se,sd)):
        filt=gate@positive@gate
        np.testing.assert_allclose(np.linalg.eigvalsh(filt),[bg["cW"],bg["fW"],1.],atol=1e-10)
        common, row=PARENT["daughter"](bg,filt)
        nu=row["neutrinos"]
        limits=bg["model"]["JOINT"]["LEPTONS"]
        charged_pulls=(np.asarray(row["MZ_leptons"])-limits[:,0])/limits[:,1]
        clocks=[nu["dm31_eV2"],nu["dm32_eV2"]]
        clock_pulls=[]
        for x,ref in zip(clocks,("juno_NO_dm31","accelerator_NO_dm32")):
            center,low,high=windows[ref]
            clock_pulls.append((x-center)/(low if x<center else high))
        angle_pass=nu["angles_in_box"]
        rows[key]=dict(tau_common=common,charged_pulls=charged_pulls,angle_readings=nu["angle_readings"],
                       angles_in_boxes=angle_pass,solar_central_atmospheric=clocks,clock_pulls=clock_pulls,
                       FC68_FC90=nu["SK_T2K_FC68_90"],delta_deg=nu["delta_deg"],
                       joint_two_widths=all(abs(charged_pulls)<=2) and all(angle_pass) and
                       max(map(abs,clock_pulls))<=2 and bool(nu["SK_T2K_FC68_90"][1]))
    assert rows["charged_sign"]["joint_two_widths"]
    return dict(weak_pole_GeV=bg["Mweak"],quark=spectator_test(down_h),lepton=spectator_test(charged),
                strange_muon_absolute_overlap=overlap,principal_angle_deg=np.degrees(angle),
                sign_operator_difference=distance,optimal_shared_reflection_max_distance=radius,
                actual_ordered_mass_row_response=rows)


def native(value):
    if isinstance(value,dict):
        return {k:native(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)):
        return [native(v) for v in value]
    if isinstance(value,(bool,np.bool_)):
        return bool(value)
    if isinstance(value,(float,np.floating)):
        return 0. if abs(value)<1e-14 else float(f"{value:.10g}")
    return value


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=native(solve())
    dest=HERE/"results.json"
    if args.check:
        assert result==json.loads(dest.read_bytes())
        print("1750 spectator signs, leakage, axis transport and physical threshold tests verified")
    else:
        with dest.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
