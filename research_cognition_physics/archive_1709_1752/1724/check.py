"""Recover a complex scalar flavor vertex and certify the aggregate-path cost."""
import argparse
from itertools import product
import json
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


CFG = read("inputs.json")


def pair(z):
    return [float(z.real), float(z.imag)]


def flavor():
    run = runpy.run_path(str(HERE / CFG["charged_run"]))
    _, rep = run["charged_background"](.1)
    y = run["unpack"](rep["down_texture_10PeV"])
    left, values, right_dag = np.linalg.svd(y)
    left, values, right_dag = left[:, ::-1], values[::-1], right_dag[::-1, :]
    p = np.outer(left[1, :].conj(), right_dag[:, 2].conj())
    np.testing.assert_allclose(left @ np.diag(values) @ right_dag, y, atol=1e-15)
    omega = np.exp(4j * np.pi * CFG["s_angle_pi"])
    z = y[1, 2]
    long = (z.imag / omega.imag) * omega
    return y, values, p, omega, long


def budget(z, omega):
    old = read(CFG["reference_cost"])
    example = old["two_relaxed_weak_examples"][0]
    tree = np.abs(old["original_named_Wilsons"][:5])
    edge_w = np.abs(example["edge_wilsons"])
    upper = max(max(edge_w[:5]), sum(edge_w[5:]))
    normal = min(abs((z * omega**i).imag) for i in range(3))
    def phase_cost(order):
        return normal * 2**(1+4*order) / (np.sqrt(3)/2)
    assert phase_cost(7) > upper
    # N=5 is the only possible remaining order. C,E are single real paths.
    allowed = [[n for n in range(0, 30, 3) if v*2**(4*n) <= upper*(1+1e-12)]
               for v in tree]
    pairs = [[c, e] for c, e in product(allowed[2], allowed[3]) if (c-e-1)%5 == 0]
    assert pairs == [[6, 0]]
    np.testing.assert_allclose(upper, tree[2]*2**24, rtol=1e-12)
    return {"total_complex_Y23":pair(z), "distance_to_nearest_CP_ray":normal,
            "N_at_least_7_edge23_l1_floor":phase_cost(7), "N5_basic_edges_n_at_or_below_cap":allowed,
            "N5_possible_C_E":pairs, "N5_upper_edge_l1":edge_w[:5].tolist()+[sum(edge_w[5:])],
            "whole_graph_sharp_minimax_edge_l1":upper}


def widths(values, p, long, c, d):
    g = c*np.diag(values) + d*long*p
    fc = abs(g)**2 + abs(g.T)**2
    diag = abs(np.diag(g))**2
    return np.array([fc[1, 2], diag[1], diag[2]]), diag, fc


def scalar(values, p, long, short, c, d):
    density = abs(p)**2 + abs(p.T)**2
    k = np.diag(p)
    A = np.column_stack([k.real, -k.imag])
    def recover(w):
        size = w[0]/(d*d*density[1, 2])
        rhs = (w[1:] - c*c*values[1:]**2 - d*d*size*abs(k[1:])**2)/(2*c*d*values[1:])
        return np.linalg.solve(A[1:], rhs)
    w, diag, fc = widths(values, p, long, c, d)
    point = recover(w)
    np.testing.assert_allclose(point, pair(long), atol=5e-15)
    reflected = np.conj(long*k[1])/k[1]
    wr, _, fcr = widths(values, p, reflected, c, d)
    np.testing.assert_allclose(wr[:2], w[:2], rtol=1e-13)
    np.testing.assert_allclose(np.triu(fc, 1), np.triu(fcr, 1), atol=1e-15)
    assert abs(wr[2]/w[2]-1)>1e-5
    singular = {}
    for i, j in ((0, 1), (0, 2), (1, 2)):
        sij = abs(np.linalg.det(A[[i,j]]))/(abs(k[i])*abs(k[j]))
        singular["dsb"[i]+"dsb"[j]] = sij
    # The inverse for signed Cartesian components is affine in three measured rates.
    extremes = [np.linalg.norm(recover(w*(1+np.array(sign)))-point)/abs(long)
                for sign in product((-1,1), repeat=3)]
    np.testing.assert_allclose(sum(k/values), 0, atol=3e-11)
    return {"y_d_s_b":values.tolist(), "complex_P_diagonal":[pair(v) for v in k],
            "P_sq_plus_reversed_for_ds_db_sb":[density[0,1],density[0,2],density[1,2]],
            "relative_two_diagonal_determinant":singular,
            "widths_sb_ss_bb":w.tolist(), "recovered_L_cartesian":point.tolist(),
            "free_checks_dd_db_ds":[diag[0],fc[0,2],fc[0,1]],
            "width_ratios_db_and_ds_to_sb":[fc[0,2]/fc[1,2],fc[0,1]/fc[1,2]],
            "reflected_L_holding_all_FC_and_ss":pair(reflected),
            "fractional_bb_change_at_reflection":wr[2]/w[2]-1,
            "fractional_bb_change_exchanging_mark_to_short":
                widths(values,p,short,c,d)[0][2]/w[2]-1,
            "worst_fractional_L_error_per_equal_fractional_width_error":max(extremes)}


def finite_width_check(values, p, long, c, d):
    """Spin-summed tree widths, including both conjugate unordered final states."""
    masses = CFG["finite_kinematics_m_b_over_scalar_mass"]*values/values[-1]
    unit = CFG["reconstruction_coordinate_unit"]
    sites = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    def rates(z):
        g = c*np.diag(values)+d*z*p
        rows = []
        for i,j in sites:
            a,b = masses[i],masses[j]
            phase_space = np.sqrt((1-(a+b)**2)*(1-(a-b)**2))
            amp = ((1-a*a-b*b)*(abs(g[i,j])**2+abs(g[j,i])**2)
                   -4*a*b*(g[i,j]*g[j,i]).real)
            rows.append(phase_space*amp/(1+(i==j)))
        return np.array(rows)
    baseline = rates(0)
    da = (rates(unit)-rates(-unit))/2
    db = (rates(1j*unit)-rates(-1j*unit))/2
    daa = (rates(unit)+rates(-unit))/2-baseline
    dbb = (rates(1j*unit)+rates(-1j*unit))/2-baseline
    dab = rates((1+1j)*unit)-baseline-da-db-daa-dbb
    data = rates(long)
    response = np.column_stack([da,db,daa,dbb,dab])/data[:,None]
    solution, _, rank, singular = np.linalg.lstsq(response,(data-baseline)/data,rcond=None)
    a,b = np.array(pair(long))/unit
    np.testing.assert_allclose(solution,[a,b,a*a,b*b,a*b],atol=2e-7)
    assert rank==5
    return {"mass_ratios_d_s_b":masses.tolist(),"rates_dd_ss_bb_ds_db_sb":data.tolist(),
            "row_normalized_singular_values":singular.tolist(),"lifted_system_rank":int(rank),
            "recovered_L":(solution[:2]*unit).tolist()}


def traces(y, long):
    old = read(CFG["named_paths"])
    patterns = {"original":[0,0,0,0,0,0,1]}
    patterns.update({f"Z{row['group']}":row["n"] for row in old["examples"]})
    rows = {}
    sites = [(0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2)]
    amps = [y[i,j] for i,j in sites[:5]]+[y[1,2]-long,long]
    for name, n in patterns.items():
        n = np.asarray(n)
        ds = np.zeros_like(y)
        dx = np.zeros_like(y)
        for site, v, ep, ex in zip(sites, amps, [0,2,0,1,2,1,1]+4*n, n%3):
            ds[site] += v*ep
            dx[site] += v*ex
        s_trace = np.trace(np.linalg.solve(y, ds))
        x_trace = np.trace(np.linalg.solve(y, dx))
        np.testing.assert_allclose(dx[1,2], long, atol=1e-15)
        assert np.count_nonzero(dx)==1
        np.testing.assert_allclose([s_trace.imag,abs(x_trace)],0,atol=1e-12)
        # Independently obtained from a finite radial displacement of a determinant.
        eps=1e-6
        def log_product(t):
            shifted=np.zeros_like(y)
            for site, v, ep in zip(sites, amps, [0,2,0,1,2,1,1]+4*n):
                shifted[site] += v*np.exp(ep*t)
            return np.linalg.slogdet(shifted)[1]
        fd=(log_product(eps)-log_product(-eps))/(2*eps)
        np.testing.assert_allclose(fd,s_trace.real,atol=2e-8)
        rows[name]={"S_trace":s_trace.real,"X_trace":x_trace.real,
                    "S_trace_over_colored_pair":s_trace.real/6}
    return rows


def main():
    y, values, p, omega, long = flavor()
    short = y[1,2]-long
    c,d = CFG["doublet_projection"], CFG["normalized_X_projection"]
    ratio = CFG["alternative_doublet_to_X_ratio"]
    alternate = scalar(values,p,long,short,ratio/np.hypot(1,ratio),1/np.hypot(1,ratio))
    return {"cost":budget(y[1,2],omega),"scalar":scalar(values,p,long,short,c,d),
            "X_dominated_example":alternate,
            "finite_mass_test":finite_width_check(values,p,long,c,d),
            "radial_tags":traces(y,long),"matching_scale_GeV":CFG["scale_GeV"]}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=main()
    if args.check:
        assert read("results.json")==results
        print("1724 aggregate cost, complex scalar reconstruction and independent radial tags verified")
    else:
        (HERE/"results.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf-8")
