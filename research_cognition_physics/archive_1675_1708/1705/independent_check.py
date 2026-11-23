"""Check the published SK+T2K 2-D confidence polygons against the protected physical column."""
import argparse
import hashlib
import json
import math
import struct
import zipfile
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TARGET = ROOT / "independent_results.json"
FILE_SHA = "25f2852ee6629a2eeb7a9dc8b9c43165dbd2feb801944f10031104681dc9680c"
ZIP_SHA = "3dcdc2f9234c923f80339a2136f78eb96de646c59af79d2e609a22e8b43ac2e8"
PREFIX = "gr2D_th23_dCP_wRC_FC_globalMO_NO_"
GRAPH_KEYS = ["bestfit", "conf68_0", "conf90_0", "conf90_1"]
# Hash of the unrounded big-endian fX and fY vectors, in GRAPH_KEYS order.
VERTEX_SHA = "59eb5327371a7482ed5ef9af307a30f2bd97bc27e756bddf84d26175a7af071e"


def released_vertices(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ZIP_SHA
    with zipfile.ZipFile(path) as archive:
        data = archive.read("FrequentistAnalysis1.root")
    assert hashlib.sha256(data).hexdigest() == FILE_SHA
    graphs = {}
    for name in GRAPH_KEYS:
        marker = (PREFIX + name).encode("ascii")
        key = data.index(b"\x06TGraph" + bytes([len(marker)]) + marker) - 26
        size, version, uncompressed, _, key_bytes, cycle, address, _ = struct.unpack_from(">ihiIhhii", data, key)
        assert (version, cycle, address) == (4, 1, key)
        payload = data[key + key_bytes:key + size]
        assert payload[:3] == b"ZL\x08" and int.from_bytes(payload[3:6], "little") == len(payload) - 9
        graph = zlib.decompress(payload[9:])
        assert len(graph) == uncompressed == int.from_bytes(payload[6:9], "little")
        assert struct.unpack_from(">H", graph, 4)[0] == 4
        # ROOT TGraph-v4 streams four byte-counted bases, then fNpoints, fX and fY.
        offset = 6
        for _ in range(4):
            count = struct.unpack_from(">I", graph, offset)[0]
            assert count >> 30 == 1
            offset += 4 + (count & 0x3fffffff)
        n = struct.unpack_from(">I", graph, offset)[0]
        offset += 4
        assert graph[offset] == graph[offset + 1 + 8*n] == 1
        x = struct.unpack_from(">" + str(n) + "d", graph, offset + 1)
        y = struct.unpack_from(">" + str(n) + "d", graph, offset + 2 + 8*n)
        graphs[name] = [list(pair) for pair in zip(x, y)]
    return graphs


def vertex_hash(graphs):
    digest = hashlib.sha256()
    for name in GRAPH_KEYS:
        points = graphs[name]
        for values in zip(*points):
            digest.update(struct.pack(">" + str(len(values)) + "d", *values))
    return digest.hexdigest()


def edges(points):
    return zip(points, points[1:] + points[:1])


def vertical_crossings(points, y):
    values = [d0 + (d1-d0)*(y-u0)/(u1-u0)
              for (u0, d0), (u1, d1) in edges(points)
              if min(u0, u1) <= y < max(u0, u1)]
    values.sort()
    assert len(values) % 2 == 0
    return values


def interior(points, y, delta):
    crossings = vertical_crossings(points, y)
    return any(a < delta < b for a, b in zip(crossings[::2], crossings[1::2]))


def slab_top(points, left, right):
    # All vertices of the polygon clipped by vertical lines are original vertices or intersections.
    values = [d for u, d in points if left <= u <= right]
    for (u0, d0), (u1, d1) in edges(points):
        if u0 == u1:
            continue
        for y in (left, right):
            if min(u0, u1) <= y <= max(u0, u1):
                values.append(d0 + (d1-d0)*(y-u0)/(u1-u0))
    return max(values) if values else None


def reactor_root(a):
    return (1-math.sqrt(1-a))/2


def delta_from_column(y, z, p):
    x = 1-p[0]/(1-z)
    denominator = 2*math.sqrt(x*(1-x)*z*y*(1-y))
    cosine = (p[1]-x*(1-y)-(1-x)*y*z)/denominator
    return x, -math.acos(cosine)


def allowed_upper_y_cells(points, boxes, x_range, z_range, y_range):
    e, m, _ = boxes
    x0 = max(x_range[0], 1-e[1]/(1-z_range[1]))
    x1 = min(x_range[1], 1-e[0]/(1-z_range[0]))
    z = z_range[1]
    allowed = []
    grid = 1400
    for i in range(grid):
        a = y_range[0] + (y_range[1]-y_range[0])*i/grid
        b = y_range[0] + (y_range[1]-y_range[0])*(i+1)/grid
        delta_top = slab_top(points, a, b)
        if delta_top is None:
            continue
        # cos is increasing on [-pi,0]; positive-part also bounds negative interference.
        cos_max = max(0, math.cos(delta_top))
        v = min(b, max(.5, a))
        bound = x1*(1-a) + (1-x0)*b*z + 2*cos_max*math.sqrt(x1*(1-x1)*z*v*(1-v))
        if m[0] <= bound + 1e-10:
            allowed.append(b)
    return max(allowed)


def calculate(graphs):
    assert vertex_hash(graphs) == VERTEX_SHA
    assert [len(graphs[k]) for k in GRAPH_KEYS] == [1, 59, 60, 14]
    assert graphs["conf68_0"][0] == graphs["conf68_0"][-1]
    assert all(abs(p[1]+math.pi)<1e-14 for p in (graphs["conf90_0"][0],graphs["conf90_0"][-1]))
    assert all(abs(p[1]-math.pi)<1e-14 for p in (graphs["conf90_1"][0],graphs["conf90_1"][-1]))
    bestfit = graphs["bestfit"][0]
    assert all(interior(graphs[k], *bestfit) for k in ("conf68_0", "conf90_0"))
    sources = json.loads((ROOT/"../1703/results.json").read_bytes())
    active = {k: row for k, row in sources["cases"].items() if row["solar_reactor_compatible"]}
    assert len(active) == 6
    x = [.2908,.3164]
    yrange = [.49,.63]
    z = [reactor_root(v) for v in (.0803,.0899)]
    low_cos = min(row["necessary_cos_delta_lower"] for row in active.values())-1e-8
    phase = math.acos(low_cos)
    phase_caps = {k:[min(a[1] for a in graphs[k]),max(a[1] for a in graphs[k])]
                  for k in GRAPH_KEYS[1:]}
    assert phase < phase_caps["conf90_1"][0]
    j_prefactor_floor = (1-z[0])*math.sqrt(z[0]*x[0]*(1-x[0])*yrange[1]*(1-yrange[1]))
    result = {
        "data_url":"https://zenodo.org/records/12702685", "official_member":"FrequentistAnalysis1.root",
        "zip_sha256": ZIP_SHA, "root_sha256": FILE_SHA, "coordinate_sha256": VERTEX_SHA,
        "official_object_prefix": PREFIX, "official_coordinates_rad":graphs,
        "printed_fc_1d_profiled_MO_683_percent_delta_rad":[-2.71,-1.03],
        "printed_fc_1d_profiled_MO_683_percent_y":[.443,.574],
        "model_abs_phase_max_rad":phase, "printed_1d_gap_rad":1.03-phase,
        "fc_NO_global_order_contours_phase_extrema_rad":phase_caps,
        "positive_branch_in_90_percent_2d":False,
        "minimum_abs_J_if_within_published_2d":{
            "conf68":j_prefactor_floor*math.sin(-phase_caps["conf68_0"][1]),
            "conf90":j_prefactor_floor*math.sin(-phase_caps["conf90_0"][1])},
        "paper_CPC_p_value_and_robustness_ceiling":[.037,.050],
        "cases":{}
    }
    for name, row in active.items():
        p = row["p_MZ"]
        y = .493
        wz = reactor_root(.0890)
        wx, delta = delta_from_column(y, wz, p)
        assert x[0] < wx < x[1]
        assert all(interior(graphs[k], y, delta) for k in ("conf68_0","conf90_0"))
        assert not any(interior(graphs[k], y, -delta) for k in GRAPH_KEYS[1:])
        m = wx*(1-y)+(1-wx)*y*wz+2*math.cos(delta)*math.sqrt(wx*(1-wx)*wz*y*(1-y))
        rec = [(1-wx)*(1-wz),m,1-(1-wx)*(1-wz)-m]
        assert max(abs(a-b) for a,b in zip(rec,p)) < 1.6e-9
        caps = {k:allowed_upper_y_cells(graphs[k],row["protected_probability_envelopes"],x,z,yrange)
                for k in ["conf68_0","conf90_0"]}
        result["cases"][name] = {"sin2_theta12":wx,"sin2_theta23":y,"sin2_2theta13":.0890,
            "delta_rad":delta,"column":rec,"official_68_and_90_2d_contains_witness":True,
            "official_68_vertical_slice_at_witness":vertical_crossings(graphs["conf68_0"],y),
            "official_90_vertical_slice_at_witness":vertical_crossings(graphs["conf90_0"],y),
            "certified_y_upper_over_whole_probability_envelopes":caps}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--release",type=Path,help="optional original Zenodo ZIP for independent extraction")
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    reference = json.loads(TARGET.read_bytes()) if TARGET.exists() else None
    graphs = released_vertices(args.release) if args.release else reference["official_coordinates_rad"]
    data = calculate(graphs)
    if args.check:
        assert data == reference
        print("1705 independent 2D contours and protected-column checks passed")
    else:
        TARGET.write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n",encoding="utf-8")
