"""Author second-algorithm checks; this is not independent-agent review."""
from pathlib import Path
from fractions import Fraction
from urllib.parse import unquote
import argparse
import ast
import hashlib
import json
import re
import numpy as np
import fusion_pairing_mass as science

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
NOTE = STAGE / "research_note_1039.md"
OUT = HERE / "research_round_1039_checks.json"
OWN = [
    "fusion_pairing_mass.py", "fusion_pairing_mass_results.json",
    "drafts/fusion_pairing_mass_derivation.md", "selection_audit.md",
    "input_dependency_update_v0_28.md", "review.md", "NEXT.md",
    "verify_round1039.py",
]
INPUTS = [
    BASE / "archive_223_230/research_note_227.md",
    BASE / "archive_531_553/research_note_531.md",
    BASE / "archive_923_934/research_note_934.md",
    BASE / "archive_923_934/934/drafts/operation_hypotheses_screening.md",
    BASE / "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    STAGE / "1009/input_dependency_ledger_v0_1.md",
    STAGE / "research_note_1013.md",
    STAGE / "research_note_1019.md",
    STAGE / "research_note_1032.md",
    STAGE / "research_note_1035.md",
    STAGE / "1035/input_dependency_update_v0_24.md",
    STAGE / "1035/bridge_ledger_v0_2.json",
    STAGE / "research_note_1036.md",
    STAGE / "research_note_1038.md",
    STAGE / "1038/parallel_integration_1036_1038.md",
    STAGE / "1038/NEXT.md",
]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def exact_rank(rows):
    """Sparse Gaussian elimination over Q; no floating rank threshold."""
    pivots = {}
    for values in rows:
        row = {j: Fraction(int(x)) for j, x in enumerate(values) if x}
        while row:
            lead = min(row)
            value = row[lead]
            if lead not in pivots:
                pivots[lead] = {j: x/value for j, x in row.items()}
                break
            for j, x in pivots[lead].items():
                new = row.get(j, Fraction(0))-value*x
                if new:
                    row[j] = new
                else:
                    row.pop(j, None)
    return len(pivots)


def independent_algorithm():
    saved = json.loads(science.OUT.read_text("utf8"))
    # Square map from the saved abstract multiplication table, not U @ U.
    indicators = {}
    for name in ["D8", "Q8"]:
        data = saved["groups"][name]
        table = data["multiplication_table"]
        char = data["character_table_on_elements"][4]
        fs = sum(Fraction(char[table[g][g]], 8) for g in range(8))
        assert fs == (1 if name == "D8" else -1)
        indicators[name] = str(fs)
    # Direct constraints on all symmetric matrices, with rational exact rank.
    gen = {
        "D8": [np.array([[0, -1], [1, 0]], complex),
               np.array([[1, 0], [0, -1]], complex)],
        "Q8": [np.array([[1j, 0], [0, -1j]], complex),
               np.array([[0, 1j], [1j, 0]], complex)],
    }
    ranks = []
    for name in ["D8", "Q8"]:
        sign = 1 if name == "D8" else -1
        for n in range(1, 7):
            size = 2*n
            columns = []
            for i in range(size):
                for j in range(i, size):
                    e = np.zeros((size, size), complex)
                    e[i, j] = e[j, i] = 1
                    pieces = []
                    for u0 in gen[name]:
                        u = np.kron(u0, np.eye(n))
                        pieces.extend((u.T @ e @ u-e).ravel())
                    columns.append(pieces)
            matrix = np.array(columns).T
            assert np.all(matrix.imag == 0)
            integer = matrix.real.astype(int)
            assert np.array_equal(integer, matrix)
            rank = exact_rank(integer.tolist())
            dimension = len(columns)-rank
            assert dimension == n*(n+sign)//2
            ranks.append(dict(group=name, n=n, unknowns=len(columns),
                              exact_constraint_rank=rank,
                              allowed_complex_dimension=dimension))
    # Closed entries of trivial projectors, avoiding group averaging.
    pd = np.array([[1, 0, 0, 1], [0, 0, 0, 0],
                   [0, 0, 0, 0], [1, 0, 0, 1]], float)/2
    pq = np.array([[0, 0, 0, 0], [0, 1, -1, 0],
                   [0, -1, 1, 0], [0, 0, 0, 0]], float)/2
    for name, projector in [("D8", pd), ("Q8", pq)]:
        data = saved["groups"][name]["trivial_pair_projector"]
        assert np.array_equal(projector, np.array(data["real"]))
        assert np.count_nonzero(data["imag"]) == 0
    assert np.trace(pd @ pq) == 0
    rational_cases = [(Fraction(1, 10), Fraction(1, 20)),
                      (Fraction(1, 5), Fraction(1, 10)),
                      (Fraction(49, 100), Fraction(1, 100)),
                      (Fraction(1, 2), Fraction(0))]
    margins = [str(1-2*eta-eps) for eta, eps in rational_cases]
    assert margins == ["3/4", "1/2", "1/100", "0"]
    return dict(exact_square_map_indicators=indicators,
                rational_constraint_ranks=ranks,
                explicit_projectors_checked=True,
                rational_readout_margins=margins)


def link_check(path):
    fence = chr(96)*3
    s = re.sub(fence+r".*?"+fence+r"|\$\$.*?\$\$", "",
               path.read_text("utf8"), flags=re.S)
    count = 0
    for target in re.findall(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)", s):
        target = unquote(target.split("#", 1)[0])
        if not target or re.match(r"^[A-Za-z]+:", target):
            continue
        assert (path.parent / target).resolve().exists(), (path, target)
        count += 1
    return count


def verify():
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text("utf8")))
    second = independent_algorithm()
    for name in ["fusion_pairing_mass.py", "verify_round1039.py"]:
        ast.parse((HERE/name).read_text("utf8"))
    docs = [NOTE]+[HERE/name for name in OWN if name.endswith(".md")]
    links = sum(link_check(p) for p in docs)
    for p in docs:
        s = p.read_text("utf8")
        assert "\ufffd" not in s
        assert not any(ord(c) < 32 and c not in "\n\r\t" for c in s)
    artifacts = [NOTE]+[HERE/name for name in OWN]
    return dict(
        round=1039, passed=True, scientific_groups_proposed=1,
        global_cumulative_count=None, new_cognitive_axioms=0,
        independent_agent_review="pending main-line assignment",
        author_second_algorithm_is_not_independent_review=True,
        exact_equality_checks=fresh["exact_equality_checks"],
        mass_space_calibrations=16, approximate_mass_calibrations=20,
        independent_algorithm=second, local_links_checked=links,
        artifact_sha256={str(p.relative_to(STAGE)): sha(p) for p in artifacts},
        input_sha256={str(p.relative_to(BASE)): sha(p) for p in INPUTS},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write:
        with OUT.open("x", encoding="utf8") as f:
            f.write("{}\n")
    result = verify()
    if args.write:
        OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n",
                       encoding="utf8")
    elif OUT.exists():
        science.compare(result, json.loads(OUT.read_text("utf8")))
    print(json.dumps(dict(round=1039, passed=True,
                          exact_mass_space_ranks=12,
                          links=result["local_links_checked"],
                          artifact_hashes=len(result["artifact_sha256"])),
                     ensure_ascii=False))
