"""Build and verify the frozen evidence inventory for construction baseline C0.

The old indexed claims retain their scope. Hash verification checks source
integrity, not the correctness of every theorem or historical test run.
Pandoc is only used to read the existing Markdown index structurally.
"""

import argparse
import copy
import hashlib
import json
import re
import unittest
from pathlib import Path

import pypandoc


ROOT = Path(__file__).resolve().parent.parent
OUTPUT = Path(__file__).with_name("inherited_results.json")
INDEX_OUTPUT = Path(__file__).with_name("INHERITED_INDEX.md")
LIMITS = {"research_information_geometry": 8, "research_cognition_physics": 216}
LIVE_NAVIGATION = {"README.md", "RESEARCH_STATE.md", "research_direction.md"}
GEOMETRY = (
    ("Static correlations do not uniquely determine interactions; dynamic responses can identify a declared model.",
     "Quantum states, interaction family and measurement access are inputs; no spacetime reconstruction."),
    ("Finite noisy dynamic observations support a scoped local identification and robustness analysis.",
     "Synthetic samples and calibrated model assumptions; not an experimental confirmation of cognitive physics."),
    ("A specified conditional-shift walk converges to a 1+1-dimensional free Dirac generator.",
     "Lattice, scale relation, unitary walk and continuum limit are declared; no gravity dynamics."),
    ("Variable propagation supports an effective geometric description with conformal ambiguity.",
     "The background is prescribed; massless data alone do not fix a unique physical metric."),
    ("A dynamical mean-field matter-background model has preparation-equivalence counterexamples.",
     "Energy conservation alone does not ensure statistical consistency or universal gravitational coupling."),
    ("A joint quantum background retains correlations missed by product mean-field models.",
     "A selected qubit-oscillator model, not quantization of all gravitational geometry."),
    ("A shared minimal metric requires compatible propagation speeds and normalized mass coefficients.",
     "Shared background variables do not themselves force universal metric coupling or field equations."),
    ("A leading soft spin-two Ward identity conditionally enforces universal coupling in connected scattering.",
     "Lorentz scattering, massless helicity two and soft-pole assumptions are inputs; emergence and nonlinear gravity remain open."),
)


def nodes(value):
    if isinstance(value, dict):
        if "t" in value:
            yield value
        for child in value.values():
            yield from nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from nodes(child)


def plain_text(value):
    if isinstance(value, list):
        return "".join(plain_text(child) for child in value)
    if not isinstance(value, dict):
        return ""
    kind = value.get("t")
    content = value.get("c")
    if kind == "Str":
        return content
    if kind in ("Space", "SoftBreak", "LineBreak"):
        return " "
    if kind in ("Code", "Math"):
        return content[1]
    if kind in ("Link", "Image"):
        return plain_text(content[1])
    return plain_text(content)


def local_links(value):
    return [node["c"][2][0] for node in nodes(value) if node.get("t") == "Link"
            and not re.match(r"^[a-zA-Z]+:", node["c"][2][0])]


def cognition_entries(document):
    rows = []
    for block in document["blocks"]:
        if block["t"] == "Table":
            table = block["c"]
            if len(table) != 6:
                raise ValueError("Expected the existing six-part Pandoc research index table.")
            rows.extend([cell[4] for cell in row[1]] for body in table[4] for row in body[2] + body[3])
        elif block["t"] == "LineBlock":
            for line in block["c"]:
                cells = [[]]
                for inline in line:
                    if inline.get("t") == "Str" and inline["c"] == "|":
                        cells.append([])
                    else:
                        cells[-1].append(inline)
                while cells and not plain_text(cells[-1]).strip():
                    cells.pop()
                rows.append(cells)
    entries = []
    for cells in rows:
        if not cells:
            continue
        note_links = [link for link in local_links(cells[0])
                      if re.fullmatch(r"research_note_\d+\.md", link)]
        if not note_links:
            continue
        if len(cells) != 4 or len(note_links) != 1:
            raise ValueError("A research row must retain four columns and one numbered note.")
        note = note_links[0]
        round_number = int(re.search(r"\d+", note).group())
        entries.append({
            "series": "research_cognition_physics", "round": round_number,
            "note": "research_cognition_physics/" + note,
            "indexed_conclusion": plain_text(cells[1]).strip(),
            "indexed_scope": plain_text(cells[2]).strip(),
            "indexed_artifacts": ["research_cognition_physics/" + link for link in local_links(cells[3])],
            "evidence_status": "inherited conditional result or audit; original assumptions retained",
        })
    expected = list(range(1, LIMITS["research_cognition_physics"] + 1))
    if [entry["round"] for entry in entries] != expected:
        raise ValueError("C0 requires each cognitive round 1..216 once and in order; review scope changes explicitly.")
    if any(not entry["indexed_conclusion"] or not entry["indexed_scope"] for entry in entries):
        raise ValueError("Every inherited result must retain both a conclusion and its scope.")
    return entries


def source_files(root):
    files = []
    for series in LIMITS:
        folder = root / series
        notes = sorted(int(match.group(1)) for path in folder.glob("research_note_*.md")
                       if (match := re.fullmatch(r"research_note_(\d+)\.md", path.name)))
        if notes != list(range(1, LIMITS[series] + 1)):
            raise ValueError("Numbered source coverage changed in " + series)
        files.extend(path for path in folder.iterdir() if path.is_file()
                     and path.suffix in {".md", ".py", ".json", ".txt"}
                     and path.name not in LIVE_NAVIGATION)
    return sorted(files)


def build_snapshot(root=ROOT):
    index = root / "research_cognition_physics" / "README.md"
    document = json.loads(pypandoc.convert_file(str(index), "json", format="markdown"))
    cognitive = cognition_entries(document)
    geometric = [{
        "series": "research_information_geometry", "round": index + 1,
        "note": "research_information_geometry/research_note_{:02d}.md".format(index + 1),
        "indexed_conclusion": conclusion, "indexed_scope": scope,
        "evidence_status": "inherited conditional physics benchmark, not a cognitive derivation",
    } for index, (conclusion, scope) in enumerate(GEOMETRY)]
    entries = geometric + cognitive
    files = source_files(root)
    artifacts = []
    for path in files:
        data = path.read_bytes()
        if path.suffix == ".json":
            json.loads(data.decode("utf-8-sig"))
        artifacts.append({"path": path.relative_to(root).as_posix(), "bytes": len(data),
                          "sha256": hashlib.sha256(data).hexdigest()})
    file_names = {item["path"] for item in artifacts}
    for entry in entries:
        note = root / entry["note"]
        entry["title"] = note.read_text(encoding="utf-8-sig").splitlines()[0].removeprefix("# ").strip()
        if entry["note"] not in file_names:
            raise ValueError("Missing note: " + entry["note"])
        for linked in entry.get("indexed_artifacts", []):
            if not (root / linked.split("#")[0]).is_file():
                raise ValueError("Missing indexed artifact: " + linked)
    return {
        "baseline": "C0", "cutoff_date": "2026-09-17", "source_round_limits": LIMITS,
        "coverage": "All 224 numbered source notes, indexed conclusions and scope, plus top-level text/code/JSON artifacts in both research libraries",
        "review_depth": "Full index consolidation and source integrity snapshot; not a new proof or rerun of every historical experiment",
        "mutable_navigation_excluded_from_hashes": sorted(LIVE_NAVIGATION),
        "old_initial_unnumbered_model": "research_cognition_physics/cognitive_loop.py",
        "all_historical_claims_adopted_as_axioms": False,
        "entries": entries, "artifacts": artifacts,
        "counts": {"numbered_notes": len(entries), "artifacts": len(artifacts),
                   "python_files": sum(path.suffix == ".py" for path in files),
                   "json_files": sum(path.suffix == ".json" for path in files)},
    }


def assert_same_snapshot(saved, current):
    if saved == current:
        return
    old_files = {item["path"]: item["sha256"] for item in saved.get("artifacts", [])}
    new_files = {item["path"]: item["sha256"] for item in current.get("artifacts", [])}
    changed = sorted(path for path in old_files.keys() | new_files.keys() if old_files.get(path) != new_files.get(path))
    raise ValueError("C0 inheritance changed; do not overwrite silently. Changed sources: " +
                     (", ".join(changed[:12]) if changed else "indexed metadata or scope"))


def render_index(snapshot):
    lines = ["# C0 逐轮继承索引", "",
             "由冻结清单生成，覆盖认知路线216轮及信息几何8轮。结论与边界一起继承，不代表全部结论已在本轮重证，也不把所有假设同时采纳。",
             "", "完整文件校验值见[inherited_results.json](inherited_results.json)；按主题整合见[FOUNDATIONS.md](FOUNDATIONS.md)。", ""]
    labels = {"research_information_geometry": "前序信息几何：8轮", "research_cognition_physics": "认知结构研究：216轮"}
    for series in LIMITS:
        lines.extend(["## " + labels[series], "", "| 记录 | 继承结论 | 适用边界 | 索引中的实现或证据 |",
                      "|:---|:---|:---|:---|"])
        for entry in snapshot["entries"]:
            if entry["series"] != series:
                continue
            conclusion = entry["indexed_conclusion"].replace("|", "\\|").replace("\n", " ")
            scope = entry["indexed_scope"].replace("|", "\\|").replace("\n", " ")
            artifacts = " / ".join("[" + Path(path).name + "](../" + path + ")"
                                   for path in entry.get("indexed_artifacts", [])) or "见原笔记及冻结清单"
            lines.append("| [" + str(entry["round"]) + "](../" + entry["note"] + ") | "
                         + conclusion + " | " + scope + " | " + artifacts + " |")
        lines.append("")
    lines.extend(["未编号的原经典闭环亦保留：[cognitive_loop.py](../research_cognition_physics/cognitive_loop.py)。",
                  "", "本清单不复制或修改旧证明；当前C0工作假设以[MODEL_CONTRACT.md](MODEL_CONTRACT.md)为准。", ""])
    return "\n".join(lines)


class InheritanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.snapshot = build_snapshot()

    def test_all_numbered_results_are_in_order_and_have_scope(self):
        self.assertEqual(self.snapshot["counts"]["numbered_notes"], 224)
        for series, limit in LIMITS.items():
            entries = [entry for entry in self.snapshot["entries"] if entry["series"] == series]
            self.assertEqual([entry["round"] for entry in entries], list(range(1, limit + 1)))
            self.assertTrue(all(entry["indexed_conclusion"] and entry["indexed_scope"] for entry in entries))

    def test_every_numbered_note_and_original_code_is_hashed(self):
        names = {item["path"] for item in self.snapshot["artifacts"]}
        self.assertTrue(all(entry["note"] in names for entry in self.snapshot["entries"]))
        self.assertIn(self.snapshot["old_initial_unnumbered_model"], names)
        self.assertTrue(all(len(item["sha256"]) == 64 for item in self.snapshot["artifacts"]))

    def test_metadata_and_scope_changes_are_detected(self):
        changed = copy.deepcopy(self.snapshot)
        changed["entries"][0]["indexed_scope"] = "Unconditional claim"
        with self.assertRaises(ValueError):
            assert_same_snapshot(self.snapshot, changed)

    def test_source_changes_are_detected(self):
        changed = copy.deepcopy(self.snapshot)
        changed["artifacts"][0]["sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            assert_same_snapshot(self.snapshot, changed)
        assert_same_snapshot(self.snapshot, copy.deepcopy(self.snapshot))

    def test_readable_index_links_every_numbered_note(self):
        rendered = render_index(self.snapshot)
        for entry in self.snapshot["entries"]:
            self.assertIn("](../" + entry["note"] + ")", rendered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write-snapshot", action="store_true")
    action.add_argument("--verify", action="store_true")
    action.add_argument("--write-index", action="store_true")
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(InheritanceTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    snapshot = InheritanceTests.snapshot
    if args.write_snapshot:
        if OUTPUT.exists():
            assert_same_snapshot(json.loads(OUTPUT.read_text(encoding="utf-8")), snapshot)
        else:
            OUTPUT.write_text(json.dumps(snapshot, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if args.verify:
        assert_same_snapshot(json.loads(OUTPUT.read_text(encoding="utf-8")), snapshot)
        if INDEX_OUTPUT.exists() and INDEX_OUTPUT.read_text(encoding="utf-8") != render_index(snapshot):
            raise ValueError("The readable inheritance index differs from the frozen inventory.")
    if args.write_index:
        assert_same_snapshot(json.loads(OUTPUT.read_text(encoding="utf-8")), snapshot)
        INDEX_OUTPUT.write_text(render_index(snapshot), encoding="utf-8")
    print(json.dumps(snapshot["counts"], indent=2))


if __name__ == "__main__":
    main()