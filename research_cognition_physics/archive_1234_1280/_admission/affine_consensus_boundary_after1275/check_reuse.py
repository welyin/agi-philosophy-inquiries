"""Check archived certificates and exact corollaries; no ODE rerun/new experiment."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
BASE = STAGE.parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen(number):
    folder = STAGE / str(number)
    for rel, expected in read(folder / "acceptance.json")["files_sha256"].items():
        assert sha(folder / rel) == expected, (number, rel)


def main():
    for number in (1250, 1272, 1274, 1275):
        frozen(number)
    old = read(STAGE / "1250/results.json")
    assert old["status"] == "PASS"
    branches = []
    for branch in old["group_same_history_null"]["branches"]:
        s = branch["screen"]
        lo, hi = map(F, s["corrected_balance"])
        assert lo <= 0 <= hi
        assert max(abs(lo), abs(hi)) < F("2.76e-19")
        gap = list(map(F, s["positive_memory_gap"]))
        margin = list(map(F, branch["read_contract"]["corrected_separation_margin"]))
        assert margin == [x - F(7, 10**9) for x in gap]
        assert margin[0] > F("1.52929e-7")
        branches.append({"alpha": branch["alpha"],
                         "old_corrected_residual": s["corrected_balance"],
                         "old_memory_gap": s["positive_memory_gap"],
                         "old_separation_margin": list(map(str, margin))})
    # 1175 eq.(14), with its original preparation and fixed ordinary inputs.
    eps, nu, omega, mass = F(3, 20), F(7, 20), F(1, 2), F(2)
    zdd = -8 * omega * eps * mass
    add = (eps - nu) * zdd
    theta_prime = 2 * add
    ricci_kk = -theta_prime
    assert (zdd, add, theta_prime, ricci_kk) == (F(-6, 5), F(6, 25), F(12, 25), F(-12, 25))
    sources = [STAGE / "1250/results.json", STAGE / "1250/proof.md",
               STAGE / "1272/source_proof.md", BASE / "archive_1161_1200/1175/proof.md",
               BASE / "archive_301_341/research_note_316.md"]
    out = {"status": "PASS", "kind": "archived_result_reuse_and_exact_corollary",
           "new_scientific_rounds": 0, "new_scientific_check_groups": 0,
           "new_simulation": False, "new_physical_instrument": False,
           "1250_uses_selected_EH_dynamics": True,
           "1250_archived_branches": branches,
           "1175_comparator_only": {"zdd0": str(zdd), "add0": str(add),
                                    "theta_prime0": str(theta_prime), "Ricci_kk0": str(ricci_kk)},
           "source_sha256": {str(p.relative_to(BASE)): sha(p) for p in sources},
           "goal_completed": False}
    (HERE / "reuse_results.json").write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "archived_branches": len(branches),
                      "new_simulation": False, "new_round": False}))


if __name__ == "__main__":
    main()
