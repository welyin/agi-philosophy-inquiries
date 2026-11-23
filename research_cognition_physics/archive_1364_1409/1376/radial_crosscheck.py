"""Retrospective common-source check using already frozen 1372 and 1369 values."""
import json
import math
from pathlib import Path

base = Path(__file__).resolve().parent
def read(path):
    return json.loads((base / path).read_text(encoding="utf-8-sig"))

radial = read("../1372/results.json")["bosonic_fourth_root_with_f_equal_v_no_neutral"]
cusp = read("results.json")["threshold_landmarks_GeV"]["cusp"][0]
masses = [radial[key]["Higgs_matching_boundary_GeV"]/14**0.25
          for key in ("box_low","center","box_high")]
min_trace = 7*math.log(cusp/max(masses))/(8*math.pi**2)
assert max(masses)<cusp and min_trace>1.5
print("single complex septet [low,center,high] (GeV):", masses)
print("lower bound on F at kink (step approximation):",min_trace)
