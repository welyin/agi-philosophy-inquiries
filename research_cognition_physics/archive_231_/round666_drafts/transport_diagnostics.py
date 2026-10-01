"""Capture first666 transport failure before selecting any repair."""
import hashlib
import json
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_reference_mass_stabilizer as model

try:
    model.transported_auxiliary_check()
except AssertionError as error:
    tb=error.__traceback__
    while tb.tb_next is not None:tb=tb.tb_next
    local=tb.tb_frame.f_locals
    result={k:local[k] for k in ('source_error','missing_error','mean_errors') if k in local}
    result['code_sha256']=hashlib.sha256((BASE/'joint_reference_mass_stabilizer.py').read_bytes()).hexdigest()
    result['failure']='First formal666 probe failed the combined source/quadrature assertion; no result published.'
    with (HERE/'transport_first_failure.json').open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    with (HERE/'joint_reference_mass_stabilizer_first_attempt.py').open('xb') as f:
        f.write((BASE/'joint_reference_mass_stabilizer.py').read_bytes())
    print(json.dumps(result))
