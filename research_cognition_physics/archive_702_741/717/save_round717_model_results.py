"""Preserve the earlier force-only draft and publish the completed model result."""
import json
import os
from pathlib import Path
import joint_record_mass_feedback as model
HERE=Path(__file__).resolve().parent
before=model.TARGET.read_bytes()
old=json.loads(before)
assert old['tests_run']==3 and old['scope']['nonthermal_counterfamily_not_Gibbs']
result=model.run()
backup=HERE/'round717_drafts/force_results_before_bounded_record.json'
with backup.open('xb') as f:f.write(before)
temp=HERE/'joint_record_mass_feedback_results.new.json'
with temp.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
assert model.TARGET.read_bytes()==before
os.replace(temp,model.TARGET)
print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
