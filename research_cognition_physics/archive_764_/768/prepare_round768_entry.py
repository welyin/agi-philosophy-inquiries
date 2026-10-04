"""Prepare a read/derive-only working-entry publication; no new round count."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
source=(HERE/'publish_round766_entry.py').read_text('utf8')
mapping={'766':'768','765':'767','3482':'3488','3691':'3714',
         'physical_hadamard_entry':'joint_source_entry'}
pat=re.compile('|'.join(map(re.escape,sorted(mapping,key=len,reverse=True))))
text=pat.sub(lambda m:mapping[m[0]],source)
text=text.replace('from round768_drafts import smooth_positivity_calibration as calibration\n','')
text=text.replace('range(584, 768)','range(584, 768)')
a,b=text.index('assert calibration.run()'),text.index('report_text = ')
text=text[:a]+"""artifacts = ("round768_drafts/joint_source_entry.md",)
"""+text[b:]
text=text.replace('ast.parse((HERE/artifacts[1]).read_text("utf8"))',
                  'text_checks = core.text_checks(HERE/artifacts[0])\nassert text_checks["display_formulas"] == 3')
a,b=text.index('summary = '),text.index('planned = ')
text=text[:a]+"""summary = "**768完整来源入口已推进，正式仍767／3488：** [工作报告]({p}round768_drafts/joint_source_entry.md)回用734—735、756和764，核成熟BRST扩展的符号映射，并分清在壳自由态与离壳来源变分。下一项接同一完整作用的辅助部门及局部来源，未宣称量子Ward完成。[入口核验]({p}round768_drafts/joint_source_entry_checks.json)。不增加完成轮次，目标保持。"
"""+text[b:]
text=text.replace('**768短距离准备入口已推进','**768完整来源入口已推进')
text=text.replace('calibration_reproduced=True,',
                  'new_scientific_tests=0, entry_text_checks=text_checks,')
text=text.replace('Hadamard_state_existence_claimed=False,',
                  'inherited_767_free_Hadamard_state_preserved=True, all_sector_Ward_claimed=False,')
with (HERE/'publish_round768_entry.py').open('x',encoding='utf8',newline='\n') as f:
    f.write(text)
print('768 working-entry publisher prepared.')
