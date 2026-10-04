"""Register an executed699 screening entry without claiming a completed round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round698_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(697|698)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('polynomial_error_entry','free_domination_entry')
s=s.replace('original-source polynomial error contract','inherited endpoint comparison for the offdiagonal lower-bound task')
s=s.replace("import verify_interaction_rounds as core","import sys\nsys.path.insert(0,str(ARCHIVE))\nimport verify_interaction_rounds as core")
s=s.replace("['display_formulas']==8","['display_formulas']==1")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [原非对角下界的逐点比较筛查]('+prefix+'round699_drafts/free_domination_entry.md)'
        '复用696精确端点，排除以足够大常数乘自由核的逐点捷径；未重复实验。'
        '下一步须处理完整积分的矩或非常数下界。B10仍待认证，正式科学轮次保持698，目标不变。')
'''+s[b:]
s=s.replace("'free_domination_entry_first.py',",'')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('699 entry publication prepared')
