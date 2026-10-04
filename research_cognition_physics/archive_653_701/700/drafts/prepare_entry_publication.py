"""Register700 scope entry without claiming an extra science round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round699_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(698|699)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('free_domination_entry','mass_time_scope_entry')
s=s.replace('inherited endpoint comparison for the offdiagonal lower-bound task','mass and shared-time scope corollary of the certified negative physical kernel')
s=s.replace("['display_formulas']==1","['display_formulas']==2")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [质量路径与保来源替代的边界]('+prefix+'round700_drafts/mass_time_scope_entry.md)'
        '将680可逆质量合同接到699负见证：固定其余基线时，仅让质量随τ变化也不能修复坏τ的全代数RP。'
        '同参数处精确保物理反射型的正扩充亦被排除；共同连续缩放仍开放。'
        '原证书已复算，正式科学轮次保持699，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('700 entry publication prepared')
