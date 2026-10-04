"""Register702 source scale entry, not a completed science round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round701_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(700|701)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('joint_scale_entry','source_scale_entry')
s=s.replace('finite-box zero-set and diagonal approximation scope audit',
            'normalized source scale corollary of699 negative certificate')
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [原负证书的观测尺度审计]('+prefix+'round702_drafts/source_scale_entry.md)'
        '699已认证数在明确的来源合同下给无量纲负裕量>0.0067；'
        '不能仅凭原始数值小而忽略，也未证明宏观来源中保留。'
        '整数证书已核，不重算积分；正式科学轮次保持701，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('702 entry publication prepared')
