"""Register701 executed joint-scale entry without claiming a science round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round700_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(699|700)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('mass_time_scope_entry','joint_scale_entry')
s=s.replace('mass and shared-time scope corollary of the certified negative physical kernel','finite-box zero-set and diagonal approximation scope audit')
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [共同尺度的辅助表示与实际物理极限]('+prefix+'round701_drafts/joint_scale_entry.md)'
        '多AP时间片的完整配置零集与热矩可接逐盒辅助收敛；对角选择不要求先加统一硬谱隙。'
        '原物理极限、非零归一和资源一致界仍未证明，不能用表示精度替代。'
        '实际AP核及有理误差合同已核，正式科学轮次保持700，目标不变。')
'''+s[b:]
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('701 entry publication prepared')
