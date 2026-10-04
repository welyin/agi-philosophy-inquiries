"""Register705 regional-compression entry, not a completed new round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round704_drafts/check_and_publish_entry.py').read_text('utf8')
s=re.sub(r'(?<!\d)(703|704)(?!\d)',lambda m:str(int(m.group())+1),s)
s=s.replace('own_gibbs_compression_entry','regional_compression_entry')
s=s.replace('own finite Gibbs initial-state corollary for original625 real histories',
            'regional compression algebra and inherited thermal tail corollary')
s=s.replace("['display_formulas']==4","['display_formulas']==6")
a=s.index('    entry=');b=s.index('    head,rest=',a)
s=s[:a]+'''    entry=(marker+' [全谱压缩与区域操作]( '+prefix+'round705_drafts/regional_compression_entry.md)'
        '精确压缩恒等式显示有限全谱截止未必保区域交换，625旧能源界控制固定热态中的误差。'
        '原读口双副本诊断已核，不将简单Kraus块冒充完整仪器；真正局部Gauss／尺度映射继续检验。'
        '正式科学轮次保持704，旧空间接口与目标不变。')
'''+s[b:]
s=s.replace('[全谱压缩与区域操作]( ','[全谱压缩与区域操作](')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
print('705 entry publication prepared')

