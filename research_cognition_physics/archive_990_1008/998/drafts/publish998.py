"""Publish the finite receiver mechanism without rewriting frozen research."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
TITLE='## 998：组织形成与有限接收者'
MARK='## 998工作审计：共同边界与候选投入'


def main():
    paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals={p:p.read_bytes() for p in paths}
    updates={}
    for p,raw in originals.items():
        s=raw.decode('utf-8-sig').replace('\r\n','\n')
        assert TITLE not in s and s.count(MARK)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(TITLE+'\n\n'
            +f'[998报告]({pre}research_note_998.md)接成熟引力热力学：明确有效支上，'
            +'有限释能、有限辐射接收者及单一热交换模稳定可以相容；接收者越大不必越稳定。'
            +f'[结果]({pre}998/finite_receiver_formation_results.json) · '
            +f'[核验]({pre}998/research_round_998_checks.json)。正式998／累计3782，整体未完成。\n\n'
            +f'[整体机制补充]({pre}998/formation_resource_adoption_v1.md)区分载体形成、信息保持和运行资源。'
            +'本例未从等离子体生成材料、未提取功或制造记录，也未推导引力。停止云模型与器件优化；'
            +f'接[999]({pre}999/drafts/STATUS.md)审这三类条件在共同环境中的相容机制。'
            +'应用目标、旧空间结论和957保持。\n\n')
        s=s.replace(MARK,block+MARK,1)
        if p==paths[0]:
            assert s.count('001—997轮共997份')==1
            s=s.replace('001—997轮共997份','001—998轮共998份',1)
        if p==paths[4]:
            before='当前正式997／累计3781，997已结项'
            assert s.count(before)==1
            s=s.replace(before,'当前正式998／累计3782，998已结项',1)
        if p==paths[5]:
            assert s.count('231—997的767份')==1
            s=s.replace('# 231—997轮阶段成果总览','# 231—998轮阶段成果总览',1)
            s=s.replace('231—997的767份','231—998的768份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至997）'
            assert s.count(before)==1
            s=s.replace(before,before.replace('997','998'),1)
            rows={
                'C19':('|C19 真空、热态、非平衡准备|1、3、6、H2；592、625、955—980、994—998|'
                       '旧共同准备条件保留；998把有限接收者的能源、熵和热响应纳入结构变化|'
                       '宏观状态方程、弱接触及准平衡支是输入；未共同实现实际热史、形成与记录|'),
                'C23':('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961、964—980、992—993、998|'
                       '992资源失配和993条件误差窗复用；998有限接收者反馈给条件性热稳定|'
                       '粗粒热熵不等于整体微观熵；载体／信息／运行窗口未共同实现，不自动产生宇宙箭头|'),
                'C26':('|C26 宏观经典、流体等|1、3、4、6；350、354、758—761、973—977、998|'
                       '977共同经典来源箱与量子内层保留；998在明示经典有效支筛选形成／散热机制|'
                       '热交换模稳定不是完整机械稳定、真实结构形成或实际记录；未签收全部流体恢复|')}
            for cid,row in rows.items():
                s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M)
                assert count==1,cid
            anchor='### 当前取舍\n'
            assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式998采用局部组织与整体内部有限接收者共同变化的机制；'
                '当前正式998／累计3782。载体、信息、运行三个条件分别保留，下一项审共同环境窗口；'
                '不继续细化云模型或热设备，C01—C27未全部关闭。以下工作审计保留其当时范围。\n',1)
            replacement=('4. 957数学假说v0.2保持；[998机制补充](../../998/formation_resource_adoption_v1.md)'
                '将结构变化、能源接收与热反馈相接；[B997](../../997/cosmological_adoption_v1.md)'
                '及宇宙边界按原范围保留。停止单项候选优化，下一项复用旧记录及热合同，判断载体、'
                '信息和运行资源在共同环境中的相容机制，不先造完整自主供能存储器。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    for p,data in updates.items():p.write_bytes(data)
    print('Published 998 to eight live navigation files; historical reports preserved.')


if __name__=='__main__':main()
