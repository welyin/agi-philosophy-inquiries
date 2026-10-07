"""Publish mechanism screening 987, preserving the earlier overlap audit."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
old=(STAGE/'986/verify_round986.py').read_text('utf-8')
prefix=old.split('    result=read(HERE/"metric_dipole_bridge_results.json")',1)[0]
for a,b in [('range(776,986)','range(776,987)'),('Delivery checks for 986','Delivery checks for 987'),
            ('research_round_986_checks.json','research_round_987_checks.json')]:prefix=prefix.replace(a,b)
prefix=prefix.replace('    for path,keys in extra.items():',
    '    extra["987/drafts/common_overlap_audit_checks.json"]=("frozen_evidence_hashes","audit_file_hashes")\n'
    '    for path,keys in extra.items():')
checks=r'''    result=read(HERE/"coalition_mechanism_screen_results.json")
    assert result["round"]==987 and result["all_scientific_checks_passed"]
    for rel,digest in result["source_hashes"].items():assert sha(ROOT/rel)==digest,rel
    import importlib.util
    from fractions import Fraction as F
    from itertools import product
    spec=importlib.util.spec_from_file_location("core987",HERE/"coalition_mechanism_screen.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.compare(mod.run(),result)
    # Independent enumeration of 4^4 training-count outcomes, rather than 2^4
    # majority labels, with a recurrence for validation instead of a CDF sum.
    def predictions(table,mask):
        return [int(sum((table>>y)&1 for y in range(4) if y&mask==x&mask)*2>
                    sum(1 for y in range(4) if y&mask==x&mask)) for x in range(4)]
    def err(truth,pred,independent=False):
        if independent:return F(1,2)
        return F(1,10)+F(1,5)*sum(pred[x]!=((truth>>x)&1) for x in range(4))
    def pick(table):
        return min(range(4),key=lambda m:(err(table,predictions(table,m))+F(m.bit_count(),8),m.bit_count(),m))
    cdf_cache={}
    def cdf(p):
        if p in cdf_cache:return cdf_cache[p]
        dist=[F(1)]+[F(0)]*32
        for _ in range(32):
            dist=[dist[j]*(1-p)+(dist[j-1]*p if j else 0) for j in range(33)]
        cdf_cache[p]=sum(dist[:6]);return cdf_cache[p]
    independent_residuals=[]
    for independent in (False,True):
        truth=0 if independent else 6
        total=F(0);commit=F(0);loss=F(0);setup=F(9,2)
        for counts in product(range(4),repeat=4):
            weight=F(1)
            for x,k in enumerate(counts):
                p=F(1,2) if independent else F(9,10) if (truth>>x)&1 else F(1,10)
                weight*=math.comb(3,k)*p**k*(1-p)**(3-k)
            table=sum(int(k>=2)<<x for x,k in enumerate(counts))
            mask=pick(table);p=err(truth,predictions(table,mask),independent)
            fallback=err(truth,predictions(table,0),independent)
            accept=cdf(p) if mask else F(0)
            total+=weight;commit+=weight*accept
            loss+=weight*(accept*(p+F(mask.bit_count(),8))+(1-accept)*fallback)
            setup+=weight*4*(mask.bit_count()+1)*bool(mask)
        expected=result["finite_trial"]["independent_summary" if independent else "xor_summary"]
        assert total==1
        assert commit==F(expected["commit_probability"]["exact"])
        assert loss==F(expected["service_loss_with_port_fees"]["exact"])
        assert setup==F(expected["expected_setup_cost"]["exact"])
        independent_residuals.append(0)
    xor=result["finite_trial"]["xor_summary"]
    assert F(xor["net_using_max_setup"]["exact"])>F(727,50)
    assert F(result["finite_trial"]["independent_summary"]["commit_probability"]["exact"])<F(1,40000)
    assert len(result["finite_trial"]["all_functions"])==16
    assert all(F(row["normalization"]["exact"])==1 for row in result["finite_trial"]["all_functions"])
    witness=result["witness"]
    assert len(witness["retained_records"])==44 and witness["remaining_credits"]==768
    assert all(row["credits_before"]-row["credits_after"]==row["ports"].bit_count()+1
               for row in witness["retained_records"])
    note=STAGE/"research_note_987.md";prose=note.read_text("utf-8")
    assert prose.count("$$")==12
    assert re.findall(r"\\tag\{(\d+)\}",prose)==[str(i) for i in range(1,7)]
    for term in ("整体目标未完成","排除G作为普遍扩展规则","保留T作为具体候选",
        "不把XOR或贪心失效称为新定理","没有构造整套控制器的同一自主Hamiltonian",
        "此费用是决策权重","不声称净任务收益等于净物理收益"):
        assert term in prose,term
    newfiles=[note,HERE/"coalition_mechanism_screen.py",HERE/"coalition_mechanism_screen_results.json",
        Path(__file__),HERE/"drafts/mechanism_priority_entry.md",HERE/"drafts/mechanism_map_v0_3.md",
        HERE/"drafts/publish987.py",STAGE/"988/drafts/STATUS.md"]
    for p in newfiles:
        if p.suffix==".py":ast.parse(p.read_text("utf-8"))
    nav=[STAGE.parent/n for n in ("README.md","research_direction.md","RESEARCH_STATE.md")]+[
        STAGE/n for n in ("README.md","文件索引.md","阶段成果总览.md","跨阶段主题索引.md",
                         "_shared/notes/unified_physics_condition_ledger_current.md")]
    mapping=nav[-1].read_text("utf-8-sig").split(
        "## 六条共同协议：全局缺口对应与检验优先级（截至987）",1)[1].split("### 当前取舍",1)[0]
    assert re.findall(r"^\|(C\d\d) ",mapping,re.M)==[f"C{i:02d}" for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix==".md"]+nav:
        content=re.sub(r"\$\$.*?\$\$","",doc.read_text("utf-8-sig"),flags=re.S)
        for link in re.findall(r"\]\(([^)]+)\)",content):
            if re.match(r"^[a-zA-Z]+://",link) or link.startswith("#"):continue
            target=(doc.parent/link.split("#")[0].strip("<>")).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob("research_note_*.md"):
        match=re.fullmatch(r"research_note_(\d+).md",p.name)
        if match and p.parent.name.startswith("archive_"):nums.append(int(match[1]))
    assert sorted(n for n in nums if n<=987)==list(range(1,988))
    assert "001—987轮共987份" in nav[0].read_text("utf-8-sig")
    assert "231—987的757份" in nav[5].read_text("utf-8-sig")
    for p in nav:assert "987：关系更新的证据机制与有限联合试探" in p.read_text("utf-8-sig")
    oldfiles=[STAGE/"986/research_round_986_checks.json",HERE/"drafts/STATUS.md",
        HERE/"drafts/common_overlap_audit_checks.json"]
    prev=read(oldfiles[0])
    out=dict(round=987,date="2026-10-07",all_delivery_checks_passed=True,
        formal_reports=987,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_986=prev["cumulative_numbered_test_groups_from_985"]+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,independent_training_count_outcomes_checked=512,
        exact_probability_residuals=independent_residuals,
        mechanism_decision="reject universal strict edgewise gains; retain bounded validated coalition trial",
        complete_autonomous_physical_controller=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for key in ("frozen_inputs","new_scientific_and_entry_files"):assert before[key]==out[key],key
    return out

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run(args.write)
    if args.write:
        with TARGET.open("x",encoding="utf-8") as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write("\n")
    print(json.dumps({k:v for k,v in out.items() if k not in
         ("frozen_inputs","new_scientific_and_entry_files")},ensure_ascii=False,indent=2))
'''
target=STAGE/'987/verify_round987.py';assert not target.exists()
target.write_text(prefix+checks,encoding='utf-8')
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
 STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                  '_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
for p,raw in original.items():
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    heading='## 987：关系更新的证据机制与有限联合试探'
    previous='## 986后合并审计：记录交叠与来源边界'
    assert heading not in s and s.count(previous)==1
    block=(heading+'\n\n'
        +f'[987报告]({pre}research_note_987.md)回到961整体机制，排除逐关系即时净改善作为普遍扩展规则，'
        +'保留有预算的联合试探＋独立验证候选。有限两位任务的全部分支精确枚举：'
        +'XOR在最坏建立费下256次服务的期望净任务收益>14.54；无关资料误接入<1/40000，仍承担试探损失。'
        +f'[结果]({pre}987/coalition_mechanism_screen_results.json) · '
        +f'[核验]({pre}987/research_round_987_checks.json)。正式987／累计3772，整体目标未完成。\n\n'
        +f'[机制图v0.3]({pre}987/drafts/mechanism_map_v0_3.md)具体化M1与M4；费用不是热价，'
        +'政策未成为同一自然Hamiltonian。P981接合暂作验证资产；停止样本与模式优化，'
        +f'接[988]({pre}988/drafts/STATUS.md)先筛选证据如何使实际关系维持或撤销。应用目标及957假说保持。\n\n')
    s=s.replace(previous,block+previous,1)
    s=s.replace('001—986轮共986份','001—987轮共987份')
    if p==paths[4]:s=s.replace('当前正式986／累计3771，986已结项','当前正式987／累计3772，987已结项')
    if p==paths[5]:
        s=s.replace('# 231—986轮阶段成果总览','# 231—987轮阶段成果总览',1)
        s=s.replace('231—986的756份','231—987的757份',1)
    if p==paths[-1]:
        s=s.replace('## 六条共同协议：全局缺口对应与检验优先级（截至986）',
                    '## 六条共同协议：全局缺口对应与检验优先级（截至987）',1)
        rows={
        'C04':'|C04 同一内部动力学与控制|1、6、H2；929、935、947、961—968、987|987用固定有限政策依据实际训练及新验证选择端口，排除严格逐关系即时获益的普遍性|可交换任务、菜单、损失及准备输入；未构造整个控制器的自主自然H|',
        'C27':'|C27 共同数据与可区分预测|1、3、5、6；945、958—986、987|987分别枚举协同任务及独立对照，验证使无关关系误接入<1/40000；全部失败计账|经典机制筛选不推出几何或物理定律，端口费用不是完整能源或热账|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M);assert count==1
        replacement='4. 957假说v0.2保持；[987机制图v0.3](../../987/drafts/mechanism_map_v0_3.md)保留有限联合试探与独立验证，排除逐关系即时获益的普遍要求。下一项先回查证据驱动的实际关系保持／撤销；P981接合暂作后续验证资产，不自动接续物理接口和误差优化。'
        s,count=re.subn(r'^4\. 957假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M);assert count==1
    if b'\r\n' in raw:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Published mechanism round 987 to eight live navigation files.')
