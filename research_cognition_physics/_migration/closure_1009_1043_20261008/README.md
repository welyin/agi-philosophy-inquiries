# 1009—1043阶段关闭与目录迁移

2026-10-08。仅将实际目录`archive_1009_`改名为`archive_1009_1043`。新增科学组0，累计3819不变。全部科学Python、JSON及冻结正文原字节保留；仅阶段README和索引加入归档标题／复算说明，其改动以可逆字符补丁保存，不保存整份导航副本。

[逐文件映射与补丁](manifest.json)登记源、目标、前后哈希及原字节数；核验入口是项目根目录的`python -B -X utf8 scripts/close_research_1009_1043.py --verify`。首次执行核验后本目录的`checks.json`保存执行证据。

冻结的`research_layout.py`、`run_research_current.py`与其测试原字节保持。新入口复用这些旧工具，仅加入此关闭阶段的路径、runpy及已知Python子进程路由。不会生成旧目录壳、junction、ZIP或第二份研究目录。历史逻辑路径只存在于只读内存视图中。

## 复算

从项目根目录运行：

```powershell
python -B -X utf8 scripts/run_research_closed.py --verify-closure
python -B -X utf8 scripts/run_research_closed.py --script archive_1009_1043/1043/verify_round1043.py
python -B -X utf8 scripts/run_research_closed.py --script archive_1009_1043/1043/legacy_1009_scope_check.py
```

旧逻辑路径亦可作为`--script`值。直接执行历史脚本仍可能按旧目录查依赖，应使用上述入口。1009原默认全目录比较失败继续保留，范围复核通过不改称原入口通过。目录迁移不补充物理证明或改变阶段结论。

只有两份导航的反向补丁用于验证原字节；它们以后若再编辑，须记录新的维护层，不能绕过旧哈希检查。新增阶段资料应保存在后继目录中。
