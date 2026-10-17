# 静态共同窗口：中央角与预算的独立交叉审阅

2026-10-17。**PASS_CONDITIONAL_WINDOW_AND_SCOPE_ONLY**。非编号审阅，不新增科学组；CT未采用、原目标不变。

审阅对象为 `research_cognition_physics/archive_1234_/_admission/static_patch_window_after1270/` 下列七文件。独立交叉审查中央窗口和主端预算；SdS正文是本审阅者先前推导的发布版，不能另计一次独立证明。

## 核对结论

1. 中央角公式正确。直接积分固定迹给 $t_N=Z^{-1}\sum_{E_j<N}d_j(e^{-\beta E_j}-e^{-\beta N})$；由 $d_N=(t_N/p_N)C_Nd$ 得 $S_N=-\tau(C_Nd\log d)/p_N+\log p_N-\log t_N$。两角分解的 Bernoulli 相对熵项符号为负。中央参考的材料边缘、能级依赖观察者上端及其关联均已保留，没有继续假称独立乘积参考。
2. 未知态量词成立：在 $\operatorname{supp}f\subset[0,B]$、$B+E_{\max}\le N$ 下，新增总能投影在整个材料—观察者输入支持及任意旁观伙伴上恒等，不只是保持能量人口。此项不认证 dressing 或实际制备。$q\le B$ 的角一般非中央，与该总能角不能混用。$p_N<1$ 的反例及小失败率不自动控制熵差的限制准确。
3. 主端中央尾项可直接由 $Z\ge d e^{-\beta E_{\max}}$ 得到。本例 $\beta>6$、$N-E_{\max}=320$，故 $1-t_N<2^{-1920}$、$-\log t_N<2^{-1919}$。40阶正项和给出的 Gamma 截断界与面积有理界正确；三项和严格小于 $3/1000$。这是任意允许工作态的统一**熵差**上界，不是全态迹距或总实验误差。所用 $\beta=2\pi\sqrt{.99}$ 与本例几何一致，未沿用旧诊断的 $\log2$。
4. 中央代码已审读：精确分段和有理人口承担有限恒等式检查，特征值及对数浮点值仅作九例诊断。主端预算此前版本的熵、面积、中央尾及相位部分已实际默认只读重放：`5b5d89`，exit 0，`READ_ONLY_MATCH`；发布版上述科学字段保持不变，本次不重复运行。新增经典稳定有理字段已读，其物理来源由单独来源审阅负责；它不改变本次熵—面积预算，也不提供量子径向证书。
5. README范围可接受。独立面积来自选后EH的SdS根；其谱函数期望、完整壳源字典、背景参考流与壳固有时的同对象身份仍须区分。有限内作用角有中心，尚未变成原局域场论的II因子或JSS物理代数；实际量子壳来源、温度仪器、独立面积读口均未获证。经典表面作用及稳定性不能自动补这些身份。停止参数优化、先核具名同源量子过程的下一步没有把完整壳工程提升为全纲领门槛。

## 发布输入 SHA256

|文件|SHA256|
|---|---|
|`README.md`|`d4dd388b528de47d94f7daf33d6f33a74e56ca2d3722ac986341068082448d38`|
|`central_energy_window_after1270.md`|`b4f087f79788e8cac3ef87305a547b3ccdf2b9d0e4209dbdd842b8480fbd383c`|
|`central_energy_window_after1270.py`|`90b6013fda5a880a833c78c29b1fb89f14f6ed3bc6900a41c4d8576e55959c68`|
|`central_energy_window_after1270.json`|`8f087c2a0ed673cecd066f2ca1fbe21cba310f36fbbacab87933ab02768a106f`|
|`sds_geometry_scope_after1270.md`|`9c552dd1a642607bf416db4b5d5fbbe5f085e609726215bd2273bb7fc31ddc84`|
|`static_window_budget_after1270.py`|`31a356083e021867d95b60e21dd61230688ab735cee1d9093cbb6127c6d2d570`|
|`static_window_budget_after1270.json`|`045b72be7565c2949f19385613e2268184aeaf9a145ab4c0cb7147cf403e7e86`|

本次发布版只读审阅回执 `1930f7`、`788ec0`，均 exit 0。未修改正文、结果、导航或任何冻结科学文件。