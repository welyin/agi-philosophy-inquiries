# 有限CT接续范围审计：独立审阅

2026-10-15。结论：**PASS_WITH_SCOPE**。本次审阅是非编号应用的核对，不新增科学轮次、正式检查组或采用条件；CT仍为待验证猜想。

## 1. 数学与共同对象

已逐段检查[范围稿](finite_ct_continuation_scope_after1245.md)、[代码](finite_ct_continuation_check_after1245.py)及[结果](finite_ct_continuation_results_after1245.json)，并独立手算以下内容。

- 从1245式(18)在原两类约束面消去总动量约束，所得向量密度$W^i$和分部积分符号正确。全部局部$H[f]$保持的必要条件是$\partial_i(N^2W^i)=0$；这不是完整Dirac算法的充分性结论。已读回正文§2.1最终限定：CMC识别沿1183的非退化近邻$\lambda\in[2/5,3/5]$，未向$\lambda=1/3$退化Legendre类别外推。
- 对$p=\sin x^1/10$，指针的空间梯度能与势能之和恒为$1/200$；再保$P_p=1$动能及三杆，得到$\mathcal H_f=131/200+3p/100$。直接收缩完整$\pi^{ij}$给$\mathcal H_g=-9-3p+3p^2/4$；记录器与钟项逐点补齐至$\mathcal H=0$。$\mathcal D_g=-dp,\mathcal D_m=dp$，两类初始约束均成立，没有删去材料反作用。
- $n_c$在$p\in[-1/10,1/10]$严格递增，端点恰为$16043/5400$、$17227/5400$；$F,M_r$及正作用量频率均为正。采用同一标签测度与其PB是必要的，不能将$n_c$作为演化中随意修改的控制。
- 固定$Q=ab$之前已给完数据，$R=0$使该初片对全部$Q>0$相同。$Q\ne1$时必要条件化为$\nabla\cdot(N^2\nabla p)=0$；紧环面上积分排除任意光滑严格正$N$。这只排除该准备的原精确接续，不能排除错误系数理论的全部初态。
- 两个积分独立核得$\langle|\nabla p|^2\rangle=1/200$与$\langle(1+\sin x^1/4)^2|\nabla p|^2\rangle=13/2560$。式(12)的测试函数$Np$在初片固定，故不增加测试函数导数项。$Q=1001/1000$给$\dot H[p]=1/200000$，符号一致。

## 2. 去重与适用范围

[1183](../../archive_1161_1200/1183/proof.md)已经处理任意小非零迹系数偏差、固定波数正源初片与正lapse障碍；本稿正确复用。新增窄诊断只是把$ab$偏差应用于1245的既有材料类别，不重报一般括号或CMC机制，因此维持非编号合适。

新空间依赖的指针准备、钟标签密度及$I=J=K=1$都不是1242原已认证准备。本文“同材料”指相同材料函数和正则类别；不表示继承原均匀实验菜单、实际制备或有限记录证书。$Q=1$只消除此初刻障碍；新初值的局部存在性未另证。原1242/1245均匀正例未被否定。

“资格条件在任意非零偏差下仍须满足”与“预测误差发生有限跳变”不同。本稿残差按$Q-1$连续变小，未给有限时间积分余项、实际读数分辨率或所有可能有效完成的反证。精确有效作用与一致父过程的受控截断已分开；不得据此要求任意截断先精确闭合，也不得把此候选的PDE修复变成研究纲领门槛。1235、1243、1245原适用结论均保留。

## 3. 实际复算与版本绑定

执行既有Python的默认只读模式（无`--write`）：

```powershell
& 'C:/Users/admin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B -X utf8 research_cognition_physics/archive_1234_/_admission/finite_ct_continuation_check_after1245.py
```

退出码0，输出`PASS_EXACT_INITIAL_CONSTRAINT_SCOPE`、`READ_ONLY_MATCH`，两个正lapse积分与上述手算一致；`new_numbered_rounds=0`、`new_formal_scientific_groups=0`。程序核的是初值多项式和精确积分，未积分ADM演化，也未重新数值验证一般HDA。

本审阅绑定下列SHA256；未修改任何被审文件、旧科学稿或导航：

| 文件 | SHA256 |
|---|---|
| `finite_ct_continuation_scope_after1245.md` | `8c46987d4c8498ce7e9e5139192189d0001837c97ca3bee82dcc81058ee64f4d` |
| `finite_ct_continuation_check_after1245.py` | `966e0789dcb243bb5e055bc44dba1a372d8f58934f30e526740491c97dbed43c` |
| `finite_ct_continuation_results_after1245.json` | `2d7e2eaae1be2f46c72ae97feef490ab50f098bbb0762e3cc514fe93c128a35a` |
