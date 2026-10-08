# 1042独立数学与原文范围审阅

2026-10-08。审阅者：`/root/next_structural_selection`。本轮报告、科学代码及结果均非本审阅者编写；只新增本文件，不修改作者冻结资产，不新增科学校准组。

**结论：通过，未发现需要修订的数学或范围问题。** 已读冻结正文、代码及结果；执行不导入作者函数的独立核算，并在作者首冻后通过默认只读核验及另写的SHA核对。

## 1. 原文与历史范围

[Harlow定理1.1、§5.1—5.3](https://arxiv.org/html/1607.03901)的对象为固定有限码、物理二分及逻辑代数；双方的精确恢复、相应全态熵式及态对相对熵关系在此合同内等价。§6.4明确不由单个分区推出任意区域的恢复或面积极值解释。1042保持了这些限制。

391的实际通道代数算法和其AB|CD互补特例、839§7的相对熵运输没有被重算为新发现。新报告在指定编码内独立核验互补性，而未将“共同中心存在”当成它的证明。

## 2. 未归一化块及两端同一中心项

逻辑限制态是直接和矩阵\(\bigoplus_\alpha w_\alpha\rho_{\alpha,a}\)，其迹为1；每个块并不分别归一。因而中心Shannon项及模算符中的\(-\log w_\alpha\)必须存在。正文式(3)—(6)和代码`logical_reduction`、`algebra_modular`均正确保留了它。表示中的b重数不应被另加为代数熵。

对固定p，区域块是\(w_\alpha\rho_{\alpha,a}\otimes\tau_\alpha\)。其对数分解后，χ在另一端的纯化满足

$$
\langle\chi_\alpha|[-\log\tau_\alpha\otimes I]|\chi_\alpha\rangle=h(p_\alpha).
$$

两边具有相同Schmidt谱，所以同一个\(L_p\)同时出现。模投影式在参考的每个未归一化限制块正定时有通常矩阵对数意义。态对相对熵的全称结论可直接按块与支持判据延伸到非忠实参考；支持不相容时为∞，数值代码并未以小数截断代替这一声明。

式(6)是固定编码及固定参考内的有限差分恒等式。其\(\Delta\langle L_p\rangle\)可非零，因为中心概率可以改变；若连p也改，原式不能静默沿用。当前正文已明确区分。

## 3. Kraus产品、矩阵单位与最大性

不用作者的SVD核算即可写出A端Kraus算符：令\(t_{\alpha0}=p_\alpha\)、\(t_{\alpha1}=1-p_\alpha\)，则

$$
K_{\alpha b k}=\sqrt{t_{\alpha k}}
\sum_a|\alpha,a,k\rangle_A\langle\alpha,a,b|_C.
$$

两个Kraus产品仅在α、k相同才非零，为

$$
K_{\alpha b k}^{\dagger}K_{\alpha b' k}
=t_{\alpha k}\sum_a|\alpha,a,b\rangle\langle\alpha,a,b'|.
$$

0<p<1保证这些权重非零。八个不重叠支撑的矩阵独立，张成\(\bigoplus_\alpha(I_a\otimes M_b)\)。其对易子为M；B端交换a、b。独立有理消元得到对易约束秩56、对易子维数8，与作者结果相同。

另将编码写成四个整数系数矩阵\(V_{\alpha k}\)乘各自Schmidt振幅。对两端全部八个逻辑矩阵单位逐项核验\(O_{\rm phys}V_{\alpha k}=V_{\alpha k}O\)，64份整数身份全部精确成立。因全部伴随矩阵单位也在菜单内，未漏掉伴随交织；结论覆盖全码态和参考，随机样本不承担该全称证明。

## 4. 独立数值与精确见证

独立核算未导入或调用`complementary_entropy.py`。使用标准库Fraction／80位Decimal和整数NumPy：

- 对易约束有理秩56，最大可纠正代数复维8；
- 64份Schmidt分量矩阵单位交织精确为零；
- 三个原编码点、双侧共24个模投影对角项及熵／相对熵／有限差分，最大80位残差\(1.3\times10^{-79}\)；
- 剩余自由见证半迹距离精确1/4，声明的总辅助能量差精确Δ/2；
- 熵差为0.1308120359411369591292…；跨编码辅助相对熵为0.1438410362258904637196…。熵差严格正等价于27>16，不依赖浮点符号。

高精度核算使用的两个独立满秩逻辑态为扇区内乘积态：ρ中心权重(1/3,2/3)，σ中心权重(3/5,2/5)；各块两侧概率如下。

|侧／扇区|ρ|σ|
|---|---|---|
|a／0|(1/4,3/4)|(1/5,4/5)|
|a／1|(3/7,4/7)|(2/3,1/3)|
|b／0|(2/5,3/5)|(3/4,1/4)|
|b／1|(1/6,5/6)|(4/5,1/5)|

这些态只作独立算术交叉，不扩大作者科学组。可复算核心如下；不导入作者资产：

```python
from fractions import Fraction as F
from decimal import Decimal as D, getcontext
getcontext().prec = 80
def dec(x): return D(x.numerator)/D(x.denominator)
def ln(x): return dec(x).ln()
def ent(xs): return -sum(dec(x)*ln(x) for x in xs)
def rel(xs,ys):
    return sum(dec(x)*(ln(x)-ln(y)) for x,y in zip(xs,ys))
wr,ws=[F(1,3),F(2,3)],[F(3,5),F(2,5)]
rr=[[[F(1,4),F(3,4)],[F(3,7),F(4,7)]],
    [[F(2,5),F(3,5)],[F(1,6),F(5,6)]]]
ss=[[[F(1,5),F(4,5)],[F(2,3),F(1,3)]],
    [[F(3,4),F(1,4)],[F(4,5),F(1,5)]]]
errors=[]
for ps in ([F(1,2),F(1,4)],[F(1,4),F(1,2)],
           [F(1,8),F(3,8)]):
    ts=[[p,1-p] for p in ps]; hs=[ent(t) for t in ts]
    for side in (0,1):
        r=[wr[a]*rr[side][a][i] for a in range(2) for i in range(2)]
        s=[ws[a]*ss[side][a][i] for a in range(2) for i in range(2)]
        rp=[wr[a]*rr[side][a][i]*ts[a][k]
            for a in range(2) for i in range(2) for k in range(2)]
        sp=[ws[a]*ss[side][a][i]*ts[a][k]
            for a in range(2) for i in range(2) for k in range(2)]
        errors.append(abs(ent(rp)-ent(r)-sum(dec(w)*h for w,h in zip(wr,hs))))
        errors.append(abs(rel(rp,sp)-rel(r,s)))
        deltaK=-sum(dec(x-y)*ln(y) for x,y in zip(r,s))
        deltaL=sum(dec(x-y)*h for x,y,h in zip(wr,ws,hs))
        errors.append(abs(ent(rp)-ent(sp)-deltaK-deltaL+rel(r,s)))
        for a in range(2):
            for i in range(2):
                proj=-sum(dec(t)*ln(ws[a]*ss[side][a][i]*t) for t in ts[a])
                errors.append(abs(proj+ln(ws[a]*ss[side][a][i])-hs[a]))
assert max(errors)<D('1e-75')
print(max(errors))  # 1.3E-79
```

## 5. 相对熵、资源与误差范围

固定编码内的辅助态在两份区域状态中相同，故从相对熵消去。跨编码的辅助态不同，不能再消去。正文保留了具体跨编码读数及正相对熵，未混淆两个量词。

全部所列恢复和相对熵数据相同，不意味着全部实际任务相同，更不意味着A7的准备／资源／来源相同。本族中的χ是编码输入；其变化不能免费解释成规范冗余。固定V后L唯一和跨V自由完全相容。

有限误差部分采用已给的全任务输出距离；部分迹收缩再接有限维熵连续界。[Audenaert原文](https://arxiv.org/abs/quant-ph/0610146)中的距离须取半迹范数。报告采用的\(\epsilon\log7+h(\epsilon)\)在\(0\le\epsilon\le7/8\)单调，故用距离上界ε合法。它没有将小恢复误差直接换成固定编码距离，也没有由熵误差给无界应力误差。

不产生面积极值、Newton系数或GR的范围陈述准确。

## 6. 冻结后只读验收

已运行`verify_round1042.py`默认入口，它含科学结果只读复算；返回通过，8份作者资产、9份历史输入、9个本地链接。另以独立SHA256循环核对收据的所有OWN和历史条目，全部一致。

审阅所见收据SHA256：

`295b7dbeeeca06c76adcc4b0f16a730f38d842e789cedf8baa7ec92e541d7ef6`

主要资产：

|文件|SHA256|
|---|---|
|research_note_1042.md|5e538633726ed3dac6c756d70b3d2887c7b213cf9a56ca04c5dbbe40af89ff5d|
|complementary_entropy.py|ab4e1e720470480c0bcc6721b8bda8ebf1430d269f9073968f111a0566c1e70b|
|complementary_entropy_results.json|4a12afee1b118f4db6cc919eecde01baedf033a520221dfe9e1e94eece6a7970|

本独立审阅不在作者OWN中，不要求修改原收据。没有待修问题；不建议继续扩张本编码的规模或噪声样例。

## 7. 主线后审：任务身份与中心置换的表示商

2026-10-08追加。主线指出了一项应明确的量词；独立核验同意。此处限定前述审阅对正文§5的解释，不改作者冻结报告、代码或结果。

原成对例p=(1/2,1/4)与q=(1/4,1/2)可由同时交换α=0、1的逻辑／物理标签对应。因此固定α=0准备及固定辅助效果的概率差，确实反驳**保持这些任务身份时**L的唯一性；但单凭该对，不能证明商去中心标签置换后仍有两个不同模型。此前“不是规范冗余”的判断必须带上述任务身份限定。

冻结结果已包含第三族r=(1/8,3/8)，无需新增样例。因h在(0,1/2)严格递增，

$$
\operatorname{spec}L_p=\{\log2,h(1/4)\},\qquad
\operatorname{spec}L_r=\{h(1/8),h(3/8)\},
\qquad \max\operatorname{spec}L_r=h(3/8)<\log2.
$$

每个谱值在码空间中的重数均为4。最后的严格差也可完全改写为整数见证：

$$
\log2-h(3/8)=\frac18\log\frac{3^3 5^5}{4^8}>0,
\qquad 84375>65536.
$$

故在保持标准代数熵约定和物理访问二分的表示商下，中心置换、块内酉或两端局部酉都不能把p族与r族的中心熵谱同化。这给商去中心标签置换后仍有编码自由的证书，所用两个谱均已在原结果中，不新增科学组。

任意混合A、B的全局主动酉不是这里的被动表示商；若另行运输访问代数、准备或读数任务，必须连同那些对象重新说明等价。这个后审仍不将χ资源变化、任意整体A7任务或物理几何视为相同。
