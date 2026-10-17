# 1282独立审阅：正常Gauss校准与共同规范边势

2026-10-17。结论：**PASS_WITH_SCOPE**。独立核验1组，新公理0。已通读[证明](proof.md)、[作者代码](check.py)、[结果](results.json)、[NEXT](NEXT.md)及[报告](../research_note_1282.md)，并回读[557](../../archive_554_584/research_note_557.md)的源Hamiltonian、群生成元归一及读口范围。没有必修问题。

## 1. 分类与正常态的解析核验

四元数交换子分类不依赖数值秩：若实线性T与全部左单位四元数乘法对易，则由线性性对任意四元数q有 $T(q)=qT(1)$，所以T恰为右乘；反之右乘确与左乘对易。右乘的虚部是反对称矩阵，实对称部分只有标量。由两端独立局部作用分别消去U依赖，作者的 $A=aI,D=dI,B(U)=B_0R(U)$ 因而成立。实不可约并未被误用为一般实交换子只有一维。

$B_0B_0^T=(b_0^2+|\mathbf b|^2)I$ 给8维二次型的成对谱与正性判据，包括退化边界。正性等价于作者式(6)，并推出

$$
q,r\ge0,\qquad
(a-d)^2+4|\mathbf b|^2\le qr.
$$

对每个U，$(X_1,X_2)\mapsto(y,z)$ 正交，作者波函数归一系数正确；归一Haar积分不改变其范数。局部描述作用使y、z共同旋转，故波函数本身在Gauss物理部门。正有限宽度保证统一Schwartz性与群光滑性，相关性、动量及群导数不会产生非法delta态。它确可作为原物理空间中的公开准备输入；这不证明自主制备过程。

独立中心Gaussian矩给

$$
\langle F\rangle_{s,t}=sq+tr.
$$

固定同一 $t>0$ 的两份不同有限s使均值相等，当且仅当q=0。由上述正性不等式，立即得到 $a=d=b_0=\kappa\ge0,\mathbf b=0$。这里选择的是557原边势的形式，系数为 $\kappa_{\rm here}=\kappa_{557}v_{557}$；整体强度、群、表示和二次有效阶仍为输入。

## 2. 有限误差与完整能源

在y、z坐标中，相对 $\kappa_*=r/4$ 的目标矩阵，差为

$$
\widetilde M-\widetilde M_*=
\begin{pmatrix}qI/2&V\\V^T&0\end{pmatrix},
\qquad
\|V\|=\frac12\sqrt{(a-d)^2+4|\mathbf b|^2}.
$$

三角不等式和正性给作者式(16)。再乘二次型的 $1/2$ 及输入二阶矩 $B_2$，得到式(17)的 $B_2/4$ 因子。界对U一致，故保物质—链路相关与任意被动伙伴；并非无界乘法势的全空间算子范数或不同演化后的通道界。实际误差必须除以两份共同宽度的平方差，且 $R_{\max}$ 有独立预算；作者均已明确。

557的群归一为 $-\sum_AJ_A^2=3I/4$。从原配置函数c的三个群导数独立计算，

$$
\left\langle\sum_A|J_Ac|^2\right\rangle=3st,\qquad
\langle-\Delta_G\rangle=\frac{3(s-t)^2}{16st}.
$$

物质动能为 $1/s+1/t$；两个概率方差1的singlet贡献 $1/2$。原四次W通过完整Gaussian四阶矩得到作者式(21)，包括非零 $L_{12}$ 的交叉项。因此两次边势均值相同，没有被误说成完整能源相同。

本独审另选正定L且 $L_{12}\ne0$ 的有理核验点。两份态边势均值均为1，而完整账分别为

$$
\frac{705977}{72000},\qquad
\frac{178775173}{3528000}.
$$

这些数值只检验一般能源公式，不替换557冻结的参数，也不是新增实际制备证书。

## 3. 独立计算与回执

[独立代码](independent_check.py)不导入作者实现。采用：

- 四元数乘法表生成左右乘矩阵，Fraction消元精确得到实交换子维数4、对称部分维数1；
- 三个有理单位四元数核对完整两端协变、y/z块变换、谱及有限误差；
- 直接逐指标Gaussian配对矩，核两个正常态的原群动能、物质动能及十变量W账，使用非零 $L_{12}$；
- 分别核删除正性、左端完整协变、二次阶、固定差别宽度的反例，并保留零作用分支。

独立有限误差例取 $q=1/50,r=199/50$，两份均值为 $607/750,316/375$；每份读数误差预算为 $1/2000$，得到 $\varepsilon=103/5000$。系数界约0.153827，实际系数差与界的比值约0.743725。浮点只用于三组有限谱诊断，最大残差约 $8.88\times10^{-16}$；全称结论由解析证明承担，未进行额外随机扫描。

|实际执行|回执|结果|
|---|---|---|
|首次排他保存独立结果|2edd4a|exit 0，CREATED_ONCE|
|异工作目录默认只读复算|909f98|exit 0，READ_ONLY_MATCH|
|最终哈希与预算核对|5402a5|exit 0|

默认模式只读复算；--write 排他创建，不覆盖结果。

## 4. 资格与停止线

作者报告和NEXT准确保留了以下限制：

- 正二次类、Gauss作用、载体和校准角色是公开建模输入，尚未成为新增认知公理。
- 待测一般F的实际读口须登记有效权限；旧557固定F的仪器不能自动授权任意未知作用辨识。准备、参考、读数和记录成本均未免费化。
- 同均值选择只作用于不含额外singlet系数的声明类别；高阶正失配项和零作用仍可通过。
- 精确选回旧边势后可按原权限继承557、559；非零有限误差的静态势均值界本身不授予它们的全部动态结论。
- 不生成规范群、维数、时间分量、链路独立动力学或完整Yang–Mills理论。

本轮因此具有明确的同载体增量：以两份有限正常Gauss态的共同幅度校准，选择过去直接指定的边势形式，并给可审计误差与实际能源账。无需扩建比较器；后继仅按NEXT检查真正未受约束的作用输入。

## 5. 最终版本绑定

以下均为本次实际审阅版本；后生成验收文件不预先绑定。

|工作区根相对文件|SHA256|
|---|---|
|research_cognition_physics/archive_1281_1287/1282/proof.md|4c26761ed867a5387dba5b4865b0043e205c3f8beee4be65ea914e15cd0eab7e|
|research_cognition_physics/archive_1281_1287/1282/check.py|8d223354ffb901d786816411cd58ff3679efd6ecfda2ae3f167eb428e915f7be|
|research_cognition_physics/archive_1281_1287/1282/results.json|b0a923f27b1acf2d2b8deaaa047f1a76e418c8cbbd025696766f465bcaa19071|
|research_cognition_physics/archive_1281_1287/1282/NEXT.md|d8d0754df3ce30d4590c0e1823319a383ff010229767e12c4b1a8887c4ebe22c|
|research_cognition_physics/archive_1281_1287/research_note_1282.md|7b1d099a5adb7ce7c386417a439574d0ca46ecf27c4e74731fb762eedc054c44|
|research_cognition_physics/archive_1281_1287/1282/independent_check.py|8f635fe9a6b9cfce8020a050b69772c2d67fa53b501fe0b97faf46688004e016|
|research_cognition_physics/archive_1281_1287/1282/independent_results.json|fe6f497b1e8567af4db1bd04dc5eafa17af170013fea0a83fdeb73581a163bcf|
