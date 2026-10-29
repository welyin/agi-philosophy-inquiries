# 1333：均匀单步关系语法的谱障碍与有限扰动

2026-10-19。[事前准入](../_admission/round1333_admission.md)。新公理0。本文给一个候选类别的必要条件排除；不交付量子形成器件，也不证明所有层级认知机制失败。

## 1. 动机、输入及判定对象

认知动机是多个角色复用同一交接步骤，两部门采用相反的接续方向。以下是其**额外候选语法**，不是现行公理的推论。令正交内部角色数为$n$，定义开链移位

$$
S=\sum_{k=0}^{n-2}|k+1\rangle\langle k|,
\qquad Y_u=aI+bS,\qquad Y_d=cI+dS^\dagger. \tag{1}
$$

四个复数任意；端点、均匀权、仅一次接续、共同左角色字典及相反方向是输入。链不是空间，也尚未证明其标签就是物理代际。两部门采用正则动能，在同一匹配尺度以$M_f=hY_f$直接解释Dirac质量，$h\ne0$是共同因子。奇异值按降序记为$\sigma_1\ge\sigma_2\ge\cdots$。规范群、表示、尺度和完整来源沿原条件菜单登记，未由(1)生成。

1316的[重中介链](../../archive_1313_1322/1316/proof.md)保留任意入口$Y_0$，其$Y=Y_0F$仍可取任意满秩矩阵；[1317](../../archive_1313_1322/1317/proof.md)审查多入口权限。这里直接限制整个味矩阵，旧定理不替代本轮判据。min–max和奇异值扰动是使用的成熟数学工具。

## 2. 所有$n\ge3$的谱障碍

先取$Y=aI+bS$，记$A=|a|^2$、$B=|b|^2$、$r=|ab|$、$P_k=|k\rangle\langle k|$。因$SS^\dagger=I-P_0$，

$$
H=YY^\dagger=(A+B)I-BP_0+a^*bS+ab^*S^\dagger. \tag{2}
$$

$n\ge3$时，子空间$L=\operatorname{span}\{|n-2\rangle,|n-1\rangle\}$不含$|0\rangle$。$H$压缩到$L$的两个特征值恰为$A+B\pm r$。降序min–max给

$$
\lambda_2(H)\ge A+B-r. \tag{3}
$$

另一方面，$\|S\|=1$，故

$$
\lambda_1(H)=\|Y\|^2\le(|a|+|b|)^2=A+B+2r. \tag{4}
$$

由$A+B\ge2r$，非零$Y$满足

$$
\frac{\sigma_2(Y)^2}{\sigma_1(Y)^2}
\ge\frac{A+B-r}{A+B+2r}\ge\frac14,
\qquad \boxed{\sigma_2(Y)\ge\tfrac12\sigma_1(Y)}. \tag{5}
$$

对$Y_d$，$S^\dagger S=I-P_{n-1}$；取前两维压缩即得同一结论。也可由链反射$RSR^\dagger=S^\dagger$运输。证明涵盖全部复相位和全部$n\ge3$，不依赖有限扫描。常数$1/2$只作统一充分下界，未声称最优或可达。

### 边界

- $a=b=0$时不定义比值，但齐次不等式(5)仍为$0\ge0$；零部门不能匹配非零质量。
- $b=0,a\ne0$时比值为1；$a=0,b\ne0$且$n\ge3$时有$n-1$个$|b|$及一个零，前两奇异值仍相等。
- $n=2$不满足避开端点的二维条件。$Y=S_2$给比值0；满秩例$a=b=1$给比值$(3-\sqrt5)/2<1/2$。不能将(5)延伸到$n=2$。
- $n=1$没有第二奇异值，不属于本判据。

## 3. 有限扰动不能以任意小代价隐藏障碍

令$Y\ne0$、$\widetilde Y=Y+\Delta Y\ne0$，并以谱范数定义相对改动

$$
\delta=\frac{\|\Delta Y\|}{\|Y\|}. \tag{6}
$$

Weyl奇异值不等式与(5)直接给

$$
\sigma_2(\widetilde Y)\ge(1/2-\delta)\|Y\|,
\qquad \sigma_1(\widetilde Y)\le(1+\delta)\|Y\|. \tag{7}
$$

若已知$\delta\le\epsilon<1/2$，则

$$
\boxed{\frac{\sigma_2(\widetilde Y)}{\sigma_1(\widetilde Y)}
\ge\frac{1/2-\epsilon}{1+\epsilon}}. \tag{8}
$$

反之，对任意相对改动，若目标要求该比值不超过$0\le\rho<1/2$，(7)仍可直接重排，不需预先假设$\delta<1/2$：

$$
1/2-\delta\le\rho(1+\delta)
\quad\Longrightarrow\quad
\boxed{\delta\ge\frac{1/2-\rho}{1+\rho}}. \tag{9}
$$

例如$\rho=0.1$需要$\delta\ge4/11\simeq0.363636$。这指相对于候选矩阵算子范数的改动，不是每个矩阵元的相对误差；下界也不是能达到目标的充分条件。任意左右酉换基不改变奇异值和谱范数，因此改变表示不能避开此必要界。将额外动能、阈值或多步作用隐入$\Delta Y$时，必须另证明其大小及来源。

## 4. 普通Dirac质量／混合中的五个连续自由

本节只取共同左酉换基、两部门独立右酉换基的普通Dirac字典。**不借此删除强CP的$\theta$或其它固定来源相位**；完整形成过程若有额外相位参考，也须另计。

在$a,b,c,d$均非零的泛型点，两部门整体相位可使$a,c$实正；共同对角角色换基将$S$乘一个相位，可再使$b$实正。只剩

$$
|a|,\ |b|,\ |c|,\ |d|,\qquad
\phi=\arg\frac{bd}{ac}. \tag{10}
$$

这首先给至多五个连续自由。要排除进一步连续冗余，令$C=|c|^2,D=|d|^2$，$H_u=Y_uY_u^\dagger,H_d=Y_dY_d^\dagger$。三角性给

$$
\det H_u=A^n,\quad\operatorname{tr}H_u=nA+(n-1)B,
\qquad
\det H_d=C^n,\quad\operatorname{tr}H_d=nC+(n-1)D. \tag{11}
$$

这四个换基不变量分别恢复$A,B,C,D$。混合不变量为

$$
\operatorname{tr}(H_uH_d)
=n(A+B)(C+D)-B(C+D)-D(A+B)
 +2(n-1)\operatorname{Re}(a^*bc^*d). \tag{12}
$$

固定四个正模后，它恢复$\cos\phi$；在$\sin\phi\ne0$处局部变化非零。因此泛型**恰有五个连续物理自由**。离散重标记、CP共轭及特殊零／退化子集不增加连续维数。与不受限两部门Dirac矩阵相比，(1)确实减少质量／混合资料，却仍不确定这些五个量或$n$。

## 5. 该语法并不强迫CP守恒

固定约定$\mathcal C=[H_u,H_d]=H_uH_d-H_dH_u$。定义

$$
u=a^*b,\quad v=cd^*,\quad
w=uv^*-u^*v=2i\operatorname{Im}(a^*bc^*d).
$$

由$[S,S^\dagger]=P_{n-1}-P_0$，直接展开得

$$
\begin{aligned}
\mathcal C={}&w(P_{n-1}-P_0)
 +B\bigl(v|1\rangle\langle0|-v^*|0\rangle\langle1|\bigr)\\
&+D\bigl(u|n-1\rangle\langle n-2|-u^*|n-2\rangle\langle n-1|\bigr).
\end{aligned} \tag{13}
$$

$n=3$时，两条边共享中点，$\operatorname{tr}\mathcal C=0$且
$\det\mathcal C=w(B^2|v|^2-D^2|u|^2)$，故$\operatorname{tr}\mathcal C^3=3\det\mathcal C$。$n\ge4$时，端点分别形成两个$2\times2$块，其余为零；用$\operatorname{tr}X^3=(\operatorname{tr}X)^3-3\det X\operatorname{tr}X$得到相同结果。因此所有$n\ge3$都有

$$
\boxed{\operatorname{tr}[H_u,H_d]^3
=6i\operatorname{Im}(a^*bc^*d)|bd|^2
\bigl(|bc|^2-|ad|^2\bigr)}. \tag{14}
$$

CP保持会使此换基不变量与其复共轭相等；但反Hermitian交换子的三次迹纯虚，因此非零(14)是CP破坏的充分证书。取$a=b=1,c=2,d=i$即得$36i$，对所有$n\ge3$成立。$n=3$且质量非简并时它还对应通常的Jarlskog判据；这里无需指定混合矩阵的符号约定。对$n>3$不把此单一不变量的零值宣称为CP守恒的完整判据。

该结果说明谱失败不能归咎于“没有复相位”。特殊实相位或其它使(14)消失的子集仍须按其自身情况判定。

## 6. 有共同尺度的历史物理示例

物理比较限定为直接三代身份$n=3$。仅取[Xing、Zhang、Zhou，2008，Table II](https://arxiv.org/html/0712.1419)的历史$\overline{\mathrm{MS}}$、$\mu=M_Z$数据作为宽裕对照：单位统一为GeV，$m_c=0.619\pm0.084$、$m_t=171.7\pm3.0$、$m_s=0.055^{+0.016}_{-0.015}$、$m_b=2.89\pm0.09$。这不是最新拟合，也不是新实验；不引用该文的Majorana或MSSM分支。

即把所列误差端点当作保守矩形范围，分别取分子上端、分母下端，也只有

$$
\frac{m_c}{m_t}\le\frac{0.703}{168.7}<0.0042,
\qquad
\frac{m_s}{m_b}\le\frac{0.071}{2.80}<0.0254. \tag{15}
$$

矩形端点计算不声称联合统计置信度；这里只需远宽于这些数的$\rho=0.1$。在声明的**直接同尺度质量身份**下，两个部门都与(5)不相容；小于$4/11$的相对矩阵改动不能使任一部门达到这一宽裕目标。对$n>3$，(5)仍是全体代际的数学约束，但不能将三代的经验排序直接当作额外未知代际的观测。

没有RG、阈值、动能及场字典匹配时，本轮不排除(1)作高尺度祖先规则，也不宣称任意后继理论失败。若要改变这些前件，必须另登记机制及其必要检验，不能事后称为同一被排除类别。

## 7. 数值核验与停止线

[check.py](check.py)用NumPy直接构造有限矩阵，检查(2)—(5)、边界、(7)—(9)、不变量(11)—(14)、换基及相位运输，并保存[results.json](results.json)。诊断含不同$n$、复相位、零系数和固定扰动；数值残差只核实现，不证明全$n$定理，不构成实验测量或浮点区间认证。

执行方式为`--write-new`以独占创建结果，`--check`只读重算比较；不写其它轮次。不新造共同父器件、不宣称已有实际形成来源。独审和主线验收另行绑定本轮文件。

**停止：** 本均匀一次语法退出直接低能三代质量候选。不继续调四标量、代际数或参考余额来修复已证明的谱障碍。后继最多先作计0多步语法／其它形成机制筛查；没有1334准入，也未选定该替代。完整认知生成目标继续开放。
