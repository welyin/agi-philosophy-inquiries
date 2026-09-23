# 第316轮：有限类光片的面积平衡与不能删除的几何项

日期：2026-09-23。接续[315轮](research_note_315.md)及[304轮](research_note_304.md)；本轮10项检查。

## 1. 问题与输入

认知动机是边界变化的代价应完整计入。本轮检验：把局部视界的一阶关系延长到有限过程，会漏掉什么？

对接[Jacobson（1995）](https://arxiv.org/abs/gr-qc/9504004)、[Amsel、Marolf与Virmani的物理过程第一定律](https://arxiv.org/abs/0708.2738)及[Eling、Guedens、Jacobson的非平衡分析](https://arxiv.org/abs/gr-qc/0602001)。本轮复算几何恒等式，不将这些既有工具记为原创，也不预先令Ricci曲率等于物质应力。

额外输入：光滑四维Lorentz度规、无旋仿射类光束、无焦散的局部片和指定末端条件。这里的类光片尚未被证明是某个整体时空的事件视界；面积尚未被解释成量子熵。

## 2. 可独立计算曲率的度规

使用Rosen形式、单位横向坐标面积：

$$
ds^2=-2\,d\lambda\,dv+a(\lambda)^2dx^2+b(\lambda)^2dy^2,
\quad a''=-t_xa,\quad b''=-t_yb,
\quad a(L)=b(L)=1,\quad a'(L)=b'(L)=0.
\tag{1}
$$

常数t_x、t_y是输入的横向潮汐系数。正t的解为cos[√t(L−λ)]，负t改为cosh[√(−t)(L−λ)]。只取a、b正的区间；达到焦散的参数被代码拒绝。

v＝常数上的k＝∂_λ为仿射类光切向量。直接从度规得到

$$
A=ab,\qquad \theta=\frac{a'}a+\frac{b'}b,
\qquad \sigma^2=\frac12\left(\frac{a'}a-\frac{b'}b\right)^2,
\qquad R_{kk}=-\frac{a''}a-\frac{b''}b=t_x+t_y.
\tag{2}
$$

代码另外从度规一阶导数构造全部Christoffel符号，再差分和收缩计算R_λλ。t_x＝0.3、t_y＝0.1时得到0.399999999999285；这一步没有直接调用t_x＋t_y作为曲率结果。

## 3. 精确积分恒等式

由A′＝Aθ和Raychaudhuri关系：

$$
\theta'=-\frac{\theta^2}{2}-\sigma^2-R_{kk},
\qquad (A\theta)'=A\left(\frac{\theta^2}{2}-\sigma^2-R_{kk}\right).
\tag{3}
$$

乘以λ，在[0,ℓ]积分一次，可得有限面积差

$$
A(\ell)-A(0)=\ell A(\ell)\theta(\ell)
+\int_0^\ell \lambda A(\lambda)
\left(R_{kk}+\sigma^2-\frac{\theta^2}{2}\right)d\lambda.
\tag{4}
$$

这里只做微积分。最后一项的负号来自对变化面积测度的精确记账，不能把整个非线性余项都命名为非负“熵产生”。在本轮未来端点约定下，剪切项为正、膨胀修正为负；与过去端点的局部Jacobson展开比较时应先统一区间与权重。

若ℓ＝L，端点项因指定条件消失；任取中间截面时它一般不消失，而且可处在一阶。

例如t_x＝0.3、t_y＝0.1、L＝1、ℓ＝0.6：

| 面积差 | Ricci积分 | 剪切积分 | 膨胀积分 | 端点项 |
|---:|---:|---:|---:|---:|
| 0.1569158907 | 0.0666751551 | 0.0014147488 | −0.0053926326 | 0.0942186194 |

完整等式残差为5.6×10⁻¹⁷。仅保留Ricci项的误差约0.09024，远大于积分误差。

## 4. 一阶近似的可控范围

先取L＝1、t_x＝t_y＝ε/2。无剪切，且

$$
\Delta A=\sin^2\sqrt{\epsilon/2}
=\frac\epsilon2-\frac{\epsilon^2}{12}+O(\epsilon^3).
\tag{5}
$$

ε＝0.4、0.2、0.1、0.05时，相对ε/2的绝对余项依次为0.01298281、0.00328920、0.000827798、0.000207640，每次减半后约缩小4倍。可用一阶平衡的理由是余项受控，而非它们精确不存在。

再取t_x＝ε、t_y＝−ε。这个给定平面波度规Ricci平坦，但横向潮汐不为零：

$$
R_{kk}=0,\qquad
\Delta A=1-\cos\sqrt\epsilon\cosh\sqrt\epsilon
=\frac{\epsilon^2}{6}+O(\epsilon^4).
\tag{6}
$$

ε＝0.4时，面积增加0.0266565085，剪切积分0.0267591953，膨胀修正−0.0001026868。ε降到0.05时，ΔA/ε²＝0.1666656746，趋近1/6。这是明确反例：有限面积变化不能只由Ricci项代表。

这不是从认知原则生成引力波；它是给定度规中的几何控制。指定R_{kk}＝0并不意味着Riemann曲率或剪切为零。

## 5. 与物理热流的边界

仅当另行建立了相应场方程、boost能流和温度识别，才可把Ricci积分替换为热量：

$$
R_{kk}=8\pi G T_{kk},\qquad
\frac1{4G}\int\lambda A R_{kk}\,d\lambda
=2\pi\int\lambda A T_{kk}\,d\lambda=\frac{\delta Q}{T_U}
\quad(\hbar=1).
\tag{7}
$$

式(7)在本轮是条件性接口，未用于证明式(4)。因此不能先用它完成能流匹配，再声称从匹配独立推出同一Einstein方程。非最小耦合下还要共同处理Wald表面项。

10项检查覆盖独立曲率、Raychaudhuri微分关系、精确积分、端点反例、两种扰动阶数、剪切／膨胀的独立作用、仿射重标度和焦散拒绝。下一轮检查量子相对熵在嵌套可访问区域中提供什么约束；几何剪切与量子相对熵暂不认作同一个对象。

[代码](null_sheet_balance_audit.py)、[结果](null_sheet_balance_audit_results.json)、[核验](research_round_316_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/null_sheet_balance_audit.py
