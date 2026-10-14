# 1213：正能连续谱、纯记录与严格保能重置的资源边界

2026-10-15。作者稿，待独立验收。新公理0。[代码](check.py) · [结果](results.json) · [下一项](NEXT.md)。本轮服务[四律目标](../_shared/goal_four_laws.md)的第三律资源来源问题，不宣称已完成第三律。

## 1. 净增、认知动机和命题范围

[1211](../1211/proof.md)已证明：原同h有限材料加零H记录，严格保裸总能时，无基态占据部门不能被冷却；纯同h基态库存可以交换进目标。这里不重证谱尾或该路由，而问：**没有同h基态库存，若加入同类正能正常钟的能源容量，纯零H记录能否承接信息并实现精确重置？**

认知动机是资源角色须由实际改变判定：名字叫“钟”“记录”不表示其只能计时或可以免费复位。既有[421](../../archive_370_428/research_note_421.md)、[457](../../archive_429_466/research_note_457.md)的未知关联/完整账，以及[1209](../1209/proof.md)的辅助归还区别直接复用。421的“材料约化是同一精确酉”前提不适用于本轮的重置通道。

**本轮证明的限定反例：** 正下界生成元、有限维材料和记录、正常连续谱辅助、有限一二能源矩、有限最坏能源转移和完整内账，仍不排除指定有限时刻材料精确落入基态。相应联合耦合在此明确新增，不是原QQ、旧实际局域钟或CC/SR几何许可。连续谱辅助并非有限维；初态Gaussian不是最大能源有界；输出位置矩甚至可能发散。因此它只排除上述**弱能源条件**推出“端点不可达”的命题，不是“所有资源均有界”、原公理整体、实际局域机器或一般第三律的反模型。

## 2. 同h材料、正能钟与纯空白资源

取1201—1206同一个四维h，其非简并谱为$e_0<e_1<e_2=1<e_3$。以下基底均为h的能源本征基，不另选一套材料谱。定义

$$
\bar H_S=h-e_0I=\sum_{i=0}^{3}\Delta_i|i\rangle\langle i|,
\quad \Delta_i=e_i-e_0,\quad \Delta_0=0,
\quad w=\Delta_3<1.094.
\tag{1}
$$

钟C取通常自由粒子Hilbert空间$L^2(\mathbb R,dp)$，$H_C=P^2/(2M)$，$M>0$。K维数4，$H_K=0$，提供独立纯空白$|0\rangle_K$。未知材料可与任意被动R相关：

$$
\rho_{SR}\otimes|\psi\rangle\langle\psi|_C\otimes|0\rangle\langle0|_K.
\tag{2}
$$

纯钟和纯K都是明确初始资源；不把其制备来源视为已经推出。没有另给同h纯基态库存。取实动量Gaussian

$$
\psi(p)=(2\pi s^2)^{-1/4}\exp[-(p-p_c)^2/(4s^2)],\qquad s>0.
\tag{3}
$$

它正常化、有有限初始位置矩，且

$$
\mu=\langle H_C\rangle=\frac{p_c^2+s^2}{2M},\qquad
\nu=\langle H_C^2\rangle=\frac{p_c^4+6p_c^2s^2+3s^4}{4M^2}.
\tag{4}
$$

此处仅保留原正能钟的Hilbert空间/自由Hamiltonian类型；**新增耦合不具有原1203—1205空间门窗的实际钟身份**。

## 3. 半线能谱上的上移是等距，不是酉平移

定义酉能谱表示$\mathcal U:L^2(\mathbb R,dp)\to L^2(\mathbb R_+,dE)\otimes\mathbb C^2$：

$$
(\mathcal U\psi)_\sigma(E)=\left(\frac{M}{2E}\right)^{1/4}
\psi(\sigma\sqrt{2ME}),\quad \sigma=\pm1.
\tag{5}
$$

Jacobian来自$|dp/dE|=\sqrt{M/(2E)}$；$E=0$是零测集。$H_C$变为乘法E，自伴定义域为$\{f:\int E^2|f(E)|^2dE<\infty\}$。

对$\Delta\ge0$定义

$$
(T_\Delta f)_\sigma(E)=\mathbf1_{E\ge\Delta}f_\sigma(E-\Delta),
\quad (T_\Delta^\dagger f)_\sigma(E)=f_\sigma(E+\Delta).
\tag{6}
$$

变量代换严格给出

$$
T_\Delta^\dagger T_\Delta=I,\qquad
T_\Delta T_\Delta^\dagger=P_{[\Delta,\infty)},\qquad
H_CT_\Delta=T_\Delta H_C+\Delta T_\Delta.
\tag{7}
$$

最后一式在$D(H_C)$成立，并有对所有t的有界强恒等式

$$
e^{-itH_C}T_\Delta=e^{-it\Delta}T_\Delta e^{-itH_C}.
\tag{8}
$$

这不是无下界$H_C=P$，也不是将半线边界藏入一个双向酉平移。以下将钟向量按式(5)识别到能谱表示，故$T_\Delta\psi$是$T_\Delta\mathcal U\psi$的简记；若完全使用动量表示，相应算符为$\mathcal U^\dagger T_\Delta\mathcal U$。

## 4. 全Hilbert空间上的自伴酉与强保能

张量顺序为S、C、K。对$i=1,2,3$定义

$$
V_i=|0\rangle\langle i|\otimes T_{\Delta_i}\otimes|i\rangle\langle0|,
\quad D_i=V_i^\dagger V_i=|i\rangle\langle i|\otimes I_C\otimes|0\rangle\langle0|,
$$
$$
R_i=V_iV_i^\dagger=|0\rangle\langle0|\otimes P_{[\Delta_i,\infty)}\otimes|i\rangle\langle i|.
\tag{9}
$$

六个$D_i,R_i$两两正交。令

$$
W=I-\sum_{i=1}^3(D_i+R_i)+\sum_{i=1}^3(V_i+V_i^\dagger).
\tag{10}
$$

每对$D_i\oplus R_i$上是由$T_{\Delta_i}$识别的交换；补空间恒等。故W有界、$W=W^\dagger$、$W^2=I$。特别是$S=0,K=0$全钟子空间保持；$S=0,K=i$但$E<\Delta_i$的子空间也保持，不需要不存在的负能源输出。

令$Q=\bar H_S+H_C$。由式(8)，$e^{-itQ}V_i=V_ie^{-itQ}$，其伴随以及式(9)投影也如此。因此

$$
[W,e^{-itQ}]=0\quad(\forall t\in\mathbb R).
\tag{11}
$$

W与Q的全部谱投影对易，保持$D(Q)$及各能源矩定义域，$W^\dagger QW=Q$；这比未处理定义域的形式交换子声明更强。全W只连接钟能源差为$0$或$\pm\Delta_i$的扇区，最坏谱能源移不超过w。

## 5. 未知SR的精确重置，信息完整转存

令$T_0=I$，对固定初始钟定义

$$
J_\psi|i\rangle=T_{\Delta_i}|\psi\rangle_C\otimes|i\rangle_K.
\tag{12}
$$

由K正交和$T_i^\dagger T_i=I$，$J_\psi^\dagger J_\psi=I_S$。式(10)对所有输入向量严格满足

$$
W(|i\rangle_S|\psi\rangle_C|0\rangle_K)
=|0\rangle_S J_\psi|i\rangle.
$$
$$
(W\otimes I_R)(\rho_{SR}\otimes|\psi,0\rangle\langle\psi,0|)(W^\dagger\otimes I_R)
=|0\rangle\langle0|_S\otimes(J_\psi\otimes I_R)\rho_{SR}(J_\psi^\dagger\otimes I_R).
\tag{13}
$$

式(13)不后选择，不要求材料平稳、已知或与R独立。R的全部相关和未知相干由同一等距进入CKR，没有信息消失。目标基态概率对全部未知SR精确为1。谱移后的钟状态之间可以不正交；记录K正是保证各输入正交性的物理承接者。

## 6. 正有界相互作用与有限端点，不是关闭后的冷产品

取指定$\tau>0$，$\hbar=1$，定义新增全局相互作用

$$
G=\frac{\pi}{2\tau}(I-W),\qquad 0\le G\le\frac\pi\tau I,
\qquad H_{\rm parent}=Q+G.
\tag{14}
$$

G有界且与Q强对易，故$H_{\rm parent}$自伴于$D(Q)$，非负。W的本征值为$\pm1$，于是

$$
e^{-i\tau H_{\rm parent}}=e^{-i\tau Q}e^{-i\tau G}
=e^{-i\tau Q}W.
\tag{15}
$$

这给有限时刻的精确材料基态。实际末等距是$J_{\psi,\tau}|i\rangle=e^{-i\tau H_C}T_{\Delta_i}\psi\otimes|i\rangle$，不能删掉自由相位后误称全通道恰为W。任意SR的相对相位仍完整保存。

G常驻，端点不满足“交互已关闭”。例如$t=2\tau$时$e^{-i2\tau G}=I$，材料回到原人口；激发输入可以再次激发。未给无代价关门、释放后持久低温、实际触发钟读出或循环操作。该有限参数时间也尚未与旧SR实际计时合同同一化。

## 7. 全能源、熵与角色账

记$p_i=\langle i|\rho_S|i\rangle$、$\bar\Delta=\sum_i p_i\Delta_i$、$\overline{\Delta^2}=\sum_i p_i\Delta_i^2$。式(13)—(15)给

$$
\Delta E_S=-\bar\Delta,\qquad \Delta E_C=\bar\Delta\in[0,w],\qquad
\langle H_C^2\rangle_{\rm out}=\nu+2\mu\bar\Delta+\overline{\Delta^2}
\le\nu+2\mu w+w^2.
\tag{16}
$$

记录裸能源始终0。初态$\langle W\rangle=p_0$，故

$$
\langle G\rangle_{\rm in}=\langle G\rangle_{\rm out}
=\frac{\pi}{2\tau}(1-p_0).
\tag{17}
$$

Q与G各自守恒；没有把交互能源设零。式(16)是能源转移，不自动命名为热或可循环提取的功。

初始CK纯且与S独立，终态$S(CK)=S(\rho_S)$；初始SR的全部联合谱在CKR中保留。对能源对角输入，K末态是$\operatorname{diag}(p_i)$，熵为$H(p)$；一般相干输入的K还含$\rho_{ij}\langle T_j\psi|T_i\psi\rangle$，不可将其擅自去相干。钟的末边缘是至多4个纯态的混合，因而其熵有限且不超过$\ln4$。对最大混合材料，K变成$I/4$，独立空白的纯度已被消耗。CK没有返回初态；反复运行必须另计恢复或新空白来源。

**Gaussian的重要限制。** 式(3)在$p=0$非零。对$\Delta>0$，把$T_\Delta\mathcal U\psi$变回动量表示，在$p>\sqrt{2M\Delta}$得到

$$
\phi_\Delta(p)=\left(\frac{p}{\sqrt{p^2-2M\Delta}}\right)^{1/2}
\psi(\sqrt{p^2-2M\Delta}),
\tag{18}
$$

另一动量分支同理。阈值处行为为$(p-\sqrt{2M\Delta})^{-1/4}$，可平方积分，但不在$D(X)$：其弱导数不是$L^2$。对有激发人口的输入，不能给予有限输出$\langle X^2\rangle$或空间局域钟资格。有限能矩、正常态、正生成元绝不意味着所有物理资源均有限。本轮不额外优化波包或将这一候选改造成实际局域制冷机。

## 8. 同一构造的删除检验

1. **删能源容量。** 只保一个同h目标和零H记录，则严格保能酉与非简并h对易，目标人口不变（1211复用）。把式(10)的谱移直接换成恒等，会在每个激发输入少掉$\Delta_i$的完整能源，不能说仍是保能重置。
2. **删指定纯空白。** 用$I_K/4$替代初始$|0\rangle$，同一W对能源对角p和Gaussian钟的末基态概率为

$$
F_{\rm mixed\ K}=\frac14+\frac{p_0}{4}\sum_{i=1}^3\Pr_\psi(H_C<\Delta_i)
\le\frac14+\frac{3p_0}{4}<1\quad(p_0<1).
\tag{19}
$$

这是同一器件的删除检验，不是所有纯钟加混合记录的普遍不能纯化定理；纯钟仍可能提供别的资源。
3. **删物理信息承接。** 若企图只保S、C，并把不同$|i\rangle\psi$都映为$|0\rangle T_i\psi$，本实Gaussian的不同谱移波函数在重叠支上同号，因此内积严格正，映射不保正交，不可能是酉。完整K可以被忽略在观察描述外，却不能被物理删除且仍声称总信息消失。
4. **增严格循环归还。** 对混合S，式(13)使CK熵增加$S(\rho_S)>0$，不能同时返回其原独立纯态；不得把本次重置称为可免费重复的循环。已有421/457资源账适用，而材料精确酉限制不适用。

## 9. 与第三律文献、认知原则和后继任务的关系

[Masanes–Oppenheim的原文“Fundamental Assumptions”及“Thermal Machines”](https://arxiv.org/html/1412.3828)除能源守恒外，还规定独立初态、正温热浴、浴态密度/体积条件，并用能量平移协变限定工作载体的角色；时间界另接机器访问体积及最坏工作随时间增长。本文纯K不是正温热浴，半线$T_\Delta$不是双边工作平移群，改变后的CK承担信息/熵。本构造不满足该不可达定理的完整合同，不构成对其反证；这些文献条件也尚未由本项目现公理全部导出。

FUCP的理想操作与CC的完整参与者/记录要求，允许提出并记清此类数学候选，不自动授予G的实际局域实现。D-P的完整能源要求在式(16)—(17)有明确账；CO/SR空间、实际钟和有限传播身份没有被本构造证明。新增耦合及初始纯资源应与已采用公理分层。

**本轮改变的下一问题：** 仅排除同h纯基态库存仍不够；纯空白、控制和钟的资源角色须同时说明。若第三律研究“从正温来源产生新纯度”，必须辨认哪些辅助可消耗纯度，哪些被认证为工作/控制载体而不得偷作未计熵汇，以及终态是否要求释放和有限几何资源。该来源与角色合同仍是待判认知猜想，不等于把“不得达到零温”加入定义，也不要求所有普通工作态处于平衡。

## 10. 复算与证明分工

代码用原h的已认证根区间和已存冷参考人口；在显式有限受支持配对子空间核交换、守恒标签、未知SR等距及完整相位。14个固定/配对扇区各带同一$L^2([0,\ell],du)$，Q在其中是能源偏置矩阵加$uI$，W保持同一u。矩阵只核扇区系数部分，公共自由相位$e^{-i\tau u}$可以提出；有限个常数包张成的空间本身并非Q的不变子空间。128维系数矩阵的其余恒等补全仅供诊断，不代替那些扇区之外的真实W。

另用有理能量包的积分恒等式核半线等距/阈值投影，检查Gaussian有限能矩、实际人口能源账和删除检验。有限矩阵既不取代式(5)—(11)的连续谱证明，也不声称有限截断钟能对任意输入完成精确降温。没有参数扫描、图像、原QQ实现或不可达速率优化。默认复算只读；$\texttt{--write}$仅允许首次创建尚不存在的结果文件，不用于覆盖已保存或冻结结果。
