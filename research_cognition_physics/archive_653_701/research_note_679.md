# 第679轮：原物理来源的反射、接触抵消与共同内积

日期：2026-10-02。接[678](research_note_678.md)、[已执行入口](679/drafts/physical_reflection_entry.md)。[代码](679/joint_physical_source_reflection.py)、[结果](679/joint_physical_source_reflection_results.json)、[核验](679/research_round_679_checks.json)、[条件账](679/unified_physics_condition_ledger_679.md)。两组新核验、十六式；主代理审查，无独立代理审查。

## 1. 本轮问题与旧成果复用

已读取README、research_direction、RESEARCH_STATE、678及679入口、已有结果；未发现运行中的Python。回查661、668—669、673—678：674只完成标量及质量来源反射；678实现全部来源的局部辅助表示，但局部表示本身没有证明物理反射正性。

本轮补上**原物理Weyl子代数全部Grassmann来源的反射共轭关系**，包含有限软调节、原精确极限、奇异质量及零权重。关键是两份观测之差为高斯方程项，且来源接触项精确抵消。完整原平均上的物理半区型由此为Hermitian。正性、非零归一及原H_F身份仍需分别证明。

依用户提醒，重新核对下表原报告；继承结果不重复计入新检查或新定理：

|旧轮次|直接复用|仍需对应或不可恢复的假设|
|---|---|---|
|[382](../archive_370_428/research_note_382.md)|完整固定迹qubit反向可分接口给三维上界|实际完整位置壳及其物理接口|
|[383](../archive_370_428/research_note_383.md)|真实连续自由对合足以承担上界|不重新要求预先给标准几何对径|
|[384](../archive_370_428/research_note_384.md)|全族可逆重定向、协变与对比完整给下界；极小作用足够|额外Lipschitz／Hausdorff输入已经消去|
|[386](../archive_370_428/research_note_386.md)|非退化强缩放连接真实邻域；有限尺度上界证书|有限认证的覆盖、误差和变化模仍有量词|
|[425](../archive_370_428/research_note_425.md)|局部紧端点群、一致双侧半幅和成本收缩给NSS／Lie坐标|这是386的替代路线；不叠加两者，也不要求群交换|
|[522](../archive_467_530/research_note_522.md)|指定热参考及全域绝对差坐标任务给维数下界|不能替换成仅局部差读或角分辨率|
|[523](../archive_467_530/research_note_523.md)|实际方向仪器与热下界有条件性共同实现|共同实现已存在；需接到当前h、s、CAR／Gauss及自主过程，二维角读反例保留|

本轮不生成新空间，不调整目标，认知系统具体设计继续后置。

|层次|本轮地位|
|---|---|
|认知动机|同一物理对象的来源、反射和内积相容|
|继承输入|原有限盒、Wilson／overlap、16通道、M、Pφ、H_b及全部辅助／规范平均|
|解析增量|方程项与接触抵消；全部物理来源反射；完整物理型Hermitian|
|数值验证|非投影、不同谱秩、奇异质量代数诊断；原非平坦背景0／2／4／6来源|
|未完成|动态规范反射正性、非零归一、原H_F时间、共同连续极限及量子GR|

反线性、反转Grassmann乘积次序的反射约定继承661，见[Kikukawa—Usui，§II式(5)—(8)及§IV](https://arxiv.org/html/1005.3751v3#S4)。该文给自由overlap及非规范Yukawa情形的正性，不能替代当前动态规范模型的证明。本轮只读文字，无图像检验。

## 2. 无质量核和原物理观测

星号为逐元素共轭，†为共轭转置，T为普通转置。使用原实手征框J±、Q±=J±J±ᵀ；L为时间格点反射乘γ₀，L实、对称、L²=I。原M按E同步反射。n为原单粒子维数，是64的倍数；r=n/2。

$$
M^{\mathsf T}=-M,\quad \bar M=M^\dagger=-M^*,\quad
M^\dagger M=I,\quad [M,Q_\pm]=0,\quad L^{\mathsf T}M^\theta L=M .
\tag{1}
$$

只需ε Hermitian及下列反射合同，不需要ε²=I。677有限有理函数为奇函数，与精确sign一样满足：

$$
D=\tfrac12(I+\gamma_5\varepsilon),\quad D^\dagger=\gamma_5D\gamma_5,\quad
T=\gamma_5L,\quad \varepsilon^\theta=-T\varepsilon T^\dagger,\quad
D^\theta=LD^\dagger L .
\tag{2}
$$

令ξ=(ψ,χ)。由P_v=Q_-+γ₅D，673／677的核与观测重写为

$$
N_0=\begin{pmatrix}MQ_+&-D^{\mathsf T}\\D&M^\dagger Q_-\end{pmatrix},
\qquad Y=\Phi\xi,\quad
\Phi=\begin{pmatrix}
J_-^{\mathsf T}(I-D)&0\\
-J_+^{\mathsf T}M(I-\tfrac12D)&J_+^{\mathsf T}
\end{pmatrix}.
\tag{3}
$$

这里没有把软P_v当成精确投影。设C=J₊ᵀLJ₋，在原框下为实置换，定义

$$
R_\xi=\begin{pmatrix}0&L^{\mathsf T}\\L&0\end{pmatrix},\qquad
R_Y=\begin{pmatrix}0&C^{\mathsf T}\\C&0\end{pmatrix},\qquad
R_\xi^{\mathsf T}N_0^\theta R_\xi=-N_0^*,\quad \det R_\xi=1 .
\tag{4}
$$

最后两式由(1)—(2)及n偶得到，只是裸无质量高斯核的反射。实际物理观测反射为

$$
\widetilde\Phi:=R_Y\Phi^\theta R_\xi
=\begin{pmatrix}
J_-^{\mathsf T}&-J_-^{\mathsf T}M(I-\tfrac12D^\dagger)\\
0&J_+^{\mathsf T}(I-D^\dagger)
\end{pmatrix}.
\tag{5}
$$

它通常不等于Φ*。661已证的裸变量时间支撑障碍仍有效。

## 3. 方程项和接触抵消的解析证明

记A=N₀*、F=Φ*，后文多项式用f、g。取不依赖D的局部矩阵

$$
S=\begin{pmatrix}
-\tfrac12J_-^{\mathsf T}M&J_-^{\mathsf T}\\
J_+^{\mathsf T}&\tfrac12J_+^{\mathsf T}M^\dagger
\end{pmatrix},\qquad \widetilde\Phi=F+SA .
\tag{6}
$$

逐块相乘可验：SA上行的ψ、χ块分别为J₋ᵀD*、−J₋ᵀM(I−D†/2)，正是(5)减Φ*的上行；下行用γ₅-Hermiticity与MM†=I得到。

接触项必须独立核对。写D的手征块为[[a,b],[c,d]]、M=diag(m₊,m₋)，两侧直接乘法给同一矩阵：

$$
SF^{\mathsf T}=\widetilde\Phi S^{\mathsf T}
=\begin{pmatrix}
-\tfrac12m_-(I-d^\dagger)&-\tfrac14m_-b^\dagger m_+^\dagger\\
-c^\dagger&-\tfrac12(I-a^\dagger)m_+^\dagger
\end{pmatrix},
\qquad SF^{\mathsf T}-FS^{\mathsf T}-SAS^{\mathsf T}=0 .
\tag{7}
$$

a、b、c、d仅为D块名，不是调节参数。末式由(6)代入第一式得到。A可逆时，Aᵀ=−A给

$$
(F+SA)A^{-1}(F+SA)^{\mathsf T}
=FA^{-1}F^{\mathsf T}
+SF^{\mathsf T}-FS^{\mathsf T}-SAS^{\mathsf T}
=FA^{-1}F^{\mathsf T}.
\tag{8}
$$

因此全部高斯收缩相同，原Pfaffian权重也未变；所有高阶矩由同一Pfaffian恒等式控制，而非由几个二点样本猜测。

## 4. 奇异核和零权重的多项式延伸

逆只在证明的稠密集使用。满足γ₅-Hermiticity的D组成实仿射空间，未归一增广Pfaffian是D实部、虚部的多项式。存在明确非奇异点：

$$
D=tI,\quad t\ne0\quad\Longrightarrow\quad \det N_0=t^{\,2n}\ne0 .
\tag{9}
$$

按手征重排后，每个2r块的行列式为t^(2r)。从任意D向I/2插值，行列式不是零多项式，故可逆点稠密；ε向0插值还保持Hermitian收缩。这是证明用的代数延伸，不宣称每个插值都是原Wilson配置。

对任意2k列物理来源Z，令P_rev反转列顺序：

$$
C_{b,Z}=\operatorname{Pf}\begin{pmatrix}
N_{0,b}&\Phi_b^{\mathsf T}Z\\-Z^{\mathsf T}\Phi_b&0
\end{pmatrix},\qquad
Z^\theta=R_Y^{\mathsf T}Z^*P_{\rm rev}.
\tag{10}
$$

反射增广矩阵经diag(Rξ,I)合同后，主块为−A。Schur Pfaffian公式及(8)给目标：2k来源块的负号贡献(−1)^k，列反转的行列式也是(−1)^k，两者抵消；主块Pf(−A)=(−1)^n Pf A=Pf A。故

$$
C_{b^\theta,Z^\theta}=\overline{C_{b,Z}},\qquad
I_{0,b^\theta}(\Theta f)=\overline{I_{0,b}(f)}
\quad\text{对全部物理多项式 }f .
\tag{11}
$$

奇数来源恒为零。两侧均为多项式，(9)将恒等式延伸至奇异N₀和零权重，不除以任何配分函数，不固定手征谱秩。sign H在零特征值可取0并保持(2)；原两时间完整积分又可直接复用673的Haar零测度结果。

## 5. 加回原质量

Pφ=diag(p,−p*)按各格点取值；λ实。原质量和反转Grassmann次序的约定给

$$
P_\phi^\theta=-R_Y^{\mathsf T}P_\phi^*R_Y,\qquad
V_\phi(Y)=\tfrac\lambda2Y^{\mathsf T}P_\phi Y,\qquad
\Theta V_\phi=V_{\phi^\theta}.
\tag{12}
$$

e^V是有限Grassmann多项式。将(11)应用到f e^V，得到

$$
N_\lambda=N_0+\lambda\Phi^{\mathsf T}P_\phi\Phi,\qquad
I_{\lambda,b^\theta}(\Theta f)=\overline{I_{\lambda,b}(f)} .
\tag{13}
$$

这包含全部物理来源、实质量插入、奇异质量和同步反射的E／φ多项式系数。不要求逆质量。所证为积分后的物理子代数身份，并非RξᵀNλθRξ=−Nλ*。实际代数样本中，错误裸核身份的最大分量偏差约0.128，裸Φ反射偏差约0.607；正确来源身份仍成立。

## 6. 完整原平均上的Hermitian型

复用673的全部双Haar、S⁹与原H_b配置热核平均dϖτ。674已经核实其在半区交换、时间接口各自取逆下不变；669、673、677保证有限物理来源、有限质量多项式绝对可积。保留全部原Gauss输送；带电列仅是协变来源，不能直接认作孤立可观测量。

$$
\mathfrak I_\lambda(\Theta f)=\overline{\mathfrak I_\lambda(f)},\qquad
Q_\lambda(f,g):=\mathfrak I_\lambda(\Theta f\,g),\qquad
\overline{Q_\lambda(g,f)}=Q_\lambda(f,g).
\tag{14}
$$

最后一步用Θ²=id及Θ(Θg f)=Θf g，没有交换奇Grassmann量。f、g取真实物理正半区代数。因此有限矩阵Q为Hermitian，将674的标量实性提升为完整物理来源合同。但是

$$
Q_\lambda=Q_\lambda^\dagger
\quad\not\Longrightarrow\quad Q_\lambda(f,f)\ge0
\quad\text{或}\quad \mathfrak I_\lambda(1)>0 .
\tag{15}
$$

动态规范正性、非零归一、原H_F和实际仪器身份没有自动完成。同一反射证明在每个有限软调节成立，完整来源极限直接复用677：

$$
a\to0,\ aL\to\infty
\quad\Longrightarrow\quad
\int|C_{a,L,Z}-C_{{\rm original},Z}|\,d\varpi_\tau\to0
\qquad(\text{固定有限盒、固定 }\tau>0).
\tag{16}
$$

这不要求统一Wilson谱隙，也不是新的支配收敛证明。不包含盒细化、τ→0、一般几何来源导数或变维辅助积分的联合极限。

## 7. 结果、复算及下一接口

第一组使用原64维M、质量和16通道，核Hermitian软ε、30／34与31／33谱秩、奇异质量及奇异无质量核；任意ε样本仅验证代数，不冒充实际Wilson背景。第二组回到677原2×2盒，保留全部原群、非平坦空间链路、独立时间链路和变化的E、φ，同时核精确sign与a=.23、L=1软调节。

- 矩阵身份最大分量残差小于1.8×10⁻¹⁵。
- 可分辨样本的0／2／4／6来源反射相对误差小于7.0×10⁻¹⁴。
- 奇异或近零权重只报告绝对残差，不除权重；一般结论由解析式(6)—(13)承担。
- 未数值抽样完整Haar或S⁹积分；完整平均结论来自解析身份与已有可积性。

既有Python／NumPy运行 joint_physical_source_reflection.py --write-results 可首次保存；省略参数只读复算。入口和代数／方程项探针保留在round679_drafts。核验复算678并校验历史证据散列，不重复全部旧实验，不覆盖历史。

本轮合并C01、C14、C16、C17、C19、C20、C22中的物理来源反射条件，没有增加认知公理。原有限H_F过程、指定连续物质、给定作用的经典几何和手征辅助候选四分支仍分别记账。

下一项[680入口](680/drafts/STATUS.md)回查668—669已有半区质量插入，检验其在完整Gauss候选上是否为保物理代数的可逆变换，能否把正性准确约化到原无质量共同泛函。只将新增的完整对象匹配计为进展，不重证一般指数插入公式，不以质量调参、有限Gram样本或第五方向正性替代动态规范正性。
