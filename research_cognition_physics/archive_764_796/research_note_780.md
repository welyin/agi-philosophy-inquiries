# 第780轮：动态标架的补偿旋转、原旋量字典与分次辅助接触

日期：2026-10-05。接[779](research_note_779.md)及[冻结入口](780/drafts/STATUS.md)。[代码](780/frame_quartet_dictionary.py)、[结果](780/frame_quartet_dictionary_results.json)、[核验](780/research_round_780_checks.json)。

## 1. 本轮增量与保留的边界

**在原背景附近，动态vierbein可以明确分成度规与六个Lorentz取向坐标；原旋量可随之重写，且坐标变换所需的、依赖动态度规的补偿旋转满足精确组合律。** 因而局部Lorentz部门在适配的BV坐标中成为六份非传播四元组。779的局部树接触处方可扩展到这些同时含偶、奇字段的四元组，而不增添物理传播模。

这补上了773“经典六方向左逆”与779“仅偶辅助处方”之间的具体连接。**适配坐标中的共同构造，不自动等于原vierbein坐标下所有量子复合插入及来源已经输送完毕。** 非线性换元的Wick接触、复合源和局部规范化继续交由781；原全N1/N2尚不结项。

|层次|本轮采用或获得的内容|
|---|---|
|认知动机|同一物理关系的内部描述变换，应与共同操作和来源相容|
|已有物理输入|753原四维背景、作用、群、全部物种与参数、F>0；背景spin结构及局部取向|
|分析范围|背景附近的实解析局部片；固定形式阶有限jet；紧支插入；不涉及远离该片的拓扑|
|解析增量|明确的标架截面、动态补偿cocycle、完整经典BV分裂、分次辅助局部接触构造|
|规范化选择|延续779：纯辅助闭环的局部接触取零；传播圈图保留|
|尚缺接口|原非线性变量的量子输送、共同绝对Wick减除及原来源核对；相互作用正态与实际记录另列|

回查[664](../archive_653_701/research_note_664.md)、[735](../archive_702_741/research_note_735.md)、[769](research_note_769.md)、[773](research_note_773.md)与775—779。本轮不重复三维或连续坐标证明，不重新引入384已消去的Lipschitz条件，也不把四维输入改称本轮产物。

## 2. 原标架的局部截面及可逆性

在原背景取coframe E₀，g₀=E₀ᵀηE₀，η=diag(−,+,+,+)。以下矩阵行是内部指标、列是坐标指标。在g₀附近定义

$$
 B(g)=(g_0^{-1}g)^{1/2},\qquad
 \widehat E(g)=E_0B(g),\qquad
 E=\Lambda(q)\widehat E(g),\qquad \Lambda(q)=e^q\in SO^+(1,3).
 \tag{1}
$$

平方根选I附近的实解析支。g₀⁻¹g相对g₀自伴，所以Bᵀg₀=g₀B，继而ÊᵀηÊ=g。任意足够接近E₀且具有同一时间/空间取向的E，都给Λ=EÊ(g)⁻¹∈SO⁺(1,3)；近I的log给六个q。全部是点态可逆变换，没有微分Green逆。

平方根的变分由Sylvester方程确定：

$$
 B\,\delta B+\delta B\,B=g_0^{-1}\delta g.
 \tag{2}
$$

左侧在B=I是2倍恒等算符，故邻域内可逆。各阶微分仍是光滑点态系数乘有限jet，与773的正则局部系数类相容。并不宣称任意Lorentz度规对都存在这一实平方根或全局光滑截面。

已有对接：[Hassan等，Metric Formulation of Ghost-Free Multivielbein Theory，§3式(3.9)—(3.10)](https://arxiv.org/pdf/1204.5202)使用度规平方根与相对标架。这里只复用其局部参数化思路并给邻域证明；没有引入该文的多度规作用或借用其动力学结论。

## 3. 同一旋量与Levi-Civita联络

在I附近选连续spin提升S(Λ)，令Ψ=S(Λ)χ。原664已固定无独立扭率的Levi-Civita路线。本轮仍取ω=ωLC(E)，不是另加一个连接变量。由无扭率Cartan方程及唯一性，

$$
 \omega_{LC}(E)=\Lambda\omega_{LC}(\widehat E)\Lambda^{-1}
                   -d\Lambda\,\Lambda^{-1},\qquad
 \Omega(E)=S\Omega(\widehat E)S^{-1}-dS\,S^{-1}.
 \tag{3}
$$

令γ满足{γᵃ,γᵇ}=2ηᵃᵇ，S⁻¹γᵃS=Λᵃᵦγᵇ。于是D_E(Sχ)=S D_Êχ。Spin变换与手征投影相容，并保持电荷共轭双线性SᵀCS=C；原Yukawa及Majorana内部矩阵未更换。所有原物种的动力学因此可在(g,χ)中表达，q不进入物理作用。内部规范群与Lorentz群的作用相互通勤。

这一步是原作用的经典可逆重写。它不产生664明确排除的独立torsion及额外四费米项，也不更换费米计数。数值程序用一般连接值检验式(3)的代数协变身份；解析论证另由LC唯一性将它限制回原无扭率模型。

## 4. 动态补偿旋转是闭合所必需的

固定背景截面Ê(g)通常不是对所有坐标变换自然的截面。对足够小的f，定义

$$
 f^*\widehat E(g)=L_f(g)\widehat E(f^*g),\qquad
 \chi\longmapsto S(L_f(g))^{-1}f^*\chi.
 \tag{4}
$$

两侧具有同一度规，所以L_f(g)是Lorentz矩阵。沿连通局部片选同一spin提升后，由pullback的结合性直接得到cocycle。常值线性坐标Jacobians的点态版本尤其明确：

$$
 L(J;g)=\widehat E(g)J\widehat E(J^TgJ)^{-1},\qquad
 L(J_1J_2;g)=L(J_1;g)L(J_2;J_1^TgJ_1).
 \tag{5}
$$

对一般f，后一因子必须取已经变换的度规，且前一矩阵随pullback评价。式(5)在代码中作有限变换检验。它保证(g,χ)上的微分同胚作用离壳闭合，毋须假装旋量是没有补偿的普通标量。

无穷小补偿为

$$
 k(c,g)=\left[\mathcal L_c\widehat E(g)
              -\widehat E'_g(\mathcal L_cg)\right]\widehat E(g)^{-1}
       \in\mathfrak{so}(1,3),\qquad
 s\chi=\mathcal L_c\chi-\rho(k(c,g))\chi+s_{\rm int}\chi.
 \tag{6}
$$

这里ρ是spin表示；ghost的符号随原BV约定共同使用。闭合证明来自式(4)的有限群作用，因此包含k随g变动的项，不是只核一个固定背景矩阵的交换子。

**负对照。** 在g₀=η取仅作用于空间1、2方向的A=diag(1,−1)、B₁₂=B₂₁=1。固定背景的反自伴投影k₀(A)=k₀(B)=0，但k₀([A,B])=[A,B]，其12、21分量为2、−2。若把补偿始终冻结在g₀，两个单步旋量作用为零，组合却需非零旋转。缺的正是度规变化引起的补偿变分。这是一个具体捷径失败，不是原标架模型的异常或反证。

## 5. 完整经典BV坐标中的六份规范四元组

在完整经典BRST系统内，不先删除原Lorentz ghost λ；改用

$$
 \vartheta^a=sq^a,\qquad sq^a=\vartheta^a,\quad s\vartheta^a=0,
 \qquad S_{\rm min}=S_{\rm base}^{BV}+\int q_a^*\vartheta^a.
 \tag{7}
$$

a=1,…,6。q处于正则坐标片，λ→sq对λ的Jacobian是点态可逆矩阵；其余部分可含微分同胚ghost及其有限导数。故这是局部可逆ghost重参数化。物理字段(g,χ,…)的BRST变换独立于q、ϑ，因为(4)已经下降到商。反场按整个变换的余切提升、含形式转置与分部积分一起变换；只有改字段却不改反场不能得到式(7)。

每方向添原非最小规范配对(\barϑ,B)，s\barϑ=B、sB=0。用背景密度配对和Ψ_fr=∫\barϑ q选代数规范。采用777括号及sF=(S,F)约定，得到

$$
 S_{\rm fr}=\int\left(q^*\vartheta-\bar\vartheta^*B+Bq
                              -\bar\vartheta\vartheta\right),\qquad
 S=S_{\rm base}^{BV}+\sum_{a=1}^{6}S_{{\rm fr},a}.
 \tag{8}
$$

配对的数值矩阵固定在背景，没有用动态度规给这个纯规范块加权。完整作用满足经典主方程。代码的单方向命名(q,r,c,h,p,t,z,a)=(q,B,ϑ,\barϑ,q*,B*,ϑ*,\barϑ*)，得到sq=c、sh=r、sp=r、st=q−a、sz=p−h、sa=c，其余为零；对321个超单项式精确核s²=0。

反场置零后的块为S_fr,0=qB+ϑ\barϑ。偶Hessian是[[0,1],[1,0]]，奇块亦可逆；都是零阶算符。因此advanced=retarded为点支撑逆，因果差为零，可延伸原W令全部四元组二点为零。它们不添传播自由度，物理自由星代数和原自由正态保持。**BV作用、时间序接触仍保留四元组，不能先令所有规范变量为零。**

## 6. 779到偶奇辅助块的扩展

令x为六份四元组字段，S_fr,0为上述固定二次作用。反场当作不收缩的中心插入参数，传播字段及原16辅助仍用778—779共同族。任意紧支偶局部V=O(t)，以左变分定义

$$
 \frac{\delta^L}{\delta u}\left[S_{\rm fr,0}(u)+V(x+u)\right]=0,
 \qquad
 \mathcal R_{\rm gr}(V;x)=S_{\rm fr,0}(u)+V(x+u),\qquad
 \mathscr S_{\rm adapted}(V)=\mathscr S_{P+e}(\mathcal R_{\rm gr}(V)).
 \tag{9}
$$

此处\mathscr S包括指数的i/ℏ因子，R本身是相互作用泛函。把不同奇插入用外部Grassmann源配成偶V，再取源系数，便定义完整分次菜单。用一份四元组(q,r,c,h)写递推为u_r=−V_q、u_q=−V_r、u_h=−V_c、u_c=+V_h，所有右侧在x+u处计算。最后一个正号来自左导数∂_h(ch)=−c，不能套用偶对称逆矩阵。

每固定t阶只用较低阶u、固定点态逆、乘法和有限微分，所以仍是局部有限jet。它生成有Koszul符号的辅助树；纯辅助闭环选局部规范化为零。树收缩成顶点后，**原传播圈图由同一\mathscr S_{P+e}继续处理**。这不是未经减除的有限维Berezin/Gaussian测度，也不宣称所有物理圈图消失。

相容性的证明与779结构相同，但变分必须分次：

$$
 \mathcal R_{\rm gr}(V)=V+O(t^2),\qquad
 \mathcal R'_{\rm gr,V}[W]=W(x+u),\qquad
 \delta_\zeta\mathcal R_{\rm gr}(V)=
          \mathcal R'_{\rm gr,V}[\delta_\zeta V].
 \tag{10}
$$

驻值方程抵消δu；奇方向采用同一左变分及源的次序。这给字段独立性和所有基本场插入身份。点态局部性给不相交扰动的Hammerstein恒等式，继承因果分解；实结构相容的二次块与实系数形式解继承实性。导数插入按同一Euler变分与分部积分延伸；反场纯插入继续因不收缩而因子化。R=V+O(t²)保持T₁。局部树不会制造新的传播波前集；固定阶导数有限，延拓类别仍是原有限标度/有限jet类别。

单方向精确收缩检查给

$$
 C=\partial_q\partial_r-\partial_h\partial_c,\qquad
 [s_0,C]=-\Delta_{BV},\qquad
 \Delta_{BV}=\partial_q\partial_p+\partial_r\partial_t
                -\partial_c\partial_z-\partial_h\partial_a.
 \tag{11}
$$

这校准偶奇接触的同一BV约定，不把有限收缩算符当成连续场论中的δ(0)。在适配变量域，原传播族和上述分次扩展可使用[Fröb，1803.10235，定理3、6—9](https://arxiv.org/pdf/1803.10235)的局部异常框架。它给一致性问题的合法共同域，并不直接给“实际反常为零”或原坐标系中的量子换元定理。

对实际不含四元组的物理作用及关系观测，u=0、R_gr=V。因而加上这一纯规范表示，不改变它们在适配变量中已定义的量子值。若还要求原E、Ψ、λ及其复合来源的同一绝对规范化，必须继续下一节。

## 7. 量子层面仍不能略过的具体项

经典余切提升保BV括号，不自动保固定背景的普通Wick乘积和时间序字段独立性。即使q的W为零，Ê(g)也是度规的非线性复合场。在一点取E₀=I、g₀=η并令g=η+h，

$$
 \widehat E(g)=I+\tfrac12\eta h-\tfrac18(\eta h)^2+O(h^3),\qquad
 \Delta\langle\widehat E\rangle_{\rm contact}
       =-\frac{\hbar}{8}\,\big\langle(\eta h)^2\big\rangle_{\Delta W}
         +O(\hbar^2).
 \tag{12}
$$

右侧用769的两Hadamard态之间**光滑差ΔW**说明状态变化接触，不能代入未减除的W(x,x)。绝对接触须用原770的同一Wick处方。这里复用769的普遍链式法则，式(12)只是把已知二阶均值修正落实到本轮具体Ê映射，不另计发现一套换帧理论。

此外，原线性源J·E应变成对Λ(q)Ê(g)的复合源；仅把J改名为新变量的线性源会换掉可观测量。[Criado—Pérez-Victoria，1811.09413，§2](https://arxiv.org/html/1811.09413)明确区分换元的Jacobian和源变换；其在特定正则化下忽略局部Jacobian的做法，不能无条件移植到本项目固定Hadamard减除。

因此780证明的是**完整经典字典及适配坐标的分次局部量子构造**。原共同来源、原全N1/N2的最终签收，仍须核非线性复合映射与共同规范化。这是已定位到具体映射的接口，不是再次要求重建传播T族或重算16个偶辅助。

## 8. 可复算结果与下一步

沿用现有Python/NumPy，运行 `780/frame_quartet_dictionary.py --check`；证据入口为 `780/verify_round780.py`。四组诊断：

1. 24组局部矩阵检验：截面、Lorentz重构、有限补偿组合最大残差低于1.4×10⁻¹⁵；Sylvester导数与中心差分差约4.0×10⁻¹¹。冻结补偿的最小有限组合缺陷约3.90×10⁻⁴，并给精确A、B反例。
2. 16组spin/Dirac检验：协变残差低于3.3×10⁻¹⁵；漏掉连接非齐次项的最小缺陷约0.372；保手征和Majorana双线性。
3. 321个超单项式精确核四元组CME、nilpotency及式(11)。这是新四元组的符号校准，不算321个物理定理。
4. 四阶分次驻值展开、50项包络/驻值身份全为有理数精确相等；去掉驻值二次项出现−qr缺陷。纯旁观物理作用保持。

有限检验不能代替上述连续局部证明；未计算原圈反常或建立相互作用正态。上一轮3522项累计分组检查加本轮4组为3526。

下一轮[781入口](781/drafts/STATUS.md)：以本轮显式非线性映射为对象，将源、Wick接触、BV插入与原770—775来源放进同一量子输送；先核真实需要的条件和已可复用的769结论。若只得到局部形式或一圈结果，就按其量词记录，不扩大成全局、UV或全物理完成。
