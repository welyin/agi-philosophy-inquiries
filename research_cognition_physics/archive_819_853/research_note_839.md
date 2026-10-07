# 第839轮：同一原区域中的记录熵—参考偏移精确分解

日期：2026-10-05。接[838](research_note_838.md)与[冻结入口](839/drafts/STATUS.md)。[工作稿](839/drafts/research_note_839_working.md)保留冻结，有限层结论在此接到原自由区域。

配套：[有限参考校准](839/gaussian_reference_entropy_split.py) · [原结果](839/gaussian_reference_entropy_split_results.json) · [嵌套校准](839/nested_region_entropy_split.py) · [新结果](839/nested_region_entropy_split_results.json) · [核验](839/research_round_839_checks.json) · [复算](839/verify_round839.py)。

## 1. 结论及共同条件

**对829—838指定的独立余部组织准备，同一个原自由费米区域中的Araki相对熵，可精确分为有限的非Gaussian记录项与共同的Gaussian参考偏移项。** 总相对熵可以为无穷；记录内容之间的相对熵始终有限。整族参考相对熵的有限性归结为同一份Gaussian偏移的有限性。

|层次|本轮地位|
|---|---|
|认知观察与候选性质|内容和共同工作条件可以分开比较；二者须在同一材料与参考下核算|
|替代解释|这是指定编码/独立余部准备的性质，不证明所有认知必须采用此准备|
|数学要求|同一区域的正常态、实际嵌套CAR限制、共同Gaussian参考；支持与极限共同检查|
|继承输入|730原固定背景自由参考、829组织族、837局部支持、838绝对二点偏移|
|解析增量|有限公式到原区域偶代数的映射；消去所有有限参考必须忠实的额外要求；统一有限性条件|
|开放|Gaussian偏移有限性、应力/模流识别、完整相互作用、面积响应与跨尺度|

对象是原自由物理费米子代数，不将其叫作含动态引力的全部区域代数。四维、群谱、作用及参数仍为输入。

## 2. 文献与历史范围

复用305的相对熵余项，830、837的有限码熵亏，以及637—645的固定图热参考与共同粗化。不把固定图能源截断搬作当前区域模式极限。

核对[Araki 1977，定理3.8(2)、3.9(2)、§5](https://www.jstage.jst.go.jp/article/kyotoms1969/13/1/13_1_173/_pdf/-char/en)：增加的有限维子代数生成目标von Neumann代数时，正常态相对熵按限制收敛，包含非忠实态。这里只用有限维子代数情形，不扩大早期定理范围。[Xu 2022，§II.A](https://doi.org/10.1063/5.0067599)也列出这一收敛性质。本轮新增是验明原模式和态满足其前提。

## 3. 同一区域内的嵌套对象

在837原正则Cauchy片选开集B，包含全部十模式的紧支持，其许可因果发展O留在原背景范围。原自由时间片性质把B内Cauchy资料与O内场代数对应。H_B为B内资料的Hilbert闭包，C为自对偶共轭，Q为20维C不变编码投影。

从B内光滑紧支C实资料的稠密列先减Q分量，再Gram正交化，跳过零向量，每次取两个实模式。减去的Q模式仍紧支于B，故全部有限步骤保支持：
$$
H_n=\operatorname{ran}Q\oplus
 \operatorname{span}_{\mathbb C}\{e_1,\ldots,e_{2n}\},
\qquad H_n\subset H_{n+1},\quad
\overline{\bigcup_nH_n}=H_B,\quad CH_n=H_n .
\tag{1}
$$

没有每级重选真空、背景或边界。CAR范数连续性保证有限模式多项式的并稠密。π为原自由表示，取偶部分：
$$
\mathfrak N_n=\pi(\operatorname{CAR}(H_n)_{\rm ev})'',
\qquad
\left(\bigcup_n\mathfrak N_n\right)''
=\mathfrak M_B:=\pi(\operatorname{CAR}(H_B)_{\rm ev})'' .
\tag{2}
$$

有限偶代数是两个宇称矩阵块的直和，嵌入保同一单位。计算相对熵必须保块概率，不能分别归一后漏掉中心的经典项。

## 4. 同一参考、两份准备及正常性

ω_bg为730原准自由参考，ν为其编码正交余部限制。σ=σ_{ε,r}沿837，0<ε<1、|r|<1，d=1024。定义偶分级乘积：
$$
\omega_r=\sigma_{\epsilon,r}\,\widehat\otimes\,\nu,\qquad
\omega_g=(I_d/d)\,\widehat\otimes\,\nu .
\tag{3}
$$

ν准自由，故ω_g是ω_r的同二点Gaussian态，协方差正是838的Q/2+(I−Q)P(I−Q)。原ω_bg的交叉相关保留，没有换成乘积参考。

正常性不依赖假设连续区域有密度矩阵。原表示的有限CAR因子为M_d；在其矩阵单位分解中，用有限个Kraus算符重置该因子、余部取原限制，得到正常态。偶准备与原偶态保证分级字典一致；限制到区域偶代数仍正常。这是829已有数学准备，不是自治控制证明。

若λmin、λmax为σ的极端本征值：
$$
\alpha\omega_g\le\omega_r\le\beta\omega_g,\qquad
\alpha=d\lambda_{\min}>0,\quad \beta=d\lambda_{\max}<\infty .
\tag{4}
$$

所以二者在所有有限限制和目标区域有同一支持。原参考可以非忠实，不添底噪改变它。

## 5. 有限分解及非忠实支持

γ_n、ρ_{r,n}、g_n分别为三态在H_n的有限CAR密度。它们均偶，完整有限CAR的迹相对熵等于偶代数直和相对熵。忠实Gaussian参考的log为常数加二次Majorana形式，同二点性给
$$
\operatorname{Tr}[(\rho_{r,n}-g_n)\log\gamma_n]=0,\qquad
\operatorname{Tr}[(\rho_{r,n}-g_n)\log g_n]=0 .
\tag{5}
$$

展开定义，复用标准Gaussian最大熵的迹论证：
$$
D(\rho_{r,n}\Vert\gamma_n)
=D(\rho_{r,n}\Vert g_n)+D(g_n\Vert\gamma_n),\qquad
D(\rho_{r,n}\Vert g_n)
=c_{\epsilon,r}:=\log d-S(\sigma_{\epsilon,r}) .
\tag{6}
$$

ν_n可非忠实：两准备均支持在完整编码空间乘suppν_n，记录项仍是c。若此支持不包含于suppγ_n，则两项与γ_n比较的相对熵同为∞，式(6)按扩展值成立。

若包含，将Gaussian γ_n作费米正规形。其纯占据模式在两准备中被同样固定：正占据算符期望0或1强制密度支持于该本征空间。去掉固定模式后，logγ_n在共同支持上由其余有限二次项表示，式(5)—(6)仍成立。g_n纯模式同理。因而无需要求所有有限限制忠实，且不作无穷相减。

## 6. 原区域精确分账

对式(1)—(2)应用增加子代数收敛，两项与原参考比较的相对熵分别趋向其区域值；c与n无关。于是
$$
D_{\mathfrak M_B}(\omega_r\Vert\omega_{\rm bg})
=c_{\epsilon,r}
+ D_{\mathfrak M_B}(\omega_g\Vert\omega_{\rm bg})
=:c_{\epsilon,r}+b_B .
\tag{7}
$$

同一极限还给
$$
D_{\mathfrak M_B}(\omega_r\Vert\omega_g)=c_{\epsilon,r}<\infty .
\tag{8}
$$

这里是Araki相对熵，未给连续区域配普通有限熵或密度矩阵。ε=.1、r=.6时c≈5.463074974 nat，但b_B尚未由原背景计算。整族有限性归结为一个共同条件：
$$
D_{\mathfrak M_B}(\omega_r\Vert\omega_{\rm bg})<\infty
\quad\Longleftrightarrow\quad b_B<\infty .
\tag{9}
$$

这减少了逐输入选择参考/验证有限性的自由度，没有证明b_B必定有限。有限秩光滑二点差和Hadamard性质，尚未在原背景上提供所需模能量界。

## 7. 内容比较与837的强化范围

同一ε的组织态共享余部，有限相对熵完全位于编码因子。取极限：
$$
D_{\mathfrak M_B}(\omega_r\Vert\omega_{r_0})
=(1-\epsilon)D(\tau_r\Vert\tau_{r_0})<\infty .
\tag{10}
$$

837对任意正常扩展只给下界；这里对共同独立余部扩展、且区域包含全部编码模式给等号，没有扩成任意扩展结论。

只有b_B有限时，式(7)才允许相减：
$$
D_{\mathfrak M_B}(\omega_r\Vert\omega_{\rm bg})
-D_{\mathfrak M_B}(\omega_{r_0}\Vert\omega_{\rm bg})
=S(\sigma_{\epsilon,r_0})-S(\sigma_{\epsilon,r}) .
\tag{11}
$$

b_B=∞时，式(10)仍比较内容，但式(11)不能用作∞−∞的定义。改用其它参考、相对扰动或更柔和准备是另一个待检验分支。

## 8. 物理边界和负对照

每个有限层的内容变化保原Gaussian二次模能量，不证明原区域模流由局部应力生成。若完整参考模算符另有高阶项K_ng，则有限比较满足
$$
\Delta D=-\Delta S+\Delta\langle K_{\rm ng}\rangle
\quad\hbox{在共同二次模响应为零时}.
\tag{12}
$$

冻结校准加入六次项.31 Z₁Z₂Z₃，κ=.6时留下.186，限定Gaussian分解范围。不是否定所有相互作用；其高阶模项须和同阶来源共同处理。

本轮没有给形式级数取实体对数，没有完成有限耦合正常网、全引力关系区域或面积响应。自由物质、记录和参考已有同一区域比较；自主准备、图连续及跨尺度仍缺。

## 9. 复算与下一项

两组首次正式签收（含冻结工作校准），累计3623：

- 四模式配对/相关参考恒等式残差小于9.5×10⁻¹⁶，保非Gaussian负对照。
- 同一个六模式态对限制到3、4、5、6模式，相对熵约.303038、.338392、.397312、.496197，均是固定记录项.192744757加逐级增加的共同Gaussian项；宇称直和保中心权重。
- 非忠实且支持相容时身份仍成立；支持不相容时两项同为∞，没有修改原参考。
- 另列趋近支持边界的参考例，说明有限模式数不给统一上界；这些参考变化不是原区域的嵌套序列。

原连续映射和极限由解析论证承担，数值不冒称原PDE解。接[840入口](840/drafts/STATUS.md)，把b_B有限性接到同一原区域协方差的模谱条件，不反复扫描小矩阵，也不将新条件直接命名为认知公理。
