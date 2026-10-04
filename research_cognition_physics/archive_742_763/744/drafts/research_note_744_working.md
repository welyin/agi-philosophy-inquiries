# 744工作报告：原相互作用诱导的未知态读口与精确选择规则

日期：2026-10-04。接[743](../../research_note_743.md)、[入口](STATUS.md)。[代码](native_instrument_selection_entry.py)、[结果](native_instrument_selection_entry_results.json)、[核验](entry_checks.json)。正式计数保持743／3424。

## 1. 从两个均值转到完整输入效果

采用743的正常Gauss玻色准备ψ，所有节点的singlet包取717的R=0对称包，链路Haar常数。编码仍用同一个节点的偶sterile空间，不增加寄存器：

$$
|0_L\rangle=\Omega,\qquad |1_L\rangle=e^{i\arg Y_s}a_v^\dagger b_v^\dagger\Omega,
\quad V|j\rangle=\psi\otimes|j_L\rangle,\quad V^\dagger V=I_2.
\tag{1}
$$

准备V是声明的输入，不认为原H已自行执行编码。任意未知编码态可与被动外部参考纠缠。用原完整H的U(t)及577既有末读

$$
E_r=\tfrac12+\tfrac r4\sin s_v,\qquad L_r=\sqrt{E_r},\quad r=\pm1,
\qquad\mathcal I_{r,t}(\rho)=L_rU(t)V\rho V^\dagger U(t)^\dagger L_r.
\tag{2}
$$

输出保留整个原物理系统，不先压回编码。该仪器CP，结果求和保迹，并保原Gauss；这些是等距编码、原U和既有L的直接性质，不算新的普适仪器定理。实际未知态效果为

$$
M_r(t)=V^\dagger U(t)^\dagger E_rU(t)V,\qquad
\tfrac14I\le M_r(t)\le\tfrac34I,\qquad M_+(t)+M_-(t)=I.
\tag{3}
$$

末读仍是操作输入。将式(2)写出，不代表已经完成自主装置、因果空间极限或稳定记忆。

## 2. 原H的联合对称性

令S将全部节点s同时翻号而不改Higgs及链路；N为全部CAR模式的总占据数。取

$$
\Theta=S\exp(i\pi N/2),\qquad\Theta^2=(-1)^N.
\tag{4}
$$

原K及目标测度在S下不变；现场势依s²，测地边只依两端同时不变的内积，全部玻色动能／电磁项也不变。Dirac及跳跃保N且质量系数在S下不变；Majorana的s因S反号，其创建或湮灭二次式因费米相位再反号。因此在原共同核上、继而在同一自伴实现上，

$$
\Theta H_F\Theta^\dagger=H_F,\qquad
\Theta E_+\Theta^\dagger=E_-,\qquad
\Theta V=VZ,\quad Z=\operatorname{diag}(1,-1).
\tag{5}
$$

最后一式使用当前对称ψ，不能扩张到任意准备。597的纯玻色符号对称不直接搬到含Majorana模型；这里用的是明确的联合变换。有限图canonical对称不自动成为相互作用连续测度中的无反常对称，后者未签收。

## 3. 任意等待时间的效果规则

由式(3)—(5)直接得到

$$
ZM_+(t)Z=I-M_+(t),\qquad
M_+(t)=\tfrac12I+b_x(t)X+b_y(t)Y,
\qquad b_x(t)^2+b_y(t)^2\le\tfrac1{16}.
\tag{6}
$$

因此没有Z分量。这是整个原H下的全时间选择规则，不是二阶近似：

$$
\langle0|M_+(t)|0\rangle=\langle1|M_+(t)|1\rangle=\tfrac12
\quad\text{对所有 }t.
\tag{7}
$$

它排除“当前对称准备＋这个固定单次奇读口”直接判别两种编码占据。未排除不同准备、其它效果、时序组合或原完整相互作用的测量能力。即使读数边缘相同，式(2)的完整后态也未被证明相同。

## 4. 有效读取的确存在，且方向由原耦合固定

复用717／743力身份。代码从原完整32模式质量中压缩得到Q_L=|Y_s|X。对称准备使单次效果的标量部分严格为1/2；真实二阶项为

$$
M_+(t)=\tfrac12I-
\frac{|Y_s|}{8w_v}\langle\sqrt F\cos s_v\rangle_\psi\,t^2X+O(t^3).
\tag{8}
$$

有限二维编码和ψ的H幂域保证这里是效果矩阵的算符范数展开，适用于未知编码输入。它不是全物理Hilbert空间上U的范数Taylor展开。原R=0包给w_v=1时X系数−0.05019934720188815；支持保证cos s>0，故解析上非零。对于X本征输入，两种概率差为−0.1003986944037763t²+O(t³)，与743由两个均值得到的旧响应一致。

原12组完整质量、目标度量、现场势与两端测地边的对称检查残差均为0；完整CAR压缩的虚部残差约1.39×10⁻¹⁷。没有数值模拟全玻色传播，式(6)—(7)依赖解析对称证明。

## 5. 为什么下一项必须保留联合历史

对于原多次末读，令K_\boldsymbol r为实际按时序相乘的L和U，则同一对称给

$$
Z\,V^\dagger K_{\boldsymbol r}^\dagger K_{\boldsymbol r}V\,Z
=V^\dagger K_{-\boldsymbol r}^\dagger K_{-\boldsymbol r}V.
\tag{9}
$$

单次规则只约束反号奇的报告。以两次结果乘积r₁r₂为例，它在全部标签反号下为偶，相应效果与Z对易，**允许**包含Z分量；这只是允许，不证明其在原H中确实非零。后续应实际检验这个联合通道，并保留中间读取的反作用，不能用两次未经扰动的Heisenberg均值代替真正instrument历史。

## 6. 当前结论与下一项

已从两个均值推进到未知编码输入的完整CP表达与严格效果选择规则；给定原读口在这一准备下读取配对相干，无法单次读取编码占据。联合历史的占据敏感性、完整后态和共同来源仍待核。复用718、623—625已有历史和域，不新增原子测量权限，不发明新的物种或免费记录器。

空间382—386、425、522—523、604及649／699保持；统一目标未改、未完成。
