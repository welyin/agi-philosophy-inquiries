# 第824轮：原标量参考的量子均值与逻辑输入产生的力

日期：2026-10-05。接[823](research_note_823.md)、[冻结入口](824/drafts/STATUS.md)及[工作稿](824/drafts/research_note_824_working.md)。

[完整原质量代码](824/scalar_force_bridge.py) · [结果](824/scalar_force_bridge_results.json) · [早期参考jet校准](824/reference_normal_jet_probe.py) · [结果](824/reference_normal_jet_probe_results.json) · [核验](824/research_round_824_checks.json) · [复算](824/verify_round824.py)。

## 1. 联合结论和未越过的边界

**822保留的原纯玻色参考资料，在当前输入族的首圈相对均值中没有额外收缩差；但参考的法向变化率仍受完整标量力约束。原场空间逆度量与全部质量导数相乘，恰好消去混合项，使singlet加速度的费米来源只取原Majorana质量分子N_s。**

在816—819已经允许的存在性家族内，小幅改变原sterile模式的粒子—孔极化，可使所选Σ_*上的这份标量力非零，同时保留旧能量非零及模式外满秩、逻辑噪声见证。没有新增物种、耦合或独立来源。

**本轮尚未把Σ_*上的非零性搬到822使用的原初片Σ₀。** 未声称旧纯粒子见证本来就有非零N_s均值；该均值在所选片实际上为零。也没有证明任意物质补偿失败，或统一模型不成立。

|层次|本轮处理|
|---|---|
|认知动机|记录内容、钟尺和来源由同一物质承担，不能要求它们彼此独立变化|
|继承输入|原四维、群谱、作用、753背景、原W和有限形式准备合同|
|复用|573参考与目标度量；769/796复合均值；802—803完整Dirac密度；816—819允许模式族；822受限补偿|
|解析增量|纯玻色收缩差消去；原加速度—Majorana算符字典；同一模式族兼有实际标量力和旧噪声；关系度规的法向剩余项|
|数值角色|原64维矩阵和原15³背景上的公式校准；不替代原连续存在证明|
|未签收|Σ₀实际标量来源非零、全部初片共同准备、π_g完整字典、有限耦合和自主装置|

上一目标轮完成821—823并保存本轮工作稿，属于有效进展。本轮先读导航、旧报告及结果，普通进程查询未见运行中的Python；保留全部冻结文件，无新应用任务及图像检查。

## 2. 原参考的收缩差确实消去

对原纯玻色局部jet函数Y(g,φ)，使用同一原单插入规范化及ε=√ℏ计数：

$$
\langle Y\rangle_\rho=Y_*+\epsilon^2
\left(Y'_B\mu_{B,\rho}+\tfrac12\operatorname{Contr}_{W_B-H_B}Y''_{BB}
+Y^{(1)}_*\right)+O(\epsilon^3),\qquad
\Delta k_Y=0 .
\tag{1}
$$

原因是当前818族和821/823首阶准备保持领先W_B；Y无经典费米腿，Y''_FF=0。原减除与量子单插入项相同，该阶后者只取背景值。比较的是同一处方中的差，没有把发散重合值当有限数。

这适用于h、s、r_h=(∇h)²、r_s=(∇s)²及有限导数，也适用于用这些参考构造的纯玻色关系字段。Higgs模长h的二阶角向项并未删除；只是其收缩在输入差中相消。不能把式(1)推广给任意显含费米双线性的canonical动量。

因此822的目标初始h/s一阶jet及四参考空间jet，在量子首圈差中没有额外k障碍。但这是目标资料的解释；823仍只给每份有限物理菜单的一份共同准备，未给整个初片的状态实现。原X在指定坐标的值也不是单独的微分同胚不变观测，物理比较须使用下面的关系字典。

## 3. 参考Jacobian的行列式不决定整个关系度规

记X=(h,s,r_h,r_s)，α_B=X_B⁻¹∘X_*，ξ=(DX_*)⁻¹δX。原关系度规的一阶变化为

$$
\delta g^{\rm rel}=\delta g-\mathcal L_\xi g,
\qquad \delta X|_\Sigma=0\ \Longrightarrow\ 
\xi|_\Sigma=0,\quad\partial_i\xi|_\Sigma=0 .
\tag{2}
$$

在共同Gaussian初始规范，822固定几何使δg_μν|Σ=0。令u为参考位移ξ的法向变化率，精确分量式为

$$
u^\alpha=(DX_*)^{-1\,\alpha A}\delta(\partial_nX^A),\qquad
\delta g^{\rm rel}_{ij}=0,\quad
\delta g^{\rm rel}_{ni}=-\gamma_{ij}u^j,\quad
\delta g^{\rm rel}_{nn}=2u^n .
\tag{3}
$$

式(1)使纯玻色关系度规的相对收缩也消去，但不能消去式(3)。特别在573原证明点，dh=(v_h,0,0,0)、ds=(v_s,0,s_y,0)，v_h s_y≠0。已保δ∂_nh=δ∂_ns=0，故u^n=u^y=0，而

$$
\begin{pmatrix}(r_h)_x&(r_h)_z\\(r_s)_x&(r_s)_z\end{pmatrix}
\begin{pmatrix}u^x\\u^z\end{pmatrix}
=\begin{pmatrix}\delta\partial_n r_h\\\delta\partial_n r_s\end{pmatrix} .
\tag{4}
$$

2×2块的可逆性复用旧连续秩证书。背景r_h/r_s的时间导数并未取零；它们在此乘u^n而消去。混合度规分量可变，原Jacobian行列式仍可不变。这是同一材料参考中的关系变化，不等于已经生成新的曲率动力学。

## 4. 原标量加速度的准确保持条件

在任一声明的Cauchy片采用共同Gaussian和内部时间规范。固定初始g/K、全部物质配置与空间jet，并保持v_h、v_s；Higgs方向取原局部内部规范。此时δv只在Higgs角向，φ·δv=0。原联络Γ^A_BC=(δ^A_Bφ_C+δ^A_Cφ_B)/(6F)，使其径向变化为零。

写a为标量的法向加速度。把完整相对方程的来源按原作用体积约定定义为K_AB δa^B+δ(其余经典项)_A+Δj_A=0，令f^A=−K^{AB}Δj_B，则

$$
\delta a_h=\frac{2}{h}v_\perp\cdot\delta v_\perp+f_h,
\qquad\delta a_s=f_s,
\qquad\delta\partial_n r_I=-2v_I\delta a_I\quad(I=h,s).
\tag{5}
$$

第一项来自模长h的二次链式法则；不能把h当成单个线性分量。空间/势项、初始K的阻尼和径向速度已固定。原规范动能系数常数，规范电动量不直接进入此刻的标量加速度；颜色不作用于Higgs或singlet。这不删除它们对后续演化的影响。

所以在这类保持合同下，若v_s f_s≠0，任何只调规范共轭动量和Higgs角向速度的方案，都不能同时保r_s的法向导数。改变几何、配置、径向速度、输入来源或扩大保持合同，属于其他候选，未排除。

## 5. 全部标量力收缩成原质量分子

原802 Einstein字段中的费米质量矩阵及573逆目标度量为

$$
\mathsf B(\phi)=\frac{\sum_A N_A\phi^A}{\sqrt F},\qquad
F=2-\phi^2/6,\qquad
\mathcal K^{AB}=F(\delta^{AB}-\phi^A\phi^B/12),\qquad
\mathcal K^{AB}\partial_B\mathsf B=\sqrt F\,N_A .
\tag{6}
$$

证明：质量导数为N_B/√F+(N·φ)φ_B/(6F^{3/2})。第一项收缩产生−√F φ^A(N·φ)/12；第二项用1−φ²/12=F/2产生相反项，恰好抵消。没有只保singlet质量、没有删原Higgs/Yukawa矩阵。

固定g和802旋量平凡化时，原费米标量依赖只在B；spin、密度及内部连接无φ导数。803的半密度变换T也只依赖g，因此此标量变分没有端点T项。原作用中的Nambu权1/2和质量负号给

$$
\Delta j_A=\tfrac12\Delta\langle\Psi^\dagger
(\partial_A\mathsf B)\Psi\rangle,
\qquad f_s=-\frac{\sqrt F}{2}\Delta\langle\Psi^\dagger N_s\Psi\rangle,
\qquad N_s=N_5 .
\tag{7}
$$

该j符号按第4节的加速度方程定义；作用导数本身带相反号。状态无关局部减除在相对差中相消。其余完整应力和规范来源仍由同一个ΔP产生，不允许把式(7)当独立外加力。

在法向标架，Ξ=γ^{1/4}Ψ为原Cauchy半密度。以Q(K)=Ξ†KΞ/2，实空间测试η给

$$
\int_\Sigma d^3x\,\eta\frac{\sqrt\gamma}{\sqrt F}f_s
=-\Delta\omega Q(\eta N_s)
=\tfrac12\operatorname{Tr}(\Delta P\,\eta N_s) .
\tag{8}
$$

这是原曲背景上的来源公式，不是辅助平直Hamiltonian的质量拟合。半密度、原F及Nambu计数已进入同一字典。

## 6. 合法的小旋转给实际标量力

819的纯sterile粒子极化s满足Γ(n)s=ςs；在该片〈s,N_s s〉=0，不能用其非零质量作用范数冒充来源均值。原非零Majorana参数给N_s²s=ν²s、ν>0，且N_s与Γ(n)反对易。定义

$$
t_\theta=\cos\theta\,s+\sin\theta\,N_ss/\nu,\qquad
\langle Ct_\theta,t_\theta\rangle=0,\quad\|t_\theta\|=1,
\qquad
\langle t_\theta,N_st_\theta\rangle=\nu\sin2\theta,\quad
\langle t_\theta,\Gamma(n)t_\theta\rangle=\varsigma\cos2\theta .
\tag{9}
$$

CAR合法性由原C N_s C=−N_s、自伴性和sterile Majorana配对的反对称性给出；N_s s是原已有孔极化，不是新物种。可在原小片平滑选s与t；原质量分子ν为固定参数。两不交包分别以t_θ替换s，仍得到正交C不变四维空间和相同形式的两模偶逻辑字典。

按819的Z⁺/Z⁻约定，ΔP=D_z=|f₁〉〈f₁|−|Cf₁〉〈Cf₁|。取η≥0支持第一包的非零区，式(8)给

$$
\int_{\Sigma_*}\eta\frac{\sqrt\gamma}{\sqrt F}f_{s,Z^+-Z^-}
=\langle f_1^\theta,\eta N_sf_1^\theta\rangle
=\nu\sin2\theta\int\eta|\chi_1|^2\ne0
\quad(0<|\theta|\ll1).
\tag{10}
$$

这次来源来自同一实际允许输入；没有把早期754诊断当成量子来源。X/Y在所选片的局部二次来源仍因包支持不交而为零。

## 7. 与旧噪声见证相容，且不移动初片

先固定816—819已取得的有限κ、时间窗及公共过程。模式依θ在所有所需光滑范数中连续；原一阶响应算符作用于这些固定有限光滑模式也连续。因旧模式外Gram严格正、完整能量对比非零，存在θ₀>0使

$$
|\theta|<\theta_0\quad\Longrightarrow\quad
G_\theta>0,\quad\Delta q_{{\rm energy},\theta}\ne0,
\qquad0<|\theta|<\theta_0\quad\Longrightarrow\quad
f_{s,\theta}\not\equiv0\ \text{于所选}\ \Sigma_* .
\tag{11}
$$

817的有限T和偶准备、818共同输入和噪声论证针对这份合法S_θ重用。因此同一过程族可同时具有原逻辑噪声、能量来源和singlet力；未重新选择玻色W或作用。

这里靠有限维严格正的开放性，不拿数值采样确定原θ₀，不交换κ→∞与θ→0。原先θ=0的冻结见证保持有效。

Σ_*依旧是816选定的片。822的显式补偿只在原Σ₀证明过；两片间虽有原Dirac传播，却不能推断非零局部双线性在另一片也非零。因此本轮只给Σ_*实际力以及任意同片保持合同的必要条件，不声称原Σ₀的完整关系度规已被迫改变。

## 8. 核验及下一项

一组联合检查，含先前未正式计数的工作probe，累计3602：

- 12组原五标量配置、全部64维质量矩阵验证式(6)，最大误差4.75×10⁻¹⁶以内；漏逆度量的秩一项产生约0.0651944误差。
- 原ν≈0.32280025；θ=0.02时〈N_s〉≈0.0129086，能量主极化≈0.9992001，CAR各向同性误差为0。θ=0的〈N_s〉确实为0。
- 原Clifford/质量系数的有限空间Gram在所检θ下保持正；这只校准连续开放性，原实际G_θ的证明沿第7节。
- 早期原15³背景的诊断保δv_h=δv_s=0，却有最大δ∂_nr_h≈0.0503300；原证明点谱插值的关系度规混合项非零。该来源仍是明确标注的754诊断，不用来证明式(10)。

所有早期结果中的working标志保留，正式签收由本报告和核验收据决定。未数值计算原曲时空传播或全部均值；未计算有限强度装置。

接[825入口](825/drafts/STATUS.md)：检验能否在**原Σ₀**上给同一合法模式家族的非零singlet力，同时保旧过程的严格噪声与能量见证。应利用完整初片的原自由资料及有限模式开放性，不能把Σ_*改名为Σ₀。全统一目标和独立输入保持；本轮仅排除满足式(5)非零条件时的过强保持合同：

$$
v_s f_s\ne0\quad\Longrightarrow\quad
\text{同片固定几何、配置及径向一阶资料，不能再固定}\ \partial_n r_s .
\tag{12}
$$

