# 第600轮：完整物质应力与Yukawa力怎样共同进入有效作用约化

日期：2026-10-01。接[598](research_note_598.md)、[599](research_note_599.md)，回查[351](../archive_342_369/research_note_351.md)—[352](../archive_342_369/research_note_352.md)、[581](../archive_554_584/research_note_581.md)—[582](../archive_554_584/research_note_582.md)、[585](research_note_585.md)和[588](research_note_588.md)。[代码](600/joint_full_matter_operator_reduction.py)、[结果](600/joint_full_matter_operator_reduction_results.json)、[核验](600/research_round_600_checks.json)、[条件账](600/unified_physics_condition_ledger_600.md)。三组复算、十二式；主代理审查，无新增独立代理审查。

## 1. 实际拼接缺口与新结果

581用纯标量领先Einstein方程约化曲率项，582把无迹规范应力加入。598之后还存在费米物质：其应力一般不再无迹，且质量直接依赖原标量。因此不能仅把旧公式中的规范应力改叫“总应力”，也不能继续使用没有Yukawa力的标量方程。

**本轮给出完整离壳恒等式：在一阶EFT范围，原已知有效项仍能约化，但必须同时生成由总应力迹、标量质量力决定的交叉项。** 系数由同一原作用及其变换固定，不能分别任意选择或省略。

这接通的是“共同物质—领先几何—有效作用表示”接口，不是认知原则到GR的完成证明。351的限定纯度规结论、352的共同形变条件与多声音锥反例继续保留，本轮不重做它们。

|层次|输入或成果|
|---|---|
|认知动机|整合新物质时，共同几何与其它部门必须一起变换|
|继承输入|原K、U、标量局部系数及598质量函数；四维Einstein领先作用、规范和spin结构|
|范围|一阶局部EFT、有限系数、有界光滑资料、F远离0、已处理边界|
|解析增量|非零应力迹分解、质量耦合的标量力补项、599曲率—质量项的物质接触项|
|数值验证|局部张量恒等式、原U及质量矩阵的真实微扰替换|
|未成立|完整费米外线反项、量子复合算符重整化、固定ℏ图极限、量子GR及物理态|

## 2. 完整领先作用与应力迹

沿581—582的Euclidean正二倍逆度规变分约定。χ统称其余规范和费米场，标量协变动能计在原标量作用中：
 
$$
I_0^{\rm all}=\int\sqrt g\,[-R/2+S/2+U]+I_m[g,\phi,\chi],
\qquad
\tau_{\mu\nu}={2\over\sqrt g}{\delta I_m\over\delta g^{\mu\nu}},
\quad t=g^{\mu\nu}\tau_{\mu\nu},\quad
J_a={1\over\sqrt g}{\delta I_m\over\delta\phi^a}.
\tag{1}
$$

S=K_ab DμφᵃD^μφᵇ，B_{μν}=K_ab DμφᵃDνφᵇ，Q=B:B。费米度规变分在共同vierbein的对称变化下定义，包含spin联络变化；不能只取质量密度当完整Hilbert应力，也不提前用费米在壳方程代替t。

τ、J来自领先I_m，不是重复加入已算的同阶圈修正。Grassmann场的偶局部多项式可作下面的形式代数；不能因此把量子平均的乘积因子化。

令W=τ−tg/2、C=B+Ug+W、r₀=S+4U−t，则
 
$$
E_{\mu\nu}=R_{\mu\nu}-C_{\mu\nu},\qquad e=R-r_0,\qquad
H_{\mu\nu}={1\over\sqrt g}{\delta I_0^{\rm all}\over\delta g^{\mu\nu}}
=-\tfrac12(E_{\mu\nu}-eg_{\mu\nu}/2).
\tag{2}
$$

完整领先度规方程给R=S+4U−t，不能继续用582的R=S+4U。Einstein关系仍来自给定领先作用，不是由原有限图H推出。

## 3. 一般曲率项的离壳分解

对常数a,b及局部系数f（可含S），平方差直接给
 
$$
a\,{\rm Ric}^2+bR^2+fR
=a\,C:C+b r_0^2+f r_0+E:D,\qquad
D_{\mu\nu}=a(R_{\mu\nu}+C_{\mu\nu})
+g_{\mu\nu}[b(R+r_0)+f].
\tag{3}
$$

四维中E=−2H+(tr H)g，所以T=−2D+(tr D)g给E:D=H:T。这是离壳恒等式，提供显式一阶度规重定义。

对580原标量系数，a=1/12、b=1/24、f=−tr A/6−S/9，A=∇²_K U。相对581纯标量曲率约化，新增
 
$$
\boxed{\Delta_\tau=
{B:\tau\over6}+{\tau:\tau\over12}
-{tS\over18}-{Ut\over2}+{t\,\operatorname{tr}A\over6}+{t^2\over24}.}
\tag{4}
$$

证明：W:W=τ:τ、B:W=B:τ−tS/2、tr W=−t。Ric²贡献B:τ/6−tS/12−Ut/6+τ²/12，R²贡献−tS/12−Ut/3+t²/24，线性R项贡献t tr A/6+tS/9，相加即得。

t=0时精确恢复582的B:τ/6+τ²/12。一般情况下，其余四项不能丢；它们不是四条新认知原则，而是完整物质变换固定产生的项。其它圈部门若改变a,b,f，应使用式(3)，不能继续套这里的标量数值。

## 4. 标量质量力必须同时保留

令jᵃ=K^{ab}J_b、𝒯为目标协变张力，M_X为580目标梯度矩阵。完整标量Euler导数及目标链式关系为
 
$$
\mathcal E_\phi=\operatorname{grad}_K U-\mathcal T+j,\qquad
\Box U=\operatorname{tr}(A M_X)
+\langle\operatorname{grad}U,\mathcal T\rangle_K.
\tag{5}
$$

所以
 
$$
-{\operatorname{tr}(A M_X)\over6}
={|\operatorname{grad}U|_K^2\over6}
+{(\operatorname{grad}_K U)^aJ_a\over6}
-{\langle\operatorname{grad}U,\mathcal E_\phi\rangle_K\over6}
-{\Box U\over6}.
\tag{6}
$$

第二项是旧纯标量和常系数规范动能分支没有的质量力耦合。只补应力而不补J仍然不相容。使用581的纯标量约化密度并保留582的规范连接项：
 
$$
b_{\rm red}^{0}={\mathcal E_4\over36}-{7S^2\over216}
+{11Q\over108}+{US\over18}
+U^2-\tfrac23U\operatorname{tr}A+\tfrac12\operatorname{tr}A^2
+{|\operatorname{grad}U|_K^2\over6},
\qquad
b_{\rm red}^{\rm all}=b_{\rm red}^{0}+\Delta_g b_4+\Delta_\tau
+{(\operatorname{grad}_K U)^aJ_a\over6}.
\tag{7}
$$

式(3)取本标量a,b,f后，完整恒等式为
 
$$
b_4^{\rm scalar,g}
=b_{\rm red}^{\rm all}+H:T
-{\langle\operatorname{grad}U,\mathcal E_\phi\rangle_K\over6}
-{\Box U\over6},\qquad
T=-2D+(\operatorname{tr}D)g.
\tag{8}
$$

左侧是已知玻色背景局部算符扇区，不是含所有费米外线的完整一圈作用。把这份已知算符嵌入完整领先物质作用再约化，右侧就产生其它物质相互作用；尚未算出的独立费米外线算符仍开放。

## 5. 共同变量替换与原质量函数

按[Criado–Pérez-Victoria原研究](https://arxiv.org/html/1811.09413)的微扰场重定义方法，采用
 
$$
g_{\rm old}^{\mu\nu}=g_{\rm new}^{\mu\nu}-\lambda T^{\mu\nu},
\qquad
\phi_{\rm old}^a=\phi_{\rm new}^a+
{\lambda\over6}(\operatorname{grad}_K U)^a,\qquad
\chi_{\rm old}=\chi_{\rm new}.
\tag{9}
$$

完整I₀的一阶变化恰消式(8)的Euler项。边界处理后，该有效扇区成为λb_red^all+O(λ²)。λ是有限匹配系数的记账量，不能设为发散极点。F有正余量、有限阶导数有界、变化足够小时有依赖这些界的Taylor余项；不是高阶PDE全局可逆定理。

源、参考、准备、边界及读口须共同拉回。598的质量M_E=N(φ)/√F中，N按实φ线性；同一标量变化给
 
$$
\partial_a M_E={N_a\over\sqrt F}
+{N(\phi)\phi_a\over6F^{3/2}},\qquad
\delta M_E={\lambda\over6}(\operatorname{grad}_K U)^a\partial_aM_E.
\tag{10}
$$

这正是式(7)J项的相应费米双线性系数。Dirac和Majorana都改变；简化标量／引力作用后若仍原样使用旧质量函数，就漏同阶耦合。本轮未设计认知装置。

## 6. 599的曲率—质量项接入总应力

599偶宇称费米UV项含f(φ)R。以全左手对称质量矩阵、Majorana半权计数时，共同b₄中的质量部分给f=−tr(M†M)/6；实际有限Wilson系数仍需匹配。式(2)给
 
$$
fR=f(S+4U-t)+2f\,\operatorname{tr}H.
\tag{11}
$$

完整约化需要−ft这一质量—总应力迹接触项，同一质量矩阵固定其关系。不能将它当无关系数任意删除。

须区分保留显式低能费米场并匹配局部算符，与先积分掉一个部门后只使用所得有效作用。后者不能把同一被积分部门的应力再作为独立经典源加回来。精确高低模划分及量子匹配尚未完成；本身份不许可重复计圈，也不证明〈τ:τ〉=〈τ〉:〈τ〉。量子复合来源还需重整化。

## 7. 守恒和引力约束的边界

已知微分同胚身份在其余完整χ方程成立时，按本应力约定为
 
$$
\nabla^\mu\tau_{\mu\nu}=-J_aD_\nu\phi^a+\mathcal F_\nu,\qquad
\mathcal F_\nu={1\over\sqrt g}{\delta I_{\rm scalar}\over\delta A_\mu^A}
F^A_{\nu\mu},\qquad
\mathcal E_\phi=0\ \Longrightarrow\
\nabla^\mu(T_{\mu\nu}^{\rm scalar}+\tau_{\mu\nu})=0.
\tag{12}
$$

其中𝓕为规范交换力；它在标量应力散度中取相反号，只有零规范曲率或相应交换消失时才可删去。此成熟身份已在303／325等工作使用，不另算新定理。须同时纳入完整规范方程及电流，不能冻结仍受力的规范源后宣布各部门独立守恒。J和𝓕都是部门交换，不能在变换时漏掉。

共同作用确定交换、度规来源及一阶表示变换；守恒仍不等于Einstein方程或局部约束闭合。590正几何H尚未映射为该连续协变作用，585 lapse歧义及588限定菜单障碍继续成立。高阶有效项下，不能不处理相空间／阶数就套351的纯度规二次动量唯一性。

## 8. 复算与接续

三组检查通过：

1. 18组离壳张量，包括无迹极限与完整领先度规关系样本，以580原系数核式(4)、(8)，最大误差2.23×10⁻¹⁶；无迹捷径一般漏非零差。样本不是已制备费米态或全局解。
2. 原574的U、K及598复Y质量矩阵：势梯度的独立差分误差3.46×10⁻¹¹；Dirac一阶系数范数0.0366274，Majorana为0.00236402。实际φ+λgrad U/6替换的二阶余项比为4.00042、4.00021、4.00010。检验的是原质量双线性系数，未求完整费米态。
3. 599同一质量块给f≈−0.20608246，一般应力资料核式(11)，最大误差1.12×10⁻¹⁶，−ft一般非零；不是复合应力重整化或实验拟合。

本轮把“曲率约化”“完整物质应力”“Yukawa标量力”“质量变换”合成一个一阶EFT映射，减少分别拼接的任意性；并未算出全部有限参数。

下一项[601有效约束与物理分支](601/drafts/STATUS.md)：在相同变量和EFT阶上，核高阶项怎样进入lapse、shift与约束传播；复用351—352、585—588，分清微扰低能分支与非扰动高阶理论。认知实现后置，统一目标保持开放。
