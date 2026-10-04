# 第320轮：Wald熵、Einstein面积与仿射参数的共同变换

日期：2026-09-23。接续[319轮](research_note_319.md)；本轮9项检查。

## 1. 假设与范围

认知动机是同一物理过程的不同数学描述应保持可比较量。本轮检验：319轮Wald加权面积的增加，是否就是另一套变量中的普通面积增加？

对接[Jacobson、Kang与Myers的共形映射／熵方法](https://arxiv.org/abs/gr-qc/9503020)及[Ford与Roman的非最小标量分析](https://arxiv.org/abs/gr-qc/0009076)。使用标准经典场重定义，复算的是同一平面波分支的光学与能量接口。

输入仍为319轮的作用量、场方程、F＞0光滑片、κ_g＝8πG＞0与未来端点。没有据此断言不同变量的完整量子测度、区域代数、外部物质或全局边界已自动等价。

## 2. 面积与类光生成元必须一起变换

在四维，令

$$
g^E_{ab}=F g^J_{ab},\qquad A_E=F A_J,
\qquad S_W^J=\frac{F A_J}{4G}=\frac{A_E}{4G}.
\tag{1}
$$

上标J与E分别指原变量和Einstein变量。光锥不因正共形因子改变，但原仿射参数通常不再是新度规的仿射参数：

$$
d\lambda_E=F\,d\lambda_J,\qquad
k_E^a=F^{-1}k_J^a,\qquad
\theta_E=\frac1{A_E}\frac{dA_E}{d\lambda_E}=\frac{\Theta_J}{F},
\qquad \sigma_E^2=\frac{\sigma_J^2}{F^2}.
\tag{2}
$$

常数归一化可另选，但λ依赖的F不能省略。代码从g_E及其导数构造联络，直接核对k_E的测地线加速度；不是先把变换式当作待测结果再相互代入。

## 3. 标量动能也改变

在作用量层面，忽略已明确分离的全导数，并在F＞0片上重定义标量ψ：

$$
I_E=\int\sqrt{-g_E}\left[\frac{R_E}{2\kappa_g}
-\frac12 Z_E(\phi)(\nabla_E\phi)^2\right],
\qquad Z_E=\frac1F+\frac{3}{2\kappa_g}\left(\frac{F_{,\phi}}F\right)^2,
\qquad \frac{d\psi}{d\phi}=\sqrt{Z_E}.
\tag{3}
$$

这是无势能、单实标量分支的标准变换。不能只把g乘F，却继续使用原来的动能系数；有其他物质时也必须一起变换其作用量。

由式(2)—(3)，期望核验的Einstein框架方程是

$$
R^E_{ab}k_E^ak_E^b=\kappa_g\left(\frac{d\psi}{d\lambda_E}\right)^2,
\qquad
\frac{d\theta_E}{d\lambda_E}
=-\frac{\theta_E^2}{2}-\sigma_E^2
-\kappa_g\left(\frac{d\psi}{d\lambda_E}\right)^2.
\tag{4}
$$

这些是同一已输入动力学的另一种表示，不是第二次独立推导出引力。

## 4. 独立的数值核对

复用319轮的同一(r,φ,b)解。程序计算完整g_E、g_E逆矩阵和联络，用五点导数与联络收缩计算R^E_λλ，再除以F²与标量能量比较。

512步结果：

| 检验量 | 最大误差／失配 |
|---|---:|
| 从新截面度规计算A_E，与F A_J比较 | 3.33×10⁻¹⁶ |
| 正确k_E的仿射测地线加速度 | 4.44×10⁻¹⁶ |
| 错把原λ继续当新仿射参数 | 1.4199558004 |
| 独立联络曲率与标量源 | 1.12×10⁻⁶ |
| Einstein框架Raychaudhuri方程 | 1.19×10⁻⁶ |
| 数值积分ψ(φ)后差分，与√Z_E φ′比较 | 1.80×10⁻⁸ |

独立曲率残差在128、256、512步依次为2.080506×10⁻⁴、1.607328×10⁻⁵、1.121051×10⁻⁶，随细化趋向四阶收敛。不能把面积的机器精度代数匹配当成整套微分动力学也已达到10⁻¹⁶精度。

Z_E沿该例从1到4.1314228002。若遗漏式(3)的第二项，类光能量源的最大失配约4.6260039231；它是物理表达式漏项，远大于差分残差。

## 5. 正性与解释边界

从式(2)和319轮得到

$$
F>0,\quad \Theta_J\ge0
\quad\Longrightarrow\quad \theta_E\ge0,
\qquad \frac{dS_W^J}{d\lambda_J}
=F\,\frac{d(A_E/4G)}{d\lambda_E}\ge0.
\tag{5}
$$

原框架面积可下降，而Einstein框架面积与原Wald熵对应增加。这消除了此经典分支的一处表示歧义，没有证明“哪个框架的面积就是所有测量读数”。有质量测量装置与其耦合尚未加入。

正的光滑F也不自动保证两个框架的全局完备性、因果边界与渐近结构相同；本轮只核验同一正F局部片和所选末端。

9项检查包括面积、仿射生成元、独立曲率、收敛、聚焦、标量重定义、动能漏项、ξ＝0恒等映射，以及F≤0时拒绝套用正F变换。下一轮审查F＝0：这一定是原几何的曲率奇点，还是正F描述本身的适用边界？

[代码](320/einstein_frame_optics_audit.py)、[结果](320/einstein_frame_optics_audit_results.json)、[核验](320/research_round_320_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/einstein_frame_optics_audit.py
