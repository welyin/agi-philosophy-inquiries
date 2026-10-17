# 独立审阅：正常态到有限球面与未来视界的同源运输

2026-10-17。结论：**PASS_WITH_SCOPE**。

审阅通过范围为同一正常 Weyl 初态、固定背景上的真实自由场传播及同一最小耦合物理来源的有限通量身份。它不认证实际热库接触、几何熵或共同反作用。非编号；不增加正式科学组；CT 仍未采用，原目标仍未完成。本审只写此文件。

## 1. 独立解析核验

已全文审读 [README](README.md)、[证明](proof.md)、[复算入口](check.py)、[结果](results.json)及[来源基线](source_baseline.json)。除默认复算外，独立从径向 KG、Hilbert 应力和视界仿射方程核对如下链条，不以复算通过代替解析证明。

- **材料身份与正常性。** 最小材料质量 $m^2=2H^2$ 是选定背景上的常数，度规变分保持 $m$ 固定。势消去 $F(m^2+F'/r)=0$ 成立，不等同采用共形改进应力。全时解 $r\varphi=f(x-t)-f(-x-t)$ 保中心光滑；单独第一项只在 $t\ge0$ 使用。紧支 Cauchy 数据及薄时间片的 $g=P(\chi\varphi)$ 给楔内局部 Weyl 实现，故同 BD 表示中正常；不由求积推断正常性。公开输入这份有限能源初态无需追加未调用的准备工厂，实际调用的外源或读取则仍须另计。

- **完整能源与有限球通量。** 独立展开得到
  $$
  Q=2\pi\int\left[f'^2+\left(f'-\frac Frf\right)^2+m^2Ff^2\right]dx
  =4\pi\int f'^2dx .
  $$
  中间的质量项、梯度项及分部边界全部保留。对 $x=L$ 的球面，
  $$
  P_R=4\pi\left[f'(L-t)^2-\frac{F(R)}R f(L-t)f'(L-t)\right].
  $$
  完整窗口的交叉项由 $[f^2]_a^b=0$ 消失；截短窗口不可删除它，瞬时功率也未被误称处处非负。所选光滑波包满足
  $$
  Q\ge16\pi e^{-2/3}>0,
  $$
  因此较弱的正文下界同样有效。

- **未来方向及成熟归一。** 出射图满足 $\Gamma^v{}_{vv}=H$，故未来仿射参数为 $\lambda=e^{Hv}/H$，且 $\zeta=H\lambda\partial_\lambda$。视界迹为 $\varphi_{\mathcal H}=Hf(-v)$，从而
  $$
  T_{\lambda\lambda}=\frac{f'(-v)^2}{\lambda^2},\qquad
  Q=H\int\lambda T_{\lambda\lambda}\,d\lambda\,dA_0 .
  $$
  复用 Fröb–Much–Papadopoulos 原定理时，其生成元为 $\zeta/H$，所以完整静态楔相对熵是 $2\pi Q/H$。没有把有限球或有限视界子段改称具有相同模 Hamiltonian 的子区域。有限球的 $\beta_R=N_R\beta_0$、$E_R=Q/N_R$ 及固有时长 $N_R(b-a)$ 只是静态几何归一；代码的乘积检查未认证实际温度计。有限仿射视界段不授有限中心固有时的到达或读取。

- **面积边界条件。** 若另行输入线性 Einstein–Raychaudhuri 响应，令 $J(\lambda)=\int T_{\lambda\lambda}dA_0$、$a=\delta A$，则
  $$
  a''=-8\pi GJ,\qquad
  \beta_0Q=\frac{[a-\lambda a']_-^+}{4G}.
  $$
  未来平稳 $a'_+=0$ 给
  $$
  \Delta a=8\pi G\int_{\lambda_-}^{\lambda_+}(\lambda-\lambda_-)J(\lambda)\,d\lambda>0.
  $$
  初始零膨胀 $a'_-=0$ 则给
  $$
  \Delta a=-8\pi G\int_{\lambda_-}^{\lambda_+}(\lambda_+-\lambda)J(\lambda)\,d\lambda<0.
  $$
  两者是不同边界问题，不能共同施加于非零正通量，也不能省略端点导数而领取有限面积熵等式。正文对此范围保留正确。

## 2. 实际复算与范围判定

默认只读运行的真实回执为 **d8c5ff，exit 0，READ_ONLY_MATCH**。得到 $Q=76.06328048638551$、相对熵定理右端 $477.9196863679372$，384 节点四种能源积分最大差 $5.684341886080801\times10^{-13}$。这些浮点差仅是积分诊断，不是完整物理误差界。

来源及目标哈希核对回执为 **383a8d，exit 0**：基线列出的 8 项来源及原 goal 均匹配。此次不另造独立数值实现、不增加核验组；独立性来自上述来源、方向、边界和量词推导。

本包确实把上一包的正常态与初片荷，接成同一场的有限完整出流。旧 1004/318、1250 与 after1270 的一般相干来源、端点工具及其他材料预算直接复用；未把它们重复算成新定理。相干位移的非零能流是纯功样激发的运输，不能仅凭 $\beta_0Q$ 或参考 KMS 身份称为实际热库热流。相对熵也不是本包已认证的几何熵。

未获得耦合物质—度规解、真空变化与非线性反作用总余项、独立面积读口或热接触仪器。也不能把 1270/1272 的原材料权限移给新最小标量。幅度正负只改变场符号，来源与相对熵仍为二次；不授双向一阶饱和。上述限制不否定有限/二阶方向，也不增加所有候选必须先完成全 QFT 或全部设备的门槛。

## 3. 最终文件绑定

以下为本次实际读取并在保存本审阅前再次核对的 SHA256；不绑定将来生成的导航或出版记录。

|文件|SHA256|
|---|---|
|[README.md](README.md)|e9bdcf25d0406864cf080edbc00391cc8f5d155d4ead0b84bdb225d4f721e9cb|
|[proof.md](proof.md)|912f62c7b8f2e3106e4df5f3dced420c21a249e820559e26112e5e03d5e1df95|
|[check.py](check.py)|f9a4c40eb2ebbdc5d50addf2a7e557315be2eefbbde9de8cfdead7a79c1f2fa5|
|[results.json](results.json)|67fd91db75086f5a09cd253e0669bc760a2e489085cd53378a4eb7ab17118c44|
|[source_baseline.json](source_baseline.json)|8dae6f6ac34cc8c22fb5d86da7ca344a12c873e2d7397c9c020233a7e9a1097b|
