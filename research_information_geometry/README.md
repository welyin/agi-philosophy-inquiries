# 信息几何与动力学重建

[返回项目目录](../README.md) · [认知物理研究主线](../research_cognition_physics/README.md)

## 研究范围

本支线考察信息与可区分性数据能否确定传播、动力学和共同几何，并检验连续极限、反作用及普适耦合需要哪些额外条件。各轮分别给出模型、推导或反例，以及数值检查。

这些结果属于指定模型或明确前提下的论证；详细条件与开放问题以各轮原报告为准。

## 逐轮入口

|轮次|研究问题|代码|结果|
|---|---|---|---|
|[01](research_note_01.md)|静态量子信息的非唯一性与局域响应重建|[experiment.py](experiment.py)|[结果](results.json)|
|[02](research_note_02.md)|制备误差、有限计数与动态重建的可识别性|[robustness.py](robustness.py)|[结果](robustness_results.json)|
|[03](research_note_03.md)|局域可逆更新与受控的一维Dirac极限|[relativistic_limit.py](relativistic_limit.py)|[结果](relativistic_results.json)|
|[04](research_note_04.md)|位置相关光锥、共形歧义与引力动力学的缺口|[variable_geometry.py](variable_geometry.py)|[结果](variable_geometry_results.json)|
|[05](research_note_05.md)|守恒反作用与量子—经典平均场的统计边界|[backreaction.py](backreaction.py)|[结果](backreaction_results.json)|
|[06](research_note_06.md)|联合量子演化、制备等价性与平均场丢失的关联|[quantum_background.py](quantum_background.py)|[结果](quantum_background_results.json)|
|[07](research_note_07.md)|两个物质探针的共同几何条件|[universal_geometry.py](universal_geometry.py)|[结果](universal_geometry_results.json)|
|[08](research_note_08.md)|软自旋2一致性对普适耦合的约束|[soft_graviton.py](soft_graviton.py)|[结果](soft_graviton_results.json)|

## 阅读与复算

先阅读对应报告的模型约定、结论范围和复算说明，再查看代码与已保存结果。[依赖说明](requirements.txt)列出运行要求。此目录沿用原支线编号，不与认知物理主线的轮次混用。
