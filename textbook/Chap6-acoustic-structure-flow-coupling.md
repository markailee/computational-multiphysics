# 第6章 流动—结构—声场耦合的建模与求解策略

流动、结构和声场可以通过界面力、界面运动、材料变形和背景流场相互影响。流体压力和黏性力驱动结构运动，结构变形改变流体域和声学边界；非定常流动为声传播提供背景介质和声源，结构振动又向周围介质辐射声波。若温度变化显著，还会通过浮力、物性、热膨胀和热应力形成流动—热—结构耦合。

这些问题可以共用一个基本观点：**每个子域先满足自身控制方程，再通过界面运动学、动力学和守恒条件组成整体问题。** 本章按“流声分解—声学人工边界—流固耦合—流固热耦合—结构声耦合—统一求解策略”的顺序组织课件内容，并补充非匹配网格数据传递、附加质量、耦合松弛、多速率时间推进和分层验证方法。

除特别说明外，$\Omega_f$、$\Omega_s$ 和 $\Omega_a$ 分别表示流体、结构和声学子域，耦合界面记为 $\Gamma$。符号 $P$ 表示不可压或水动力压力，$p'$ 表示声学扰动压力，$p$ 在可压缩流动中表示总压力。不同子域的时间步可不同，但所有数据交换都必须明确对应的时间层。

## 6.1 耦合系统的统一描述

### 6.1.1 三类耦合路径

流—固—声系统包含以下主要作用路径。

- **流动与结构**：流体向结构施加压力和黏性面力，结构向流体提供界面位移和速度。
- **流动与声场**：背景流动决定声波对流、折射和部分声源项，声学扰动在弱可压分解中通常不显著反作用于平均流动。
- **结构与声场**：声压向结构施加扰动载荷，结构法向加速度向声场提供振动边界。
- **热效应**：流固界面交换温度和热流，结构温度产生热应变，流体温度反过来影响密度、黏度和声速。

并非所有问题都需要同时保留全部路径。刚性结构的流致噪声可省略结构变形；小振幅振动与静止介质可只解结构声耦合；温差很小且热膨胀可忽略时，可在流固耦合中省略热方程。模型简化的依据应是耦合强度和时间尺度，而不是求解器是否提供了相应模块。

### 6.1.2 子域方程与界面残差

将流体、结构、声学和热场的离散未知量分别记为 $\mathbf x_f$、$\mathbf x_s$、$\mathbf x_a$ 和 $\mathbf x_T$。整体系统可抽象为

$$
\mathbf R_m(\mathbf x_m,\mathbf x_{\setminus m})=\mathbf 0,
\qquad m\in\{f,s,a,T\},
\tag{6.1}
$$

其中 $\mathbf x_{\setminus m}$ 表示对子域 $m$ 有输入的其他场。分区法依次求解各个 $\mathbf R_m$，通过外迭代使界面条件收敛；整体法则把所有残差和未知量放入同一个非线性系统。

界面条件通常包括：

1. **运动学连续**：位移、速度或法向加速度相容；
2. **动力学平衡**：面力或热流满足作用与反作用；
3. **几何守恒**：运动网格体积变化与网格通量一致；
4. **功率一致**：非匹配网格映射不凭空产生或耗散界面功。

以结构外法线 $\mathbf n_s$ 指向声学域为例，小振幅结构声耦合中的界面功率可写为

$$
\mathcal P_{a\leftarrow s}
=\int_\Gamma p'\,\dot{\mathbf u}_s\cdot\mathbf n_s\,\mathrm d\Gamma.
\tag{6.2}
$$

结构从声场获得的功率应为 $-\mathcal P_{a\leftarrow s}$。离散压力映射和加速度映射若不能保持这一共轭关系，即使总力误差很小，也可能在长期计算中产生虚假能量。

### 6.1.3 单向、双向和整体耦合

单向耦合只保留一个作用方向，例如先把流场压力传递给结构，而不把结构变形反馈给流场。它适合反馈影响很小的问题，但不能通过增加外迭代“变成”双向耦合。

双向分区耦合分为：

- **松耦合或显式分区**：每个时间步只交换一次或少数几次数据，不在当前时间层消除全部界面残差；
- **紧耦合或隐式分区**：每个时间步内反复求解子域并交换数据，直至界面残差满足容差；
- **整体耦合**：所有变量在同一个矩阵或统一非线性框架中同步更新。

“紧耦合”描述时间层内的一致性，“整体耦合”描述代数系统组织方式，两者不是同义词。紧耦合分区法同样需要求解大量代数子系统，只是没有显式构造整体矩阵。

## 6.2 流声分解与声学 SIMPLE 算法

### 6.2.1 可压缩方程与热力学关系

黏性可压缩流动的连续性和动量方程为

$$
\frac{\partial\rho}{\partial t}
+\frac{\partial(\rho u_i)}{\partial x_i}=0,
\tag{6.3}
$$

$$
\frac{\partial(\rho u_i)}{\partial t}
+\frac{\partial(\rho u_iu_j)}{\partial x_j}
=-\frac{\partial p}{\partial x_i}
+\frac{\partial\tau_{ij}}{\partial x_j}+\rho b_i.
\tag{6.4}
$$

对牛顿流体，$\tau_{ij}$ 为黏性应力张量。状态方程写成 $p=p(\rho,S)$，其中 $S$ 为单位质量熵，则

$$
\mathrm dp=c^2\mathrm d\rho
+\left(\frac{\partial p}{\partial S}\right)_\rho\mathrm dS,
\qquad
c^2=\left(\frac{\partial p}{\partial\rho}\right)_S.
\tag{6.5}
$$

声速是等熵压缩性决定的热力学量。理想气体满足 $c=\sqrt{\gamma p/\rho}$；液体则可采用适当的状态方程或经验声速关系。不能在所有流体中都用理想气体公式。

### 6.2.2 不可压背景场与声学扰动

低马赫数流声分解把瞬时变量写为

$$
u_i=U_i+u_i',
\qquad
p=P+p',
\qquad
\rho=\rho_0+\rho'.
\tag{6.6}
$$

其中 $\mathbf U$、$P$ 和 $\rho_0$ 是不可压背景场的速度、压力和密度，$u_i'$、$p'$ 和 $\rho'$ 是声学扰动。$P$ 可以随时间变化，并可进一步分成时间平均值 $\overline P$ 和水动力脉动 $P^*$；它不等同于声压 $p'$。涡脱落引起的近场压力脉动通常主要属于 $P$，传播到远场的可压缩压力扰动才主要属于 $p'$。

定义声学动量扰动矢量

$$
f_i=\rho'(U_i+u_i')+\rho_0u_i'
=\rho u_i-\rho_0U_i.
\tag{6.7}
$$

将式（6.6）代入可压缩连续性方程，并利用 $\partial U_i/\partial x_i=0$，得到

$$
\frac{\partial\rho'}{\partial t}
+\frac{\partial f_i}{\partial x_i}=0.
\tag{6.8}
$$

若快速压力扰动近似等熵，而熵变化主要归入背景场，则压力与密度扰动满足时间关系

$$
\frac{\partial p'}{\partial t}
+\frac{\partial P}{\partial t}
=c^2\frac{\partial\rho'}{\partial t}.
\tag{6.9}
$$

因此声学连续性方程也可写为

$$
\boxed{
\frac1{c^2}\frac{\partial p'}{\partial t}
+\frac{\partial f_i}{\partial x_i}
=-\frac1{c^2}\frac{\partial P}{\partial t}
}.
\tag{6.10}
$$

右端是不可压压力时间变化产生的源项。若 $P$ 随时间变化，不能简单地在所有时刻令 $p'=c^2\rho'$；应根据式（6.9）或一致的状态方程更新 $\rho'$。

忽略黏性项后，课件中的声学动量方程可整理为

$$
\frac{\partial f_i}{\partial t}
+\frac{\partial[f_i(U_j+u_j')]}{\partial x_j}
=-\frac{\partial(\rho_0U_i u_j')}{\partial x_j}
-\frac{\partial p'}{\partial x_i}.
\tag{6.11}
$$

完整形式还可保留黏性扩散、背景速度梯度以及密度扰动相关项。实现时应固定采用完整或线性化方程，不能在一套推导中随意丢弃非线性项，却在另一套系数中重新加入。

在静止均匀介质、小振幅且 $P$ 定常的极限下，$f_i\approx\rho_0u_i'$，式（6.10）和式（6.11）退化为线性声学方程，并进一步得到

$$
\frac1{c^2}\frac{\partial^2p'}{\partial t^2}
-\nabla^2p'=0.
\tag{6.12}
$$

这一极限是检查流声分解程序的重要基准。对于强剪切、大振幅、激波或显著熵波，简单低马赫分解不一定足够，应使用可压缩方程或更完整的声学扰动方程。

### 6.2.3 有限体积离散与时间格式

对控制体 $V_P$ 积分式（6.10），有

$$
\frac{V_P}{c_P^2}\frac{\mathrm Dp'_P}{\mathrm Dt}
+\sum_{f\in\mathcal F_P}\mathbf f_f\cdot\mathbf A_f
=-\frac{V_P}{c_P^2}\frac{\mathrm DP_P}{\mathrm Dt}.
\tag{6.13}
$$

一阶隐式欧拉具有较强数值耗散。二阶 backward differentiation formula（BDF2）写为

$$
\left.\frac{\partial\phi}{\partial t}\right|^{n+1}
\approx\frac{3\phi^{n+1}-4\phi^n+\phi^{n-1}}{2\Delta t},
\tag{6.14}
$$

可用于声压和声学动量时间项。BDF2 起步时需要一阶格式或其他二阶起步方法；变时间步时不能直接沿用等步长系数。

对流项可采用迎风、高阶 TVD、中心或低耗散高阶格式。课件采用“一阶迎风加限制器修正”的 TVD 形式。其优点是抑制间断附近振荡，但限制器也可能在光滑声波极值处生效并增加耗散。长距离声传播更关注色散和耗散，应通过平面波数值色散试验选择格式，而不能只以“高阶”或“TVD”判断声学精度。

### 6.2.4 声学动量预测与 Rhie–Chow 插值

将式（6.11）的通用对流—扩散形式离散，使用当前声压估计 $p'^*$ 求解声学动量预测方程：

$$
A_P\mathbf f_P^*
=\sum_{NB}A_{NB}\mathbf f_{NB}^*
-V_P\nabla p_P'^*+\mathbf S_P.
\tag{6.15}
$$

系数包含时间项、声对流、扩散和松弛影响，$\mathbf S_P$ 包含历史项、延迟对流修正、背景场源项及非正交修正。

同位网格上，界面动量不能直接只做单元中心算术平均，否则可能产生声压棋盘振荡。Rhie–Chow 形式可写为

$$
\mathbf f_f^*
=\overline{\mathbf f^*}_f
+D_f\left[
\overline{\nabla p'^*}_f
-(\nabla p'^*)_f^{1-\delta}
\right],
\qquad
D_f=\overline{\left(\frac{V}{A_P}\right)}_f,
\tag{6.16}
$$

其中 $(\nabla p'^*)_f^{1-\delta}$ 由界面两侧单元中心声压直接构造。用于连续性方程的界面通量必须采用这一动量插值结果。

### 6.2.5 声压修正与变量更新

定义修正量

$$
\mathbf f=\mathbf f^*+\mathbf f^c,
\qquad
p'=p'^*+p^c.
\tag{6.17}
$$

修正后的动量也满足同一层次的离散方程。与流场 SIMPLE 类似，忽略邻点动量修正后得到

$$
A_P\mathbf f_P^c=-V_P\nabla p_P^c,
\qquad
\mathbf f_P^c=-\frac{V_P}{A_P}\nabla p_P^c.
\tag{6.18}
$$

代入离散连续性方程，声压修正方程为

$$
\boxed{
\frac{V_P}{c_P^2}\mathrm D_tp_P^c
-\sum_{f\in\mathcal F_P}D_f(\nabla p^c)_f\cdot\mathbf A_f
=-\sum_{f\in\mathcal F_P}\mathbf f_f^*\cdot\mathbf A_f
-\frac{V_P}{c_P^2}\mathrm D_t(P_P+p_P'^*)
}.
\tag{6.19}
$$

右端是预测场的负声学质量残差。求解后可按

$$
p'\leftarrow p'^*+\alpha_ap^c,
\qquad
\mathbf f\leftarrow\mathbf f^*+\mathbf f^c
$$

更新，其中 $\alpha_a$ 为声压亚松弛因子。随后由式（6.9）更新 $\rho'$，并由式（6.7）恢复质点振速：

$$
\mathbf u'
=\frac{\mathbf f-\rho'\mathbf U}{\rho_0+\rho'}.
\tag{6.20}
$$

一个时间步内的声学 SIMPLE 迭代为：

```text
由不可压流场更新 U、P、rho_0、c 和声学源项
用当前 p'* 预测声学动量 f*
用 Rhie–Chow 插值构造界面声学动量通量
组装并求解声压修正方程，得到 p^c
更新 p'、f、rho'、u'、界面通量和延迟项
检查声学连续性、动量、状态关系和耦合界面残差
未收敛则继续当前时间步迭代
```

不可压背景流场和声场可以采用不同的内迭代次数。若背景流场尚未收敛就把其压力源项当作最终输入，声场可能收敛到错误的源场。

## 6.3 声学人工边界

### 6.3.1 特征波与出流边界

有限计算域不能简单把外边界设成固定声压或零法向梯度，因为这两种处理都可能产生明显反射。在线性化 Euler 方程中，沿边界法向的特征速度通常为

$$
\lambda_1=U_n,
\qquad
\lambda_{2,3}=U_n\pm c.
\tag{6.21}
$$

$U_n$ 对应随平均流动输运的熵波和涡波，$U_n\pm c$ 对应两个声波特征。亚声速出流边界上，部分特征向域外传播，另一部分可能由边界条件决定；边界处理必须与波的类型一致，而不是对所有变量施加同一种条件。

### 6.3.2 Tam 辐射与出流条件

Tam 型辐射条件利用远场声波的渐近形式。设声源参考位置为 $\mathbf x_s$，边界点为 $\mathbf x_b$，$r=|\mathbf x_b-\mathbf x_s|$，则对声学量 $\phi$ 可采用

$$
\left[
\frac1{V(\theta)}\frac{\partial}{\partial t}
+\frac{\partial}{\partial r}
+\frac{d-1}{2r}
\right]\phi=0,
\qquad d=1,2,3.
\tag{6.22}
$$

$(d-1)/(2r)$ 表示柱面波或球面波的几何扩散。对 $d=1,2,3$，该系数分别为 $0$、$1/(2r)$ 和 $1/r$。在仅含声波的辐射边界，可对 $\rho'$、质点振速和 $p'$ 使用相应关系；在同时存在涡波和熵波的出流边界，压力条件与速度、密度条件需要区别处理。

一种平均流下的传播速度估计为

$$
V(\theta)=\mathbf U\cdot\mathbf e_r
+c\sqrt{1-M^2\sin^2\theta},
\qquad
M=\frac{|\mathbf U|}{c},
\tag{6.23}
$$

其中 $\theta$ 是平均流动方向与 $\mathbf e_r$ 的夹角。该关系来自特定的均匀平均流和远场假设，不能机械用于强非均匀近场或超声速边界。

有限体积实现时，可把 $\partial/\partial r$ 转为坐标梯度，用单元面积向量重构单元中心梯度，再把边界条件整理成边界单元的中心系数、邻点系数和源项。Tam–Webb 的 DRP 空间格式与 Tam 辐射/出流边界是相关但不同的概念：前者降低数值色散，后者控制域边界的波反射。

### 6.3.3 PML 的坐标伸缩与辅助变量

完美匹配层（PML）在物理计算域外增加吸收层。频域坐标伸缩可概括为

$$
\frac{\partial}{\partial x_i}
\longrightarrow
\frac1{s_i(\omega)}\frac{\partial}{\partial x_i},
\qquad
s_i(\omega)=1+\frac{\sigma_i}{s_\omega\omega},
\tag{6.24}
$$

其中采用 $e^{-\mathrm i\omega t}$ 时取 $s_\omega=\mathrm i$，采用相反的 Fourier 约定时取 $s_\omega=-\mathrm i$。$\sigma_i=0$ 时，PML 方程必须与内部声学方程一致。

直接求解频域方程不便与内部时域求解器耦合，因此引入辅助变量 $q_i$，把伸缩算子恢复为时域松弛方程。一种典型形式为

$$
\frac{\partial q_i}{\partial t}+\sigma_iq_i
=\frac{\partial G'}{\partial x_i},
\tag{6.25}
$$

再用 $q_i$ 的时间导数替换原方程中的空间导数组合。若还采用时间坐标变换，右端会增加相应的 $\beta$ 项。辅助变量与主变量必须在同一时间层更新，PML 界面处变量和通量保持连续。

### 6.3.4 吸收函数、网格与时间推进

常用多项式吸收函数为

$$
\sigma_i(x_i)=\sigma_{\max}
\left(\frac{x_i-x_{i,0}}{L_i}\right)^m,
\qquad x_i>x_{i,0},
\tag{6.26}
$$

其中 $m$ 常取2或3。若外层采用反射壁面，可按目标往返反射系数 $R$ 估算

$$
\sigma_{\max}
\approx-\frac{(m+1)c\ln R}{2L_i}.
\tag{6.27}
$$

这只是初始设计值，最佳吸收强度还受网格拉伸、斜入射、离散格式和 Evanescent 波影响。吸收系数突变或过强会产生离散界面反射；过弱则吸收不足。PML 起始网格应与内部网格平滑衔接，吸收函数在界面处从零连续增加。

课件采用二阶段 SSP Runge–Kutta 显式推进。对 $\dot G'=R(G')$，其形式为

$$
G'^*=G'^n+\Delta tR(G'^n),
\qquad
G'^{n+1}=\frac12\left[G'^n+G'^*+\Delta tR(G'^*)\right].
\tag{6.28}
$$

显式 PML 的时间步除声学 CFL 外，还受最大 $\sigma_i$ 限制；强吸收层可能成为稳定性控制区域。

### 6.3.5 人工边界的选择

| 项目 | Tam 辐射/出流边界 | PML |
|---|---|---|
| 主要思想 | 按特征波和远场渐近形式抑制反射 | 在附加吸收层内衰减传播波 |
| 附加网格 | 通常不需要 | 需要一定厚度的吸收层 |
| 斜入射与复杂波场 | 依赖远场和平均流假设 | 通常更宽容，但仍有离散反射 |
| 实现成本 | 边界单元公式较集中 | 辅助变量、参数和网格更多 |
| 涡波、熵波 | 必须区别处理 | 需使用与完整扰动方程相容的形式 |
| 主要验证 | 脉冲出域反射和远场声压 | 层厚、吸收强度、网格拉伸和入射角扫描 |

验证时应记录域内监测点在脉冲出域前后的波形，并与更大计算域参考解比较。只观察计算“没有发散”不能证明人工边界无反射。

## 6.4 流固耦合与 ALE 动网格

### 6.4.1 ALE 流体方程

在任意 Lagrange–Euler（ALE）框架中，网格以速度 $\mathbf U_g$ 运动。积分形式的连续性方程为

$$
\frac{\mathrm d}{\mathrm dt}\int_{V(t)}\rho\,\mathrm dV
+\int_{\partial V(t)}\rho(\mathbf U-\mathbf U_g)\cdot\mathbf n\,\mathrm dS=0.
\tag{6.29}
$$

动量方程为

$$
\begin{aligned}
&\frac{\mathrm d}{\mathrm dt}\int_{V(t)}\rho\mathbf U\,\mathrm dV
+\int_{\partial V(t)}\rho\mathbf U(\mathbf U-\mathbf U_g)\cdot\mathbf n\,\mathrm dS\\
&\quad=\int_{\partial V(t)}(-p\mathbf I+\boldsymbol\tau)\mathbf n\,\mathrm dS
+\int_{V(t)}\rho\mathbf b\,\mathrm dV.
\end{aligned}
\tag{6.30}
$$

网格速度只改变相对对流通量，不应被误加到物理速度边界条件中。控制体体积随时间变化，时间项中的 $V_P^{n+1}$ 和 $V_P^n$ 也必须保留。

### 6.4.2 几何守恒定律

均匀流在任意网格运动下应保持不变，这要求离散几何守恒定律（GCL）成立：

$$
\frac{V_P^{n+1}-V_P^n}{\Delta t}
=\sum_{f\in\mathcal F_P}\mathbf U_{g,f}\cdot\mathbf A_f
=\sum_{f\in\mathcal F_P}\dot V_f.
\tag{6.31}
$$

右端更适合由面在一个时间步内扫过的体积 $\dot V_f$ 计算，而不是仅用某个任意时刻的网格速度乘面积。若网格位移更新与时间积分不同步，即使结构位移很小，也可能产生虚假质量源和压力振荡。自由流保持试验是动网格程序的必要检查。

### 6.4.3 结构方程与流固界面条件

结构离散方程沿用第5章的形式：

$$
\mathbf M_s\ddot{\mathbf d}
+\mathbf C_s\dot{\mathbf d}
+\mathbf K_s\mathbf d
=\mathbf F_s^{\mathrm{other}}+\mathbf F_{f\rightarrow s}.
\tag{6.32}
$$

非线性结构应把 $\mathbf K_s\mathbf d$ 替换为内力矢量 $\mathbf f_{\mathrm{int}}(\mathbf d)$。

无滑移流固界面满足位移和速度连续：

$$
\mathbf x_g=\mathbf x_{g,0}+\mathbf d_\Gamma,
\qquad
\mathbf U_f=\dot{\mathbf d}_\Gamma
\qquad \text{在 }\Gamma_{fs}\text{ 上}.
\tag{6.33}
$$

无滑移条件还包括切向速度；无黏或滑移模型通常只约束法向运动。动力学条件要求两侧面力作用与反作用。令 $\mathbf n_f$ 和 $\mathbf n_s$ 分别为流体和结构外法线，则 $\mathbf n_s=-\mathbf n_f$，且

$$
\boldsymbol\sigma_f\mathbf n_f
+\boldsymbol\sigma_s\mathbf n_s=\mathbf 0,
\qquad
\boldsymbol\sigma_f=-p\mathbf I+\boldsymbol\tau.
\tag{6.34}
$$

结构载荷由流体总 Cauchy 应力积分得到，不能只保留压力而遗漏重要黏性面力。若流声分解同时存在，则总法向压力应由 $P+p'$ 组成，但平均预载和声学扰动载荷不能重复计算。

### 6.4.4 非匹配网格上的数据传递

非匹配流固网格需要空间映射。位移映射应满足界面几何相容，力映射则应满足

$$
\sum_a\mathbf F_{s,a}=\int_\Gamma\mathbf t\,\mathrm d\Gamma,
\qquad
\sum_a(\mathbf x_a\times\mathbf F_{s,a})
=\int_\Gamma\mathbf x\times\mathbf t\,\mathrm d\Gamma.
\tag{6.35}
$$

若位移映射矩阵为 $\mathbf H_{s\leftarrow f}$，虚功共轭要求力映射使用相应的转置加权形式，而不是另行采用不一致的最近节点平均。还应检查常面力、线性面力及常虚功场的映射误差。

对于大滑动、接触或浸没边界，不能始终假设节点一一对应。需要采用投影、切割面元或重叠网格守恒插值，并处理投影失败、空隙和穿透。

### 6.4.5 动网格方法

| 方法 | 基本处理 | 优点 | 主要限制 |
|---|---|---|---|
| 弹簧类比 | 网格边视为弹簧，解节点力平衡 | 实现简单、速度快 | 大变形易畸变，弹簧刚度影响大 |
| 扩散光顺 | 求解位移扩散方程 | 稳健、可控制边界层 | 极大变形仍可能失败 |
| 网格重构 | 删除畸变网格并重新生成 | 可处理较大变形 | 成本高，场插值引入误差 |
| 层铺法 | 增减或压缩边界层网格 | 法向活塞运动效率高 | 适用运动形式有限 |
| 浸没边界 | 固定背景网格，用边界力或切割单元表示结构 | 适合大位移 | 近壁精度与守恒处理复杂 |
| 重叠网格 | 多套网格独立运动并插值 | 适合多体相对运动 | 洞映射、插值和守恒实现复杂 |

弹簧类比可组装为 $\mathbf K_g\Delta\mathbf X=\mathbf B_g$；扩散光顺则常求解 $\nabla\cdot(D_g\nabla\mathbf X_g)=0$ 并施加界面位移。近壁网格可增大 $D_g$ 或弹簧刚度，以尽量保持边界层厚度和单元正交性。每次更新后必须检查负体积、扭曲和 GCL，不能只检查界面位置。

## 6.5 流固热耦合与结构热应力

### 6.5.1 流体和固体能量方程

不可压或低速流体的焓形式能量方程可写为

$$
\frac{\partial(\rho h)}{\partial t}
+\nabla\cdot[\rho(\mathbf U-\mathbf U_g)h]
=\nabla\cdot(k_f\nabla T)+\dot q_f.
\tag{6.36}
$$

固体热传导方程为

$$
\rho_sc_{p,s}\frac{\partial T_s}{\partial t}
=\nabla\cdot(k_s\nabla T_s)+Q_s.
\tag{6.37}
$$

密度、黏度、导热系数和声速若随温度变化，应及时更新。可压缩流应使用一致的总能量或总焓形式，并保留压力功、黏性耗散和 ALE 相对通量。

### 6.5.2 共轭传热界面

理想热接触界面满足温度连续和热流守恒：

$$
T_f=T_s,
\qquad
-k_f\frac{\partial T_f}{\partial n}
=-k_s\frac{\partial T_s}{\partial n}
\qquad \text{在 }\Gamma_{fs}\text{ 上},
\tag{6.38}
$$

这里两侧使用同一个界面法向和法向导数方向。若分别使用各子域外法线，则应写成 $\mathbf q_f\cdot\mathbf n_f+\mathbf q_s\cdot\mathbf n_s=0$，避免热流符号混淆。

采用湍流壁面函数时，流体侧热流可由壁面函数计算，但它必须与固体侧导热通量在界面收敛时一致。存在接触热阻 $R_c$ 时，温度不再连续，而应满足 $q=(T_f-T_s)/R_c$；辐射、蒸发或冷凝也应加入同一界面能量平衡。

共轭传热分区迭代可采用 Dirichlet–Neumann 交换：一侧接收温度，另一侧接收热流，下一轮交换角色或保持固定方向。高温差、强非线性或热物性变化显著时，Robin–Robin 条件和适当松弛通常更有利于收敛。

### 6.5.3 热应变与热应力

各向同性线弹性体的热应变为

$$
\boldsymbol\varepsilon_{\mathrm{th}}
=\alpha_T(T-T_{\mathrm{ref}})\mathbf I,
\tag{6.39}
$$

其中 $\alpha_T$ 为线膨胀系数，$T_{\mathrm{ref}}$ 为无热应力参考温度。总机械本构关系变为

$$
\boldsymbol\sigma
=\mathsf C:(\boldsymbol\varepsilon
-\boldsymbol\varepsilon_{\mathrm{th}}
-\boldsymbol\varepsilon_0),
\tag{6.40}
$$

其中 $\boldsymbol\varepsilon_0$ 表示其他初始或本征应变。热膨胀只在受约束或温度梯度不相容时产生热应力；自由均匀膨胀不应产生虚假应力。

结构方程因此写为

$$
\mathbf M_s\ddot{\mathbf d}
+\mathbf C_s\dot{\mathbf d}
+\mathbf K_s\mathbf d
=\mathbf F_f+\mathbf F_a+\mathbf F_{\mathrm{th}}
+\mathbf F_s^{\mathrm{other}}.
\tag{6.41}
$$

在常用位移有限元记号下，等效热载荷可写为

$$
\mathbf F_{\mathrm{th}}
=\int_{\Omega_s}\mathbf B^{\mathrm T}
\mathsf C\boldsymbol\varepsilon_{\mathrm{th}}\,\mathrm d\Omega.
\tag{6.42}
$$

采用第5章的控制体面通量法时，应直接把式（6.40）代入面力通量，或推导出与式（6.42）等价的离散热载荷，不能重复施加两项。

## 6.6 结构声耦合

### 6.6.1 线性声学方程

静止、均匀、无黏介质中的小振幅声学扰动满足

$$
\frac{\partial\rho'}{\partial t}
+\rho_0\nabla\cdot\mathbf u'=q_m,
\tag{6.43}
$$

$$
\rho_0\frac{\partial\mathbf u'}{\partial t}
=-\nabla p'.
\tag{6.44}
$$

利用 $p'=c^2\rho'$ 并消去速度，得到带源波动方程

$$
\frac1{c^2}\frac{\partial^2p'}{\partial t^2}
-\nabla^2p'=q_a.
\tag{6.45}
$$

若介质密度和声速随空间变化，应保留守恒形式 $\nabla\cdot(\rho_0^{-1}\nabla p')$，不能直接用 $\nabla^2p'$ 替代。黏热损耗、平均流和频散材料需要相应扩展模型。

### 6.6.2 界面力平衡和加速度连续

取结构外法线 $\mathbf n_s$ 指向声学域。小变形结构声界面上，声压对结构的面力为 $-p'\mathbf n_s$，因此

$$
\boldsymbol\sigma_s\mathbf n_s=-p'\mathbf n_s.
\tag{6.46}
$$

由声学 Euler 方程和界面运动连续，有

$$
\frac{\partial p'}{\partial n_s}
=-\rho_0\ddot{\mathbf u}_s\cdot\mathbf n_s.
\tag{6.47}
$$

式（6.46）是结构的力边界，式（6.47）是声场的 Neumann 边界。二者只约束法向加速度，理想声学介质不对结构施加切向速度约束。

若结构还承受不可压背景压力 $P$，则总压力载荷应包含 $P+p'$；若平均压力已形成独立预应力状态，则振动方程中只保留扰动部分，并考虑预应力对刚度的影响。

### 6.6.3 界面离散与耦合矩阵

设结构和声学界面形函数分别为 $\mathbf N_s$、$\mathbf N_a$，定义耦合面积分矩阵

$$
\mathbf G_{sa}
=\int_\Gamma\mathbf N_s^{\mathrm T}
\mathbf n_s\mathbf N_a\,\mathrm d\Gamma.
\tag{6.48}
$$

结构受到的声学节点力为 $-\mathbf G_{sa}\mathbf p$。采用一致的弱形式或与其相容的有限体积离散，可得耦合系统

$$
\mathbf M_s\ddot{\mathbf d}
+\mathbf C_s\dot{\mathbf d}
+\mathbf K_s\mathbf d
+\mathbf G_{sa}\mathbf p
=\mathbf F_s,
\tag{6.49}
$$

$$
\mathbf M_a\ddot{\mathbf p}
+\mathbf C_a\dot{\mathbf p}
+\mathbf K_a\mathbf p
+\rho_0\mathbf G_{sa}^{\mathrm T}\ddot{\mathbf d}
=\mathbf Q_a.
\tag{6.50}
$$

转置耦合项来自界面虚功共轭。不同文献可能把 $\rho_0$、负号或声速因子并入 $\mathbf G_{sa}$，实现时必须从自己的界面公式推导，不能只按矩阵外观判断。

在频域中，令 $\mathbf d=\widehat{\mathbf d}e^{\mathrm i\omega t}$、$\mathbf p=\widehat{\mathbf p}e^{\mathrm i\omega t}$，得到

$$
\begin{bmatrix}
\mathbf K_s+\mathrm i\omega\mathbf C_s-\omega^2\mathbf M_s&
\mathbf G_{sa}\\
-\omega^2\rho_0\mathbf G_{sa}^{\mathrm T}&
\mathbf K_a+\mathrm i\omega\mathbf C_a-\omega^2\mathbf M_a
\end{bmatrix}
\begin{bmatrix}\widehat{\mathbf d}\\\widehat{\mathbf p}\end{bmatrix}
=
\begin{bmatrix}\widehat{\mathbf F}_s\\\widehat{\mathbf Q}_a\end{bmatrix}.
\tag{6.51}
$$

该系统通常非对称且可能病态。可通过变量缩放改善条件数，并采用合适的稀疏直接法、迭代法和块预条件。耦合固有频率不能仅由结构真空模态或声学刚性壁模态单独判断。

### 6.6.4 显式中心差分与稳定时间步

结构采用速度半时间层的中心差分时，可写为

$$
\begin{aligned}
\mathbf a_s^n&=\mathbf M_s^{-1}
(\mathbf F_s^n-\mathbf G_{sa}\mathbf p^n
-\mathbf C_s\mathbf v_s^{n-1/2}-\mathbf K_s\mathbf d_s^n),\\
\mathbf v_s^{n+1/2}&=\mathbf v_s^{n-1/2}
+\Delta t_s\mathbf a_s^n,\\
\mathbf d_s^{n+1}&=\mathbf d_s^n+\Delta t_s\mathbf v_s^{n+1/2}.
\end{aligned}
\tag{6.52}
$$

无阻尼声学波动方程的显式中心差分为

$$
\mathbf p^{n+1}
=2\mathbf p^n-\mathbf p^{n-1}
+\Delta t_a^2\mathbf M_a^{-1}
(\mathbf Q_a^n-\mathbf K_a\mathbf p^n
-\rho_0\mathbf G_{sa}^{\mathrm T}\mathbf a_s^n).
\tag{6.53}
$$

若有 $\mathbf C_a\dot{\mathbf p}$ 项，一阶时间导数可用 $(\mathbf p^{n+1}-\mathbf p^{n-1})/(2\Delta t_a)$ 中心差分，并与二阶导数联立求解新时间层。

显式推进受 CFL 条件限制：

$$
Co=\frac{c\Delta t}{h_{\mathrm{eff}}}\leq Co_{\max}.
\tag{6.54}
$$

课件给出的特定二维格式参考值为：声学三角形单元 $Co_{\max}=\sqrt{2/3}$、四边形单元为1；结构纵波对应值分别为 $1/\sqrt2$ 和1。实际稳定上限依赖网格类型、维度、质量集中、时间格式和耦合方式，应通过数值稳定分析或试验确认，并留安全裕量。

稳定性不等于精度。低阶单元通常每个波长至少需要约10个节点作为起点，长距离传播或要求低相位误差时可能需要更高空间阶数和更多每波节点数。结构纵波速度 $c_p$ 和声速 $c$ 不同，故两个子域常有不同的最佳时间步。

## 6.7 统一分区求解与多速率策略

### 6.7.1 求解策略比较

| 策略 | 时间层内数据交换 | 优点 | 主要风险或成本 |
|---|---|---|---|
| 单向耦合 | 只沿一个方向传递一次 | 成本低、稳定、实现简单 | 无法描述反馈作用 |
| 双向松耦合 | 每步一次或少量交换 | 可复用独立求解器，易于并行 | 时间滞后、能量误差、强耦合失稳 |
| 双向紧耦合 | 每步多次 Picard 或 quasi-Newton 迭代 | 界面一致性好，适合强耦合 | 每步成本高，需要松弛和收敛控制 |
| 整体耦合 | 所有变量同步求解 | 可避免分区外迭代，耦合一致 | 矩阵大、非对称、预条件和软件开发困难 |

整体耦合不保证总计算时间一定更短。它能减少或消除分区外迭代，但可能增加单次线性求解和内存成本。选择策略时应比较达到同一界面残差和物理精度所需的总成本。

### 6.7.2 附加质量与 Picard 松弛

轻质结构与稠密流体耦合时，结构微小加速度会带动显著流体惯性，形成很强的附加质量效应。传统串行显式分区即使各子域单独稳定，也可能在耦合后失稳；简单增大结构质量或减小时间步只是数值处理，不能替代物理模型和耦合算法改进。

定义结构界面位移的定点残差

$$
\mathbf r^{(k)}
=\mathcal H(\mathcal F(\mathbf d_\Gamma^{(k)}))
-\mathbf d_\Gamma^{(k)},
\tag{6.55}
$$

其中 $\mathcal F$ 表示结构位移经流场求解产生的载荷，$\mathcal H$ 表示结构求解返回的新位移。固定点松弛更新为

$$
\mathbf d_\Gamma^{(k+1)}
=\mathbf d_\Gamma^{(k)}+\omega_k\mathbf r^{(k)}.
\tag{6.56}
$$

Aitken 动态松弛可写为

$$
\omega_k=-\omega_{k-1}
\frac{(\mathbf r^{(k-1)})^{\mathrm T}
(\mathbf r^{(k)}-\mathbf r^{(k-1)})}
{\|\mathbf r^{(k)}-\mathbf r^{(k-1)}\|_2^2},
\tag{6.57}
$$

并应对 $\omega_k$ 设置合理上下限、处理分母过小，必要时退回固定松弛。强耦合问题还可采用接口 quasi-Newton、Anderson 加速或整体 Newton 方法。松弛因子越小不一定越好：它能改善稳定性，却会显著增加耦合迭代成本。

### 6.7.3 多速率时间推进

设声学同步时间步与结构子步满足

$$
\Delta t_a=m\Delta t_s,
\qquad m\in\mathbb N.
\tag{6.58}
$$

在 $[t_0,t_0+\Delta t_a]$ 内，结构执行 $m$ 个子步，声场在同步时刻推进一次。结构子步所需声压可用已知声学时间层插值：

$$
p'(t)\approx\sum_{\ell}p'_\ell L_\ell(t),
\tag{6.59}
$$

其中 $L_\ell$ 为时间插值基函数。声场所需的结构加速度可在同步区间累积、平均或用相容的高阶时间重构获得。

多速率实现应遵守：

- 每个子域的实际步长都满足自身稳定性和精度要求；
- 所有子域在同步时刻对齐，避免一个求解器使用另一求解器尚未计算的未来值；
- 时间插值阶数与子域时间格式相容，谨慎使用高阶外推；
- 压力和加速度的时间映射尽量保持界面功一致；
- 多速率结果与统一小时间步结果进行时间收敛比较。

若声学步长反而比结构步长更小，可交换子步角色或选择共同同步周期。简单取所有子域最小时间步最稳健但成本最高；盲目使用超过局部稳定上限的“最优步长”则不可接受。

### 6.7.4 界面与耦合收敛判据

耦合收敛不能只看某个子域的内部残差。可组合使用界面运动、力、温度和热流残差：

$$
\eta_\Gamma=\max\left\{
\frac{\|\mathbf r_d\|}{d_{\mathrm{ref}}},
\frac{\|\mathbf r_t\|}{F_{\mathrm{ref}}},
\frac{\|\mathbf r_T\|}{T_{\mathrm{ref}}},
\frac{\|\mathbf r_q\|}{q_{\mathrm{ref}}}
\right\},
\tag{6.60}
$$

并同时检查各场位移、声压、速度和温度在相邻耦合迭代间的变化。参考尺度应与当前问题一致，避免在零载荷或零位移时使用不稳定的相对误差。

每个时间步结束还应检查：

- 流体质量、动量和能量残差；
- 结构动力平衡和约束反力；
- 声学连续性、状态关系和人工边界残差；
- 界面总力、总力矩、热流和功率；
- GCL、网格质量和多速率同步误差。

## 6.8 综合求解流程与验证

### 6.8.1 流—固—声统一时间步

一种可扩展的紧耦合分区流程如下。没有热耦合或结构反馈时，可跳过相应步骤。

```text
初始化流场、结构位移/速度、声场、温度场和网格历史
进入同步时间步 [t_n, t_(n+1)]
    由上一时间层或上一轮迭代预测界面位移和温度
    开始耦合外迭代 k
        更新 ALE 网格、网格速度和扫过体积，检查 GCL
        求解背景流动与流体能量方程，得到 U、P、T_f
        更新流声源项和背景声速
        求解声学 SIMPLE 迭代和人工边界，得到 p'、u'、rho'
        传递流体温度和热流，求解固体热传导
        合成流动面力、声学面力和热载荷，求解结构运动
        将结构界面速度分别反馈给 ALE 流场和声场
        松弛界面位移、温度或其他强耦合量
        检查各子域残差、界面残差和守恒误差
    直到耦合收敛或达到最大迭代数
    接受当前时间层，保存历史变量和统计数据
进入下一同步时间步
```

结构载荷合成时应明确

$$
\mathbf F_{f\rightarrow s}
=\int_\Gamma[-(P+p')\mathbf n_s
+\boldsymbol\tau\mathbf n_s]\,\mathrm d\Gamma
$$

或其他与法向约定等价的形式。若流场求解器已经把 $p'$ 包含在总压力中，则不能再额外添加一次声学压力。声学子域与流固子域也可采用不同界面网格，但映射后必须重新检查总力和功率。

### 6.8.2 分层验证算例

复杂耦合计算应分层验证，而不是直接比较一张最终云图。

1. **背景流场**：验证网格、湍流模型、动网格自由流保持和 GCL。
2. **声场**：用静止介质平面波、球面波或单极子验证声速、波长、数值色散和 PML/Tam 反射。
3. **流声耦合**：检查均匀平均流下声波对流和折射，以及不可压压力源项的时间层一致性。
4. **流固耦合**：用活塞、弹性腔体或标准弹性构件验证界面位移、总力、附加质量和能量。
5. **共轭传热**：用一维双层导热和给定对流系数验证温度连续、热流守恒和接触热阻。
6. **热应力**：用自由热膨胀和完全约束热应力检查参考温度及符号。
7. **结构声耦合**：用振动活塞辐射、耦合腔体模态或解析阻抗验证式（6.46）和式（6.47）。
8. **多速率与完整系统**：与统一小时间步、整体耦合或高保真参考解比较，并执行网格和时间步收敛。

声压后处理可使用

$$
p_{\mathrm{rms}}=\sqrt{\frac1{T_s}\int_{t_1}^{t_1+T_s}p'^2\,\mathrm dt},
\qquad
SPL=20\log_{10}\left(\frac{p_{\mathrm{rms}}}{p_{\mathrm{ref}}}\right),
\tag{6.61}
$$

其中空气声学常取 $p_{\mathrm{ref}}=20\,\mu\mathrm{Pa}$，水声学常取 $1\,\mu\mathrm{Pa}$。采样时间必须覆盖足够周期，频谱分析还应说明窗函数、频率分辨率和泄漏控制。

### 6.8.3 能量与守恒检查

线性声场的扰动能量可近似写为

$$
E_a=\frac12\int_{\Omega_a}
\left(\frac{p'^2}{\rho_0c^2}+\rho_0|\mathbf u'|^2\right)\mathrm d\Omega.
\tag{6.62}
$$

结构机械能包括动能和应变能。对于无外加功、无阻尼的理想结构声系统，两者总和应只在界面间转移；若采用显式松耦合，能量漂移往往比位移误差更早暴露时间层不一致。

流固热耦合还需分别检查热量和机械功，不能把热能与声能不加区分地合并。人工边界会带走物理上离开计算域的能量，PML 的吸收属于预期出流，不应误判为数值不守恒；但内部区域与 PML 界面的能通量仍应与吸收层衰减相符。

### 6.8.4 常见错误

1. 把水动力压力 $P$ 与声压 $p'$ 混为一谈，或在结构载荷中重复加入总压力；
2. 流声分解中忽略 $\partial P/\partial t$，却仍直接使用 $p'=c^2\rho'$；
3. 用一阶迎风长时间传播声波，导致幅值和相位严重衰减；
4. 把零梯度或固定压力边界当成普适无反射边界；
5. PML 吸收系数在内部界面不为零，或辅助变量与主变量时间层不同步；
6. 流固界面只匹配压力大小，不检查黏性力、总力矩和虚功；
7. 网格速度未满足 GCL，引起自由流压力振荡；
8. 非匹配网格位移和力采用互不共轭的映射，产生虚假界面能量；
9. 共轭传热中两侧使用不同法向却仍直接比较带符号热流；
10. 自由均匀热膨胀产生非零热应力，通常是参考温度或边界条件错误；
11. 结构声耦合只施加声压载荷，遗漏声学侧的法向加速度边界；
12. 紧耦合只检查子域残差，不检查界面残差；
13. 强流固耦合出现发散时只盲目减小结构时间步，不分析附加质量和分区稳定性；
14. 多速率插值使用未来值或不一致外推，造成非因果更新和能量漂移；
15. 只展示收敛后的位移或声压云图，没有网格、时间步和守恒验证。

## 本章小结

流—固—声耦合的核心不是简单串联多个求解器，而是在每个时间层同时满足子域方程、界面运动学、界面动力学和守恒条件。流声分解把变量分为不可压背景量和声学扰动量，声学 SIMPLE 通过动量预测、声压修正和动量修正满足声学连续性。Tam 边界利用特征波与远场辐射关系，PML 则通过坐标伸缩和辅助变量构造吸收层。

流固耦合采用 ALE 相对对流速度描述运动网格，并必须满足 GCL。界面位移、速度和面力在非匹配网格上还要满足总力、总力矩和虚功共轭。共轭传热通过温度与热流条件连接流体和固体能量方程，热应变再进入结构本构和等效载荷。结构声耦合则通过声压力平衡和法向加速度连续形成双向系统，其耦合矩阵具有明确的转置关系。

单向、松耦合、紧耦合和整体耦合各有适用条件。轻质结构—稠密流体等强耦合问题需要隐式分区、Aitken 或 quasi-Newton 加速，必要时采用整体求解。多速率方法可以匹配不同波速和网格尺度，但必须保证同步、时间插值和界面功率一致。最终可信度来自分层验证、界面残差、能量平衡以及网格和时间步收敛，而不是某个单一收敛指标。

## 习题

1. 说明流—固—声系统中至少四条耦合路径，并给出一个可以安全采用单向耦合的例子。
2. 从可压缩连续性方程和式（6.6）推导式（6.8），说明使用了哪项不可压背景场条件。
3. 推导式（6.10），并解释为什么 $P$ 随时间变化时不能始终令 $p'=c^2\rho'$。
4. 在静止均匀介质和小振幅条件下，从流声分解方程推导经典波动方程。
5. 从式（6.15）—（6.18）推导声压修正方程式（6.19），指出 SIMPLE 近似忽略了哪一项。
6. 比较一阶迎风、TVD 和低耗散高阶中心格式用于长距离声传播时的色散与耗散。
7. 根据特征速度式（6.21）说明亚声速出流边界为什么不能对所有变量施加同一种条件。
8. 推导二维和三维辐射条件中的几何扩散系数，并讨论声源中心选择对边界误差的影响。
9. 从式（6.26）和式（6.27）估算给定层厚与目标反射系数下的 $\sigma_{\max}$，并说明还需进行哪些数值试验。
10. 从均匀流保持要求推导 GCL 式（6.31），说明面扫过体积比任意时刻 $\mathbf U_g\cdot\mathbf A$ 更可靠的原因。
11. 分别写出无滑移和理想滑移流固界面的运动学条件，并推导两侧面力的符号关系。
12. 证明若位移和力映射满足虚功共轭，则常界面面力下可保持总功；讨论如何构造相应数值试验。
13. 推导共轭传热界面在“共同法向”和“各子域外法向”两种写法下的热流守恒条件。
14. 说明自由均匀热膨胀为何不应产生热应力，并分析完全约束均匀加热时的应力方向。
15. 从线性声学 Euler 方程推导结构声界面条件式（6.47）。
16. 从界面虚功推导式（6.49）和式（6.50）中的 $\mathbf G_{sa}$ 与 $\rho_0\mathbf G_{sa}^{\mathrm T}$，说明转置关系的意义。
17. 解释附加质量为何会使串行显式流固耦合失稳，并比较减小步长、Picard 松弛和整体求解的作用。
18. 设计 $\Delta t_a=4\Delta t_s$ 的多速率交换时间表，标出每次数据交换对应的时间层。
19. 为一维结构—声学活塞问题设计能量守恒检查，列出结构功率、声学辐射功率和人工边界能通量。
20. 编写完整流—固—热—声计算的验证矩阵，至少包含自由流保持、平面波、CHT、热膨胀、流固耦合和结构声耦合六类算例。

## 延伸阅读

1. Patankar, S. V. *Numerical Heat Transfer and Fluid Flow*. Hemisphere, 1980.
2. Tam, C. K. W., and Webb, J. C. Dispersion-relation-preserving finite difference schemes for computational acoustics. *Journal of Computational Physics*, 1993, 107: 262–281.
3. Tam, C. K. W. Computational aeroacoustics: issues and methods. *AIAA Journal*, 1995, 33: 1788–1796.
4. Berenger, J. P. A perfectly matched layer for the absorption of electromagnetic waves. *Journal of Computational Physics*, 1994, 114: 185–200.
5. Hu, F. Q. On absorbing boundary conditions for linearized Euler equations by a perfectly matched layer. *Journal of Computational Physics*, 1996, 129: 201–219.
6. Donea, J., Huerta, A., Ponthot, J.-Ph., and Rodríguez-Ferran, A. Arbitrary Lagrangian–Eulerian methods. In *Encyclopedia of Computational Mechanics*, Wiley, 2004.
7. Felippa, C. A., Park, K. C., and Farhat, C. Partitioned analysis of coupled mechanical systems. *Computer Methods in Applied Mechanics and Engineering*, 2001, 190: 3247–3270.
8. Badia, S., Nobile, F., and Vergara, C. Fluid–structure partitioned procedures based on Robin transmission conditions. *Journal of Computational Physics*, 2008, 227: 7027–7051.
9. Zienkiewicz, O. C., Taylor, R. L., and Zhu, J. Z. *The Finite Element Method: Its Basis and Fundamentals*. 7th ed., Butterworth-Heinemann, 2013.
10. Bathe, K. J. *Finite Element Procedures*. 2nd ed., Prentice Hall, 1996.
11. Fahy, F. J., and Gardonio, P. *Sound and Structural Vibration: Radiation, Transmission and Response*. 2nd ed., Academic Press, 2007.
