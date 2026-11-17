# CPU 路径积分中光源采样、BSDF 采样与发光表面命中的 MIS 分析

> 分析对象：`src/pbrt/cpu/integrators.cpp` 中的 `SimplePathIntegrator`、`PathIntegrator`、
> `VolPathIntegrator` 与 `BDPTIntegrator`，配合 `bsdf.h`、`bxdfs.h`、`lights.h`、
> `util/sampling.h` 的实现。
>
> 核心问题：同一条辐射度传输路径可以由 **光源采样（NEE, next event estimation）**、
> **BSDF 采样随机游走命中**、以及 **射线直接命中发光表面/环境** 三种方式产生贡献，
> 本文说明 pbrt-v4 如何保证这些贡献不被重复计数，delta（Dirac）事件下采样概率与
> 方向密度的边界如何处理，并区分无偏成立的条件与当前实现并未保证无偏的场景。

---

## 1. 三类贡献在代码中的位置

以 `PathIntegrator::Li()`（`integrators.cpp:628-762`）为主线，每个路径顶点上有三个
潜在的辐射度来源：

| 贡献来源 | 代码位置 | 估计量形式 |
|---|---|---|
| 光源采样（NEE） | `SampleLd()`，`integrators.cpp:764-805`，调用点 `:724-730` | `β · f(wo,wi)·cosθ · L_i / p_l`，再乘 MIS 权重 `w_l` |
| BSDF 采样延续路径 | `integrators.cpp:732-748` | `β *= f·cosθ / p_b`，下一段射线自行寻找交点 |
| 发光表面命中（面积光） | `integrators.cpp:663-676` | `β · Le`，乘 MIS 权重 `w_b`（或全权重 1） |
| 环境/无限远光命中（逃逸射线） | `integrators.cpp:644-660` | `β · Le(ray)`，乘 MIS 权重 `w_b`（或全权重 1） |

其中：

- `p_l = lightSampler.PMF(ctx, light) · light.PDF_Li(ctx, wi, true)`
  （光源选择概率 × 光源上点/方向的采样密度），见 `integrators.cpp:797`、`:670-671`。
- `p_b` 是 BSDF 对该方向的立体角密度；当 `pdfIsProportional` 为真时重新求值精确 PDF
  （`integrators.cpp:740`），保证进入 MIS 的密度是真实密度而非成比例的近似。
- 功率启发式权重（`util/sampling.h:72-77`）：

```cpp
Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
    Float f = nf * fPdf, g = ng * gPdf;
    if (IsInf(Sqr(f))) return 1;
    return Sqr(f) / (Sqr(f) + Sqr(g));
}
```

---

## 2. 避免重复计数的基本机制：MIS 分区恒等式

### 2.1 非 delta 情形：两个策略共享同一测度

考虑一个非镜面顶点 `x` 看向面积光源上点 `y` 的贡献 `f(x)·Le(y)·G(x,y)`。它有两条
生成途径：

1. **光源策略**：在 `x` 处用 NEE 显式采样 `y`，密度为 `p_l`（面积测度经
   `LightLiSample` 转换后与 BSDF 密度同处于立体角/面积测度体系），权重
   `w_l = PH(1, p_l, 1, p_b)`（`integrators.cpp:801-803`）。
2. **BSDF 策略**：在 `x` 处按 `p_b` 采样方向，射线恰好命中发光面 `y`，权重
   `w_b = PH(1, p_b, 1, p_l)`（`integrators.cpp:668-674`）。

功率启发式满足代数恒等式

```
w_l + w_b = p_l²/(p_l²+p_b²) + p_b²/(p_b²+p_l²) ≡ 1   （当 p_l, p_b > 0）
```

即同一份传输贡献被两个策略按权重分割，**期望中恰好计数一次**。环境光的处理完全同构
（`integrators.cpp:651-656`），只是对 `infiniteLights` 逐个光源分别计算 `w_b`；
NEE 侧每个光源的 `p_l` 内含 `lightSampler.PMF`，两侧在"按光源分别"的意义上仍满足
分区恒等式。每顶点只采样一个光源不破坏无偏性：BSDF 策略"选择光源"是通过几何命中
隐式完成的，其方向密度 `p_b` 已天然包含对任意目标光源的密度。

### 2.2 门控：什么时候发光命中取全权重 1

`integrators.cpp:665` 与 `:648` 的条件 `depth == 0 || specularBounce`：

- **`depth == 0`（相机射线直接命中光源）**：不存在"前一个着色顶点"可以执行 NEE，
  光源策略对该贡献的采样密度为 0，唯一能生成它的是相机/BSDF 途径 → 权重必须为 1。
- **`specularBounce == true`（上一跳是镜面反射/折射）**：见第 3 节，delta BSDF
  不存在普通方向密度，NEE 在数学上无法沿该 δ 方向建立连接 → 只有 BSDF 采样能生成
  这条贡献，权重为 1。

`SimplePathIntegrator`（无 MIS 的教学版，`integrators.cpp:389-482`）用更粗暴但等价
的方式实现同一分区：`:411-412` 处 `if (!sampleLights || specularBounce)` 才累加发光
命中——当 NEE 开启且上一跳非镜面时，发光面贡献**整体丢弃**（权重 0/1 的硬分区），
完全交给 NEE；镜面链之后则只能靠命中获取。这仍然无偏，但方差显著高于 MIS 版本。

### 2.3 delta 光源：NEE 取全权重 1

`SampleLd()` 中（`integrators.cpp:798-799`）：

```cpp
if (IsDeltaLight(light.Type()))
    return ls->L * f / p_l;        // 无 MIS，权重 1
else {
    Float p_b = bsdf->PDF(wo, wi);
    Float w_l = PowerHeuristic(1, p_l, 1, p_b);
    return w_l * ls->L * f / p_l;
}
```

`PointLight`（`lights.h:231-233`）与 `DistantLight`（`lights.h:258-260`）的
`PDF_Li()` 恒返回 0：点光源在面积测度下是 δ 分布，方向光在方向测度下是 δ 分布，
BSDF 随机游走以连续方向密度**不可能**以正概率命中它们，因此 BSDF 策略对该贡献的
密度为 0，`PH` 中 gPdf=0 时 `w_l` 本就应为 1——代码用显式分支避免了无意义的
PDF 求值与数值问题。对称地，delta 光源没有可被射线相交的面积，发光表面命中路径
（`:663`）对它们天然不触发。

### 2.4 NEE 只在非镜面顶点执行

`integrators.cpp:724`：`if (IsNonSpecular(bsdf.Flags()))`。纯 specular BSDF 的
`f(wo,wi)` 对任意给定方向都是 0（δ 函数只能在自己的采样例程中通过测度变换给出
有限的 `f·cos/pdf`），NEE 即使采样光源也会得到 `f=0`，跳过是正确且必要的。
对**混合 BSDF**（如 `LayeredBxDF`，flags 同时含 `Specular|Diffuse`，见
`bxdfs.h:456-474`）：`IsNonSpecular` 为真，NEE 正常执行，而 `bsdf->f()` 与
`bsdf->PDF()` 只计非 specular 分量，NEE 与命中两侧的密度口径一致，分区恒等式仍然成立。

---

## 3. Delta 事件下采样概率与方向密度的边界

### 3.1 数学边界

delta BxDF 的散射分布含 Dirac δ：`f(wo,wi) = ρ · δ(wi - wi*(wo)) / |cosθ_i|`。
它**不存在普通的方向概率密度函数**——对任意连续方向求值 `PDF(wo,wi)=0`，
`f(wo,wi)=0`；δ 只能在采样例程内部通过测度变换消去，返回离散的分支概率与
有限的加权吞吐量。delta 光源同理。因此 MIS 的标准前提（两策略在同一连续测度上
互有正密度）在 delta 事件处失效，必须把含 δ 的策略隔离为"权重 1 的唯一策略"。

### 3.2 代码中的具体边界行为

**Specular BxDF（`bxdfs.h`）**

- `ThinDielectricBxDF::f()` 恒返回 0（`:217-219`），`PDF()` 恒返回 0（`:256-259`）。
- `Sample_f()` 返回离散分支概率 `pr/(pr+pt)` 与 `R/|cosθ|` 形式的采样值
  （`:241-252`）；光滑导体同理（`ConductorBxDF`，`:301-305`，pdf=1）。
- 光滑微面元的判定阈值是 `TrowbridgeReitzDistribution::EffectivelySmooth()`：
  粗糙度低于阈值时 **glossy 被当作 specular δ 处理**——这是一个离散化的建模边界，
  意味着"接近镜面"的材质在该阈值处发生估计量结构的突变。

**积分器侧的隔离机制**

- `specularBounce = bs->IsSpecular()`（`integrators.cpp:742`），一路传递到下一次
  发光命中的门控判断（`:665`、`:648`）。
- 介质边界（`GetBSDF` 返回空对象）被显式当作 specular 处理
  （`integrators.cpp:681-685`），注释明确写明 "disable MIS if the indirect ray
  hits a light"：边界本身不发生散射，发光归属应由另一侧顶点的 NEE/命中逻辑承担。

**BDPT 中的 δ 处理（`integrators.cpp:2097-2207`）**

- 镜面顶点标记 `vertex.delta = true` 且 `pdfRev = pdfFwd = 0`（`:2105-2108`）。
- `remap0(f) = f != 0 ? f : 1`（`:2136`）：比值中遇到 0 密度时替换为 1，防止
  0/0；配合 `delta` 标记，**凡是涉及 delta 顶点的连接策略都不进入 `sumRi` 累加**
  （`:2190`、`:2198-2201`），等价于该策略权重取 1、其他竞争策略取 0。
- `Vertex::IsConnectible()` 对 specular 表面顶点返回 false（`:1653-1654`），
  禁止在 δ 顶点上做显式连接（连接方向匹配 δ 支撑集的概率为零）。
- `PowerHeuristic` 自身带 `IsInf` 保护（`sampling.h:74-75`），密度比溢出时退化为
  权重 1。

**Delta 光源**

- `LightType::DeltaPosition`（点光、聚光等）与 `DeltaDirection`（远距光）由
  `IsDeltaLight()` 统一判定（`lights.h:37-39`），NEE 侧走全权重分支
  （`integrators.cpp:798`），体积积分器中同样分支（`:1386-1387`）。

---

## 4. 无偏所需条件 vs. 当前实现未保证的场景

### 4.1 无偏成立所依赖的条件

1. **测度一致与分区恒等式**：在所有 `p_l>0 且 p_b>0` 的方向上 `w_l+w_b=1`；
   密度必须在同一测度下表达（面积光的 `PDF_Li` 内含 `cosθ/r²` 转换，
   `Vertex::ConvertDensity`，`integrators.cpp:1731-1743`）。
2. **δ 事件隔离**：delta BxDF / delta 光只由唯一能采样它的策略承担，权重 1，
   竞争策略密度视为 0（第 2、3 节的门控与 BDPT `delta` 机制）。
3. **真实密度进入 MIS**：`pdfIsProportional` 时重新调用 `bsdf.PDF()` 求精确值
   （`integrators.cpp:740`），不能把成比例近似密度直接塞进权重。
4. **Russian roulette 补偿**：存活路径除以 `1-q`（`:752-757`，
   VolPath 版 `:1263-1268`），终止期望不改变估计量均值。
5. **采样与评估口径一致**：混合 BSDF 的 NEE 只与非 specular 分量配对（2.4 节）；
   光源选择 PMF 计入 `p_l`。
6. **不完整 PDF 的补偿闭环**：BVH 光源采样器以 `allowIncompletePDF=true` 工作时，
   采不到的光源方向 `PDF_Li=0`（如 `UniformInfiniteLight`，`lights.cpp:959-960`、
   `970-971`），对应贡献改由发光命中侧以全权重 1 收取——两策略支撑集不重叠时
   权重 1 仍是合法的 MIS 分区，期望无偏；`ImageInfiniteLight` 的
   `compensatedDistribution`（`lights.cpp:1034-1040`、`1048-1049`）即按此补偿
   关系构造。

在上述条件全部满足时，`PathIntegrator` 对无限深度、无限采样的极限是无偏的。

### 4.2 当前实现未保证无偏（有意为之或固有截断）的场景

| 场景 | 代码位置 | 性质 |
|---|---|---|
| **最大深度硬截断** | `integrators.cpp:720-721`（`if (depth++ == maxDepth) break;`） | 超出 `maxDepth` 的路径贡献直接丢弃，无 RR 式补偿 → **有偏**。这是全局光照通用的偏差/方差折中 |
| **BSDF 正则化 `regularize`** | `integrators.cpp:712-715`，`bsdf.Regularize()` | 故意把光泽/镜面微面元分布展宽以抑制萤火虫 → 估计的是被修改后的 BSDF，**有意有偏** |
| **`SimplePathIntegrator` 无 MIS** | `integrators.cpp:411-412`、`:428-444` | 估计量本身仍无偏（硬分区 0/1 也是合法 MIS），但小光源、间接光泽反射等情形方差极大；它不保证收敛效率，不保证均值之外的任何东西 |
| **`SimpleVolPathIntegrator`** | `integrators.cpp:822-935` | 教学型 delta tracking，无 MIS、不支持表面散射（`:928-931` 直接报错），构造函数拒绝 delta 光源（`:828-831`）；适用范围之外不提供正确性保证 |
| **体积 MIS 的数值边界** | `integrators.cpp:1059`（`ClampZero(sigma_maj - σ_a - σ_s)`）、`:1348-1353` | `r_u/r_l` 重加权依赖 majorant 不小于真实消光系数；约束被违反时负项被钳为 0，理论闭环被数值修补替代。单通道（`[0]`）决策 RR 与事件模式（`:1004-1011`、`:1259-1265`）在光谱强变化时引入近似 |
| **`EffectivelySmooth()` 阈值** | `bxdfs.h:291` 等 | 跨阈值时 glossy→specular 的估计量结构突变；阈值之下连续密度被强制为 0，属于离散建模近似 |
| **几何/数值量** | `ShadowEpsilon`（`integrators.h:52-54`）、有限浮点、RR 在极小吞吐量下的截断 | 自相交规避与浮点精度带来量级很小的系统性偏差；`beta` 归零即终止路径（`:1257-1258`） |
| **有限 spp 的期望本身** | 全文件 | 单次估计噪声不是偏差，但 `maxDepth`、正则化等使"无限 spp 极限"也不等于渲染方程真解的场景，才是严格意义上的有偏 |

另需注意一个**无偏但易误判**的点：`PathIntegrator` 命中发光面时若
`lightSampler.PMF` 对该光源返回 0（该光源不在采样器的可采样集合中），则 `p_l=0`，
`w_b = PH(p_b, 0) = 1`，贡献全由命中侧收取——这正是 4.1 第 6 条的补偿闭环，
不是 bug；只有当同一贡献在两侧**同时**以正权重被收取时才构成重复计数，而
`w_l+w_b≡1` 从代数上排除了这种可能。

---

## 5. 结论

1. **重复计数的避免**依靠 MIS 分区恒等式：非 delta 传输在 NEE（`w_l`）与 BSDF
   命中（`w_b`）之间按功率启发式分割，两者之和恒为 1；delta BxDF、delta 光源、
   相机首命中、介质边界等不存在共同连续测度的情形，通过 `specularBounce` 门控、
   `IsDeltaLight` 分支和 BDPT 的 `vertex.delta`/`remap0` 机制隔离给唯一可行策略，
   权重取 1。
2. **delta 事件的密度边界**是：δ 分布没有普通 PDF（求值为 0），只在自身采样例程
   中以离散分支概率 + 测度变换给出有限吞吐量；任何 MIS 权重都不能跨 δ/连续两种
   测度直接构造，代码对此的全部特判（门控、跳过连接、0→1 重映射、Inf 保护）都是
   这一数学事实的直接后果。
3. **无偏性是条件性的**：测度一致、δ 隔离、真实密度、RR 补偿、不完整 PDF 闭环等
   条件满足时，无截断的路径追踪无偏；而 `maxDepth` 硬截断与 `regularize` 是当前
   实现明确接受的偏差来源，教学版积分器与体积重加权的数值钳制则只在限定模型与
   majorant 条件下保证正确性。
