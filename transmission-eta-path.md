# pbrt-v4 透射折射率因子（η）与路径权重、轮盘终止的相互作用

> 分析对象：本仓库 pbrt-v4（分支 `task/73`，提交 `20e41fd`）。
> 所有结论均标注源码位置；文末单独区分“源码可证明的数值行为”与“缺少实验依据的物理解释”。

---

## 1. 结论速览

1. pbrt-v4 中折射率对路径权重有**两条独立的作用通道**：
   - **BxDF 层**：透射 BTDF 的返回值 `f` 在**辐亮度（Radiance）模式**下额外除以相对折射率平方 `etap²`；在**重要性（Importance）模式**下不做该除法（`src/pbrt/bxdfs.cpp:108-114, 160-167, 203-205`）。
   - **积分器层**：相机路径（辐亮度传播）维护标量 `etaScale`，每次透射事件乘以 `Sqr(bs->eta)`，且**只用于俄罗斯轮盘（RR）的终止判据，从不参与最终辐射度 `L` 的累加**（`src/pbrt/cpu/integrators.cpp:744-758`）。
2. 两条通道在数值上恰好抵消：**辐亮度模式下 `beta * etaScale` 等于同一条路径以重要性模式传播时的 `beta`**（第 4 节给出数值验证）。因此 RR 实际依据的是“介质无关”的重要性通量，而相机路径的无偏估计量 `L += beta * ...` 保持不变。
3. 进入光密介质（η>1）时 `etaScale` 增大，使 RR 终止概率 `q` **降低**（路径更不易被杀）；离开光密介质时 `etaScale` 回落，`q` 回升。往返穿过同一界面后 `etaScale` 精确回到 1。
4. 全内反射（TIR）时透射采样直接失败返回空，路径在该分支终止；此时 `FrDielectric` 恒为 1，只能走反射分支，`etaScale` 不更新。
5. 镜面极限（`EffectivelySmooth()`）下 BTDF 是 Dirac δ，`f` 与 `pdf` 都含 `1/|cosθ|` 奇异项，积分器的 `f·|cos|/pdf` 化简为有限的 Fresnel 透射率（再按模式乘除 η²）；粗糙界面走微表面分支，公式不同但 η² 修正位置完全一致。

---

## 2. 涉及的源码位置

| 主题 | 文件与行号 |
|---|---|
| `TransportMode { Radiance, Importance }` 枚举 | `src/pbrt/base/bxdf.h:109-116` |
| `BSDFSample.eta` 字段（默认 1） | `src/pbrt/base/bxdf.h:125-151` |
| `Refract()`：Snell 定律、TIR 判定、返回相对 η | `src/pbrt/util/scattering.h:22-47` |
| `FrDielectric()`：Fresnel 透射/反射率，TIR 时返回 1 | `src/pbrt/util/scattering.h:61-79` |
| 光滑/粗糙电介质 BTDF 的 `ft /= Sqr(etap)`（仅 Radiance） | `src/pbrt/bxdfs.cpp:108-114, 160-167, 199-207` |
| BSSRDF 出口 BxDF 的反向 η² 修正（Radiance 下 `*= Sqr(eta)`） | `src/pbrt/bxdfs.h:1116-1128` |
| 薄电介质：透射样本 `eta` 用默认值 1（不更新 etaScale） | `src/pbrt/bxdfs.h:247-252` |
| 分层 BxDF 出口样本固定 `eta=1` | `src/pbrt/bxdfs.h:762-768` |
| CPU `PathIntegrator`：`etaScale` 更新与 RR | `src/pbrt/cpu/integrators.cpp:635, 739-758` |
| CPU `VolPathIntegrator`：含 `r_u` 归一化的 RR | `src/pbrt/cpu/integrators.cpp:960, 1183-1185, 1256-1268` |
| CPU SPPM 相机路径：带 0.05 下限的 RR | `src/pbrt/cpu/integrators.cpp:2893, 2978-2990` |
| GPU wavefront 表面/介质/次表面 RR | `src/pbrt/wavefront/surfscatter.cpp:205-222`；`media.cpp:326-339`；`subsurface.cpp:99-112` |
| 重要性传播（光路径、光子、BDPT）不维护 etaScale | `integrators.cpp:556-596`（LightPath）、`3156-3182`（SPPM 光子）、`1965-2113`（BDPT RandomWalk） |
| 射线跨介质后 medium 的切换 | `src/pbrt/interaction.h:117-121`；`src/pbrt/interaction.cpp:99-102` |

---

## 3. 机制逐层分析

### 3.1 BxDF 层：相对折射率 `etap` 如何进入 `f`

电介质 BxDF（`DielectricBxDF`）持有构造时给定的标量 `eta`（`src/pbrt/bxdfs.h:204`，材质默认 1.5，见 `src/pbrt/materials.cpp:60, 92`）。真正传给积分器的是 `Refract()` 计算出的**沿射线方向的相对折射率 `etap`**：

```cpp
// src/pbrt/util/scattering.h:22-47（节选）
Float cosTheta_i = Dot(n, wi);
if (cosTheta_i < 0) { eta = 1 / eta; cosTheta_i = -cosTheta_i; n = -n; }
Float sin2Theta_t = sin2Theta_i / Sqr(eta);
if (sin2Theta_t >= 1) return false;          // 全内反射
...
if (etap) *etap = eta;                       // 进入光密体 etap=η；出射 etap=1/η
```

即：**同一块玻璃，从空气进入时 `etap = 1.5`，从内部出射时 `etap = 1/1.5 ≈ 0.667`**，由入射点的几何朝向自动翻转，无需积分器知道当前介质。

光滑（镜面）透射分支（`src/pbrt/bxdfs.cpp:98-115`）：

```cpp
SampledSpectrum ft(T / AbsCosTheta(wi));          // δ-BTDF 的 1/|cosθ_t|
if (mode == TransportMode::Radiance)
    ft /= Sqr(etap);                               // 仅辐亮度模式
return BSDFSample(ft, wi, pt / (pr + pt),
                  BxDFFlags::SpecularTransmission, etap);
```

粗糙微表面透射分支（`bxdfs.cpp:145-168`）与求值路径 `DielectricBxDF::f`（`bxdfs.cpp:199-207`）在相同位置做同样的 `ft /= Sqr(etap)`。反射分支的 `fr` 不含任何 η² 因子（`bxdfs.cpp:92-96, 132-143`）。

注意分支选择概率 `pr = R, pt = T` 本身**不含 η² 修正**，两种传播模式使用同一套采样概率。

### 3.2 积分器层：`etaScale` 只服务于轮盘判据

CPU `PathIntegrator::Li` 的核心三行（`src/pbrt/cpu/integrators.cpp:739-758`）：

```cpp
beta *= bs->f * AbsDot(bs->wi, isect.shading.n) / bs->pdf;   // 路径通量
...
if (bs->IsTransmission())
    etaScale *= Sqr(bs->eta);                                // 累积 η² 乘积
...
SampledSpectrum rrBeta = beta * etaScale;                    // RR 专用量
if (rrBeta.MaxComponentValue() < 1 && depth > 1) {
    Float q = std::max<Float>(0, 1 - rrBeta.MaxComponentValue());
    if (sampler.Get1D() < q) break;                          // 以概率 q 终止
    beta /= 1 - q;                                           // 存活则补偿
}
```

可由源码直接确认的事实：

- `etaScale` 初值为 1（`integrators.cpp:635`），**只在 `IsTransmission()` 时更新**，反射事件不改变它。
- `etaScale` **从不出现于任何 `L += ...` 语句**；辐射度估计只用 `beta`。因此它对无偏估计量的期望值没有影响，只影响路径存活分布与方差。
- RR 仅在 `depth > 1` 后启用；`q = max(0, 1 - max_c(rrBeta[c]))`，取光谱各通道最大分量（`SampledSpectrum::MaxComponentValue`，`util/spectrum.h:249`）。
- 存活路径执行 `beta /= 1-q`，这是标准 RR 无偏补偿；注意补偿的是 `beta` 而非 `rrBeta`，`etaScale` 不被缩放，继续沿路径累积。

GPU wavefront 路径与体积积分器逻辑相同，只是因为采用 hero-wavelength + 重采样（ratio tracking），判据多除以路径概率均值 `r_u.Average()`：

```cpp
// src/pbrt/wavefront/surfscatter.cpp:206-222
Float etaScale = w.etaScale;
if (bsdfSample->IsTransmission()) etaScale *= Sqr(bsdfSample->eta);
SampledSpectrum rrBeta = beta * etaScale / r_u.Average();
if (rrBeta.MaxComponentValue() < 1 && w.depth >= 1) { ... }
```

`VolPathIntegrator`（`integrators.cpp:1259`）、wavefront 介质散射（`media.cpp:330`）、次表面出口（`subsurface.cpp:104`）公式一致。SPPM 相机路径额外给 `q` 加了 0.05 的下限（`integrators.cpp:2986`，`max(.05f, 1 - rrBeta...)`）。

### 3.3 进出介质时权重与终止概率的关联（数值验证）

以下数值由独立脚本按仓库公式（`FrDielectric`/`Refract`/`ft`）复算，η=1.5、`|cosθ_i|=0.9`（进入与出射取相同入射角模值，仅符号不同）：

| 量 | 进入玻璃（etap=1.5） | 离开玻璃（etap=1/1.5） |
|---|---|---|
| Fresnel R / T | 0.040792 / 0.959208 | 0.046333 / 0.953667 |
| 折射角 \|cosθ_t\| | 0.956847 | 0.756637 |
| Radiance 模式 `ft`（含 `/etap²`） | 0.445541 | 2.835905 |
| Importance 模式 `ft` | 1.002468 | 1.260402 |
| `etaScale` 乘子 | ×2.25 | ×0.4444 |
| `ft_radiance × etap²` | **1.002468 = ft_importance** | **1.260402 = ft_importance** |

恒等式（每个透射事件上成立，故整条路径上也成立）：

```
beta_radiance · etaScale  ≡  beta_importance
```

法向入射往返（空气→玻璃→空气）的完整复算：

```
after enter: beta_rad=0.426667  etaScale=2.25  beta_rad*etaScale=0.960000 = beta_imp(0.96)
after exit:  beta_rad=0.921600  etaScale=1.00  beta_rad*etaScale=0.921600 = beta_imp(0.9216)
```

由此得到终止概率的关联（设入射前 `beta=b0=0.5`，cos=0.9，进入玻璃）：

- **进入光密介质**：`ft` 被 `/η²` 压低，裸 `beta` 变小，朴素 RR 会给出偏大的 `q`；`etaScale ×η²` 把判据拉回到重要性通量水平。该例中透射后 `beta=0.222771`，`q` 从朴素的 0.777 降到 0.499——即**相机路径进入玻璃后被轮盘误杀的概率显著降低**；且 0.499 恰好等于重要性模式同路径的 RR 概率。
- **离开光密介质**：`etap²<1`，`etaScale` 回落，`ft` 中 `1/etap²` 的放大被抵消，`q` 回升。
- **往返闭合**后 `etaScale=1`，介质的存在对 RR 判据不留痕迹。

### 3.4 辐亮度传播 vs. 重要性传播

仓库中两条传播方向的处理是对称设计的：

| | 相机路径（Radiance） | 光源路径 / 光子 / BDPT 光子子路径（Importance） |
|---|---|---|
| BxDF 透射修正 | `ft /= etap²` | 不除（`bxdfs.cpp:110, 164, 204` 的 mode 判断） |
| 积分器 etaScale | 维护，仅用于 RR | **不维护**（`LightPathIntegrator` 无该变量，`integrators.cpp:553-596`；BDPT `RandomWalk` 同样没有，`1965-2113`） |
| RR 判据 | `beta · etaScale`（体积版再 `/r_u`） | SPPM 光子用局部比值 `bnew/beta`（`integrators.cpp:3176-3180`）；LightPath/BDPT 不做 RR，靠 `maxDepth` 截断 |
| 显式 mode 传参 | 默认 `TransportMode::Radiance`（`base/bxdf.h:182`、`bsdf.h:59`） | `Sample_f(..., TransportMode::Importance)`（`integrators.cpp:591, 2093, 3169`） |

BDPT 连接时两个方向的 BxDF 分别用相反 mode 求值（`integrators.cpp:2352, 2399, 2414-2415`），保证连接项的量纲一致。

特例——BSSRDF 出口的 `NormalizedFresnelBxDF` 符号相反：它在 **Radiance 下 `f *= Sqr(eta)`**（`bxdfs.h:1123-1125`），因为它建模的是“从介质内部出射”的项，η 修正的方向与界面 BTDF 相反。

### 3.5 全内反射（TIR）

两条代码路径共同保证 TIR 时不存在透射样本：

1. `FrDielectric` 在 `sin²θ_t ≥ 1` 时直接返回 1（`scattering.h:72-73`），于是 `T=0, pt=0`；
2. `Refract` 在同一条件下返回 `false`（`scattering.h:36-37`），光滑分支 `if (!valid) return {};`（`bxdfs.cpp:104-106`），粗糙分支 `if (... || tir) return {};`（`bxdfs.cpp:149-152`）。

数值核对（玻璃→空气，η=1.5，临界角 arcsin(1/1.5)=41.81°）：45°、60° 入射时 `FrDielectric=1` 且 `Refract` 失败；40° 时 R≈0.245、折射有效。TIR 发生时：

- 透射 BSDF 样本为空 → 积分器 `if (!bs) break;`（`integrators.cpp:736-737`），**该路径终止且不更新 `etaScale`、不补偿**（这是 Fresnel 物理的确定性结果，不是随机终止，不引入偏差）；
- 若采样器先选到反射分支（`pr=R=1` 时必然如此），则按普通镜面反射继续，`etaScale` 不变。

### 3.6 镜面极限与粗糙界面

- **镜面极限**：`eta==1` 或 `mfDistrib.EffectivelySmooth()` 走 δ 分支（`bxdfs.cpp:80`）。`fr=R/|cosθ_i|`、`ft=T/|cosθ_t|` 与离散概率 `pr/(pr+pt)` 配套；积分器乘上 `AbsDot(wi,n)/pdf` 后 `1/|cos|` 奇异项抵消，得到有限权重 `R` 或 `T/etap²`（Radiance）。`eta==1` 时 `Flags()` 只报 Transmission（`bxdfs.h:177`），但 `Refract` 给出 `etap=1`，`etaScale` 实际不变。
- **粗糙界面**：微表面半矢量 `wm` 分支，`ft` 含 Jacobian 项 `|wi·wm · wo·wm|/(cosθ_i cosθ_o (wi·wm + wo·wm/etap)²)`（`bxdfs.cpp:154-162`），η² 修正仍在最后一行；返回的 `etap` 同样喂给积分器的 `etaScale`。
- **薄电介质 `ThinDielectricBxDF`**：透射方向直接取 `-wo`（`bxdfs.h:249`），`BSDFSample` 不传 eta → 字段默认 1（`base/bxdf.h:150`），所以虽然 `IsTransmission()` 为真，`etaScale *= 1²` 不变。这与“薄板两侧介质相同、无净折射率跃迁”的建模一致。
- **分层 BxDF `LayeredBxDF`**：内部随机游走自带独立 RR（`bxdfs.h:543-549, 697-703`，阈值 0.25、`depth>3`，**只用局部 beta，无 etaScale**）；对积分器而言出口样本固定 `eta=1`（`bxdfs.h:768`），外层 `etaScale` 同样不更新——界面内部的 η 效应被封装在 BxDF 求值内部。

---

## 4. 源码可证明的数值行为（清单）

以下各项均可由引用的代码与按代码复算的脚本直接验证，不依赖渲染实验：

1. `etaScale` 是各次透射 `etap²` 的连乘积，初值 1，反射不改变它，且**只出现在 RR 判据中，不进入辐射度累加**（`integrators.cpp:744-758`）。
2. Radiance 模式透射 `f` 带 `/etap²`，Importance 模式不带；两种模式的 Fresnel `R/T` 与采样概率 `pr/pt` 完全相同（`bxdfs.cpp:82-115`）。
3. 恒等式 `beta_rad · etaScale = beta_imp` 在逐事件层面成立（cos=0.9 入射脚本验证：1.002468=1.002468；法向往返后 `beta_rad·etaScale=0.9216=beta_imp`，`etaScale` 回到 1）。
4. RR 公式 `q=max(0,1-max_c(rrBeta[c]))`，门限 `rrBeta<1` 且 `depth>1`（wavefront 为 `depth>=1`，注释说明是为对齐 VolPath 的 depth 自增时机，`surfscatter.cpp:213-214`）。
5. 进入 η=1.5 介质、cos=0.9、`b0=0.5` 的例子中，etaScale 把 RR 终止概率从 0.777 修正到 0.499，与重要性模式的 RR 概率一致。
6. TIR 时 `FrDielectric` 精确返回 1、`Refract` 精确返回 false，透射样本为空，路径确定性终止，`etaScale` 不更新。
7. 临界角由 `sin²θ_i/η² ≥ 1` 判定；η=1.5 时为 41.8103°，脚本在 45°/60° 复现 TIR、40° 复现折射。
8. `ThinDielectricBxDF`、`LayeredBxDF` 出口、`DiffuseTransmissionBxDF`（`bxdfs.h:125`，eta 默认 1）均不改变 `etaScale`；只有 `DielectricBxDF` 的真实折射分支返回非 1 的 `etap`。
9. 重要性传播的三处实现（LightPath、BDPT RandomWalk、SPPM 光子）均无 `etaScale` 变量；其中前两者无 RR，SPPM 光子 RR 用相邻两代 beta 的比值而非绝对通量。
10. RR 存活补偿 `beta /= 1-q` 只改 `beta`；体积/wavefront 版判据中的 `/r_u.Average()` 是 hero 波长重采样概率的归一化，与 η 因子互相独立、乘法可交换。

## 5. 尚缺实验依据的物理解释（需谨慎对待）

以下是对“为什么这样设计”的物理/数学解释，方向与 pbrt 教材（PBRTv4 第 4 章、第 13 章）的标准论述一致，但**本仓库内没有测试或实验数据直接证实**，`bsdfs_test.cpp`/`integrators_test.cpp` 也没有针对 η-RR 耦合的断言：

1. **η² 的来源是立体角测度在两种介质间的变换**：Snell 定律下 `dω_t = (η_t²/η_i²)(cosθ_t/cosθ_i) ... ` 的度量换算使 BTDF 天然非互易；Radiance 传输除以 `etap²`、而伴随（重要性）传输不除，是 Veach 框架中 radiance/adjoint 量纲约定的结果。代码注释只写了 “Account for non-symmetry with transmission to different medium”（`bxdfs.cpp:109, 163, 203`），未给出推导。
2. **RR 依据重要性通量可降低方差**：解释是重要性通量对应“伴随粒子数”式的守恒量，在折射率界面处不因测度换算而人为收缩/膨胀，用它做存活判据可避免进入玻璃的路径被系统性过早终止。但仓库没有 A/B 实验（例如把 `rrBeta` 改成裸 `beta`）来量化方差或偏差差异；无偏性本身由 RR 补偿保证，与判据选择无关——这一点是数学结论，不是实验结论。
3. **`MaxComponentValue` 取最大通道**的物理动机（保证任一光谱通道都不被过度终止）属于常见工程解释，源码无注释说明。
4. **临界角附近的数值行为**（`sin2θ_t` 逼近 1 时 `cosθ_t→0` 导致 `ft=T/|cosθ_t|` 发散）在理论上由 δ-BSDF 的测度处理吸收，`CHECK_RARE(1e-5f, !valid)`（`bxdfs.cpp:104, 150`）只做罕见事件统计；仓库没有数值稳定性专项测试证明极限处无有限精度问题。
5. BSSRDF 中 `NormalizedFresnelBxDF` 的反向 η² 符号（乘而非除）的物理推导依赖扩散近似边界条件，仓库仅留 “account for adjoint light transport” 一句注释（`bxdfs.h:1123`），无独立验证。

> 建议的验证手段（仓库当前不包含）：构造单侧玻璃平板场景，分别启用/禁用 `etaScale` 项对比 RR 存活率与收敛曲线；在临界角附近扫描入射角统计 `Refract` 失败率；用 BDPT 与 PathIntegrator 互验含嵌套介质场景的像素均值。
