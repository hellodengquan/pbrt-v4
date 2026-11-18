# 纹理足迹与过滤：pbrt-v4 实现分析

本文分析 pbrt-v4 中"斜视表面上的纹理足迹"如何从**光线微分（ray differentials）**经**表面参数导数（∂uv/∂screen）**传递到 **MIP/EWA 过滤器**，并解释三条退化路径：UV 接缝、缺失微分、以及过滤被禁用时的行为。所有结论均来自源码阅读，不依赖 GPU 跑分或外部标定数据。

## 1. 足迹传递链：从像素到纹素

一次带过滤的纹理查询经过四级传递，任何一级断裂都会改变足迹形状：

```
相机光线微分 (rxOrigin/rxDirection, ryOrigin/ryDirection)
   │  CameraBase::GenerateRayDifferential   cameras.cpp:116
   ▼
交点处屏幕空间导数 dpdx, dpdy
   │  SurfaceInteraction::ComputeDifferentials   interaction.cpp:41
   ▼
纹理参数导数 dudx, dudy, dvdx, dvdy（最小二乘解）
   │  interaction.cpp:67-88（CPU 积分器）
   │  wavefront/surfscatter.cpp:74-103（wavefront 积分器，内联同一份推导）
   ▼
映射后的 (s,t) 足迹 dst0={dsdx,dtdx}, dst1={dsdy,dtdy}
   │  UVMapping::Map / SphericalMapping::Map / ...   textures.h:95-167
   ▼
MIPMap::Filter(st, dst0, dst1) → Point / Bilinear / Trilinear / EWA
      util/mipmap.cpp:229
```

### 1.1 光线微分的产生与缩放

- 相机微分通过对胶片样本偏移 `±0.05` 像素做前向/后向差分得到（`cameras.cpp:124-147`）；两个方向都成功才置 `hasDifferentials`（`cameras.cpp:150`）。
- 积分器按采样率缩小足迹：`rayDiffScale = max(0.125, 1/sqrt(spp))`（`cpu/integrators.cpp:252-254`、`2520`、`2889`）。含义：spp>1 时每个样本只覆盖亚像素区域，足迹应收窄；`0.125` 的下限防止高 spp 下足迹缩到零而重新走样。
- 镜面反弹时由 `SurfaceInteraction::SpawnRay`（`interaction.cpp:99-157`）按反射/折射定律传播微分（含曲率项 `dndx/dndy`）；若微分长度平方超过 `1e16` 则直接**丢弃**（`interaction.cpp:150-154`），后续命中只能走 1.2 的近似回退。

### 1.2 交点处的两条求导路径

`ComputeDifferentials`（`interaction.cpp:41-89`）分两种情况：

1. **有光线微分**：把辅助光线与交点切平面求交得 `px, py`，`dpdx = px - p`。这是平面近似——对斜平面准确，对强曲率表面会低估/高估足迹。
2. **无光线微分**（`hasDifferentials == false` 或辅助光线与切平面平行）：回退到 `Camera::Approximate_dp_dxy`（`cameras.h:154-183`），用**全胶片最小微分**（`FindMinimumDifferentials` 预计算）在相机空间构造替代光线再与切平面求交，并乘 `sppScale = max(0.125, 1/sqrt(spp))`。

随后用 `dpdu, dpdv` 对 `dpdx, dpdy` 做最小二乘（法方程 `AᵀA`），解出 `∂u/∂x` 等四个导数（`interaction.cpp:67-82`），并 clamp 到 `±1e8`、非有限值归零（`interaction.cpp:84-88`）。

注意 **wavefront 积分器不走光线微分**：`wavefront/surfscatter.cpp:77-103` 无条件使用 `Approximate_dp_dxy` 回退路径，即 GPU/wavefront 模式下每个命中点的足迹都来自"相机最小微分 + 切平面投影"的近似，而非逐光线传播的真实微分。这是 CPU 与 wavefront 路径在足迹精度上的结构性差异。

### 1.3 映射层的解析雅可比

`UVMapping` 只做缩放平移，微分线性传递（`textures.h:97-101`）。`SphericalMapping`/`CylindricalMapping`（`textures.h:118-167`）则计算**解析的** `ds/dp`、`dt/dp` 再与 `dpdx` 点积——这是接缝处不退化的关键设计，见 §2。

## 2. UV 接缝与缺失微分时的足迹退化

### 2.1 UV 接缝

pbrt-v4 的 `∂uv/∂screen` 是在**命中点局部**用 `dpdu/dpdv` 解出来的，不是跨像素对 `uv` 做差分，因此微分本身在 UV 图集（atlas）边界两侧各自成立，不会因为相邻像素落在另一张 chart 上而产生 `±1` 量级的假导数。但足迹退化仍有三处来源：

1. **足迹跨接缝的采样语义由 WrapMode 决定**（`util/image.h:96-146` 的 `RemapPixelCoords`）：
   - `Repeat`：纹素坐标取模，EWA/bilerp 的核跨 `u=1→0` 边界取到正确的环绕纹素，过滤在接缝处连续；
   - `Clamp`：核被边界截断，接缝处过滤宽度实际变窄（偏清晰但能量不守恒）；
   - `Black`：越界纹素返回 0，接缝附近整体变暗；
   - 对**图集纹理**（多张 chart 打包在一张图里、本不该环绕），跨边界的核会采到相邻 chart 的纹素——经典的 atlas bleeding。pbrt 不在渲染时处理，只能靠资产侧留 padding；MIP 层级越深，box 下采样（`util/image.cpp:353-363` 的 2×2 平均）把相邻 chart 平均得越厉害，所以高层级（远处、斜视角）渗色最严重。
2. **金字塔构造不跨接缝混合**：`GeneratePyramid` 先把非 2 幂分辨率**向上**取整（`FloatResizeUp`，`image.cpp:319-322`），之后每级只做局部 2×2 box，不在环绕方向闭环。对 `Repeat` 纹理，深层的接缝两侧纹素从未被预平均，跨接缝过滤完全依赖采样时核的环绕读取——这在 EWA 下是正确的，但意味着深层 MIP 在接缝处的"预过滤"假设（核内纹素已含邻域均值）在接缝法向上不成立。
3. **球面/柱面映射的坐标跳变**：`SphericalMapping` 的 `s = atan2` 在 ±π 处从 1 跳回 0（`textures.h:137`）。由于微分是解析雅可比（`dsdp` 在 3D 空间连续，`textures.h:124-133`），足迹**宽度**不退化；配合 `Repeat` 包裹，取值也连续。若资产把 `wrap` 误设为 `clamp`，跳变处会出现一条被边缘色拉伸的接缝线。

### 2.2 缺失微分

微分在以下情况缺失或被清零，足迹按确定性方式退化：

| 触发点 | 后果 |
|---|---|
| 相机微分生成失败（`cameras.cpp:150`） | `hasDifferentials=false` → 回退 `Approximate_dp_dxy` |
| 镜面微分爆炸被 squash（`interaction.cpp:150-154`） | 同上，后续所有反弹失去真实微分 |
| `AᵀA` 行列式退化（球极点处 `dpdu→0`），`invDet` 非有限 → 置 0（`interaction.cpp:71-72`） | `dudx=dudy=dvdx=dvdy=0` → **零足迹** |
| `--disable-texture-filtering` / 场景 `Option "disabletexturefiltering"`（`interaction.cpp:43-47`、`scene.cpp:505`） | 强制零足迹 |

**零足迹在 MIPMap 内的行为**（`mipmap.cpp:229-284`）：

- 非 EWA：`width=0` → `level = nLevels-1 + Log2(1e-8)` 为很大的负数 → `iLevel = max(0, floor(level)) = 0` → 在**第 0 层（全分辨率）**做点采样或双线性。即退化为逐点查询：无模糊，但在 minification 区域产生混叠（moiré、闪烁）。
- EWA：`shorterVecLength == 0` → 直接 `Bilerp(0, st)`（`mipmap.cpp:276-277`），同样退化为全分辨率双线性。

也就是说"缺失微分"不会得到错误的颜色均值，而是得到**未过滤的全分辨率样本**——方差以走样形式回到图像里，由像素采样（spp）去吸收。这解释了为什么 `disableTextureFiltering` 主要用于调试走样来源，而不是加速手段。

## 3. 纹理分辨率（MIP 层级）选择与各向异性限制的取舍

### 3.1 层级选择公式

- 非 EWA（`mipmap.cpp:232-239`）：`width = 2·max(|dst0|, |dst1|)`（取四个微分分量的最大绝对值，各向同性近似），`level = nLevels-1 + Log2(max(width, 1e-8))`。`width` 是 [0,1] 纹理坐标下的足迹宽，`Log2(width)≤0` 从顶层向下数。`level ≥ Levels()-1` 时直接返回 1×1 顶层的平均色（`mipmap.cpp:237-238`）——超宽足迹收敛到纹理均值。
- EWA（`mipmap.cpp:264-283`）：先交换使 `dst0` 为长轴；**用短轴长度选层**：`lod = max(0, Levels()-1 + Log2(shorterVecLength))`，保证短轴方向约 1 纹素/核，长轴方向的各向异性由椭圆核在同一层内覆盖，最后在相邻两层 EWA 结果间 lerp（三线性式过渡，`mipmap.cpp:282-283`）。

分辨率取舍的含义：

- 基础图被向上取整到 2 幂（`image.cpp:319-322`），非 2 幂资产会付出最多 4× 的内存（金字塔总量约为取整后基图的 4/3），且取整时的重采样已经引入一次低通——基图有效分辨率略低于原图。
- 层级公式按 `[0,1]` 归一化坐标计算，因此**提高纹理分辨率不会改变选中层的相对位置**，只是把同一足迹映射到更深的绝对层；混叠/模糊的拐点始终由"足迹宽 ≈ 1 纹素"决定。
- `Trilinear` 与 `EWA` 都在层间 lerp，避免层切换处的条带；`Point`/`Bilinear` 用 `floor(level)`，跨层边界有可见的不连续。

### 3.2 各向异性限制（maxAnisotropy）

EWA 前对椭圆偏心率的钳制（`mipmap.cpp:269-275`）：

```cpp
if (shorterVecLength * maxAnisotropy < longerVecLength && shorterVecLength > 0) {
    Float scale = longerVecLength / (shorterVecLength * maxAnisotropy);
    dst1 *= scale;            // 放大短轴
    shorterVecLength *= scale;
}
```

- 默认 `maxAnisotropy = 8`（`MIPMapFilterOptions`，`mipmap.h:41`；场景文件 `"maxanisotropy"` 参数默认同为 8，`textures.cpp:418/452`）。
- 钳制方式是**放大短轴**而不是截断长轴：代价是短轴方向过度模糊（斜视地面上垂直于视线方向变糊），收益是 (a) 选层 `lod` 随短轴变长而上移，EWA 扫描的椭圆面积（`mipmap.cpp:334-347` 的 O(长短轴乘积) 双重循环）被压到 `maxAnisotropy` 倍纹素量级，过滤成本有硬上界；(b) 过瘦椭圆在 LUT 半径截断下采样不足的问题被缓解。
- 调大 `maxanisotropy`：斜视角度下沿各向异性方向更锐利，代价是 EWA 内层循环纹素数近似线性增长，且短轴选层更深、更接近走样边界。调小：性能稳定但斜视纹理发灰。
- **GPU 路径上限不同**：`texDesc.maxAnisotropy = Clamp(maxAniso, 1, 16)`（`textures.cpp:1298`），由硬件各向异性过滤执行，且 `"ewa"` 在 GPU 上只是映射为 `cudaFilterModeLinear` 的 MIP 过渡 + 硬件 aniso（`textures.cpp:1301-1304`），并非真正的椭圆加权平均——CPU 的 EWA 与 GPU 的 aniso-trilinear 在斜视足迹上保真度不同，同一场景两种后端不保证逐像素一致。

### 3.3 EWA 核本身

- 椭圆系数 `A, B, C` 由两个微分向量的外积形式给出并加 `+1` 正则化（`mipmap.cpp:314-320`），保证即使微分退化椭圆也不会奇异；
- 权重是 α=2 的截断高斯 LUT（128 项，`mipmap.cpp:59-191`），`r²<1` 截断；
- 扫描椭圆包围盒、逐纹素查 LUT 加权（`mipmap.cpp:322-348`）——成本随足迹面积增长，这正是需要 `maxAnisotropy` 和 MIP 选层共同封顶的原因。

## 4. 过滤保真 / 逐点查询 / 禁用过滤：三条路径的对照

| 路径 | 触发方式 | 足迹处理 | 结果特征 |
|---|---|---|---|
| **完整过滤** | 默认；`"filter" "ewa"/"trilinear"/"bilinear"` | 微分 → 最小二乘 → 映射 → MIP/EWA | minification 区域无走样；EWA 保各向异性，trilinear 各向同性近似会在斜视时整体过糊 |
| **逐点查询** | `"filter" "point"`（`mipmap.cpp:241-246`） | **仍用足迹宽选 MIP 层**，但在选中层取最近纹素 | 无混叠（层级正确）但无核内平均：纹理内部细节呈块状，动画下层级间跳变（无层间 lerp） |
| **禁用纹理过滤** | `--disable-texture-filtering` 或场景 Option（`interaction.cpp:43-47`、`wavefront/surfscatter.cpp:77`） | 微分全部清零 → 第 0 层双线性/点采样 | 纹理仍读取、颜色正确，但 minification 走样全部暴露；用于诊断与对照 |
| **禁用图像纹理** | `--disable-image-textures`（`mipmap.cpp:199-203`） | 金字塔被替换为只剩顶层 1×1 | 每张纹理退化为单一平均色，几何/光照调试用；注意它在 **MIPMap 构造期**生效，与逐查询的 `disableTextureFiltering` 是正交的两个开关 |

另外两个容易混淆的"部分保真"路径：

- **`FilterFunction::Bilinear`**：选层逻辑与 trilinear 相同（`mipmap.cpp:236` 的 `level` 照常计算），只是不做层间 lerp——它是"单层的正确层级查询"，不是逐点查询。
- **Ptex 纹理**（`textures.cpp:644-683`）：同一组 `dudx…dvdy` 微分直接交给 Ptex 的 bspline 过滤器，CPU 限定（GPU 下 `LOG_FATAL`，`textures.cpp:699-701`）；Ptex 按 face 过滤，天然没有 UV 接缝问题，但微分缺失时同样退化为零足迹。

## 5. 要点回顾

1. 足迹精度的瓶颈在**微分传播**而非过滤器本身：wavefront 路径恒用相机最小微分近似，镜面微分超界即被丢弃，这两处决定了斜视表面足迹的下限精度。
2. UV 接缝处微分不爆炸（局部最小二乘 + 解析雅可比），退化只发生在**采样语义**层：wrap 模式决定核跨接缝的行为，图集渗色随 MIP 层级加深而加剧。
3. 零微分一律退化为第 0 层采样——走样而非错误均值；`disableTextureFiltering` 与 `disableImageTextures` 分别是"查询期清零足迹"与"构造期塌缩金字塔"两个不同机制。
4. 分辨率选择由"短轴 ≈ 1 纹素"驱动，`maxAnisotropy` 通过放大短轴同时封顶 EWA 成本与模糊程度；GPU 后端以硬件 aniso（≤16）近似 EWA，与 CPU 的 LUT 高斯椭圆滤波不保真一致。
