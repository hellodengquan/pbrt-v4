# pbrt-v4 几何法线与着色法线：法线贴图、次级光线与介质边界分析

> 分析对象：`src/pbrt`（pbrt-v4，提交 `261a042` 工作区），CPU 与 GPU/wavefront 两条路径。
> 约定：**几何法线** = `Interaction::n`（下文记 `ng`）；**着色法线** = `SurfaceInteraction::shading.n`（下文记 `ns`）。
> 本文只做静态源码分析。凡标注【未验证】的条目均为基于代码的视觉/数值推断，未通过实际渲染或单元测试确认。

---

## 0. 结论速览：每类约束依赖哪条法线

| 光传输约束 | 依赖法线 | 直接源码位置 |
|---|---|---|
| 自相交偏移（ray origin offset） | **几何法线 `ng`** | `interaction.h:92` → `ray.h:75-92` |
| 透射/反射后新光线的介质归属（inside/outside） | **几何法线 `ng`** | `interaction.h:117-121`；GPU: `wavefront/surfscatter.cpp:236-237,312-313` |
| BSDF 局部坐标系 / 半球方向（z 轴） | **着色法线 `ns`** | `bsdf.h:24-25` |
| 散射权重的余弦项 `AbsDot(wi, n)` | **着色法线 `ns`** | `cpu/integrators.cpp:439,453,473,739,2098` |
| 法线贴图 / bump 贴图 | 只改写 **`ns`**，`ng` 不变 | `interaction.cpp:175-188`；`materials.h:86-138` |
| 朝向翻转（reverseOrientation） | 同时翻转 `ng` 与初始 `ns`，之后二者保持同半球 | `interaction.h:179-182`；`shapes.h:936-938` |
| 遮挡/阴影求交（几何可见性） | 与法线**无关**，只看光线-三角形求交结果 | `cpu/integrators.h:52-53`；`cpu/primitive.cpp:50-77` |
| 均匀半球采样时的半球侧选择（SimplePathIntegrator） | **几何法线 `ng`** | `cpu/integrators.cpp:468-471` |
| 光源采样上下文 | 同时携带两者（`n` 与 `ns`），各取所需 | `lights.h:83-101`；`shapes.h:35-47` |

---

## 1. 两条法线的建立与分离

### 1.1 数据结构（源码证据）

`SurfaceInteraction` 同时保存两套几何（`interaction.h:251-258`）：

```cpp
Vector3f dpdu, dpdv;
Normal3f dndu, dndv;
struct {
    Normal3f n;
    Vector3f dpdu, dpdv;
    Normal3f dndu, dndv;
} shading;
```

基类 `Interaction::n`（`interaction.h:134`）即几何法线；`shading.n` 即着色法线。

### 1.2 几何法线来自真实面（源码证据）

三角形命中时（`shapes.h:936-938`）：

```cpp
isect.n = isect.shading.n = Normal3f(Normalize(Cross(dp02, dp12)));
if (mesh->reverseOrientation ^ mesh->transformSwapsHandedness)
    isect.n = isect.shading.n = -isect.n;
```

- `ng` 是由三角形两条边叉积得到的**面法线**，代表该表面片元真实的局部平面；
- 若网格/变换要求翻转朝向，`ng` 与 `ns` 一起取反。

### 1.3 平滑着色只改 `ns`（源码证据）

存在顶点法线 `mesh->n` 时，重心插值得到 `ns`，再调用
`isect.SetShadingGeometry(ns, ss, ts, dndu, dndv, true)`（`shapes.h:944-947,1006`）。

`SetShadingGeometry`（`interaction.h:194-214`）：

```cpp
shading.n = ns;
if (orientationIsAuthoritative)
    n = FaceForward(n, shading.n);      // 三角形平滑着色走这一支
else
    shading.n = FaceForward(shading.n, n);
```

- 三角形传入 `orientationIsAuthoritative = true`：允许调整 `ng` 的**符号**去贴合 `ns` 所在半球，但 `ng` 的方向向量本质不变（仅可能整体取反）；
- 法线贴图路径传入 `false`（见 §2），此时反而是 `ns` 被折回 `ng` 所在半球。
- `SetIntersectionProperties` 中有不变量断言 `CHECK_GE(Dot(n, shading.n), 0.)`（`interaction.h:223`）：**两条法线永远位于同一半球**，不会出现 `ns` 与 `ng` 夹角 ≥ 90° 的状态。

### 1.4 反向朝向（backfacing / reverseOrientation）如何处理（源码证据）

1. 场景层：`Shape` 构造记录 `reverseOrientation` 与 `transformSwapsHandedness`（如 `shapes.h:118-121`）；
2. 求交层：`flipNormal = reverseOrientation ^ transformSwapsHandedness`（`shapes.h:276,496,727,926`），在 `SurfaceInteraction` 构造函数中统一翻转（`interaction.h:179-182`）：

```cpp
if (flipNormal) {
    n *= -1;
    shading.n *= -1;
}
```

3. 工具层：`FaceForward(n, n2) = Dot(n,n2) < 0 ? -n : n`（`util/vecmath.h:1106-1108`），所有朝向对齐都通过点积符号完成，不依赖任何材质信息。

**含义**：反向朝向是几何/场景语义，在材质与法线贴图介入之前就确定了，且同时作用于 `ng`、`ns`。因此介质内外侧的判断不会因为模型被标记为反向而错乱——它读到的 `ng` 已经是翻转后的结果（见 §4）。

---

## 2. 法线贴图与 bump 贴图：只影响着色法线

### 2.1 调用点（源码证据）

`SurfaceInteraction::GetBSDF`（`interaction.cpp:175-188`）：

```cpp
const Image *normalMap = material.GetNormalMap();
if (displacement || normalMap) {
    Vector3f dpdu, dpdv;
    if (normalMap)
        NormalMap(*normalMap, *this, &dpdu, &dpdv);
    else
        BumpMap(..., displacement, *this, &dpdu, &dpdv);
    Normal3f ns(Normalize(Cross(dpdu, dpdv)));
    SetShadingGeometry(ns, dpdu, dpdv, shading.dndu, shading.dndv, false);
}
```

注意最后一个参数是 **`false`**：法线贴图不是朝向权威，`shading.n = FaceForward(shading.n, n)`（`interaction.h:201-202`），即贴图法线被强制折回几何法线所在半球。`ng` 在整个过程中不被修改。

### 2.2 NormalMap 的数学过程（源码证据，`materials.h:86-105`）

1. 采样 RGB，从 `[0,1]` 映射到切线空间 `[-1,1]` 并归一化（v 坐标做 `1-v` 翻转）；
2. 用**当前着色几何**（`shading.dpdu`、`shading.n`）构建切线空间 `Frame::FromXZ`，把贴图法线变换到渲染空间；
3. Gram–Schmidt 正交化重算 `dpdu/dpdv`，使 `Cross(dpdu,dpdv)` 恰好给出新的 `ns`。

BumpMap（`materials.h:108-138`）则通过对位移纹理做有限差分，扰动 `shading.dpdu/dpdv` 间接得到 `ns`，同样不触碰 `ng`。

### 2.3 GPU 路径的等价实现（源码证据）

`wavefront/surfscatter.cpp:111-130` 做同样的事，并显式 `ns = FaceForward(ns, w.n)`（`w.n` 即几何法线，由 `wavefront/intersect.h:73,132` 从 `intr.n` 填入）。CPU/GPU 两条路径在"贴图法线不得越过几何法线半球"这一点上一致。

### 2.4 法线贴图最终如何进入光传输（源码证据）

材质求值上下文只接收 `ns`（`materials.h:32-33`）：

```cpp
MaterialEvalContext(const SurfaceInteraction &si)
    : ..., ns(si.shading.n), dpdus(si.shading.dpdu) {}
```

BSDF 构造（`bsdf.h:24-25`）：

```cpp
BSDF(Normal3f ns, Vector3f dpdus, BxDF bxdf)
    : bxdf(bxdf), shadingFrame(Frame::FromXZ(Normalize(dpdus), Vector3f(ns))) {}
```

`shadingFrame.z == ns`。BxDF 的所有反射/折射计算都在这个局部坐标系中进行，例如镜面透射 `Refract(wo, Normal3f(0,0,1), eta, ...)`（`bxdfs.cpp:103`），粗糙透射以微面法线 `wm` 为折射轴（`bxdfs.cpp:149`）。也就是说，**法线贴图改变的是 BSDF 半球的朝向与余弦权重，而不是表面的几何位置或几何侧别**。

---

## 3. BSDF 半球与次级光线生成

### 3.1 半球定义与采样（源码证据）

- 局部坐标 z 轴 = `ns`：`bsdf.h:24-25`；
- `BSDF::Sample_f` 先把世界空间的 `wo` 变换进 `shadingFrame`，采样后再把 `wi` 变回世界空间（`bsdf.h:57-77`）；
- 反射/透射的判别在局部坐标中以 z 分量符号给出，例如 `SameHemisphere(wo, wi)`（`bxdfs.cpp:135,151`）、`CosTheta` 系列（`bxdfs.cpp:177-190`）；
- `DielectricBxDF` 同时声明反射与透射标志（`bxdfs.h:176-181`），采样结果带 `SpecularReflection` / `SpecularTransmission` 等标志（`bxdfs.cpp:96,113,143,167`）。

### 3.2 次级光线的生成分两步，分别用不同法线（源码证据）

**第一步：BSDF 决定方向 `wi`（用 `ns`）。**
路径权重的余弦项一律是 `AbsDot(bs->wi, isect.shading.n)`：

- `cpu/integrators.cpp:439`（直接光照）、`:453`（SimplePath BSDF 采样）、`:473`（均匀方向兜底）、`:739`（VolPath）、`:2098`（BDPT 随机游走）；
- GPU: `wavefront/surfscatter.cpp:191`。

**第二步：把方向变成光线（用 `ng`）。**

`Interaction::SpawnRay`（`interaction.h:98-100`）：

```cpp
RayDifferential SpawnRay(Vector3f d) const {
    return RayDifferential(OffsetRayOrigin(d), d, time, GetMedium(d));
}
```

两个关键决策都在基类 `Interaction` 中、用 `this->n`（即 `ng`）完成：

1. **起点偏移**：`OffsetRayOrigin(pi, n, w)`（`ray.h:75-92`）——
   - 偏移量 `d = Dot(Abs(n), pi.Error())`，沿 `n` 方向；
   - `Dot(w,n) < 0` 时反向偏移，保证起点沿出射方向位于表面误差盒之外；
   - 再用 `NextFloatUp/Down` 把浮点坐标向远离表面的方向舍入。
   - `ShapeSampleContext::OffsetRayOrigin`（`shapes.h:60-78`）是同一份算法，同样用几何法线。
2. **介质归属**：`GetMedium(d)`（`interaction.h:117-121`），见 §4。

阴影/连接光线 `SpawnRayTo` 同理：`ray.h:103-108` 对两端各用各自的**几何法线** `nFrom/nTo` 做偏移。`Integrator::Unoccluded` 发射的就是这种光线（`cpu/integrators.h:52-53`）。

### 3.3 一个容易混淆的例外（源码证据）

`SimplePathIntegrator` 不采样 BSDF、而在全局均匀半球上抽方向时，用**几何法线**判定该翻到哪一侧（`cpu/integrators.cpp:468-471`）：

```cpp
if (IsReflective(flags) && Dot(wo, isect.n) * Dot(wi, isect.n) < 0) wi = -wi;
else if (IsTransmissive(flags) && Dot(wo, isect.n) * Dot(wi, isect.n) > 0) wi = -wi;
```

但同段的权重余弦仍用 `isect.shading.n`（`:473`）。这再次体现分工：**几何侧别用 `ng`，着色响应量用 `ns`**。

### 3.4 光线微分（源码证据）

镜面反射/透射的光线微分计算使用 `shading.n` 及其导数 `shading.dndu/dndv`（`interaction.cpp:99-157`，透射处在 `:133-137` 还会按 `Dot(wo,n)<0` 翻转法线）。这只影响纹理过滤足迹，不影响求交与介质。

---

## 4. 介质边界：inside/outside 完全由几何法线裁决

### 4.1 介质接口的绑定（源码证据）

- 场景解析中 `MediumInterface(inside, outside)` 挂在图元上（`scene.cpp:1323-1324`，另见 `:1420,:1479`）；
- 结构定义（`base/medium.h:103-118`）：`IsMediumTransition() = inside != outside`；
- 命中后 `GeometricPrimitive::Intersect` 调用
  `SetIntersectionProperties(material, areaLight, &mediumInterface, r.medium)`（`cpu/primitive.cpp:74`）；
- `SetIntersectionProperties`（`interaction.h:218-229`）：仅当图元声明了介质转换时才保存 `mediumInterface`，否则新交点继承光线当前介质 `rayMedium`。

### 4.2 内外侧判断规则（源码证据）

`Interaction::GetMedium`（`interaction.h:117-121`）：

```cpp
Medium GetMedium(Vector3f w) const {
    if (mediumInterface)
        return Dot(w, n) > 0 ? mediumInterface->outside : mediumInterface->inside;
    return medium;
}
```

判断式是 **`Dot(光线方向, 几何法线 ng)` 的符号**：

- `Dot(w, ng) > 0`：光线与 `ng` 同向，朝"外侧"走 → `outside` 介质；
- `Dot(w, ng) < 0`：逆着 `ng`，朝"内侧"走 → `inside` 介质；
- 与 `ns`、法线贴图、bump 完全无关。

GPU/wavefront 路径把同一规则内联写出（`surfscatter.cpp:236-237` 间接光线、`:312-313` 阴影光线；`subsurface.cpp:118-120` 次表面出射）：

```cpp
ray.medium = Dot(ray.d, w.n) > 0 ? w.mediumInterface.outside
                                 : w.mediumInterface.inside;
```

### 4.3 透射时为什么不会判错侧别（源码证据 + 推理）

1. BxDF 在 `ns` 局部坐标系中用 Snell 定律生成透射方向，`Refract` 内部会按入射余弦符号自行翻转界面法线并返回相对 IOR `etap`（`util/scattering.h:22-47`）；全反射时返回 `false`，采样作废（`bxdfs.cpp:103-106,149-152`）。
2. 但透射方向 `wi` 变回世界空间后，光线进入哪一侧介质并不再询问 BxDF，而是统一走 `GetMedium(wi)` 的 `Dot(wi, ng)` 符号判断。由于 §1.3 的不变量 `Dot(ng, ns) ≥ 0`，BxDF 按 `ns` 判出的"穿过界面"方向，相对 `ng` 也必然是穿到另一侧——两条法线同半球保证了两套判定不会互相矛盾。
3. 反向朝向的表面：`ng` 在求交时已随 `flipNormal` 取反（§1.4），`inside/outside` 的语义随几何朝向一起翻转，因此 `GetMedium` 仍然给出与建模意图一致的结果。
4. 相机初始光线的介质由相机参数直接给定（`cameras.cpp:940` `ray.medium = medium`；场景侧 `scene.cpp:780`），不依赖任何法线；此后每次命中都按上述规则切换或继承。

### 4.4 透射的辐射度修正（源码证据，与法线问题相邻但独立）

`DielectricBxDF` 在 `TransportMode::Radiance` 下对透射结果除以 `etap²`（`bxdfs.cpp:110-111,164-165`），积分器侧再用 `etaScale *= Sqr(bs->eta)` 配合 RR/MIS（`integrators.cpp:744-745`；GPU `surfscatter.cpp:207-208`）。这是跨介质非对称性修正，使用的是 BxDF 返回的标量 IOR，与法线选择是两条独立机制。

---

## 5. 为什么平滑着色/法线贴图不能改变几何遮挡关系

### 5.1 遮挡由求交决定，而求交发生在着色之前（源码证据）

1. 光线求交只与真实三角形平面做数学相交：`Primitive::Intersect → shape.Intersect`（`cpu/primitive.cpp:50-52`）；遮挡测试 `IntersectP` 同样如此（`:79-83`）。
2. 阴影光线 `Unoccluded`（`cpu/integrators.h:52-53`）只消费求交布尔结果，整条路径不读取任何法线：

```cpp
bool Unoccluded(const Interaction &p0, const Interaction &p1) const {
    return !IntersectP(p0.SpawnRayTo(p1), 1 - ShadowEpsilon);
}
```

3. 法线贴图/bump 在 `GetBSDF` 中才执行（`interaction.cpp:159-188`），时间上晚于求交；且它只写 `shading.*`，从不修改 `pi`、`dpdu/dpdv`（几何那组）或网格数据。
4. 自相交防护偏移用的是 `ng`（§3.2），法线贴图对 `ns` 的任何修改都不会放大或缩小偏移方向——贴图无法把一条光线"推出"几何阴影区，也无法让被遮挡的光线越过实体。
5. alpha 测试是唯一能在求交后"忽略"某个命中的机制，但它由 alpha 纹理的随机判定驱动并继续追踪后续几何（`cpu/primitive.cpp:57-70`），与法线无关。

### 5.2 概念区分（基于源码的推理）

- `ng` 描述表面**真实的分片线性几何**：决定光线从哪一侧离开、介质如何切换、起点往哪个方向偏移、阴影光线是否被挡住；
- `ns` 是一个**着色代理量**：只决定 BSDF 局部坐标系、反射/折射方向的表观朝向、余弦权重。
- 平滑着色（顶点法线插值）与法线贴图都只生产 `ns`。因此它们可以让一个硬立方体看起来有圆润的高光，但立方体轮廓处的剪影、两个相邻片元之间的互相遮挡、阴影边缘仍然是硬的几何结果。

### 5.3 已知的副作用边界（源码证据 + 推理）

法线贴图虽然不能改遮挡，但通过 `ns` 仍会影响：
- BSDF 采样方向分布与 `f·cosθ` 权重（`bsdf.h:57-77`，`integrators.cpp:739`）；
- 光源采样的上下文 `LightSampleContext` 同时带 `n` 和 `ns`（`lights.h:88-89`），面积光把两者传给 `ShapeSampleContext`（`lights.cpp:745,766`），形状的朝向加权 PDF 可使用 `ns`；
- 镜面反射/透射光线微分（纹理足迹），见 §3.4。

【未验证】当 `ns` 与 `ng` 偏差很大时（极端法线贴图或低多边形平滑着色），BSDF 可能采样出指向几何表面"错误一侧"的方向。代码层面的保障仅有：(a) 两法线同半球断言（`interaction.h:223`）；(b) 起点按 `ng` 偏移使光线不会立即自交。但光线仍可能在短距离内命中同一表面的相邻几何，产生黑边/光泄漏类 artifact——pbrt-v4 源码中没有针对这种情况的额外钳制逻辑，这一视觉现象未在本次分析中渲染验证。

---

## 6. 自相交偏移机制汇总（源码证据）

`ray.h:75-92` 的 `OffsetRayOrigin` 是全系统唯一的几何偏移实现（`Interaction` 与 `ShapeSampleContext` 各自封装同算法）：

1. 取交点浮点误差界 `pi.Error()`，在几何法线方向上投影得偏移距离；
2. 按出射方向与 `ng` 的点积符号决定偏移正负；
3. 逐分量 `NextFloatUp/Down` 强制把新起点舍入到误差盒之外；
4. `SpawnRayTo` 对线段两端都做偏移（`ray.h:103-108`），阴影连接时双端防自交；
5. 阴影光线另有 `tMax = 1 - ShadowEpsilon`（`ShadowEpsilon = 1e-4`，`util/math.h:30,42`）作为参数化距离上的兜底。

整个机制的输入法线只有 `Interaction::n`（`ng`）。`ns` 不参与。

---

## 7. 尚未验证的视觉/数值效果清单

以下条目在源码中找到机制性依据，但本次未运行渲染或测试，**不作为已证实结论**：

1. **法线贴图裂缝/阴影渗漏**：§5.3 所述 `ns` 与 `ng` 大角度偏差时的 artifact 倾向，需构造对比场景（强法线贴图球体 + 面积光阴影）实拍确认。
2. **平滑着色的剪影保持**：代码上可断定遮挡关系不变，但视觉上"高光圆润而轮廓硬边"的具体表现未渲染验证。
3. **反向朝向 + 介质转换的实际正确性**：`flipNormal` 与 `GetMedium` 的组合在逻辑上自洽，但未用 inside/outside 不同介质（如玻璃球内嵌参与介质）的场景验证折射率与介质采样结果。
4. **偏移量在极端尺度场景下的充分性**：`pi.Error()` 依赖 `gamma(n)` 误差界（如 `shapes.h:929` 的 `gamma(7)`），超大世界坐标或超薄几何下是否仍杜绝自交，属于数值鲁棒性问题，未做压力测试。
5. **CPU 与 GPU 路径的像素级一致性**：两条路径法线/介质规则的源码逻辑一致（对比 `interaction.h:117-121` 与 `surfscatter.cpp:236-237`），但未实际比对输出图像；GPU 代码中 `subsurface.cpp:117` 还留有 `// TODO: should always just take outside in this case?` 注释，说明次表面出射的介质选择在作者处也存疑。
6. **法线贴图 v 坐标翻转（`1-v`，`materials.h:91`）与外部资产约定的匹配**：仅记录该约定存在，未验证具体资产的贴图朝向是否符合预期。
7. **BDPT/MMLT 中 `ns` 进入几何项 G**：`Vertex::ns()` 返回着色法线（`integrators.cpp:1620-1625`），`G()` 用 `AbsDot(v.ns(), d)`（`:2122-2125`），即连接策略的几何项采用着色法线；这对能量守恒/偏差的实际影响未做数值实验验证。

---

## 附：关键源码索引

| 主题 | 文件:行 |
|---|---|
| 两套法线字段 | `src/pbrt/interaction.h:134,251-258` |
| 朝向翻转构造 | `src/pbrt/interaction.h:179-182` |
| SetShadingGeometry 与同半球不变量 | `src/pbrt/interaction.h:194-214,223` |
| SpawnRay / GetMedium / OffsetRayOrigin 入口 | `src/pbrt/interaction.h:92-121` |
| 偏移与 SpawnRayTo 实现 | `src/pbrt/ray.h:75-108` |
| 三角形几何法线与顶点法线插值 | `src/pbrt/shapes.h:936-1006` |
| NormalMap / BumpMap | `src/pbrt/materials.h:86-138` |
| 法线贴图调用点 | `src/pbrt/interaction.cpp:175-188` |
| BSDF shadingFrame | `src/pbrt/bsdf.h:24-25,57-77` |
| Dielectric 反射/透射采样 | `src/pbrt/bxdfs.cpp:77-170` |
| Refract（界面法线自动翻转） | `src/pbrt/util/scattering.h:22-47` |
| FaceForward | `src/pbrt/util/vecmath.h:1101-1118` |
| 介质接口与转换判定 | `src/pbrt/base/medium.h:103-118` |
| 命中时介质属性装配 | `src/pbrt/cpu/primitive.cpp:74`；`src/pbrt/interaction.h:218-229` |
| 路径余弦权重（ns）与 SpawnRay | `src/pbrt/cpu/integrators.cpp:439,453,468-475,739,748` |
| 遮挡测试 | `src/pbrt/cpu/integrators.h:52-53` |
| GPU 法线贴图/介质/偏移 | `src/pbrt/wavefront/surfscatter.cpp:111-130,233-237,309-313` |
| GPU work item 法线来源 | `src/pbrt/wavefront/intersect.h:73-79,132-137` |
