# 从波长采样到胶片颜色：光谱 → 传感器 → 胶片的责任链分析

> 分析对象：pbrt-v4 本仓库中的 `src/pbrt/util/spectrum.{h,cpp}`、`util/sampling.h`、
> `film.{h,cpp}`、`base/film.h`、`filters.{h,cpp}`、`base/sampler.h`、`samplers.h`、
> `cpu/integrators.cpp`、`wavefront/{camera,film}.cpp`、`util/color.{h,cpp}`、
> `util/colorspace.{h,cpp}`、`util/image.cpp`、`scene.cpp`。
>
> 核心问题：一条相机路径在路径追踪阶段始终携带 **4 个波长上的光谱辐射度**
> （`SampledSpectrum`，`NSpectrumSamples = 4`，`util/spectrum.h:36`），它如何经过
> 波长 PDF 反除、RGB 传感器响应、曝光成像比、像素滤波权重累积，最终变换为输出色彩空间
> 中的线性 RGB。本文同时区分 **蒙特卡洛噪声（方差）**、**确定性的能量偏差（clamp 等）**
> 与 **色彩映射（传感器矩阵/白平衡/色彩空间）** 三类责任的归属，并说明默认传感器、
> `VisibleSurface` 与采样权重的边界。所有结论只引用本仓库代码与仓库内已有的测试、注释，
> 不引入外部标定数据，也不包含未执行的实验。

---

## 0. 责任链总览

```
                          每像素每个采样（per camera sample）
 Sampler.Get1D() ──► Film.SampleWavelengths(u)
                      │  RGBFilm/GBufferFilm: SampleVisible  (按可见度重要性采样)
                      │  SpectralFilm:        SampleUniform  (360–830 nm 均匀)
                      ▼
            SampledWavelengths { lambda[4], pdf[4] }
                      │
 Camera.GenerateRay ──┴──► 路径积分 Li(ray, lambda) 全程光谱运算
                      │      （README: RGB 只用于场景描述与最终输出，传输一律光谱）
                      ▼
            SampledSpectrum L[4]  × cameraRay.weight
                      │
 Film.AddSample(pPixel, L, lambda, visibleSurface, filterWeight)
                      │
   ┌──────────────────┴───────────────────────────────────────┐
   │ PixelSensor.ToSensorRGB:                                 │
   │   L ← L / pdf                       （光谱 MC 估计量反 PDF）│
   │   rgb_c = imagingRatio · mean( r̄_c(λ)·L )  c∈{R,G,B}    │
   │   可选 maxComponentValue 钳制（确定性偏差）               │
   └──────────────────┬───────────────────────────────────────┘
                      ▼  传感器 RGB（仍为线性、场景参照量）
   pixel.rgbSum[c] += filterWeight · rgb[c]
   pixel.weightSum  += filterWeight
                      │  渲染结束 GetPixelRGB / GetImage
                      ▼
   rgb = rgbSum / weightSum                         （像素内归一化，去滤波权重）
   rgb += splatScale · rgbSplat / filterIntegral    （SPPM/光路径 splat 项）
   rgb = outputRGBFromSensorRGB · rgb               （传感器 RGB → XYZ → 输出 RGB）
                      │
                      ▼  线性 HDR（默认 Half/EXR），无色调映射、无 OETF
                 Image.Write(filename, metadata)
```

整条链路由四个职责清晰的阶段组成：

| 阶段 | 入口代码 | 拥有的状态/责任 | 不负责的事 |
|---|---|---|---|
| ① 波长采样 | `FilmBase::SampleWavelengths` `film.h:211-214`；`SampledWavelengths::SampleVisible/SampleUniform` `spectrum.h:283-343` | 决定 `lambda[4]` 与每个波长的 `pdf[4]` | 不决定辐射度大小，不做颜色变换 |
| ② 光谱传输 | `RayIntegrator::EvaluatePixelSample` `cpu/integrators.cpp:227-289`；wavefront 对应 `wavefront/camera.cpp:53-79` | 在 4 个波长上做路径积分，产出 `SampledSpectrum L` | 不把光谱折叠成 RGB |
| ③ 传感器响应 + 曝光 | `PixelSensor::ToSensorRGB` `film.h:94-100`；`PixelSensor::Create` `film.cpp:213-253` | 反 PDF、与 `r̄/ḡ/b̄` 响应曲线内积、乘成像比、可选钳制 | 不做像素间加权平均，不做输出色彩空间变换 |
| ④ 滤波累积 + 色彩空间输出 | `RGBFilm::AddSample` `film.h:238-255`、`GetPixelRGB` `film.h:257-274`、`GetImage` `film.cpp:533-565` | 滤波权重累积、权重和归一化、splat、传感器→输出 RGB 矩阵、写文件 | 不改变样本的光谱估计本身（除显式钳制） |

---

## 1. 光谱路径权重如何进入 RGB 传感器响应，再到曝光、滤波与色彩空间输出

### 1.1 波长采样：两种 PDF，决定后续“反 PDF”的分母

每个相机样本先取一个一维随机数 `lu`（`cpu/integrators.cpp:230-233`；
`--disable-wavelength-jitter` 时固定为 0.5，`integrators.cpp:231-232`），交给胶片：

```cpp
// film.h:211-214（FilmBase，RGBFilm / GBufferFilm 继承此默认实现）
SampledWavelengths SampleWavelengths(Float u) const {
    return SampledWavelengths::SampleVisible(u);
}
```

**(a) 可见度重要性采样（RGB/GBuffer 胶片默认）。**
`SampledWavelengths::SampleVisible`（`spectrum.h:331-343`）对 `i = 0..3` 取
`up = frac(u + i/4)`，调用 `SampleVisibleWavelengths(up)` 与 `VisibleWavelengthsPDF`
（`util/sampling.h:163-171`）：

```cpp
Float VisibleWavelengthsPDF(Float lambda) {          // sampling.h:163
    if (lambda < 360 || lambda > 830) return 0;
    return 0.0039398042f / Sqr(std::cosh(0.0072f * (lambda - 538)));
}
Float SampleVisibleWavelengths(Float u) {            // sampling.h:169
    return 538 - 138.888889f * std::atanh(0.85691062f - 1.82750197f * u);
}
```

这是一个以 538 nm 为中心、按人眼可见度（与 CIE 匹配曲线形状相关的 sech² 形式）偏置的
密度，支持域为 `[360, 830]` nm（即常量 `Lambda_min..Lambda_max`，`spectrum.h:34`）。
4 个波长是对同一密度的分层样本（stratified），各自携带自己的 `pdf[i]`。
仓库测试 `spectrum_test.cpp:123-159`（`SamplingPdfY`、`SamplingPdfXYZ`）验证了用
`Y(λ)/pdf(λ)` 与 `(X+Y+Z)/pdf` 做 MC 积分能恢复 `CIE_Y_integral`，即该 PDF 与其逆函数
在仓库内是自洽且归一的。

**(b) 均匀采样（SpectralFilm 覆盖默认实现，`film.h:407-410`）。**
`SampleUniform`（`spectrum.h:283-303`）把首个波长在 `[lambda_min, lambda_max]` 上等距
分层铺开，4 个波长等间距、环绕取模，PDF 恒为 `1/(lambda_max-lambda_min)`
（默认区间同样是 360–830 nm）。SpectralFilm 的光谱桶正是依赖“均匀 PDF + 等宽桶”这一
性质，才可以不反 PDF 而直接按桶平均（见 `film.h:443-467` 的注释与实现）。

**(c) 次级波长终止。** 参与色散/介质 null-scattering 的积分器可调用
`TerminateSecondary()`（`spectrum.h:312-320`）：把 `pdf[1..3]` 置 0、`pdf[0]` 除以
`NSpectrumSamples`。此后只有 `lambda[0]` 的辐射度有效（例如
`SimpleVolPathIntegrator::Li` 在随机游走开始处调用，`integrators.cpp:842-843`）。
这是对**估计量权重的重新归一化**，不是丢弃能量——SpectralFilm 在桶平均时会给
`lambda[0]` 补回 `NSpectrumSamples` 倍权重（`film.h:453-460`），并有专门测试
`film_test.cpp:97-154` 覆盖。

### 1.2 路径传输：光谱量在 4 个波长上独立传播

相机射线由 `Camera::GenerateRay(Differential)` 生成，`lambda` 以引用传入并贯穿整条路径
（`base/camera.h:46-51`）。`CameraRay` 自带一个光谱权重 `weight`，默认构造为
`SampledSpectrum(1)`（`cameras.h:116-126`）；普通 `PerspectiveCamera::GenerateRay`
直接返回默认权重 1（`cameras.cpp:404-427`），而 `RealisticCamera` 会把镜头系统的
透过率、`cos⁴θ` 与出瞳采样 PDF 折算成每个波长的权重（`cameras.cpp:934-948`）。

积分器返回的 `L` 是 4 波长辐射度，CPU 端再乘相机权重：

```cpp
// cpu/integrators.cpp:259-260
L = cameraRay->weight * Li(cameraRay->ray, lambda, sampler, scratchBuffer, ...);
```

wavefront 端等价地在更新胶片时相乘（`wavefront/film.cpp:22-23`，
`Lw = L · cameraRayWeight`）。README“Features/Spectral rendering”明确：传输计算始终
使用点采样光谱，RGB 只出现在场景描述（纹理等）与最终图像输出两端。因此**路径阶段不存在
RGB 噪声或 RGB 偏差，只有 4 个波长通道各自的 MC 估计误差**。

### 1.3 进入传感器：先反 PDF，再与响应曲线内积

`Film::AddSample` 分发到具体胶片。RGB 胶片的第一步是把光谱样本交给像素传感器
（`film.h:238-242`）：

```cpp
// film.h:94-100  PixelSensor::ToSensorRGB
RGB ToSensorRGB(SampledSpectrum L, const SampledWavelengths &lambda) const {
    L = SafeDiv(L, lambda.PDF());                                   // (1) 反波长 PDF
    return imagingRatio * RGB((r_bar.Sample(lambda) * L).Average(),// (2) 响应曲线内积
                              (g_bar.Sample(lambda) * L).Average(),
                              (b_bar.Sample(lambda) * L).Average());
}
```

两个数学动作：

1. **反波长采样 PDF**：`L_i ← L_i / pdf_i`。这是重要性采样估计量的标准分母，使得无论
   用可见度偏置还是均匀采样，期望值都对应光谱积分而不是“按采样密度加权”的积分。
   `SafeDiv`（`spectrum.h:631-636`）在 `pdf_i = 0`（次级波长被终止）时返回 0，与
   “该波长无有效估计”一致。
2. **与传感器三条响应曲线 `r̄, ḡ, b̄`（`DenselySampledSpectrum`，`film.h:112`）做
   蒙特卡洛内积并对 4 个样本取平均**：`mean(r̄(λ_i)·L_i)`。这是对
   `∫ r̄(λ)·L(λ) dλ` 的 4 样本估计（响应曲线本身在 `Spectra::Init` 中由 CIE 数据表
   建成，`spectrum.cpp:2582-2590`）。

作为对照，`SampledSpectrum::ToXYZ`（`spectrum.cpp:205-216`）用完全相同的
`SafeDiv(X·L, pdf).Average() / CIE_Y_integral` 模式计算 CIE XYZ，`ToRGB` 再经
`RGBColorSpace::ToRGB` 转换（`spectrum.cpp:224-228`）。这条 `ToXYZ/ToRGB` 路径用于
反照率等辅助量（见 1.6 与第 3 节），**相机主路径不经过它，而是经过 `PixelSensor`**。

### 1.4 曝光：成像比 imagingRatio 是唯一的曝光乘子

曝光不在积分器或相机里，而在传感器创建时被烘焙为一个标量 `imagingRatio`
（`film.cpp:216-233`）：

```cpp
Float ISO = parameters.GetOneFloat("iso", 100.);           // film.cpp:220
Float whiteBalanceTemp = parameters.GetOneFloat("whitebalance", 0); // :221
std::string sensorName = parameters.GetOneString("sensor", "cie1931"); // :223
...
Float imagingRatio = exposureTime * ISO / 100;            // film.cpp:233
```

- `exposureTime` 来自相机快门：`shutterclose - shutteropen`，默认 `1 - 0 = 1`
  （`scene.cpp:756-757`，快门默认值 `cameras.cpp:274-275`）。
- 默认 `ISO = 100`，因此**默认成像比恰为 1**。源码注释称这是“passthrough”设定，保证
  未按物理相机编写的场景渲染结果与预期一致（`film.cpp:217-219, 230-232`）。
- `imagingRatio` 在 `ToSensorRGB` 返回值处整体乘到三个通道（`film.h:97-99`），即
  曝光是对线性传感器信号的**等比缩放**，不改变色度，只改变亮度量级。

### 1.5 像素滤波权重：样本空间的重建，不是颜色处理

每个相机样本在胶片上的子像素位置与滤波权重由 `GetCameraSample` 生成
（`samplers.h:797-814`）：`Filter::Sample` 返回 `FilterSample{p, weight}`，
其中 `weight = f(x)/pdf(x)`（`filters.h:33-38` 的通用 `FilterSampler`；Box 滤镜
退化为常数 1，`filters.h:67-70`）。默认滤镜是 Gaussian（`scene.cpp:94`，半径 1.5、
σ=0.5，`filters.cpp:40-47`）。

RGB 胶片把“传感器 RGB × 滤波权重”累加，并单独累加权重和（`film.h:249-254`）：

```cpp
for (int c = 0; c < 3; ++c)
    pixel.rgbSum[c] += weight * rgb[c];   // weight = cameraSample.filterWeight
pixel.weightSum += weight;
```

渲染结束时 `GetPixelRGB` 用权重和归一化（`film.h:257-264`）：

```cpp
RGB rgb(pixel.rgbSum[0], pixel.rgbSum[1], pixel.rgbSum[2]);
if (weightSum != 0) rgb /= weightSum;
```

因此滤波权重在分式中相互抵消，只负责把样本**重建**到像素网格（抗锯齿/重建模糊），
不承担曝光或色彩职责。注意 `AddSample` 一次只写入样本所在的单个像素
（`pFilm` 已由滤镜采样决定落点）；跨像素的 `filter.Evaluate` 加权只出现在 splat 路径
（`RGBFilm::AddSplat`，`film.cpp:499-525`，用于 SPPM/光路径），splat 项在最终归一化后
以 `splatScale · rgbSplat / filterIntegral` 加入（`film.h:266-268`），
`filterIntegral = filter.Integral()`（`film.cpp:491`）。

### 1.6 色彩空间变换：两个 3×3 矩阵的乘积，线性、无色调映射

传感器 RGB 与输出 RGB 之间隔着一个在胶片构造时就算好的矩阵（`film.cpp:495-496`，
GBuffer/Spectral 胶片相同，见 `film.cpp:656`、`:861`）：

```cpp
outputRGBFromSensorRGB = colorSpace->RGBFromXYZ * sensor->XYZFromSensorRGB;
```

- `sensor->XYZFromSensorRGB` 把传感器 RGB 映射到 CIE XYZ；
- `colorSpace->RGBFromXYZ` 再把 XYZ 映射到输出色彩空间（默认 sRGB，
  `scene.cpp:99`、`parser.cpp:1116` 等处默认 `RGBColorSpace::sRGB`；可选 DCI-P3、
  Rec2020、ACES2065-1，见 `util/colorspace.h:35`）。

该矩阵在像素归一化之后才应用（`film.h:271`；`ToOutputRGB` 提供不经像素累积的同款变换，
`film.h:292-296`）。**整条链路没有色调映射（tone map）、没有 OETF/伽马编码**：
`GetImage` 产出的是线性 HDR 像素，默认存储为 fp16 Half（`film.cpp:536-537`，
`savefp16` 默认 true），色彩空间写进 `ImageMetadata`（`film.cpp:560-562`）。
EXR 原样写线性值（`image.cpp:1173` 的 `WriteEXR`）；只有写 PNG/PFM/QOI 等非 EXR 格式
且色彩空间不是 sRGB 时，才用 `ConvertRGBColorSpace` 做一次线性矩阵转换到 sRGB
（`image.cpp:989-1006`），仍非线性传递函数编码。fp16 超出 65504 的值会被钳制并发出
警告（`film.cpp:543-558`）——这是存储格式限制，属第 2 节的偏差范畴。

`GBufferFilm` 中另有一条**独立于主路径**的反照率颜色流：可见表面的光谱反照率先乘输出
色彩空间照明体，再走 `SampledSpectrum::ToRGB`（即 1.3 的 XYZ 路径），累加到
`rgbAlbedoSum`（`film.cpp:631-635`），与辐射度主通道分开存储。

---

## 2. 噪声（方差）、能量偏差与色彩映射的责任区分

### 2.1 噪声/方差：蒙特卡洛采样的责任，胶片只做无偏归一化

下列环节决定像素最终的随机误差，它们都属于**采样与估计**而非胶片的确定性处理：

| 方差来源 | 代码位置 | 性质 |
|---|---|---|
| 波长随机采样（4 波长/样本） | `spectrum.h:283-343`，随机数入口 `integrators.cpp:230` | 光谱积分估计误差；可见度重要性采样降低 CIE 匹配函数积分方差（测试见 `spectrum_test.cpp:123-159`） |
| 路径传输的 BSDF/光源/NEE 采样 | `cpu/integrators.cpp` 各 `Li()` | 4 个波长通道各自的路径估计方差，与颜色管线无关 |
| 像素子像素位置 + 滤波权重 | `samplers.h:797-814`、`filters.h:33-38` | 空间重建；`weight=f/pdf` 保证权重归一化后期望无偏 |
| 像素内 `rgbSum/weightSum` 归一化 | `film.h:257-264` | 比值估计量；`weightSum==0` 时输出 0（未命中任何样本的像素） |
| 显式方差统计（仅 GBuffer） | `VarianceEstimator` `util/sampling.h:484-520`，调用 `film.cpp:601-602`，输出 `Variance/RelativeVariance` 通道 `film.cpp:785-792` | 只**记录**样本方差供降噪/训练，不回改像素值 |

要点：`rgbSum += filterWeight·rgb` 与 `weightSum += filterWeight` 是加权平均的标准
无偏形式；`TerminateSecondary` 调整的也是估计量权重（PDF），其正确性由
`film_test.cpp:97-154` 锁定。胶片本身不注入随机性，也不做任何时域降噪。

### 2.2 能量偏差：代码中显式、有界、可识别的确定性修改

下列操作会改变能量的期望值，属于**偏差**而非噪声，全部位于传感器/胶片/输出层，路径
传输的 MC 估计在进入这些操作之前保持无偏：

1. **`maxComponentValue` 钳制（传感器 RGB，默认关闭）。**
   `RGBFilm::AddSample` 在转成传感器 RGB 后，若最大通道超过 `maxComponentValue`
   （默认 `Infinity`，即不钳制；`film.cpp:576`），则整体等比缩小
   （`film.h:244-247`）。这是对高光/火灾flies 的硬性削峰，方向是系统性压低能量。
   GBuffer（`film.cpp:592-594`）、Spectral（`film.h:421-424` 对 RGB，`:437-441`
   对光谱量）同理。
2. **fp16 存储钳制（输出格式限制，默认开启）。** 超过 65504 的线性值被截断并计数告警
   （`film.cpp:543-558`、`:761-769`；Spectral `:987-994, 1009-1012`）。把
   `savefp16` 设为 false 即用 32 位 float 输出可消除该项。
3. **NaN/Inf 防护（积分器入口）。** `Li` 返回 NaN 或 Inf 时整条样本辐射度置 0 并记
   错误日志（`integrators.cpp:263-273`）。这是故障兜底，正常场景不应触发。
4. **曝光成像比（线性缩放，默认 1）。** `imagingRatio = exposureTime·ISO/100`
   （`film.cpp:233`）。默认参数下不改变信号；非默认时是确定性的物理量纲换算
   （辐射度 → 传感器曝光量），不改变相对色彩关系。
5. **splat 与样本的尺度对齐。** splat 项除以 `filterIntegral`
   （`film.h:266-268`），SpectralFilm 的光谱 splat 还有 `CIE_Y_integral·nBuckets /
   ((lambdaMax-lambdaMin)·filterIntegral)` 的换算（`film.cpp:1002-1008`）。
   `film_test.cpp:21-89`（`SplatMatchesSample`）验证了样本与 splat 混合后 RGB 与光谱
   桶结果一致，说明该尺度设计在仓库测试内自洽。

### 2.3 色彩映射：传感器响应、白平衡与色彩空间矩阵的责任

“光谱长什么样”与“最终写成哪三个数”之间的色度责任全部集中在 `PixelSensor` 与
`outputRGBFromSensorRGB`：

- **默认传感器 `cie1931`（`film.cpp:223, 239-241`）。** 走 `PixelSensor` 的 XYZ 构造
  函数（`film.h:80-92`）：三条响应曲线直接绑定 CIE `X/Y/Z`
  （`r_bar=&Spectra::X(), g_bar=&Spectra::Y(), b_bar=&Spectra::Z()`，`film.h:82-84`）。
  此时 `XYZFromSensorRGB` 仅在指定白平衡温度时才是 Bradford 白平衡矩阵
  （`WhiteBalance`，`util/color.h:556-565`）；默认 `whitebalance = 0` 且传感器为
  cie1931 时 `sensorIllum = nullptr`（`film.cpp:235-237`），构造函数不施加白平衡，
  矩阵保持中性。也就是说**默认配置下传感器 RGB 即 CIE XYZ 的线性对应量，再由输出色彩
  空间的 `RGBFromXYZ` 转成 sRGB 线性值**。
- **真实传感器（非 cie1931）。** 仓库内置了多组相机响应曲线数据
  （如 `canon_eos_*_r/g/b`、`sony_ilce_9_*`，数据表在 `util/spectrum.cpp:1440`
  起，经 `GetNamedSpectrum` 按 `<name>_r/g/b` 取用，`film.cpp:243-251`、
  `spectrum.cpp:2825-2830`）。此时走 `film.h:45-78` 的构造函数：用 24 个
  ColorChecker 反射率样本（`swatchReflectances`，`film.cpp:259-263` 起，数据来源在
  源码注释中注明为 BabelColor/Pascale 测量数据，属仓库内自带材料）分别在传感器照明下
  投影出传感器 RGB、在输出照明体下投影出 XYZ，再用线性最小二乘
  （`LinearLeastSquares`，`util/math.h:698-702`）解出 `XYZFromSensorRGB`。
  非 cie1931 时白平衡温度默认 6500 K（`film.cpp:226-228`），照明体由
  `Spectra::D(T)` 按 CIE D 系列基函数生成（`spectrum.cpp:2555-2566`）。
  这是一个**基于仓库内 24 个反射率样本的最小二乘色彩映射**，其残差与适用范围由该训练
  集决定；仓库没有提供额外的相机标定输入机制。
- **输出色彩空间。** `RGBFromXYZ` 由 `RGBColorSpace` 按其原色与白点构造
  （`util/colorspace.h:34, 56`），默认 sRGB。矩阵只做线性三刺激值空间转换，不做
  色域压缩或色调曲线；非 EXR 格式的 sRGB 强制转换也是线性矩阵（`image.cpp:989-1006`）。

### 2.4 三类责任的边界一句话总结

- **噪声**回答“这次估计离期望有多远”——由波长/路径/像素采样决定，胶片用权重和归一化
  保持无偏，GBuffer 只额外记录方差。
- **偏差**回答“期望本身被改动了多少”——仅来自显式钳制、fp16 截断、NaN 置零与曝光
  缩放，默认配置下其中三项关闭或为 1，唯一默认可达的偏差是 fp16 上限钳制。
- **色彩映射**回答“光谱三刺激值如何对应到输出 RGB 三个通道”——由传感器响应曲线、
  （可选）白平衡与两个 3×3 矩阵负责，与采样方差完全正交。

---

## 3. 默认传感器、可见表面与采样权重的边界

### 3.1 默认传感器设置的边界

默认链路（`scene.cpp:92-99` 的实体默认值 + 各 `Create` 的参数默认值）：

| 组件 | 默认值 | 代码依据 |
|---|---|---|
| 胶片类型 | `rgb`（RGBFilm） | `scene.cpp:98`；分发 `film.cpp:1069-1083` |
| 输出色彩空间 | sRGB（线性，HDR） | `scene.cpp:99`、`parser.cpp:1116` |
| 传感器 | `cie1931`，响应曲线 = CIE X/Y/Z | `film.cpp:223`、`film.h:80-92` |
| ISO / 白平衡 | 100 / 0（cie1931 下不做白平衡） | `film.cpp:220-221, 226-237` |
| 曝光时间 | `shutterclose - shutteropen = 1` | `scene.cpp:756-757`、`cameras.cpp:274-275` |
| 成像比 | `1 × 100/100 = 1`（passthrough） | `film.cpp:233` |
| 波长采样 | 可见度重要性采样，360–830 nm，4 分层波长 | `film.h:211-214`、`sampling.h:163-171` |
| 重建滤镜 | Gaussian，半径 1.5，σ 0.5 | `scene.cpp:94`、`filters.cpp:40-47` |
| 像素钳制 / 像素格式 | 不钳制（∞）/ fp16 Half | `film.cpp:576-577`、`film.cpp:536` |
| 采样器 | zsobol | `scene.cpp:93` |

边界含义：

1. 默认传感器**不是**某个具体相机型号，而是 CIE 1931 标准观察者的 XYZ 匹配曲线；
   它感知的是整个 360–830 nm 光谱（经 4 个随机波长 MC 积分），不模拟 Bayer 马赛克、
   快门闪光同步、读出噪声等器件效应——仓库里没有这些机制。
2. 默认成像比为 1，意味着输出数值的尺度直接对应积分器算出的线性场景参照辐射度
   （物理单位场景下配合 PhysLight 光源单位；README 第 67 行说明可在绝对物理单位下
   配合真实相机建模）。胶片层不做任何“好看”的调色。
3. 切换到真实传感器或 SpectralFilm 才会改变上述边界：前者替换响应曲线并引入最小二乘
   矩阵与（默认 6500 K 的）白平衡；后者把波长采样改成均匀采样并额外输出光谱桶
   （只能写 EXR，`film.cpp:1048-1050`）。

### 3.2 VisibleSurface 的边界：仅 GBufferFilm 使用，且与辐射度通道分离

- `Film::UsesVisibleSurface()` 是开关：RGBFilm 与 SpectralFilm 返回 false
  （`film.h:236`、`:405`），GBufferFilm 返回 true（`film.h:346`）。
- 积分器仅在该开关为真时才构造 `VisibleSurface`
  （`integrators.cpp:258-260`），其内容是几何与材质属性：位置、几何/着色法线、UV、
  光线微分、**光谱反照率**（`film.h:133-156`，构造 `film.cpp:191-204`）。
- 在 `GBufferFilm::AddSample` 中，`visibleSurface` 只驱动 G-buffer 通道
  （位置/法线/UV/反照率/方差，`film.cpp:597-636`）；**辐射度主通道
  `rgbSum/weightSum` 的更新与 `visibleSurface` 完全无关**
  （`film.cpp:638-640`），与 RGBFilm 的累积式一致。
- 因此 `VisibleSurface` 不参与传感器响应、曝光或色彩矩阵；它是旁路的调试/降噪数据流。
  RGBFilm 下 `AddSample` 收到的该指针被显式忽略（参数无名，`film.h:239-240`）。

### 3.3 采样权重的边界：三个权重各有归属，互不越界

| 权重 | 持有者 | 作用点 | 边界 |
|---|---|---|---|
| 波长 `pdf[4]` | `SampledWavelengths` | 仅在 `ToSensorRGB`（`film.h:96`）与 `ToXYZ/y`（`spectrum.cpp:212-221`）反除一次 | 不参与像素滤波；被终止的次级波长 pdf=0，经 `SafeDiv` 安全归零 |
| 相机/路径 `cameraRay.weight`（光谱） | `CameraRay` | 积分器出口乘到 `L`（`integrators.cpp:259`；wavefront `film.cpp:22-23`） | 普通相机恒为 1；仅 RealisticCamera 非 1（镜头透过率/余弦/出瞳 PDF，`cameras.cpp:944-948`） |
| 滤波 `filterWeight`（标量） | `CameraSample`（`base/sampler.h:19-25`） | `AddSample` 的 `weight` 参数，进 `rgbSum` 与 `weightSum`（`film.h:252-254`） | 三通道同权，只做空间重建；归一化后对颜色无方向性影响 |

`AddSample` 的契约（`base/film.h:28-31`）正是这三条权重的汇合点：传入光谱辐射度 `L`、
波长 `lambda`（含 PDF）、可选 `VisibleSurface` 与标量 `weight`（滤波权重）。传感器
反 PDF 与内积在函数内部完成，调用方（CPU 与 wavefront 两条渲染路径，
`integrators.cpp:287-288` 与 `wavefront/film.cpp:31-37`）不需要知道响应曲线或色彩空间
的存在——这就是责任链在 API 层面的边界。

---

## 4. 结论

1. **光谱路径权重进入 RGB 的入口是 `PixelSensor::ToSensorRGB`（`film.h:94-100`）**：
   先按波长采样 PDF 反除得到光谱 MC 估计量，再与三条响应曲线（默认为 CIE X/Y/Z）
   内积、对 4 个波长平均，乘曝光成像比（默认 1）得到线性传感器 RGB。随后由
   `filterWeight` 做像素加权累积与权重和归一化，最后经
   `RGBFromXYZ · XYZFromSensorRGB` 矩阵变换到输出色彩空间（默认线性 sRGB）并以
   fp16 EXR 写出。全链无色调映射、无传递函数编码。
2. **噪声、偏差、色彩映射分属不同责任层**：方差来自波长/路径/像素采样，胶片的加权
   平均保持无偏（GBuffer 的 `VarianceEstimator` 只记录不改值）；确定性偏差只来自
   `maxComponentValue` 钳制（默认关）、fp16 上限（默认开，有警告）、NaN 置零与曝光
   缩放（默认 1）；色度由传感器响应曲线、可选 Bradford 白平衡与两个线性 3×3 矩阵负责，
   与方差正交。
3. **默认配置的边界明确**：cie1931 传感器等价于 CIE 标准观察者、无白平衡、成像比 1、
   可见度重要性波长采样、Gaussian 重建滤镜、线性 HDR 输出；`VisibleSurface` 仅
   GBufferFilm 使用且只影响旁路 G-buffer 通道；波长 PDF、相机权重、滤波权重三个权重
   在 `Film::AddSample` 汇合，各自只在自己的数学环节生效。
4. 以上全部结论均可在本仓库源码与仓库自带测试（`spectrum_test.cpp`、`film_test.cpp`）、
   源码注释及 README 中核对；仓库内置的相机响应曲线与 ColorChecker 数据是其自带的唯一
   标定材料，本文未引入任何外部标定数据，也未执行或声称任何未运行的实验结果。
