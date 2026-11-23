# Qhull 投影面积初步实现

该原型从原版 `develop` 的 `82c9417` 创建分支 `feature/projected-area-qhull`，
独立工作树位于 `/Volumes/ORICO/LAMMPS_projects/projected_area_qhull`。
它不包含 `feature/realtime-shape` 的 OBB 实现，也不调用该实现的取点流程。

## 使用

在原来的 `fix abrade ...` 后加入：

```lammps
compute ap owned projected/area 1 0 0
dump projected all local 100 projected_area.dump c_ap[1] c_ap[2]
dump_modify projected sort 1 format float %.17g
```

`owned` 应包含每颗目标颗粒的完整表面 atom，排除中心 owning atom。
这个组名沿用现有水槽输入的定义，不是新增命令自动创建的组。
`1 0 0` 表示沿世界坐标 x 方向投影到 yz 平面；方向可以更换为任意有限、非零常量向量。
输出两列是 molecule ID 和投影凸包面积，`units cgs` 时面积单位为 cm²。
没有新增有效性列；退化点集或 Qhull 失败写 NaN，并首次警告。

复制 `in.projected_area` 和 `data.Box` 到空目录，用下面的程序运行，即可试验一个
平移、旋转、穿越周期边界的长方体：

```sh
/Volumes/ORICO/LAMMPS_projects/projected_area_qhull_build/lammps/lmp -in in.projected_area
```

这份示例的初始长方体尺寸为 2×3×4，因此初始沿 x 的投影面积应为 12。
MPI 运行使用本机已有的 `/Volumes/ORICO/opt/anaconda3/bin/mpiexec`。

## 修改位置

生产实现仅涉及以下三个文件：

| 文件 | 改动 | 职责 |
| --- | --- | --- |
| `cmake/CMakeLists.txt` | 原文件增加 5 行 | RIGID 启用时查找并链接 Qhull C++ 与可重入静态库 |
| `src/RIGID/compute_projected_area.h` | 新增 | 注册 `compute projected/area` 和声明类 |
| `src/RIGID/compute_projected_area.cpp` | 新增 | 构造投影平面、直接读取 atom 坐标、MPI 汇集、按 molecule 调用 Qhull、输出两列 |

另新增本说明、命令文档 `doc/src/compute_projected_area.rst`、示例输入、解析长方体数据及
`validate_projected_area.py`。`fix_rigid_abrade.cpp/.h` 和
`compute_rigid_local_abrade.cpp/.h` 与原版基线相同。

核心源码入口：`projected_hull_area()` 是 Qhull 适配；构造函数处理固定投影方向；
`compute_local()` 读取当前 `atom->x/image`，通过 `domain->unmap()` 展开周期坐标，
只收集本 rank 的 local atom。两个 `MPI_Gatherv` 分别传送整数 molecule ID 和二维坐标；
rank 0 按 molecule 分组计算，每颗颗粒输出一行。

没有自行实现凸包算法。`std::sort/unique` 只做点排序和完全重复坐标去重；
凸包和面积交给 Qhull 的官方 C++ 接口 `runQhull(..., 2, ..., "Qt")` 和 `volume()`。
二维的 `volume()` 是面积，`area()` 是周长。未开启随机扰动 `QJ`，也没有其他算法回退。

## 依赖与构建

使用 [Qhull 2020.2](https://github.com/qhull/qhull/tree/2020.2)（库版本 8.0.2），
官方 C++ 封装的底层是可重入的 `libqhullstatic_r`。
官方建议新代码使用[可重入接口](https://www.qhull.org/)。
依赖源码、编译目录、安装目录均放在独立的
`/Volumes/ORICO/LAMMPS_projects/projected_area_qhull_build`，没有安装到系统目录。
源码与许可证保存在 `deps/qhull-2020.2`，来源归档及 SHA256 记录在该构建目录的
`BUILD_INFO.json`。编译 Qhull 时启用了静态库和位置无关代码。

当前可复用的 Qhull 安装前缀：

```text
/Volumes/ORICO/LAMMPS_projects/projected_area_qhull_build/qhull-install
```

从本工作树重新配置和编译：

```sh
cmake -S cmake -B ../projected_area_qhull_build/lammps \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -DBUILD_MPI=ON \
  -DMPI_CXX_COMPILER=/Volumes/ORICO/opt/anaconda3/bin/mpicxx \
  -DMPI_CXX_SKIP_MPICXX=ON \
  -DPKG_GRANULAR=ON -DPKG_MOLECULE=ON -DPKG_RIGID=ON \
  -DDOWNLOAD_POTENTIALS=OFF \
  -DCMAKE_PREFIX_PATH=/Volumes/ORICO/LAMMPS_projects/projected_area_qhull_build/qhull-install
cmake --build ../projected_area_qhull_build/lammps --parallel 8
```

本原型分支启用 RIGID 时会要求 Qhull；传统 make 构建尚未接入。

## 验证结果

已完成 20 次隔离运行、19 组检查，全部符合预期；其中 2 次非法输入按预期报错。
详细结果在 `projected_area_qhull_build/validation/validation_results.json`。

- 长方体 x/y/z/斜向投影与解析面积一致，方向缩放和反向不改变面积。
- 1/2 MPI rank 均验证了旋转、跨周期边界和逐 molecule 独立输出。
- 旋转核对使用实际三条棱的解析投影公式；原 fix 会由角动量重算角速度，不能直接将
  初始节点速度中的系数当成运行时刚体角速度。没有为此修改原 fix。
- 空组输出 0 行；共线点集经 Qhull 失败后输出 NaN；零方向和非法 molecule ID 正确报错。
- 手动移动和删除节点后，新输出面积跟随当前点集变化。
- 初始 `data.Gravel` 沿 x 的面积为 `2.68643409438313 cm²`，与已知凸网格的独立面片公式一致。
- 各自相同 MPI 核数下，原版、新版未启用测量、新版启用测量的 atom 坐标、速度、受力输出
  逐字节一致。没有以 1/2 核轨迹逐字节一致作为要求。

复现（Python 仅使用标准库）：

```sh
python3 examples/rigid/Fix_Rigid_Abrade/Projected_Area/validate_projected_area.py \
  --lmp ../projected_area_qhull_build/lammps/lmp \
  --mpiexec /Volumes/ORICO/opt/anaconda3/bin/mpiexec \
  --baseline /private/tmp/shape-baseline-20261116/build/lmp \
  --gravel /Volumes/ORICO/LAMMPS_projects/Linear_Flume_Hardness_Test/data.Gravel \
  --out /Volumes/ORICO/LAMMPS_projects/projected_area_qhull_build/validation
```

没有可用的原版程序时可以省略 `--baseline`；这会跳过原版轨迹对照。
MPI 需要本机通信端口，受限沙箱可能禁止其启动。

## 初版边界

- 面积来自 atom **中心**点的二维凸包；不包含球半径，也不恢复 `offset` 前的网格表面。
- 输出按用户提供的固定方向，按请求频率计算；不自动使用水流或相对速度，不反馈到阻力。
- 凹轮廓和孔洞会被凸包填平；不同颗粒间的遮挡不在这个量中。
- MPI 采用 rank 0 汇集，优先保持实现简单；尚未测量生产水槽的长期运行开销。
- 本次检查包括节点移动、删除，但没有声称完成真实磨损重网格、offset、restart 的完整验证。
- 当前改动保留在新分支工作树中，尚未 commit 或 push。
