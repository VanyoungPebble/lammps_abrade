# fix hydro/drag 最简示例

一颗可磨损的球形颗粒（`data.Sphere`，387 个原子）从静止开始，被沿 x 方向的均匀水流带动。
没有重力、没有碰撞，只有水流阻力。`units cgs`。

## 运行

复制 `in.hydro_drag` 和 `data.Sphere` 到空目录，然后：

```sh
lmp -in in.hydro_drag
```

`lmp` 换成编译好的可执行文件路径（需要 RIGID 包和 qhull）。
**只能单进程运行**：`fix hydro/drag` 第一版在多于 1 个 MPI 进程时会报错退出。
在 MacBook 上约 2 秒跑完 50000 步（0.5 s 模拟时间）。

## 输入里的关键几行

```lammps
fix water all hydro/drag rho 1.0 nu 0.01 flow 50.0 0.0 0.0
fix abrade all rigid/abrade ${H} 0.25 ${density} molecule normals drag water
```

- 第一行定义流体参数：密度 `rho`、运动黏度 `nu`、均匀流速 `flow`。`water` 是给这个 fix 起的名字，可以随便取。
- 第二行末尾的 `drag water` 让 rigid/abrade 每一步把名为 `water` 的 fix 算出的阻力加到每颗颗粒的质心合力上。
  阻力不加到原子上，所以**不会引起磨损**，也不产生力矩。
- 要加重力和浮力：定义 `fix gravity`，在 rigid/abrade 里写 `gravity <ID>`，并在 hydro/drag 里写 `buoyancy yes`
  （输入文件里有注释掉的示例）。写了 `buoyancy yes` 却没有 `gravity` 关键字会报错。

## 输出

- `thermo`：`f_water[1-3]` 为所有颗粒的总阻力，`f_water[4-6]` 为总浮力（本例未开浮力，为 0）。
- `drag.dump`（`dump local`，每颗颗粒一行，10 列）：
  1 颗粒 ID（质心原子的 ID）、2 相对速度大小 |u_rel|、3 Re、4 C_D、5 球形度 Φ、6 横向球形度 Φ⊥、
  7 投影面积 A⊥、8–10 阻力 Fx Fy Fz。
- `body.dump`：`compute rigid/local_abrade` 输出的颗粒 ID 和质心速度。

注意：hydro/drag 的逐颗粒输出从第 1 步才开始有数据（第 0 步 0 行），
所以它要单独用一个 `dump local`，不能和 `compute rigid/local_abrade` 这类第 0 步就有数据的量放在同一个 dump 里，
否则 LAMMPS 报 "Dump local count is not consistent across input fields"。

## 预期结果

颗粒速度 vx 从 0 逐渐增大并趋近 50 cm/s：0.05 s 时约 16.6，0.5 s 时约 41.8 cm/s；
阻力随相对速度减小而减小（0.05 s 约 1377 dyne，0.5 s 约 87 dyne）。
因为高 Re 时阻力约与相对速度的平方成正比，越接近流速加速越慢。

这个网格球的 Φ = 0.9965 而不是 1。Holzer 公式在 Φ 接近 1 时变化很陡，
所以它的阻力系数比理想球偏大（详见方案文档中的说明）。本例只演示用法，尚未做数值验证。

## 其他说明

- 本例没有碰撞，所以用了比通常 DEM 大得多的时间步长 `1e-5 s`；有碰撞时请用正常的 DEM 时间步长。
- 投影面积 A⊥ 每步用 qhull 求二维凸包，再加原子半径修正 A + P·r + π·r²（始终开启，与 rigid/abrade 的 `offset` 无关）。
- 每颗颗粒每步算一次凸包（约 34 µs/颗），小算例会明显变慢。
