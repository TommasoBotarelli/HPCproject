# High Performance Scientific Computing: Shallow Water Model Simulation

**Authors**: Tommaso Botarelli, Francesco Moglia
**Institution**: LIÈGE université

Presentation: [link](https://github.com/TommasoBotarelli/HPCproject/blob/master/Presentazione.pdf)

## 📌 Project Overview
This project focuses on the numerical simulation of the **Shallow Water Model**. The model solves for free-surface elevation ($\eta$) and depth-averaged velocity ($u, v$) over a variable bottom topography ($h$) with a constant reference level. The primary goal is to optimize a baseline serial implementation and scale it using High-Performance Computing (HPC) techniques such as MPI, OpenMP, and GPU parallelization. 

## 🧮 Numerical Methods & Stability
* **Bilinear Interpolation**: Implemented to replace the original nearest-neighbor algorithm, enhancing the accuracy of converting raw bathymetry data into computational datapoints.
* **Stability Analysis**: Ensures reliable and physically meaningful results. We identified a linear relation between space ($dx$, $dy$) and time ($dt$) discretizations, determining the threshold where solutions become unstable (yielding NaN or unexpectedly large/small values).

## 📈 Arithmetic Intensity & Serial Bottlenecks

### Roofline Model & Memory Boundness
* The code operates on 3 main data structures, requiring 27 Flops and achieving 17.28 MFlops/cycle. 
* The calculated Arithmetic Intensity is **1.125 Flops/byte**.
* Due to this low intensity, the algorithm is severely **memory bandwidth bound**. Even on high-end hardware like the AMD EPYC 7542 processor, the application scales linearly with memory bandwidth and cannot reach the processor's peak compute capabilities.

### Serial Optimizations & Results
Multiple bottlenecks in the initial serial code were fixed:
1. **Boundary Conditions**: Split useless nested cycles into independent cycles.
2. **Data Locality**: Reordered nested loops for the $\eta$, $u$, and $v$ updates to improve spatial locality, which significantly increased data cache hits.
3. **Loop Invariants**: Moved constant multipliers outside of the main computational loops to reduce the number of operations.

**Optimization Results (Execution Time for dx=5, dy=5, dt=0.05, max_time=200):**
* **Baseline (Unoptimized):** ~41.032s 
* **Optimization 1 (Memory Access / Loop Reordering):** ~7.929s 
* **Optimization 2 (+ Inner Loop Constants Precomputed):** ~7.865s 
* **Optimization 3 (+ Outer Loop Constants Precomputed):** ~7.848s

*Note: Data collected via **Likwid** confirmed a massive reduction in data cache misses (thanks to memory access optimizations) and a linear reduction in data cache requests (thanks to optimized constant calculations).*

## 🚀 Parallel Implementations & Scalability Results

### 1. MPI (Distributed Memory)
* The domain was divided into sub-domains and assigned to distinct processing elements using a **Cartesian Topology** with **Non-Blocking Communications**.
* **Strong Scaling Result**: Exhibited an initial superlinear speedup due to better local cache utilization. However, efficiency rapidly declines with more ranks due to growing inter-process communication overhead.
* **Weak Scaling Result**: Efficiency diminishes starting from 2 ranks, primarily driven by cache performance. Spreading processes beyond 8 ranks across multiple NUMA cores caused bandwidth contention and unbalanced `waitAll()` waits, worsening execution time.

<table style="width: 100%;">
  <tr>
    <td align="center" width="33%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/mpi_strong_scaling.png" style="width: 100%;"/><br />
    </td>
    <td align="center" width="33%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/mpi_weak_scaling.png" style="width: 100%;"/><br />
    </td>
    <td align="center" width="33%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/mpi_weak_scaling_spread.png" style="width: 100%;"/><br />
    </td>
  </tr>
</table>

### 2. OpenMP (Shared Memory)
* Parallelized loops using threads accessing shared memory. We attempted to use the `collapse` clause to flatten nested loops, but this destroyed cache optimizations and yielded no performance benefit.
* **Strong Scaling Result**: Like MPI, efficiency heavily correlates with cache performance.
* **Weak Scaling Result**: Yielded **no speedup**. This aligns with the arithmetic intensity analysis: since the algorithm is strictly memory bandwidth-bound, adding more threads does not overcome the bandwidth limit.

<table style="width: 80%;">
  <tr>
    <td align="center" width="50%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/openmp_strong.png" style="width: 100%;"/><br />
    </td>
    <td align="center" width="50%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/openmp_weak.png" style="width: 100%;"/><br />
    </td>
  </tr>
</table>

### 3. Hybrid (MPI + OpenMP)
* Merged the two approaches to maximize node performance on large clusters.
* **Strong Scaling Result**: Superlinear speedup occurs up to 2–4 threads due to optimal cache utilization, dropping as resources increase further.
* **Weak Scaling Result**: Increasing the number of threads successfully mitigates the steep efficiency drop that occurs when scaling up MPI ranks.

<div align="center">
  <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/hybrid_strong_scaling_1.png" width="100%">
</div>

<div align="center">
  <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/results/images/hybrid_scaling.png" width="100%">
</div>

### 4. GPU Parallelization
* Implemented GPU offloading via OpenMP.
* The main bottleneck was device-to-host communication, which was resolved by employing **structured data regions** to persist data on the GPU.
* **Result**: Achieved an overall **10x speedup**.

## 🌊 Real Bathymetry Data Testing
To validate the model, complex bathymetry data of the **Tyrrhenian Sea (Pianosa Island)** was sourced from the EMODnet portal mapped to a 150x150 grid. The raw latitude/longitude coordinates were converted into strict meter discretization steps ($dx$, $dy$) using `opendem`. (Missing island values were handled as NaNs).

<div align="center">
  <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/pianosa_simulation/pianosa.png" width="50%">
</div>

**Simulation Results:**
The physical presence of the island correctly reflected the simulated waves. Testing configurations included:
1. A sinusoidal wave originating from the left side.
2. A sinusoidal wave originating from the upper-left corner.
3. Two opposing sinusoidal waves originating simultaneously from the upper-left and lower-right corners.

<table style="width: 100%;">
  <tr>
    <td align="center" width="33%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/pianosa_simulation/video-1.gif" style="width: 100%;"/><br />
      <sub><b>1</b></sub>
    </td>
    <td align="center" width="33%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/pianosa_simulation/video-2.gif" style="width: 100%;"/><br />
      <sub><b>2</b></sub>
    </td>
    <td align="center" width="33%">
      <img src="https://github.com/TommasoBotarelli/HPCproject/blob/master/pianosa_simulation/video-3.gif" style="width: 100%;"/><br />
      <sub><b>3</b></sub>
    </td>
  </tr>
</table>

## 🔮 Future Developments
* Experimenting with Dynamic Loop Scheduling in OpenMP to better balance thread workloads.
* Further optimization to increase the GPU speedup.
* Implementing dynamic boundaries on island coasts and transparent domain boundaries.
* Transitioning to an implicit integration scheme to expose more potential parallelism.
