# constellation-orbital-revisit-time-and-ground-access
# Walker Satellite Constellation Analysis Toolkit

This repository contains MATLAB scripts for optimizing and analyzing Walker-style satellite constellations with respect to **global revisit time**, with a focus on **air traffic monitoring**.  

The toolkit includes two main scripts:  

1. **`Opt_toolbox_pso_rev_time_walker_V3.m`** – PSO-based constellation optimization.  
2. **`analyze_best_solution_walker_V3_2.m`** – Revisit time analysis of optimized constellations.

---

## Objective

The goal of this toolkit is to provide a **flexible, easy-to-use framework** for:  

- Designing Walker satellite constellations for given mission constraints.  
- Computing global coverage and revisit metrics for a given constellation.  
- Performing preliminary assessments for air traffic monitoring, considering aircraft at typical cruising altitudes.  

This toolkit can be used for **research, teaching, or preliminary mission design**.

---

## Features

### WalkerCircularPSO_Grid.m

- **Particle Swarm Optimization (PSO)** to find optimal constellation parameters:
  - Number of planes (Np)  
  - Satellites per plane (Ns)  
  - Phasing factor (F)  
  - Orbital semi-major axis (a)  
  - Inclination (inc)  
- **Equal-area Earth grid** generation for evaluating global coverage.  
- **Revisit time calculation** for each grid point:
  - Worst-case revisit  
  - Average revisit  
- Considers **aircraft signals at 6 km altitude** with a local elevation mask.  
- Penalizes infeasible constellations automatically (e.g., exceeding maximum satellites).  
- Parallelized computation using MATLAB’s `parpool`.  
- Saves the best solution for later analysis.

### WalkerCircularPSO_Analysis.m

- **Load previously optimized constellations** or manually define constellation parameters.  
- Computes **revisit metrics** (worst, average, best) for all Earth grid points.  
- Focuses on latitudes relevant for air traffic monitoring (|lat| ≤ 75°).  
- Generates **visualizations** of revisit time distributions:
  - Global maps of worst, average, and best revisit times.  
  - Maps highlighting only areas of interest.  
- Accounts for **aircraft signal cone geometry**, though it is negligible for typical satellite FOVs.  
- Saves heatmap data for further analysis or plotting.

---

## How to Use

1. **Optimization**  
   Run `WalkerCircularPSO_Grid.m` to generate an optimized constellation.  
   ```matlab
   run('WalkerCircularPSO_Grid.m');
