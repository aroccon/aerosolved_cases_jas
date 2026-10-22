This repository contains the simulation setups and configuration files for the test cases presented in the paper:

**"Comparison of Eulerian and Lagrangian Frameworks for Aerosol Simulation"** *(in preparation)*

**Authors:** A. Roccon, F. Lucci, and C. Marchioli

The test cases cover different flow configurations and particle distributions, providing a basis for comparing Eulerian and Lagrangian approaches to aerosol transport.

## Carrier-Phase Properties

All simulations consider a Newtonian carrier fluid with the following constant properties:

| Property | Symbol | Value |
|---|---|---|
| Density | $\rho$ | 1 kg/m³ |
| Dynamic viscosity | $\mu$ | $1 \times 10^{-5}$ Pa·s |

## Simulation Cases

The repository includes the following configurations:

| Case | Flow configuration | Reynolds number | Particle distribution | Status |
|---|---|---|---|---|
| `cavity_1000` | Lid-driven cavity | 1,000 | 10,000 particles, 10 µm diameter | Available |
| `cavity_3200` | Lid-driven cavity | 3,200 | 10,000 particles, 10 µm diameter | Available |
| `taylor_green` | Periodic Taylor–Green vortex | 1,600 | 9,000 particles, 1–10 µm diameter (9 sections) | Available |
| `hit_frozen` | Frozen homogeneous isotropic turbulence | 57 (Taylor) | 9,000 particles, 1–10 µm diameter (9 sections) | Available |
| `hit_decaying` | Decaying homogeneous isotropic turbulence | 57 (Taylor) | — | Work in progress |

### 1. Lid-Driven Cavity (`cavity_1000`, `cavity_3200`)

These cases consider particle transport in a lid-driven cavity flow.

**Flow parameters:**
- Lid velocity: $U = 1$ m/s
- Reynolds numbers: $Re = 1000$ and $Re = 3200$
- Characteristic length $L_0$: adjusted to obtain the target Reynolds number.

**Particle parameters:**
- Number of particles: 10,000
- Particle diameter: $d_p = 10$ µm
- Monodisperse particle distribution.

### 2. Taylor–Green Vortex (`taylor_green`)

This case considers aerosol transport in a periodic Taylor–Green (TG) vortex.

**Flow parameters:**
- Periodic flow configuration
- Reynolds number: $Re = 1600$

**Particle parameters:**
- Number of particles: 9,000
- Particle diameter range: 1–10 µm
- Polydisperse distribution represented by 9 discrete size sections.

### 3. Frozen Homogeneous Isotropic Turbulence (`hit_frozen`)

This case considers aerosol transport in a homogeneous isotropic turbulence box (frozen).

**Flow parameters:**
- 3D Periodic flow configuration
- Taylor Reynolds number: $Re_\lambda = 57$

**Particle parameters:**
- Number of particles: 50,000
- Particle diameter range: 1–10 µm
- Polydisperse distribution represented by 9 discrete size sections.

### 4. Decaying Homogeneous Isotropic Turbulence (`hit_decaying`)

**Status:** Work in progress (WIP).


