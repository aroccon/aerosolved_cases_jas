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
| `hit_frozen` | Frozen homogeneous isotropic turbulence | 57 (Taylor) | 50,000 particles, 1–10 µm diameter (9 sections) | Available |
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
- 3D periodic box of side $L = 0.01$ m, $150^3$ cells (grid study: other resolutions with the same field)
- Taylor Reynolds number: $Re_\lambda = 57$
- Frozen flow (`frozenFlow true`): momentum and pressure are not solved, the velocity field stays as in `0/U`

| Quantity | Value |
|---|---|
| $u'_{rms}$ | 0.744 m/s |
| Dissipation $\varepsilon$ | 140.5 m²/s³ |
| Taylor microscale $\lambda$ | 0.768 mm |
| Kolmogorov length $\eta$ | 51.7 µm |
| Kolmogorov time $\tau_\eta$ | $2.67 \times 10^{-4}$ s |
| $\Delta x/\eta$ ($150^3$) | 1.29 |

**Particle parameters:**
- Number of particles: 50,000, diameters uniformly distributed in 1–10 µm (one `manualInjection` injector, positions in `constant/kinematicCloudPosition`)
- Polydisperse distribution represented by 9 discrete size sections (EE)
- Stokes number based on the section diameters (1.5–9.5 µm): $St_\eta = 0.047$–$1.88$
- Particles initially at rest in both frameworks (EL: `U0 (0 0 0)`; EE: `0/V` $= -$`0/U`)

**Carrier velocity field.** The initial (frozen) field is a snapshot of a forced HIT simulation performed with MHIT36 ($256^3$, $\nu = 0.015$, ABC forcing at $k = 2$; see `input_mhit36.inp`), rescaled to the physical box and viscosity (all Reynolds numbers preserved). `0/U` and `0/V` are not stored here because of their size (about 160 MB each); they are generated with the scripts in `hit_frozen/scripts`:

1. `hit_stats.m` — statistics of the MHIT36 field ($\varepsilon$, $\eta$, $\lambda$, $Re_\lambda$, spectrum) and of the rescaled field.
2. `hit_to_openfoam3d.m` — maps the staggered MHIT36 field onto the OpenFOAM cell centres (Fourier interpolation, exact phase shift of the staggered components, projection with OpenFOAM's discrete divergence so that `div(phi)` is zero to round-off) and writes `0/U`. Set `no` (cells per direction) and `L`.
3. `make_particles.py` — writes `constant/kinematicCloudPosition` and `0/V = -U`:
   `python3 make_particles.py hit_frozen --N 50000 --L 0.01`

**Case-specific settings:**
- `system/fvOptions`: the temperature is held at 293.15 K. With a frozen flow the energy equation is still solved and, without this constraint, the pressure-work term makes $T$ (and therefore $\rho$) drift unphysically.
- `0/M.1` … `0/M.9`: the section fields are provided explicitly. Otherwise each section is initialised from the patch `upperBoundary` (`initFromPatch`), which fails in parallel when the cyclic patches have no faces on any processor.
- `system/controlDict`: `maxCo 0.25` (`maxCo 0.5` from the start makes the $M_i$ equations diverge during the initial slip transient), output every 1 ms, binary format, `endTime 0.02` (about 75 $\tau_\eta$).
- `constant/kinematicCloudProperties`: one-way coupling, drag only, `maxCo 0.1`.

**Solver.** AeroSolved, branch `feature-electrostatic` (commit `7b4c48d70ee`), OpenFOAM v2312, with the patch in `aerosolved_patch/frozenFlow_correctDriftFlux.patch`. In the original code the frozen-flow branch of `aerosolEulerFoam` does not call `aerosol->correctDriftFlux()`, so the drift velocities are not updated and the convection schemes used in `YEqn.H` are not created.

### 4. Decaying Homogeneous Isotropic Turbulence (`hit_decaying`)

**Status:** Work in progress (WIP).


