# Balanced Iso-Spring Oscillator (BISO)

A MATLAB simulation and analysis toolkit for a 2D Balanced Iso-Spring Oscillator mechanism using parallel pantograph architecture.

## Overview

The Balanced Iso-Spring Oscillator (BISO) is a mechanical oscillator design that achieves isotropic (direction-independent) spring behavior in two dimensions. This is accomplished through a symmetric arrangement of four pantograph mechanisms oriented in the cardinal directions (North, South, East, West).

This repository contains MATLAB code for analyzing:
- **Equivalent Spring Stiffness**: How the effective spring constant varies across the oscillation plane
- **Reduced Mass**: The effective inertia of the moving system
- **Mass Isotropy**: Evaluation of the center-of-mass stability
- **Power Requirements**: Energy dissipation and power reserve calculations for watch mechanisms

## Repository Structure

```
balanced-isospring-oscillator/
├── src/
│   ├── equivalent_stiffness.m    # K_eq analysis and visualization
│   ├── reduced_mass.m            # Reduced mass and isotropy analysis
│   └── power_analysis.m          # Power dissipation calculations
├── lib/
│   ├── polarplot3d.m             # 3D polar plotting utility (third-party)
│   └── polarplot3d_LICENSE.txt   # BSD license for polarplot3d
├── figures/                       # Generated plot images
├── docs/                          # Technical documentation
├── LICENSE
└── README.md
```

## Requirements

- MATLAB R2018b or later
- No additional toolboxes required

## Usage

### Equivalent Stiffness Analysis

Computes and visualizes the equivalent spring stiffness as a function of position in the oscillation plane:

```matlab
cd src
equivalent_stiffness
```

This generates a 3D polar surface plot showing K_eq = f(x, y).

### Reduced Mass Analysis

Calculates the equivalent reduced mass and evaluates mass isotropy:

```matlab
cd src
reduced_mass
```

Generates multiple plots:
- 3D surface: Reduced mass vs frequency and radius
- 2D curves: Reduced mass vs radius (parametric in frequency)
- 2D curves: Reduced mass vs frequency (parametric in radius)
- Isotropy error analysis plots

### Power Analysis

Analyzes power requirements and reserve margins for watch applications:

```matlab
cd src
power_analysis
```

Evaluates the relationship between:
- Number of gear train stages
- Mainspring barrel capacity
- Power reserve margins

## System Parameters

### Geometric Constants

| Parameter | Value | Description |
|-----------|-------|-------------|
| L | 25.213 mm | Reference length |
| l | 4 mm | Half platform length |
| a | 15 mm | Pantograph arm length |

### Spring Constants

| Parameter | Value | Description |
|-----------|-------|-------------|
| k_gamma | 7.85 mN·m/rad | Spring constant at gamma pivot |
| k_beta | 6.44 mN·m/rad | Spring constant at beta pivot |

### Operating Range

| Parameter | Min | Max | Description |
|-----------|-----|-----|-------------|
| rho | 0.3 mm | 0.4 mm | Orbit radius |
| f | 1 Hz | 15 Hz | Oscillation frequency |

## Theory

The BISO mechanism achieves 2D isotropy through:

1. **Symmetric Architecture**: Four identical pantograph units arranged at 90° intervals
2. **Balanced Mass Distribution**: Counter-masses on each pantograph arm maintain center-of-mass stability
3. **Geometric Compensation**: The parallel pantograph geometry provides consistent restoring force regardless of oscillation direction

The equivalent stiffness K_eq is derived from the elastic potential energy stored in the pivot springs:

```
K_eq = (1/r²) × Σ[2·k_γ·(Δγ₁² + Δγ₂²) + 2·k_β·Δβ²]
```

where Δγ and Δβ represent angular deviations from equilibrium.

## Documentation

Detailed technical documentation is available in the `docs/` folder:
- `Pantographe Parallele BISO Balanced-IsoSpring-Oscillator.pdf` - Complete theoretical analysis

## License

Copyright (c) 2022 Daniel Abraham Elmaleh. All rights reserved.

See [LICENSE](LICENSE) for details.

### Third-Party Components

- `polarplot3d.m` by Ken Garrard is included under BSD license. See `lib/polarplot3d_LICENSE.txt`.

## Author

**Daniel Abraham Elmaleh**

Project developed as part of mechanical design coursework (CdM2 - Groupe 33).

## Citation

If you use this code in academic work, please cite:

```
@software{elmaleh2022biso,
  author = {Elmaleh, Daniel Abraham},
  title = {Balanced Iso-Spring Oscillator (BISO) - MATLAB Analysis Toolkit},
  year = {2022},
  url = {https://github.com/[username]/balanced-isospring-oscillator}
}
```
# balanced-isospring-oscillator
