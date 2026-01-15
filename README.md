# Balanced Iso-Spring Oscillator (BISO)

A MATLAB simulation and analysis toolkit for a 2D Balanced Iso-Spring Oscillator mechanism using parallel pantograph architecture.

## Overview

The Balanced Iso-Spring Oscillator (BISO) is a mechanical oscillator that achieves isotropic (direction-independent) spring behavior in two dimensions. This is accomplished through a symmetric arrangement of four pantograph mechanisms oriented in the cardinal directions (North, South, East, West).

The BISO design is particularly relevant to watchmaking, where achieving consistent timekeeping requires oscillators that perform uniformly regardless of orientation.

### Key Features

- **Isotropic Spring Behavior**: Uniform restoring force in all directions
- **Balanced Mass Distribution**: Counter-masses maintain center-of-mass stability
- **Geometric Compensation**: Parallel pantograph geometry provides consistent behavior

## Repository Structure

```
balanced-isospring-oscillator/
├── src/
│   ├── equivalent_stiffness.m    # K_eq analysis and visualization
│   ├── reduced_mass.m            # Reduced mass and isotropy analysis
│   └── power_analysis.m          # Power dissipation calculations
├── lib/
│   ├── polarplot3d.m             # 3D polar plotting utility
│   └── polarplot3d_LICENSE.txt   # BSD license for polarplot3d
├── figures/                      # Generated plot images
├── docs/                         # Technical documentation (PDF)
├── LICENSE
└── README.md
```

## Requirements

- MATLAB R2018b or later
- No additional toolboxes required

## Usage

### Equivalent Stiffness Analysis

Computes the equivalent spring stiffness as a function of position in the oscillation plane:

```matlab
cd src
equivalent_stiffness
```

Generates a 3D polar surface plot showing K_eq = f(x, y).

### Reduced Mass Analysis

Calculates the equivalent reduced mass and evaluates mass isotropy:

```matlab
cd src
reduced_mass
```

Generates:
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

Evaluates relationships between gear train stages, mainspring capacity, and power reserve margins.

## System Parameters

### Geometry

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

| Parameter | Min | Max |
|-----------|-----|-----|
| Orbit radius (rho) | 0.3 mm | 0.4 mm |
| Frequency (f) | 1 Hz | 15 Hz |

## Theory

The BISO mechanism achieves 2D isotropy through:

1. **Symmetric Architecture**: Four identical pantograph units arranged at 90° intervals
2. **Balanced Mass Distribution**: Counter-masses on each arm maintain center-of-mass stability
3. **Geometric Compensation**: Parallel pantograph geometry provides consistent restoring force

The equivalent stiffness K_eq is derived from the elastic potential energy:

```
K_eq = (1/r²) × Σ[2·k_γ·(Δγ₁² + Δγ₂²) + 2·k_β·Δβ²]
```

where Δγ and Δβ represent angular deviations from equilibrium.

## Documentation

Technical documentation is available in `docs/`:
- `Pantographe Parallele BISO Balanced-IsoSpring-Oscillator.pdf` - Complete theoretical analysis

## License

MIT License - Copyright (c) 2022 Daniel Abraham Elmaleh

See [LICENSE](LICENSE) for details.

**Third-Party**: `polarplot3d.m` by Ken Garrard (BSD license) - see `lib/polarplot3d_LICENSE.txt`

## Author

**Daniel Abraham Elmaleh**

Developed as part of mechanical design coursework (CdM2 - Groupe 33).
