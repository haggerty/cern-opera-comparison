# sPHENIX Solenoid Field Map Comparison: CERN Measurement vs. OPERA Calculation

Comparison of the measured sPHENIX solenoid magnetic field map (from CERN acceptance testing) against the field map used in sPHENIX offline tracking analysis.  The major physics question addressed is: **is the measured field consistent with the calculated field, up to a small rigid-body translation and/or rotation of the solenoid?**

## Data Sources

### Measured map (CERN, 2018)
- Repository: [`haggerty/sphenix-cernfinal-map`](https://github.com/haggerty/sphenix-cernfinal-map) (git submodule)
- Raw data: two CSV files in the surveyor coordinate system
  - Fine map: dense sampling near the solenoid bore
  - Rough map: extended z coverage outside the fine-map range
- C++ class `sPHENIXFieldMap` (in the submodule) reads the CSVs, applies the surveyor→sPHENIX coordinate transformation, azimuthally averages onto a regular (r, z) grid (NR=37, 25 mm steps; NZ=241, 20 mm steps), and enforces ∇·B = 0 to obtain a Maxwell-consistent Br.
- Coordinate transform: surveyor (x, y, z) → sPHENIX (x, −z, y), so the sPHENIX beam axis is the surveyor −y axis.

### Analysis map (OPERA calculation)
- File: `8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root`
  (18 MB ROOT file, not tracked in git — place in `fieldmap-used-in-analysis/`)
- Format: ROOT `TNtuple` named `fieldmap` with branches `x, y, z` (cm) and `bx, by, bz` (T)
- Uniform Cartesian grid: x, y, z ∈ [−110, +110] cm, 2 cm steps (111³ = 1,367,631 points)
- Sort order: x slowest, y middle, z fastest

## Comparison Methodology

The macro `compareFieldMaps.C` is a self-contained ROOT macro compiled with ACLiC.  It `#include`s `sphenix-cernfinal-map/sPHENIXFieldMap.cxx` so that both the field map class and the comparison code are compiled as a single translation unit.

**Comparison grid:** cylindrical (r, φ, z)

| Dimension | Range | Step | Points |
|-----------|-------|------|--------|
| r | 0–75 cm | 2.5 cm | 31 |
| φ | 0°–350° | 10° | 36 |
| z | −110–+110 cm | 2 cm | 111 |

**Step 1 — Interpolation**

At each (r, φ, z) grid point the analysis map is queried by trilinear interpolation on the Cartesian grid after converting (r, φ, z) → (x, y, z) = (r cos φ, r sin φ, z).  The measured map is queried via `sPHENIXFieldMap::GetField()` (bilinear interpolation on the (r, z) grid; the measured map has been azimuthally averaged, so it has no φ dependence).

Cylindrical components are obtained from the Cartesian interpolation as:

```
Br   =  Bx cos φ + By sin φ
Bφ   = −Bx sin φ + By cos φ
```

**Step 2 — Azimuthal Fourier decomposition**

At each (r, z) the 36-point φ series is decomposed into Fourier modes:

```
B(φ) = m0 + a1 cos φ + b1 sin φ + …
m0  = (1/N) Σ B(φ_i)                   [azimuthal average]
a1  = (2/N) Σ B(φ_i) cos φ_i
b1  = (2/N) Σ B(φ_i) sin φ_i
A1  = √(a1² + b1²)                      [m=1 amplitude]
φ1  = atan2(b1, a1)                     [m=1 phase]
```

The m=0 average of the analysis map is then compared directly with the measured (azimuthally averaged) map.  The m=1 component of the analysis map is used as a diagnostic for solenoid misalignment.

**Step 3 — Maxwell residual diagnostics**

For each map the macro evaluates ∇·B and ∇×B numerically on the comparison grid.

*Analysis map*: central finite differences with step h = 2 cm (one grid cell) evaluated at each (r, φ, z) point and averaged over φ:

```
∇·B  = ∂Bx/∂x + ∂By/∂y + ∂Bz/∂z
∇×B  = (∂Bz/∂y−∂By/∂z, ∂Bx/∂z−∂Bz/∂x, ∂By/∂x−∂Bx/∂y)
```

*Measured map*: cylindrical identities evaluated on the internal (r, z) grid:

```
∇·B      = (1/r)∂(rBr)/∂r + ∂Bz/∂z
(∇×B)_φ  = ∂Br/∂z − ∂Bz/∂r
```

**Step 4 — Misalignment diagnostics**

A solenoid tilted by angle θ around the x-axis produces a first-harmonic modulation in Bz with:

```
Bz,m=1(r, z) ≈ θ r ∂²Bz/∂z²
```

which is proportional to r (A1/r ≈ constant) and has a fixed phase across the (r, z) plane equal to the tilt direction.  A pure transverse translation produces a Bφ m=1 modulation that is r-independent and a Br m=1 proportional to ∂Bz/∂z.

## Key Findings

### Z-axis offset
The peak on-axis Bz occurs at **z = −24 cm** in the measured map and **z = +4 cm** in the analysis map — a **28 cm relative offset**.  This is the dominant discrepancy in the comparison.  The sPHENIX offline reconstruction does not currently apply any correction for this offset.

### Field amplitude
- Measured peak: **1.3975 T** (on-axis)
- Analysis peak: **1.3848 T** (~0.9% lower)

### Central-region agreement
Within |z| < 60 cm (the tracking volume), the maximum azimuthally-averaged difference |ΔBz| = 35.5 mT — good agreement once the z-offset is accounted for qualitatively.

### Bφ
Essentially zero in both maps (< 0.1 mT in the tracking volume center), consistent with azimuthal symmetry.

### Effective tilt in the analysis map (m=1 diagnostic)
The analysis map's Bz m=1 amplitude scales linearly with r: **A1/r ≈ 0.0095 mT/cm** (constant across all r at a given z), and the m=1 phase is constant at **−86°** across the (r, z) plane.

This is the signature of a rigid-body tilt.  The implied effective tilt is **~4.1 mrad** toward **−y** (phase −86°, i.e., nearly the −y direction).

For comparison, the physical sPHENIX solenoid was surveyed to be tilted **2.39 mrad** toward **+y** (phase +71°) — roughly opposite in direction.  The calculated (OPERA) field does not appear to incorporate the survey-measured tilt.

### Maxwell residuals

The macro evaluates ∇·B and ∇×B numerically for each map at all (r > 0) points in the tracking volume.

| Map | |∇·B| max | |∇·B| RMS | |∇×B| max | |∇×B| RMS |
|-----|-----------|-----------|-----------|-----------|
| Measured | 0.017 mT/cm | 0.003 mT/cm | 0.302 mT/cm | 0.065 mT/cm |
| Analysis (OPERA) | 334 mT/cm | 43 mT/cm | 35 mT/cm | 2.6 mT/cm |

**Measured map**: ∇·B ≈ 0 to better than 0.02 mT/cm — consistent with machine precision — confirming that `EnforceMaxwell()` enforces the divergence condition correctly on the discrete grid.  The non-zero curl (max 0.3 mT/cm, RMS 0.065 mT/cm) is real: the measured Bz contains field structure that is not curl-free at the ~0.05% level, from a combination of measurement noise and genuine field non-idealities.

**Analysis map**: The large divergence (max 334 mT/cm, RMS 43 mT/cm) is a direct consequence of storing Bx, By, Bz as three **independent** trilinear grids.  Trilinear interpolation of each component separately imposes no constraint coupling them, so ∇·B ≠ 0 wherever the field curvature is significant.  As a cross-check, the curl is much smaller (max 35 mT/cm, RMS 2.6 mT/cm), consistent with the OPERA FEM solution respecting ∇×B = 0 before discretisation.  The divergence violation means the stored analysis map is not a valid Maxwell field; it violates ∇·B = 0 at the ~3% level relative to the natural scale ∂Bz/∂z ~ 1.4 T / 2 cm ≈ 700 mT/cm.

The practical implication is that tracking through the analysis field map accumulates a systematic error proportional to the divergence residual.  A divergence-free representation (e.g., storing a vector potential, or post-processing the stored grid to enforce ∇·B = 0 via the same integral method used for the measured map) would reduce this error.

## How to Run

### Prerequisites
- ROOT (≥ 6.20) with ACLiC
- The `sphenix-cernfinal-map` submodule initialized:
  ```
  git submodule update --init
  ```
- The analysis ROOT file placed at:
  ```
  fieldmap-used-in-analysis/8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root
  ```

### Running the comparison
```bash
mkdir -p plots
root -l -b -q 'compareFieldMaps.C+'
```

Output: eleven PDF plots in `plots/` and `plots/comparison_histograms.root`.

### Output plots

| File | Contents |
|------|----------|
| `01_onaxis_Bz.pdf` | On-axis Bz vs z: measured and analysis overlaid |
| `02_Bz_2d_maps.pdf` | 2D (r, z) maps of Bz for each source |
| `03_dBz_m0.pdf` | ΔBz (azimuthal average) over the (r, z) plane |
| `04_Br_comparison.pdf` | Br comparison (measured Maxwell-consistent vs analysis m=0) |
| `05_m1_amplitudes.pdf` | m=1 amplitudes A1 for Bz, Br, Bφ |
| `06_m1_amplitude_and_phase.pdf` | m=1 amplitude and phase for Bz and Bφ |
| `07_Bphi_overview.pdf` | Bφ m=0 and m=1 overview |
| `08_Bz_radial_profiles.pdf` | Bz radial profiles at selected z slices |
| `09_phi_dependence_r30_z0.pdf` | φ dependence at r=30 cm, z=0 |
| `10_dBz_phi0_vs_phi180.pdf` | ΔBz at φ=0° vs φ=180° (asymmetry check) |
| `11_maxwell_residuals.pdf` | ∇·B and ∇×B residuals for both maps |

## Repository Structure

```
cern-opera-comparison/
├── compareFieldMaps.C              # Main comparison macro
├── sphenix-cernfinal-map/          # git submodule (measured field map class + CSVs)
├── fieldmap-used-in-analysis/      # (not tracked) place analysis ROOT file here
└── plots/                          # (not tracked) output PDFs and histograms
```

## Dependencies

- [ROOT](https://root.cern) — data containers, interpolation, histogramming
- [`haggerty/sphenix-cernfinal-map`](https://github.com/haggerty/sphenix-cernfinal-map) — `sPHENIXFieldMap` C++ class and raw CSV data
- sPHENIX tracking field map: `8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root` (not publicly distributed; contact sPHENIX offline group)
