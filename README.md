# sPHENIX Solenoid Field Map Comparison: CERN Measurement vs. OPERA Calculation

Comparison of the measured sPHENIX solenoid magnetic field map (from CERN acceptance testing) against the OPERA field map used in sPHENIX offline tracking.  The major physics question addressed is: **is the measured field consistent with the OPERA field, up to a small rigid-body translation and/or rotation of the solenoid?**

## Data Sources

### Measured map (CERN, 2022)
- Repository: [`haggerty/sphenix-cernfinal-map`](https://github.com/haggerty/sphenix-cernfinal-map) (git submodule)
- Raw data: two CSV files in the surveyor coordinate system
  - Fine map: dense sampling near the solenoid bore
  - Rough map: extended z coverage outside the fine-map range
- C++ class `sPHENIXFieldMap` (in the submodule) reads the CSVs, applies the surveyor→sPHENIX coordinate transformation, azimuthally averages onto a regular (r, z) grid (NR=37, 25 mm steps; NZ=241, 20 mm steps), and enforces ∇·B = 0 to obtain a Maxwell-consistent Br.
- Coordinate transform: surveyor (x, y, z) → sPHENIX (x, −z, y), so the sPHENIX beam axis is the surveyor −y axis.

### OPERA map
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

At each (r, φ, z) grid point the OPERA map is queried by trilinear interpolation on the Cartesian grid after converting (r, φ, z) → (x, y, z) = (r cos φ, r sin φ, z).  The measured map is queried via `sPHENIXFieldMap::GetField()` (bilinear interpolation on the (r, z) grid; the measured map has been azimuthally averaged, so it has no φ dependence).

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

The m=0 average of the OPERA map is compared directly with the measured (azimuthally averaged) map.  The m=1 component of the OPERA map is used as a diagnostic for solenoid misalignment.

**Step 3 — Maxwell residual diagnostics**

For each map the macro evaluates ∇·B and ∇×B numerically on the comparison grid.

*OPERA map*: central finite differences with step h = 2 cm (one grid cell) evaluated at each (r, φ, z) point and averaged over φ:

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
The peak on-axis Bz occurs at **z = −24 cm** in the measured map and **z = +4 cm** in the OPERA map — a **28 cm relative offset**.  This is the dominant discrepancy in the comparison.  The sPHENIX offline reconstruction does not currently apply any correction for this offset.

### Field amplitude
- Measured peak: **1.3975 T** (on-axis)
- OPERA peak: **1.3848 T** (~0.9% lower)

### Central-region agreement
Within |z| < 60 cm (the tracking volume), the maximum azimuthally-averaged difference |ΔBz| = 35.5 mT — good agreement once the z-offset is accounted for qualitatively.

### Br sign change at the Bz peak
In both maps, Br changes sign at the same z as the Bz peak, to within a few µT.

| Map | Bz peak | Br zero-crossing | Coincidence |
|-----|---------|-----------------|-------------|
| Measured | z = −24 cm, 1.39754 T | z = −24 cm | within 1 µT |
| OPERA | z = +4 cm, 1.38486 T | z = +4 cm | within 4 µT |

For the measured map this is required by construction: `EnforceMaxwell()` derives Br from the integral of ∂Bz/∂z, so Br = 0 wherever ∂Bz/∂z = 0 — i.e., at the Bz maximum.  At r = 25 mm: Br = −0.000032 T at z = −26 cm, ≈ 0 at z = −24 cm, and +0.000029 T at z = −22 cm.

For the OPERA map the same coincidence holds (Br = −0.000026 T at z = +2 cm, ≈ 0 at z = +4 cm, +0.000033 T at z = +6 cm) even though OPERA violates ∇·B = 0 at the ~40 mT/cm level overall.  This is because the ∇·B violations arise from curvature terms (second derivatives), while the Br zero-crossing at the Bz peak is a first-derivative condition — Br ∝ ∫∂Bz/∂z dr changes sign where ∂Bz/∂z = 0.  That relationship is preserved even when the radial balance of ∂(rBr)/∂r and ∂Bz/∂z is imperfect.

### Bφ
Essentially zero in both maps (< 0.1 mT in the tracking volume center), consistent with azimuthal symmetry.

### Solenoid axis tilt

**Method (`checkTilt.C`):** reads the raw measurement CSVs directly (before azimuthal averaging) and fits Bz(φ) = A₀ + A₁ cos φ + B₁ sin φ per (r, z) bin by least squares.  The fine map provides 18–31 azimuthal measurements per bin — sufficient for a clean m=1 extraction.

**Tilt signature:** for a solenoid axis tilted by angle α toward azimuth φ₀ the first-order m=1 perturbation is:
```
ΔBz ≈ −(∂Bz/∂r) · (z − z_center) · α · cos(φ − φ₀)
```
Key discriminants vs. a pure axis translation: (1) the amplitude reverses sign at z = z_center (the field maximum, where ∂Bz/∂r → 0); (2) it grows with |z − z_center|.  A translation would give constant phase and amplitude with no sign reversal.

**Measured data result at r = 300 mm:**

| z (cm) | A₀ (T) | A₁ (mT) | A₁/A₀ (%) | phase (°) |
|--------|--------|---------|-----------|----------|
| −100 | 1.354 | 2.38 | 0.18 | −67 |
| −50  | 1.397 | 0.25 | 0.02 | −65 |
| −24  | 1.401 | 0.40 | 0.03 | +44 (near zero-crossing) |
|   0  | 1.398 | 0.68 | 0.05 | +64 |
| +50  | 1.359 | 2.42 | 0.18 | +76 |
| +100 | 1.216 | 4.77 | 0.39 | +84 |
| +140 | 1.006 | 4.82 | 0.48 | +101 |

The amplitude minimum is at z ≈ −24 cm (the Bz peak), as expected.  For z > 0 the phase is consistently **+75° ± 15°**.  For z < −50 cm the phase is **−70° ± 10°** — reversed by ~145°, close to the 180° flip expected for a tilt.  Both the sign reversal and the growing amplitude with |z − z_center| confirm a rigid-body tilt rather than a translation.

**OPERA map result:** the OPERA Bz m=1 amplitude scales linearly with r — A1/r ≈ **0.0095 mT/cm**, constant phase **−86°** — also a tilt signature, implying ~4.1 mrad toward −y.

**Comparison:**

| Source | m=1 phase | Implied direction | Amplitude |
|--------|-----------|-------------------|-----------|
| Physical survey | +71° | +y (upward) | 2.39 mrad |
| Raw measured data (z > 0) | **+75°** | +y ✓ | consistent |
| OPERA calculation | **−86°** | −y ✗ | ~4.1 mrad spurious |

The OPERA dipole points **157° away** from the survey direction and does not reproduce the physical tilt.  The sPHENIXFieldMap discards φ information during loading, so the measured Cartesian map (and tracking using it) is azimuthally symmetric — the tilt-induced ~0.2–0.5% Bz asymmetry is not preserved.

### Maxwell residuals

The macro evaluates ∇·B and ∇×B numerically for each map at all (r > 0) points in the tracking volume.

| Map | |∇·B| max | |∇·B| RMS | |∇×B| max | |∇×B| RMS |
|-----|-----------|-----------|-----------|-----------|
| Measured | 0.017 mT/cm | 0.003 mT/cm | 0.302 mT/cm | 0.065 mT/cm |
| OPERA | 334 mT/cm | 43 mT/cm | 35 mT/cm | 2.6 mT/cm |

**Measured map**: ∇·B ≈ 0 to better than 0.02 mT/cm — consistent with machine precision — confirming that `EnforceMaxwell()` enforces the divergence condition correctly on the discrete grid.  The non-zero curl (max 0.3 mT/cm, RMS 0.065 mT/cm) is real: the measured Bz contains field structure that is not curl-free at the ~0.05% level, from a combination of measurement noise and genuine field non-idealities.

**OPERA map**: The large divergence (max 334 mT/cm, RMS 43 mT/cm) is a direct consequence of storing Bx, By, Bz as three **independent** trilinear grids.  Trilinear interpolation of each component separately imposes no constraint coupling them, so ∇·B ≠ 0 wherever the field curvature is significant.  As a cross-check, the curl is much smaller (max 35 mT/cm, RMS 2.6 mT/cm), consistent with the OPERA FEM solution respecting ∇×B = 0 before discretisation.  The divergence violation means the stored OPERA map is not a valid Maxwell field; it violates ∇·B = 0 at the ~3% level relative to the natural scale ∂Bz/∂z ~ 1.4 T / 2 cm ≈ 700 mT/cm.

The practical implication is that tracking through the OPERA map accumulates a systematic error proportional to the divergence residual.

## Measured Map in OPERA Format

The macro `makeMeasuredCartesianMap.C` writes the measured map onto the same 111³ Cartesian grid as the OPERA file, producing a drop-in replacement for use in PHField3DCartesian.

```bash
mkdir -p output
root -l -b -q 'makeMeasuredCartesianMap.C+'
```

Output: `output/sphenix_measured_fieldmap_cartesian.root` — same `fieldmap` TNtuple, same branches (`x:y:z:bx:by:bz:hz`), same units (cm, T), same grid layout.

At each grid node `(x, y, z)` the macro calls `sPHENIXFieldMap::GetFieldXYZ(x·10, y·10, z·10, Bx, By, Bz)` (converting cm to mm), so the output inherits the physical properties of the measured map:
- Ground-truth hardware measurement from CERN acceptance testing
- Correct axial position (peak at z = −24 cm) and field amplitude (1.3975 T)
- Azimuthally symmetric: Bφ = 0
- No spurious tilt (m=1 amplitudes < 0.01 mT throughout the tracking volume)

Of 1,367,631 grid points, 661,560 (48%) lie at r > 900 mm and receive zero field — these are corners of the Cartesian cube well outside the tracking volume (r < ~70 cm).

### Self-consistency check

`compareFieldMaps.C` was run with the measured Cartesian map in place of OPERA to verify the conversion and characterise the round-trip interpolation error (cylindrical bilinear → Cartesian grid → trilinear).  Results in `plots_meas_cartesian/`.

| Quantity | OPERA vs Measured | Meas. Cartesian vs Measured |
|---|---|---|
| On-axis ΔBz at z = 0 | −9.1 mT | 0.0 mT |
| Max \|ΔBz\| in tracking volume | 105.5 mT | **0.3 mT** |
| Max \|ΔBr\| | 74.0 mT | **0.0 mT** |
| Bz m=1 amplitude | 0.93 mT | **0.00 mT** |
| \|∇·B\| max | 334 mT/cm | 349 mT/cm |
| \|∇×B\| max | 35 mT/cm | 49 mT/cm |

The 0.3 mT max ΔBz is the round-trip interpolation error; all m=1 amplitudes are consistent with noise (phase spread 105°, only 68 cells above threshold), confirming there is no spurious azimuthal structure in the converted map.

The ∇·B residual of the measured Cartesian map (349 mT/cm) is essentially the same as OPERA (334 mT/cm).  This is expected: storing any smooth field on independent Cartesian trilinear grids introduces ∇·B violations at this level regardless of whether the source was Maxwell-consistent in cylindrical coordinates.  The cylindrical divergence-free property does not survive resampling onto a Cartesian grid with independent component interpolation.  A truly divergence-free Cartesian representation would require a different storage scheme (staggered grid, vector potential, or post-processing to enforce ∇·B = 0 in the Cartesian basis).

## How to Run

### Prerequisites
- ROOT (≥ 6.20) with ACLiC
- The `sphenix-cernfinal-map` submodule initialized:
  ```
  git submodule update --init
  ```
- The OPERA ROOT file placed at:
  ```
  fieldmap-used-in-analysis/8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root
  ```

### Running the comparison

`compareFieldMaps()` accepts three optional arguments:

| Argument | Default | Description |
|---|---|---|
| `calcRoot` | OPERA file in `fieldmap-used-in-analysis/` | Path to the Cartesian TNtuple to compare against the measured map |
| `calcLabel` | `"OPERA"` | Label used in histogram titles and legend entries |
| `outDir` | `"plots"` | Output directory for PDFs and ROOT file |

**OPERA vs measured** (default):
```bash
root -l -b -q 'compareFieldMaps.C+'
```

**Measured Cartesian vs measured** (self-consistency check):
```bash
root -l -b -q 'compareFieldMaps.C+("output/sphenix_measured_fieldmap_cartesian.root","Meas. Cartesian","plots_meas_cartesian")'
```

Output: eleven PDF plots and `comparison_histograms.root` in the specified directory.

### Output plots

| File | Contents |
|------|----------|
| `01_onaxis_Bz.pdf` | On-axis Bz vs z: measured and OPERA overlaid |
| `02_Bz_2d_maps.pdf` | 2D (r, z) Bz maps: measured and OPERA |
| `03_dBz_m0.pdf` | ΔBz = OPERA − measured (azimuthal average) |
| `04_Br_comparison.pdf` | Br: measured, OPERA (φ-avg), and difference |
| `05_m1_amplitudes.pdf` | OPERA m=1 amplitudes A1 for Bz, Br, Bφ |
| `06_m1_amplitude_and_phase.pdf` | OPERA m=1 amplitude and phase for Bz and Bφ |
| `07_Bphi_overview.pdf` | OPERA Bφ: azimuthal mean and m=1 amplitude |
| `08_Bz_radial_profiles.pdf` | Bz radial profiles at selected z slices |
| `09_phi_dependence_r30_z0.pdf` | OPERA φ dependence at r=30 cm, z=0 |
| `10_dBz_phi0_vs_phi180.pdf` | ΔBz at φ=0° vs φ=180° (asymmetry check) |
| `11_maxwell_residuals.pdf` | ∇·B and ∇×B residuals for OPERA and measured maps |

### Running the tilt analysis

```bash
root -l -b -q 'checkTilt.C+'
```

Output: six PDFs in `plots_tilt/`.  No OPERA file required — reads only the raw measurement CSVs.

| File | Contents |
|------|----------|
| `tilt_A_nphi_per_bin.pdf` | N(φ) measurements per (r, z) bin — shows azimuthal coverage |
| `tilt_B_m1_frac_map.pdf` | A₁/A₀ in the (z, r) plane |
| `tilt_C_m1_abs_map.pdf` | A₁ in mT in the (z, r) plane |
| `tilt_D_phase_map.pdf` | m=1 phase in the (z, r) plane |
| `tilt_E_phase_vs_z.pdf` | Phase vs z profiles at r = 50, 100, 200, 300, 450 mm with survey and OPERA reference lines |
| `tilt_F_amplitude_vs_z.pdf` | A₁/A₀ vs z profiles at the same r values |

### Producing the measured Cartesian map
```bash
mkdir -p output
root -l -b -q 'makeMeasuredCartesianMap.C+'
```

Output: `output/sphenix_measured_fieldmap_cartesian.root`

## Repository Structure

```
cern-opera-comparison/
├── compareFieldMaps.C              # Comparison macro (parameterised; see How to Run)
├── makeMeasuredCartesianMap.C      # Converts measured map to OPERA Cartesian format
├── checkTilt.C                     # Azimuthal Fourier analysis of raw CSVs for tilt
├── sphenix-cernfinal-map/          # git submodule (sPHENIXFieldMap class + CSV data)
├── plots/                          # OPERA vs measured: PDFs and comparison_histograms.root
├── plots_meas_cartesian/           # Meas. Cartesian vs measured: self-consistency check
├── plots_tilt/                     # Tilt analysis: m=1 amplitude and phase maps
├── fieldmap-used-in-analysis/      # (not tracked) place OPERA ROOT file here
└── output/                         # (not tracked) output of makeMeasuredCartesianMap.C
```

## Dependencies

- [ROOT](https://root.cern) — data containers, interpolation, histogramming
- [`haggerty/sphenix-cernfinal-map`](https://github.com/haggerty/sphenix-cernfinal-map) — `sPHENIXFieldMap` C++ class and raw CSV data
- OPERA field map: `8e4d6c3b1660540a658da3a275af2bde_sphenix3dtrackingmapxyz.root` (not publicly distributed; contact sPHENIX offline group)

## Next Steps

The natural follow-on is to measure the effect on tracking by switching from the OPERA map to the measured field map.  The drop-in replacement (`output/sphenix_measured_fieldmap_cartesian.root`) uses the identical TNtuple format, so it can be pointed at by `PHField3DCartesian` with a one-line config change.  Three differences are expected to matter:

1. **28 cm z-offset** — the measured Bz peak is at z = −24 cm vs z = +4 cm in OPERA.  This is the dominant discrepancy and will shift the sagitta correction for all tracks.
2. **∇·B residuals** — the measured Cartesian map has the same ~350 mT/cm divergence residual as OPERA (inherent to independent-component trilinear storage), so this is a wash unless the field integrator is changed.
3. **Spurious tilt removed** — OPERA contains a ~4.1 mrad dipole toward −y that is absent from the measured map.  Removing it could reduce the φ-dependence of tracking residuals.
