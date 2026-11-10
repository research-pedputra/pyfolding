
# Quantum ESPRESSO Band-Unfolding and Spectral-Function Plotter

## Overview

This script reads a Quantum ESPRESSO (QE) band-structure calculation performed for a supercell, projects its supercell Bloch states onto the reciprocal-space basis of a chosen primitive cell, and constructs an unfolded spectral function along a primitive-cell k-path.

The resulting intensity map displays where supercell eigenstates retain primitive-cell Bloch character. Each eigenstate contributes with a projection weight that indicates how strongly it belongs to a particular primitive-cell wavevector.

The script aligns the energy axis to the Fermi level, applies Lorentzian broadening, plots the spectral function, and saves the numerical data used to make the plot.

## 1. Physical Background: Why Band Unfolding Is Needed

A primitive crystal has translational symmetry under its primitive lattice vectors. A supercell calculation is often used to model defects, vacancies, dopants, disorder approximants, or magnetic ordering.

The supercell is larger than the primitive cell, so its Brillouin zone is correspondingly smaller. The supercell calculation consequently contains many bands that appear folded into the reduced Brillouin zone.

A conventional supercell band plot shows all these folded bands with equal visual status. It does not directly distinguish states that preserve a well-defined primitive-cell wavevector from states whose primitive translational character is weak or mixed.

Band unfolding addresses this by assigning a spectral weight to each supercell eigenstate at each primitive-cell wavevector.

## 2. Lattice and Reciprocal-Space Mapping

Let the primitive lattice vectors be the rows of the matrix:

\[
A_p =
\begin{pmatrix}
\mathbf a_1\\
\mathbf a_2\\
\mathbf a_3
\end{pmatrix}
\]

The supercell is defined by an integer transformation matrix:

\[
A_s = M A_p
\]

The script accepts the nine integer elements of \(M\) through `--matrix` and reshapes them into a \(3\times3\) matrix. The matrix must be nonsingular.

Its determinant gives the ratio of supercell to primitive-cell volumes:

\[
N = |\det M|
\]

For row-vector fractional reciprocal coordinates, the reciprocal-space mapping used by the script is:

\[
\mathbf K_{\mathrm{raw}} = \mathbf k_p M^T
\]

where \(\mathbf{k}_p\) is a primitive-cell k-point and \(\mathbf{K}_{\mathrm{raw}}\) is its corresponding supercell-coordinate representation before reduction by a reciprocal lattice vector.

The script compares this target with the k-point stored in each QE wavefunction file. Their difference should be an integer reciprocal-lattice translation:

\[
\mathbf T =
\operatorname{round}
\left(
\mathbf K_{\mathrm{raw}}-\mathbf k_{\mathrm{QE}}
\right)
\]

A tolerance check verifies that the difference is sufficiently close to an integer vector.

The reciprocal transformation used to identify primitive reciprocal vectors is:

\[
M^{-T} = (M^{-1})^T
\]

For each supercell plane-wave reciprocal index \(\mathbf G_s\), read from the QE `mill` array, the script evaluates:

\[
\mathbf G_p =
(\mathbf G_s-\mathbf T)M^{-T}
\]

A plane wave is selected when all three components of \(\mathbf G_p\) are sufficiently close to integers.

## 3. Unfolding Spectral Weight

At a fixed supercell k-point, a Kohn–Sham Bloch state can be expanded in plane waves:

\[
|\Psi_{n\mathbf K}\rangle =
\sum_{\mathbf G_s}
C_{n\mathbf K}(\mathbf G_s)
|\mathbf K+\mathbf G_s\rangle
\]

Here, \(n\) is the band index, \(\mathbf K\) is the supercell wavevector, and \(C_{n\mathbf K}(\mathbf G_s)\) is the complex plane-wave coefficient read from the QE portable wavefunction file.

The primitive-cell projection weight used in this script is:

\[
W_n(\mathbf k_p) =
\frac{
\displaystyle
\sum_{\mathbf G_s\in\mathcal S(\mathbf k_p)}
|C_{n\mathbf K}(\mathbf G_s)|^2
}{
\displaystyle
\sum_{\mathbf G_s}
|C_{n\mathbf K}(\mathbf G_s)|^2
}
\]

where \(\mathcal S(\mathbf k_p)\) is the set of plane waves selected by the integer reciprocal-space condition above.

The denominator normalizes by the total plane-wave coefficient norm for that band.

The weight is dimensionless and is expected to lie between zero and one for a correctly normalized wavefunction and a valid selected subset, subject to numerical and file-format conventions.

A weight near one means that most of the state's plane-wave content belongs to the selected primitive-cell translational sector. A small weight means that little of the state's plane-wave content projects onto that primitive wavevector.

This implementation calculates the weight directly from the squared magnitudes of the complex plane-wave coefficients:

\[
|C|^2=C^*C
\]

## 4. Spectral Function and Lorentzian Broadening

The ideal unfolded spectrum is a set of weighted delta functions:

\[
A(\mathbf k,E) =
\sum_n W_n(\mathbf k)
\delta(E-\epsilon_{n\mathbf K})
\]

where \(\epsilon_{n\mathbf K}\) is the eigenvalue of band \(n\) at the corresponding supercell k-point.

The script replaces each delta function with a Lorentzian:

\[
L(E-\epsilon_n;\eta) =
\frac{1}{\pi}
\frac{\eta}
{(E-\epsilon_n)^2+\eta^2}
\]

The broadened spectral function is:

\[
\boxed{
A(\mathbf k,E) =
\sum_n W_n(\mathbf k)
\frac{1}{\pi}
\frac{\eta}
{(E-\epsilon_{n\mathbf K})^2+\eta^2}
}
\]

The Lorentzian half-width at half-maximum is \(\eta\), and its full width at half-maximum is:

\[
\mathrm{FWHM}=2\eta
\]

The command-line option `--width-ev` specifies the FWHM in eV. The code converts it to Hartree and sets:

\[
\eta=\frac{\mathrm{FWHM}}{2}
\]

The default is \(0.002\) Hartree, equivalent to approximately \(0.0544\) eV FWHM.

A smaller width produces sharper spectral features; a larger width produces smoother, broader features. Broadening is a visualization parameter and should not automatically be interpreted as a physical scattering lifetime.

### Energy Units and Fermi-Level Alignment

QE eigenvalues and Fermi energies are read from `data-file-schema.xml` in Rydberg.

The conversion is:

\[
1\,\mathrm{Ha}=2\,\mathrm{Ry}
\]

and:

\[
1\,\mathrm{Ha}=27.211386245988\,\mathrm{eV}
\]

The plotted energy is shifted relative to the Fermi energy:

\[
\boxed{
E_{\mathrm{plot}} =
(E-E_F)\times27.211386245988
}
\]

Consequently, the Fermi level is at zero eV on the vertical axis.

The plotted range is controlled by `--emin` and `--emax`; these limits crop the displayed energy window and do not change the eigenvalues or the spectral-function calculation itself.

The Lorentzian is evaluated in Hartree, so the stored spectral function has units of inverse Hartree under this convention.

## 5. Primitive K-Path and Horizontal Axis

The input path file is a CSV containing three fractional primitive reciprocal coordinates per non-comment row:

\[
(k_1,k_2,k_3)
\]

The script ignores blank rows and lines whose first field begins with `#`. Each non-comment row must contain exactly three numeric values.

To display the path as a continuous horizontal coordinate, fractional k-points are converted to Cartesian reciprocal coordinates. The distance between consecutive points is:

\[
\Delta s_i =
\left\|
\mathbf k_{i+1}^{\mathrm{cart}}
-
\mathbf k_i^{\mathrm{cart}}
\right\|
\]

The cumulative path coordinate is:

\[
s_0=0,\qquad
s_i=\sum_{j=0}^{i-1}\Delta s_j
\]

This cumulative distance is the horizontal coordinate in the spectral-function plot. It is a reciprocal-space path length, not a real-space distance.

High-symmetry-point labels can be added using `--tick-indices` and `--tick-labels`. The indices refer to rows in the input path file, starting at zero.

## 6. What the Code Does, Step by Step

1. **Read the primitive path.** Reads the three fractional reciprocal coordinates from the CSV path file.
2. **Read QE output files.** Locates `data-file-schema.xml` and portable wavefunction files (`wfc*.dat`).
3. **Read eigenvalues and Fermi energy.** Extracts the Kohn–Sham eigenvalues and Fermi energy from the XML file.
4. **Read plane-wave coefficients.** Reads the reciprocal vectors and complex coefficients from each WFC file.
5. **Map the k-points.** Maps primitive-cell k-points to corresponding supercell k-points using the supercell transformation matrix.
6. **Select compatible plane waves.** Identifies reciprocal vectors belonging to the primitive-cell translational sector.
7. **Calculate unfolding weights.** Calculates the normalized squared coefficient sum for each band.
8. **Construct the spectral function.** Broadens each eigenvalue with a Lorentzian and multiplies it by its unfolding weight.
9. **Align the energy axis.** Subtracts the Fermi energy and converts the result to eV.
10. **Plot and save.** Saves the spectral intensity map and an NPZ archive containing the numerical data.

## 7. Input Files and Directory Structure

Expected QE output structure:

```text
<outdir>/
└── <prefix>.save/
    ├── data-file-schema.xml
    ├── wfc1.dat
    ├── wfc2.dat
    ├── wfc3.dat
    └── ...
```

The primitive path is supplied separately as a CSV file.

Example:

```csv
# primitive fractional reciprocal coordinates
0.000000,0.000000,0.000000
0.125000,0.000000,0.125000
0.250000,0.000000,0.250000
0.375000,0.000000,0.375000
```

## 8. Running the Script

Save the Python code as:

```bash
unfold_qe.py
```

A representative command is:

```bash
python3 unfold_qe.py \
    --prefix aiida \
    --outdir . \
    --path primitive_path.csv \
    --matrix 2 0 0 0 2 0 0 0 2 \
    --spin up \
    --emin -5 \
    --emax 4 \
    --output unfolded_QE_up.png
```

For a \(2\times2\times2\) diagonal supercell, the matrix is:

\[
M=
\begin{pmatrix}
2&0&0\\
0&2&0\\
0&0&2
\end{pmatrix}
\]

For the spin-down channel, use:

```bash
--spin down
```

To provide high-symmetry labels, add indices and labels corresponding to the appropriate rows in the path file:

```bash
--tick-indices 0 20 40 60 80 \
--tick-labels G M R X X1
```

The values above are illustrative and must be replaced with the actual indices and labels for the path being used.

### Main Command-Line Options

| Option | Meaning | Default |
|---|---|---|
| `--prefix` | QE calculation prefix | Required |
| `--outdir` | QE output directory | Required |
| `--path` | CSV file of primitive fractional k-points | Required |
| `--matrix` | Nine integer elements of the \(3\times3\) supercell matrix | Required |
| `--spin` | Spin channel: `up` or `down` | `up` |
| `--fermi-ry` | Optional Fermi energy override in Ry | Read from XML |
| `--width-ev` | Lorentzian FWHM in eV | \(0.002\) Ha converted to eV |
| `--nenergy` | Number of energy-grid points | 10000 |
| `--emin` | Lower displayed energy limit relative to \(E_F\), in eV | -5.0 |
| `--emax` | Upper displayed energy limit relative to \(E_F\), in eV | 4.5 |
| `--tick-indices` | Zero-based indices of high-symmetry points | None |
| `--tick-labels` | Labels corresponding to `--tick-indices` | None |
| `--output` | Output PNG filename | `unfolded_QE_<spin>_GPAW_comparison.png` |

## 9. Output Files

The script produces two files.

### Unfolded Spectral-Function Plot

A PNG image containing the spectral function \(A(\mathbf k,E)\). The horizontal axis is cumulative reciprocal-space path distance, and the vertical axis is energy relative to the Fermi level in eV. The colormap represents spectral intensity. A dashed horizontal line marks the Fermi level at zero eV.

### Numerical NPZ Archive

The script also saves a compressed NumPy archive alongside the plot. It contains:

- `k_primitive`: primitive fractional k-points.
- `path_distance`: cumulative reciprocal-space path coordinate.
- `energy_ev`: energy grid relative to the Fermi level, in eV.
- `energy_ha`: absolute energy grid used in the Lorentzian calculation, in Hartree.
- `spectral_function`: unfolded spectral function array.
- `weights`: unfolding weights for each k-point and band.
- `eigenvalues_ha`: QE eigenvalues in Hartree.
- `fermi_ha`: Fermi energy in Hartree.
- `width_ha`: Lorentzian FWHM in Hartree.
- `width_ev`: Lorentzian FWHM in eV.

The archive can be loaded later without rerunning QE:

```python
import numpy as np

data = np.load("unfolded_QE_up_data.npz")

k_distance = data["path_distance"]
energy = data["energy_ev"]
spectral_function = data["spectral_function"]
weights = data["weights"]
```

## 10. Requirements and Implementation Assumptions

The script requires Python, NumPy, Matplotlib, and the QE output files. It uses standard-library modules for argument parsing, CSV reading, file matching, XML parsing, and path handling.

Important assumptions and limitations:

- QE must provide portable plane-wave wavefunction files named `wfc*.dat`.
- The binary reader assumes a specific portable WFC record layout. QE versions, compilation options, endianness, or file-format changes may require changes to the reader.
- The implementation explicitly rejects spinor wavefunctions (`npol != 1`) and Gamma-only wavefunction files.
- The path points, supercell matrix, QE k-point sampling, and wavefunction-file ordering must describe the same calculation.
- The script expects the XML eigenvalue list and WFC band count to match at every k-point.
- The current weight calculation uses plane-wave coefficient norms. It does not include a separate real-space projection or an explicit orbital-character decomposition.
- The energy broadening is Lorentzian and is a visualization parameter.
- `--emin` and `--emax` control the plotted window; the energy grid itself is generated from the calculated eigenvalue range and broadening width.
- Comparisons with another code require consistent structural models, lattice transformations, k-paths, spin channels, exchange-correlation functionals, Hubbard parameters, pseudopotentials, Fermi-level references, band ranges, and broadening conventions.

## 11. Interpreting the Unfolded Plot

Bright regions indicate energies and primitive k-points with relatively large unfolded spectral weight.

A strong, narrow feature that follows a smooth dispersion is characteristic of a state with well-defined primitive-cell Bloch character. Weak or diffuse intensity can indicate that a supercell eigenstate has only a small projection onto the selected primitive wavevector.

In doped or defective systems, additional features may arise from impurity-derived states, band splitting, or mixing between primitive-cell wavevectors.

The plot is a spectral-weight representation, not an ordinary band plot in which every eigenvalue is drawn as an equally weighted line. Therefore, a missing or weak feature should not automatically be interpreted as the absence of an electronic state: it may have low projection weight for the chosen primitive k-point.

## 12. Comparison with GPAW

The script includes defaults intended to make the energy broadening and spectral-function display comparable to a GPAW unfolding workflow:

- Default Lorentzian FWHM: \(0.002\) Hartree.
- Default energy grid: 10000 points.
- Fermi level placed at zero.
- `viridis` colormap with `PowerNorm(gamma=0.4)`.

These matching plotting conventions do not by themselves guarantee numerically identical QE and GPAW spectra.

A meaningful comparison additionally requires consistent crystal structures, supercell definitions, k-point paths, spin treatment, exchange-correlation and Hubbard settings, pseudopotential/PAW datasets, eigenvalue windows, and unfolding-weight definitions.

Differences in basis sets, wavefunction normalization, projection conventions, and implementation details can remain even when the plotting settings match.

## Summary

The script converts QE supercell eigenvalues and plane-wave coefficients into a primitive-cell unfolded spectral function. It identifies plane-wave components compatible with each primitive k-point, sums their normalized squared amplitudes into an unfolding weight, broadens the weighted eigenvalues with Lorentzians, aligns the result to the Fermi level, and saves both a spectral map and its underlying numerical data.
