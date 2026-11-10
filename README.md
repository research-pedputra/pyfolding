
# Quantum ESPRESSO Band-Unfolding and Spectral-Function Plotter

## Overview

This script reads a Quantum ESPRESSO (QE) band-structure calculation performed for a supercell, projects its supercell Bloch states onto the reciprocal-space basis of a chosen primitive cell, and constructs an unfolded spectral function along a primitive-cell k-path.

The resulting intensity map displays where supercell eigenstates retain primitive-cell Bloch character. Each eigenstate contributes with a projection weight that indicates how strongly it belongs to a particular primitive-cell wavevector.

The script aligns the energy axis to the Fermi level, applies Lorentzian broadening, plots the spectral function, and saves the numerical data.

## 1. Physical Background

A supercell calculation is often used to model defects, vacancies, dopants, disorder, or magnetic ordering. Because the supercell is larger than the primitive cell, its Brillouin zone is smaller, causing the calculated bands to fold into the reduced Brillouin zone.

Band unfolding assigns a spectral weight to each supercell eigenstate at each primitive-cell wavevector. This allows the electronic structure to be interpreted in terms of the primitive-cell Brillouin zone.

## 2. Lattice and Reciprocal-Space Mapping

The primitive lattice vectors are represented by:

$$
A_p =
\begin{pmatrix}
\mathbf a_1\\
\mathbf a_2\\
\mathbf a_3
\end{pmatrix}
$$

The supercell is defined by an integer transformation matrix:

$$
A_s = M A_p
$$

The determinant of the transformation matrix gives the supercell-to-primitive-cell volume ratio:

$$
N = |\det M|
$$

For row-vector fractional reciprocal coordinates, the mapping used in the script is:

$$
\mathbf K_{\mathrm{raw}} = \mathbf k_p M^T
$$

where $\mathbf k_p$ is a primitive-cell k-point and $\mathbf K_{\mathrm{raw}}$ is its corresponding supercell-coordinate representation.

The reciprocal lattice translation is determined by:

$$
\mathbf{T} = \mathrm{round}\left(\mathbf{K}_{\mathrm{raw}} - \mathbf{k}_{\mathrm{QE}}\right)
$$

The inverse-transpose matrix is:

$$
M^{-T} = (M^{-1})^T
$$

For each supercell reciprocal vector $\mathbf G_s$, the corresponding primitive reciprocal vector is:

$$
\mathbf G_p =
(\mathbf G_s-\mathbf T)M^{-T}
$$

A plane wave is selected when all three components of $\mathbf G_p$ are sufficiently close to integers.

## 3. Unfolding Spectral Weight

A supercell Bloch state can be expanded in plane waves:

$$
|\Psi_{n\mathbf K}\rangle =
\sum_{\mathbf G_s}
C_{n\mathbf K}(\mathbf G_s)
|\mathbf K+\mathbf G_s\rangle
$$

where:

- $n$ is the band index.
- $\mathbf K$ is the supercell wavevector.
- $\mathbf G_s$ is a supercell reciprocal lattice vector.
- $C_{n\mathbf K}(\mathbf G_s)$ is the complex plane-wave coefficient.

The primitive-cell projection weight is:

$$
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
$$

Here, $\mathcal S(\mathbf k_p)$ represents the selected set of plane waves compatible with the primitive-cell reciprocal-space condition.

The denominator normalizes the selected plane-wave contribution by the total plane-wave coefficient norm.

The coefficient magnitude is calculated as:

$$
|C|^2 = C^*C
$$

The unfolding weight is dimensionless. A large weight indicates strong primitive-cell Bloch character at the selected k-point, whereas a small weight indicates weak projection onto that primitive reciprocal-space sector.

## 4. Spectral Function and Lorentzian Broadening

The ideal unfolded spectral function is:

$$
A(\mathbf k,E) =
\sum_n W_n(\mathbf k)
\delta(E-\epsilon_{n\mathbf K})
$$

where $\epsilon_{n\mathbf K}$ is the eigenvalue of band $n$ at the corresponding supercell k-point.

To obtain a continuous intensity map, the delta function is replaced by a Lorentzian:

$$
L(E-\epsilon_n;\eta) =
\frac{1}{\pi}
\frac{\eta}
{(E-\epsilon_n)^2+\eta^2}
$$

The resulting broadened spectral function is:

$$
\boxed{
A(\mathbf k,E) =
\sum_n W_n(\mathbf k)
\frac{\eta}
{\pi\left[(E-\epsilon_{n\mathbf K})^2+\eta^2\right]}
}
$$

The Lorentzian full width at half maximum is:

$$
\mathrm{FWHM}=2\eta
$$

Therefore:

$$
\eta=\frac{\mathrm{FWHM}}{2}
$$

The script accepts the broadening through `--width-ev`, specified in eV, and converts it to Hartree.

The default broadening is $0.002$ Hartree, equivalent to approximately $0.0544$ eV FWHM.

A smaller broadening produces sharper spectral features, whereas a larger broadening produces smoother and wider features. This is a numerical visualization parameter and should not automatically be interpreted as a physical carrier lifetime.

### Energy Conversion and Fermi-Level Alignment

QE eigenvalues and Fermi energies are read from `data-file-schema.xml` in Rydberg.

The conversion factors are:

$$
1\ \mathrm{Ha}=2\ \mathrm{Ry}
$$

$$
1\ \mathrm{Ha}=27.211386245988\ \mathrm{eV}
$$

The plotted energy is shifted relative to the Fermi energy:

$$
\boxed{
E_{\mathrm{plot}} =
(E-E_F)\times27.211386245988
}
$$

Consequently, the Fermi level is located at zero eV.

The displayed energy range is controlled by `--emin` and `--emax`.

The Lorentzian calculation is performed in Hartree, so the spectral function is expressed in inverse Hartree under this convention.

## 5. Primitive K-Path and Horizontal Axis

The input path file contains three fractional primitive reciprocal coordinates per row:

$$
\mathbf k_i=(k_{i,1},k_{i,2},k_{i,3})
$$

The fractional k-points are converted into Cartesian reciprocal coordinates. The distance between consecutive points is:

$$
\Delta s_i =
\left\|
\mathbf k_{i+1}^{\mathrm{cart}}
-
\mathbf k_i^{\mathrm{cart}}
\right\|
$$

The cumulative path coordinate is:

$$
s_0=0
$$

$$
s_i=\sum_{j=0}^{i-1}\Delta s_j
$$

This cumulative distance is used as the horizontal coordinate of the spectral-function plot.

High-symmetry-point labels can be added using `--tick-indices` and `--tick-labels`. The indices refer to rows in the input path file, starting from zero.

## 6. Computational Workflow

1. Read the primitive reciprocal-space path from the CSV file.
2. Locate the QE XML and portable wavefunction files.
3. Extract eigenvalues and Fermi energy from the XML file.
4. Read plane-wave coefficients and reciprocal vectors from the WFC files.
5. Map primitive-cell k-points to corresponding supercell k-points.
6. Identify plane waves compatible with the primitive reciprocal lattice.
7. Calculate the unfolding weight for every band.
8. Construct the Lorentzian-broadened spectral function.
9. Shift the energy axis relative to the Fermi energy.
10. Generate the spectral-function plot and save the numerical data.

## 7. Input Files

Expected directory structure:

```text
<outdir>/
└── <prefix>.save/
    ├── data-file-schema.xml
    ├── wfc1.dat
    ├── wfc2.dat
    ├── wfc3.dat
    └── ...
```

The primitive reciprocal-space path is provided separately as a CSV file.

Example:

```csv
# Primitive fractional reciprocal coordinates
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

Example command:

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

For a $2\times2\times2$ diagonal supercell, the transformation matrix is:

$$
M =
\begin{pmatrix}
2&0&0\\
0&2&0\\
0&0&2
\end{pmatrix}
$$

For the spin-down channel, use:

```bash
--spin down
```

High-symmetry labels can be specified using:

```bash
--tick-indices 0 20 40 60 80 \
--tick-labels G M R X X1
```

The indices and labels must correspond to the actual rows in the primitive path file.

### Command-Line Options

| Option | Description | Default |
|---|---|---|
| `--prefix` | QE calculation prefix | Required |
| `--outdir` | QE output directory | Required |
| `--path` | CSV primitive k-path | Required |
| `--matrix` | Nine integer elements of the supercell matrix | Required |
| `--spin` | Spin channel: up or down | up |
| `--fermi-ry` | Fermi energy override in Ry | Read from XML |
| `--width-ev` | Lorentzian FWHM in eV | 0.002 Ha converted to eV |
| `--nenergy` | Number of energy-grid points | 10000 |
| `--emin` | Lower displayed energy limit in eV | -5.0 |
| `--emax` | Upper displayed energy limit in eV | 4.5 |
| `--tick-indices` | High-symmetry-point indices | None |
| `--tick-labels` | High-symmetry-point labels | None |
| `--output` | Output PNG filename | Automatic |

## 9. Output Files

The script generates two main outputs.

### Spectral-Function Plot

The PNG image displays:

- Horizontal axis: cumulative reciprocal-space path distance.
- Vertical axis: energy relative to the Fermi level in eV.
- Color intensity: unfolded spectral function.
- Dashed horizontal line: Fermi level at zero eV.

### Numerical NPZ Archive

The compressed NumPy archive contains:

| Variable | Description |
|---|---|
| `k_primitive` | Primitive fractional k-points |
| `path_distance` | Cumulative reciprocal-space path coordinate |
| `energy_ev` | Energy grid relative to the Fermi level |
| `energy_ha` | Absolute energy grid in Hartree |
| `spectral_function` | Unfolded spectral function |
| `weights` | Unfolding weights |
| `eigenvalues_ha` | QE eigenvalues in Hartree |
| `fermi_ha` | Fermi energy in Hartree |
| `width_ha` | Lorentzian FWHM in Hartree |
| `width_ev` | Lorentzian FWHM in eV |

The archive can be loaded using:

```python
import numpy as np

data = np.load("unfolded_QE_up_data.npz")

k_distance = data["path_distance"]
energy = data["energy_ev"]
spectral_function = data["spectral_function"]
weights = data["weights"]
```

## 10. Requirements and Limitations

The script requires:

- Python
- NumPy
- Matplotlib
- Quantum ESPRESSO output files

Important implementation assumptions:

- QE must provide portable plane-wave wavefunction files named `wfc*.dat`.
- The binary reader assumes a specific portable WFC record layout.
- The implementation rejects spinor wavefunctions (`npol != 1`) and Gamma-only wavefunction files.
- The primitive k-path, supercell matrix, QE k-point sampling, and WFC ordering must describe the same calculation.
- The XML eigenvalue count must match the WFC band count.
- The current unfolding calculation uses plane-wave coefficient norms.
- The broadening is a numerical visualization parameter.
- `--emin` and `--emax` control the displayed energy range, not the underlying eigenvalues.

## 11. Interpreting the Unfolded Plot

Bright regions indicate energies and primitive k-points with relatively large unfolded spectral weight.

Strong, narrow features that follow smooth dispersions indicate states with well-defined primitive-cell Bloch character.

Weak or diffuse intensity can indicate that a supercell eigenstate has only a small projection onto the selected primitive wavevector.

In doped or defective systems, additional features may arise from impurity-derived states, band splitting, or mixing between primitive-cell wavevectors.

A weak spectral feature should not automatically be interpreted as the absence of an electronic state. It may simply have a low projection weight at the selected primitive k-point.

## 12. Comparison with GPAW

The script uses plotting defaults intended to facilitate comparison with a GPAW unfolding workflow:

- Lorentzian FWHM: $0.002$ Hartree.
- Energy grid: 10000 points.
- Fermi level aligned to zero.
- `viridis` colormap.
- `PowerNorm(gamma=0.4)` intensity normalization.

These matching plotting conventions do not guarantee numerically identical QE and GPAW spectra.

A meaningful comparison requires consistent:

- Crystal structures and supercell definitions.
- Primitive reciprocal-space paths.
- Spin treatment.
- Exchange-correlation functionals.
- Hubbard U parameters.
- Pseudopotential or PAW datasets.
- Fermi-level references.
- Energy windows and broadening conventions.
- Unfolding-weight definitions.

Differences in basis sets, wavefunction normalization, projection conventions, and implementation details can remain even when plotting settings match.

## Summary

The script converts QE supercell eigenvalues and plane-wave coefficients into a primitive-cell unfolded spectral function. It identifies plane-wave components compatible with each primitive k-point, calculates normalized unfolding weights, broadens the weighted eigenvalues using Lorentzians, aligns the energy axis to the Fermi level, and saves both the spectral map and its numerical data.
