# Quantum ESPRESSO Band-Unfolding and Spectral-Function Plotter

## Overview

This script reads a Quantum ESPRESSO (QE) band-structure calculation
performed for a supercell, projects its supercell Bloch states onto the
reciprocal-space basis of a chosen primitive cell, and constructs an
unfolded spectral function along a primitive-cell k-path. The resulting
intensity map displays where supercell eigenstates retain primitive-cell
Bloch character. It is therefore more informative than simply folding
the supercell eigenvalues into a smaller Brillouin zone: each eigenstate
contributes with a projection weight that indicates how strongly it
belongs to a particular primitive-cell wavevector.

The script also aligns the energy axis to the Fermi level, applies
Lorentzian broadening, plots the spectral function, and saves the
numerical data used to make the plot. It is designed for non-spinor,
non-Gamma-only QE portable wavefunction files and treats the spin-up and
spin-down channels separately.

## 1. Physical background: why band unfolding is needed

A primitive crystal has translational symmetry under its primitive
lattice vectors. A supercell calculation is often used to model defects,
vacancies, dopants, disorder approximants, or magnetic ordering. The
supercell is larger than the primitive cell, so its Brillouin zone is
correspondingly smaller. The supercell calculation consequently contains
many bands that appear folded into the reduced Brillouin zone.

A conventional supercell band plot shows all these folded bands with
equal visual status. It does not directly distinguish states that
preserve a well-defined primitive-cell wavevector from states whose
primitive translational character is weak or mixed. Band unfolding
addresses this by assigning a spectral weight to each supercell
eigenstate at each primitive-cell wavevector.

For a perfect crystal represented by a commensurate supercell, the
unfolded spectrum recovers the primitive-cell band structure, with the
appropriate spectral weights. In a doped or defective supercell, the
weights can be distributed among several states and energies. This
distribution can reveal band splitting, impurity-derived features,
disorder-induced spectral broadening, and loss of primitive
translational coherence.

## 2. Lattice and reciprocal-space mapping

Let the primitive lattice vectors be the rows of the matrix

\[ A_p =
```{=tex}
\begin{pmatrix}
\mathbf a_1\\
\mathbf a_2\\
\mathbf a_3
\end{pmatrix}
```
. \]

The supercell is defined by an integer transformation matrix (M):

\[ A_s = M A_p. \]

The script accepts the nine integer elements of (M) through `--matrix`
and reshapes them into a (3`\times3`{=tex}) matrix. The matrix must be
nonsingular. Its determinant gives the ratio of supercell to
primitive-cell volumes in the commensurate case:

\[ N = \|`\det `{=tex}M\|. \]

For row-vector fractional reciprocal coordinates, the reciprocal-space
mapping used by the script is

\[ `\mathbf `{=tex}K\_{`\mathrm{raw}`{=tex}} = `\mathbf `{=tex}k_p M\^T,
\]

where (`\mathbf `{=tex}k_p) is a primitive-cell k-point and
(`\mathbf `{=tex}K\_{`\mathrm{raw}`{=tex}}) is its corresponding
supercell-coordinate representation before reduction by a reciprocal
lattice vector. The script compares this target with the k-point stored
in each QE wavefunction file. Their difference should be an integer
reciprocal-lattice translation:

\[ `\mathbf `{=tex}T = `\operatorname{round}`{=tex} `\left`{=tex}(
`\mathbf `{=tex}K\_{`\mathrm{raw}`{=tex}}-`\mathbf `{=tex}k\_{`\mathrm{QE}`{=tex}}
`\right`{=tex}). \]

A tolerance check verifies that the difference is sufficiently close to
an integer vector. This is an important consistency check: it helps
detect a mismatch between the supplied primitive path, the supercell
matrix, and the k-points used in the QE calculation.

The reciprocal transformation used to identify primitive reciprocal
vectors is

\[ M\^{-T} = (M^{-1})^T. \]

For each supercell plane-wave reciprocal index (`\mathbf `{=tex}G_s),
read from the QE `mill` array, the script evaluates

\[ `\mathbf `{=tex}G_p = (`\mathbf `{=tex}G_s-`\mathbf `{=tex}T)M\^{-T}.
\]

A plane wave is selected when all three components of
(`\mathbf `{=tex}G_p) are sufficiently close to integers. In code, this
is implemented using `np.isclose` with an absolute tolerance of
(10\^{-6}). These selected plane waves form the subset compatible with
the chosen primitive-cell Bloch wavevector.

## 3. Unfolding spectral weight

At a fixed supercell k-point, a Kohn--Sham Bloch state can be expanded
in plane waves:

\[
`\left`{=tex}\|`\Psi`{=tex}*{n`\mathbf `{=tex}K}`\right`{=tex}`\rangle`{=tex}
= `\sum`{=tex}*{`\mathbf `{=tex}G_s}
C\_{n`\mathbf `{=tex}K}(`\mathbf `{=tex}G_s)
`\left`{=tex}\|`\mathbf `{=tex}K+`\mathbf `{=tex}G_s`\right`{=tex}`\rangle`{=tex}.
\]

Here, (n) is the band index, (`\mathbf `{=tex}K) is the supercell
wavevector, and (C\_{n`\mathbf `{=tex}K}(`\mathbf `{=tex}G_s)) is the
complex plane-wave coefficient read from the QE portable wavefunction
file.

The primitive-cell projection weight used in this script is

\[ W\_{n}(`\mathbf `{=tex}k_p) = `\frac{
\displaystyle\sum_{\mathbf G_s\in\mathcal S(\mathbf k_p)}
\left|C_{n\mathbf K}(\mathbf G_s)\right|^2
}{
\displaystyle\sum_{\mathbf G_s}
\left|C_{n\mathbf K}(\mathbf G_s)\right|^2
}`{=tex}, \]

where (`\mathcal `{=tex}S(`\mathbf `{=tex}k_p)) is the set of plane
waves selected by the integer reciprocal-space condition above. The
denominator normalizes by the total plane-wave coefficient norm for that
band. A small numerical floor is used to avoid division by zero.

The weight is dimensionless and, for a correctly normalized wavefunction
and a valid selected subset, is expected to lie between zero and one,
subject to numerical and file-format conventions. A weight near one
means that most of the state's plane-wave content belongs to the
selected primitive-cell translational sector. A small weight means that
little of the state's plane-wave content projects onto that primitive
wavevector.

This implementation calculates the weight directly from the squared
magnitudes of the complex plane-wave coefficients:

\[ \|C\|\^2 = C\^\*C. \]

It does not use the phase of the coefficients in the final weight,
because the projection probability is formed from the squared
coefficient magnitudes.

## 4. Spectral function and Lorentzian broadening

The ideal unfolded spectrum is a set of weighted delta functions:

\[ A(`\mathbf `{=tex}k,E) = `\sum`{=tex}*n W_n(`\mathbf `{=tex}k),
`\delta`{=tex}(E-`\epsilon`{=tex}*{n`\mathbf `{=tex}K}), \]

where (`\epsilon`{=tex}\_{n`\mathbf `{=tex}K}) is the eigenvalue of band
(n) at the corresponding supercell k-point.

A delta function is not directly suitable for a continuous intensity
plot. The script replaces it with a Lorentzian:

\[ L(E-`\epsilon`{=tex}\_n;`\eta`{=tex}) = `\frac{1}{\pi}`{=tex}
`\frac{\eta}`{=tex} {(E-`\epsilon`{=tex}\_n)^2+`\eta`{=tex}^2}. \]

The broadened spectral function is therefore

\[ A(`\mathbf `{=tex}k,E) = `\sum`{=tex}*n W_n(`\mathbf `{=tex}k)
`\frac{1}{\pi}`{=tex} `\frac{\eta}`{=tex}
{(E-`\epsilon`{=tex}*{n`\mathbf `{=tex}K})^2+`\eta`{=tex}^2}. \]

The Lorentzian half-width at half-maximum is (`\eta`{=tex}), and its
full width at half-maximum is

\[ `\mathrm{FWHM}`{=tex}=2`\eta`{=tex}. \]

The command-line option `--width-ev` specifies the FWHM in eV. The code
converts it to Hartree and sets

\[ `\eta `{=tex}= `\frac{\mathrm{FWHM}}{2}`{=tex}. \]

The default is (0.002) Hartree, equivalent to approximately (0.0544) eV
FWHM. A smaller width produces sharper spectral features; a larger width
produces smoother, broader features. Broadening is a plotting/numerical
representation of the discrete spectrum and should not be interpreted
automatically as a physical scattering lifetime.

The energy grid is generated from the minimum and maximum eigenvalues
found in the calculation, extended by five Lorentzian widths at each
end. The number of grid points is controlled by `--nenergy` (default:
10000). The code evaluates the Lorentzian contribution of every band at
every energy-grid point and sums the contributions weighted by the
unfolding weights.

### Energy units and Fermi-level alignment

QE eigenvalues and Fermi energies are read from `data-file-schema.xml`
in Rydberg. The script converts Rydberg to Hartree using

\[ 1 `\mathrm{Ha}`{=tex}=2 `\mathrm{Ry}`{=tex}, \]

and Hartree to eV using

\[ 1 `\mathrm{Ha}`{=tex}=27.211386245988 `\mathrm{eV}`{=tex}. \]

The plotted energy is shifted relative to the Fermi energy:

\[ E\_{`\mathrm{plot}`{=tex}} = (E-E_F)`\times 27.211386245988`{=tex}.
\]

Consequently, the Fermi level is at zero eV on the vertical axis. The
plotted range is controlled by `--emin` and `--emax`; these limits crop
the displayed energy window and do not change the eigenvalues or the
spectral-function calculation itself.

The Lorentzian is evaluated in Hartree, so the stored spectral function
has units of inverse Hartree under this convention. The script labels
the colorbar accordingly.

## 5. Primitive k-path and horizontal axis

The input path file is a CSV containing three fractional primitive
reciprocal coordinates per non-comment row:

\[ (k_1,k_2,k_3). \]

The script ignores blank rows and lines whose first field begins with
`#`. Each non-comment row must contain exactly three numeric values.

To display the path as a continuous horizontal coordinate, fractional
k-points are converted to Cartesian reciprocal coordinates using the
reciprocal-lattice row matrix supplied by the QE wavefunction header and
the supercell transformation used in the script. The distance between
consecutive points is calculated as

\[ `\Delta `{=tex}s_i = `\left`{=tex}\|
`\mathbf `{=tex}k\_{i+1}\^{`\mathrm{cart}`{=tex}} -
`\mathbf `{=tex}k_i\^{`\mathrm{cart}`{=tex}} `\right`{=tex}\|, \]

and the cumulative path coordinate is

\[ s_0=0,`\qquad`{=tex}
s_i=`\sum`{=tex}\_{j=0}\^{i-1}`\Delta `{=tex}s_j. \]

This cumulative distance is the horizontal coordinate in the
spectral-function plot. It is a reciprocal-space path length, not a
real-space distance.

High-symmetry-point labels can be added using `--tick-indices` and
`--tick-labels`. The indices refer to rows in the input path file,
starting at zero. The number of indices and labels must match.

## 6. What the code does, step by step

1.  **Read the primitive path.** The script reads the three fractional
    reciprocal coordinates from the CSV path file and validates the
    input.
2.  **Read QE output files.** It locates `data-file-schema.xml` and the
    portable wavefunction files (`wfc*.dat`) inside
    `<outdir>/<prefix>.save/`.
3.  **Read eigenvalues and Fermi energy.** The XML parser extracts the
    Kohn--Sham eigenvalues at each k-point and the Fermi energy. XML
    namespaces are handled by comparing local tag names.
4.  **Read plane-wave coefficients.** The binary wavefunction reader
    extracts the k-point index, k-point coordinates, spin index,
    reciprocal lattice, `mill` vectors, band count, and complex
    plane-wave coefficients.
5.  **Validate the calculation.** The script checks that the number of
    XML k-points, path points, and wavefunction files is consistent;
    checks the spin channel; checks that the reciprocal lattice is
    consistent; and verifies the k-point mapping.
6.  **Select primitive-compatible plane waves.** The integer
    transformation matrix and reciprocal translation are used to
    identify the plane-wave subset associated with each primitive
    k-point.
7.  **Calculate unfolding weights.** For each band, the squared
    coefficient magnitudes in the selected subset are divided by the
    total squared coefficient norm.
8.  **Construct the spectral function.** Each eigenvalue is broadened
    with a Lorentzian and multiplied by its unfolding weight.
    Contributions from all bands are summed.
9.  **Align the energy axis.** The Fermi energy is subtracted, and the
    energy is converted to eV.
10. **Plot and save.** The script plots the spectral function as an
    intensity map, optionally adds high-symmetry labels and vertical
    guides, saves a PNG, and writes an NPZ archive containing the path,
    energy grid, spectral function, weights, eigenvalues, Fermi energy,
    and broadening information.

## 7. Input files and directory structure

The expected QE output structure is:

``` text
<outdir>/
└── <prefix>.save/
    ├── data-file-schema.xml
    ├── wfc1.dat
    ├── wfc2.dat
    ├── wfc3.dat
    └── ...
```

The primitive path is supplied separately as a CSV file, for example:

``` text
primitive_path.csv
```

Example path file:

``` csv
# primitive fractional reciprocal coordinates
0.000000,0.000000,0.000000
0.125000,0.000000,0.125000
0.250000,0.000000,0.250000
0.375000,0.000000,0.375000
```

The path must correspond to the primitive cell associated with the
supercell matrix. The QE calculation must have been performed at the
corresponding mapped supercell k-points and must contain the
wavefunction files needed by the script.

## 8. Running the script

Save the Python code as, for example:

``` bash
unfold_qe.py
```

A representative command is:

``` bash
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

For a (2`\times2`{=tex}`\times2`{=tex}) diagonal supercell, the matrix
is

\[ M=
```{=tex}
\begin{pmatrix}
2&0&0\\
0&2&0\\
0&0&2
\end{pmatrix}
```
. \]

For the spin-down channel, use:

``` bash
--spin down
```

To provide high-symmetry labels, add indices and labels that correspond
to the appropriate rows in the path file. For example:

``` bash
--tick-indices 0 20 40 60 80 \
--tick-labels G M R X X1
```

The values above are illustrative; they must be replaced with the actual
indices and labels for the path being used.

### Main command-line options

  ------------------------------------------------------------------------------------------
  Option                  Meaning                 Default
  ----------------------- ----------------------- ------------------------------------------
  `--prefix`              QE calculation prefix;  Required
                          used to locate          
                          `<prefix>.save`         

  `--outdir`              QE output directory     Required
                          containing              
                          `<prefix>.save`         

  `--path`                CSV file of primitive   Required
                          fractional k-points     

  `--matrix`              Nine integer elements   Required
                          of the                  
                          (3`\times3`{=tex})      
                          supercell matrix        

  `--spin`                Spin channel: `up` or   `up`
                          `down`                  

  `--fermi-ry`            Optional Fermi energy   Read from XML
                          override in Ry          

  `--width-ev`            Lorentzian FWHM in eV   (0.002) Ha converted to eV

  `--nenergy`             Number of energy-grid   10000
                          points                  

  `--emin`                Lower displayed energy  -5.0
                          limit relative to       
                          (E_F), in eV            

  `--emax`                Upper displayed energy  4.5
                          limit relative to       
                          (E_F), in eV            

  `--tick-indices`        Zero-based indices of   None
                          high-symmetry points in 
                          the path                

  `--tick-labels`         Labels corresponding to None
                          `--tick-indices`        

  `--output`              Output PNG filename     `unfolded_QE_<spin>_GPAW_comparison.png`
  ------------------------------------------------------------------------------------------

## 9. Output files

The script produces two files:

### Unfolded spectral-function plot

A PNG image containing the spectral function (A(`\mathbf `{=tex}k,E)).
The horizontal axis is cumulative reciprocal-space path distance, and
the vertical axis is energy relative to the Fermi level in eV. The
colormap represents spectral intensity. A dashed horizontal line marks
the Fermi level at zero eV.

### Numerical NPZ archive

The script also saves a compressed NumPy archive alongside the plot. It
contains:

-   `k_primitive`: primitive fractional k-points.
-   `path_distance`: cumulative reciprocal-space path coordinate.
-   `energy_ev`: energy grid relative to the Fermi level, in eV.
-   `energy_ha`: absolute energy grid used in the Lorentzian
    calculation, in Hartree.
-   `spectral_function`: unfolded spectral function array.
-   `weights`: unfolding weights for each k-point and band.
-   `eigenvalues_ha`: QE eigenvalues in Hartree.
-   `fermi_ha`: Fermi energy in Hartree.
-   `width_ha`: Lorentzian FWHM in Hartree.
-   `width_ev`: Lorentzian FWHM in eV.

The numerical archive can be loaded later without rerunning QE:

``` python
import numpy as np

data = np.load("unfolded_QE_up_data.npz")

k_distance = data["path_distance"]
energy = data["energy_ev"]
spectral_function = data["spectral_function"]
weights = data["weights"]
```

## 10. Requirements and implementation assumptions

The script requires Python, NumPy, Matplotlib, and the QE output files.
It uses only standard-library modules for argument parsing, CSV reading,
file matching, XML parsing, and path handling.

Important assumptions and limitations:

-   The QE calculation must provide portable plane-wave wavefunction
    files named `wfc*.dat` in the expected `<prefix>.save` directory.
-   The binary reader is written for the specific portable WFC record
    layout assumed by the script. QE versions, compilation options,
    endianness, or file-format changes may require changes to the
    reader.
-   The implementation explicitly rejects spinor wavefunctions
    (`npol != 1`) and Gamma-only wavefunction files. It handles the `up`
    and `down` spin channels as separate collinear-spin calculations.
-   The path points, supercell matrix, QE k-point sampling, and
    wavefunction-file ordering must describe the same calculation. The
    script performs several consistency checks, but the user must still
    verify that the path and matrix are physically correct.
-   The script expects the XML eigenvalue list and WFC band count to
    match at every k-point.
-   The current weight calculation uses the plane-wave coefficient norms
    in the QE files. It does not include a separate real-space
    projection or an explicit orbital-character decomposition.
-   The energy broadening is Lorentzian and is a visualization
    parameter. It should be held consistent when comparing different
    samples or spin channels.
-   The `--emin` and `--emax` parameters control the plotted window; the
    energy grid itself is generated from the calculated eigenvalue range
    and the broadening width.
-   The script's output should be compared with another code only when
    the structural model, lattice transformation, k-path, spin channel,
    exchange-correlation functional, Hubbard parameters,
    pseudopotentials, Fermi-level reference, band range, and broadening
    convention are made consistent.

## 11. Interpreting the unfolded plot

Bright regions indicate energies and primitive k-points with relatively
large unfolded spectral weight. A strong, narrow feature that follows a
smooth dispersion is characteristic of a state with well-defined
primitive-cell Bloch character. Weak or diffuse intensity can indicate
that a supercell eigenstate has only a small projection onto the
selected primitive wavevector. In doped or defective systems, additional
features may arise from impurity-derived states, band splitting, or
mixing between primitive-cell wavevectors.

The plot is a spectral-weight representation, not an ordinary band plot
in which every eigenvalue is drawn as an equally weighted line. The
intensity depends on the calculated wavefunctions and on the unfolding
projection. Therefore, a missing or weak feature should not
automatically be interpreted as the absence of an electronic state: it
may have low projection weight for the chosen primitive k-point.

## 12. Comparison with GPAW

The script includes defaults intended to make the energy broadening and
spectral-function display comparable to a GPAW unfolding workflow: the
default Lorentzian FWHM is (0.002) Hartree, the default energy grid
contains 10000 points, the Fermi level is placed at zero, and the image
uses a `viridis` colormap with `PowerNorm(gamma=0.4)`.

These matching plotting conventions do not by themselves guarantee
numerically identical QE and GPAW spectra. A meaningful comparison
additionally requires consistent crystal structures, supercell
definitions, k-point paths, spin treatment, exchange-correlation and
Hubbard settings, pseudopotential/PAW datasets, eigenvalue windows, and
unfolding-weight definitions. Differences in basis sets, wavefunction
normalization, projection conventions, and implementation details can
remain even when the plotting settings match.

------------------------------------------------------------------------

**In short:** the script converts QE supercell eigenvalues and
plane-wave coefficients into a primitive-cell unfolded spectral
function. It does so by identifying the plane-wave components compatible
with each primitive k-point, summing their normalized squared amplitudes
into an unfolding weight, broadening the weighted eigenvalues with
Lorentzians, aligning the result to the Fermi level, and saving both a
spectral map and its underlying numerical data.
