
#!/usr/bin/env python3

import argparse
import csv
import glob
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import PowerNorm

# ============================================================
# CONSTANTS
# ============================================================
RY_TO_EV = 13.605693122994
HA_TO_EV = 27.211386245988
RY_TO_HA = 0.5

# ============================================================
# READ PRIMITIVE K-PATH
# ============================================================
def read_path(filename):
    rows = []

    with open(filename, newline="", encoding="utf-8") as stream:
        for row in csv.reader(stream):
            if not row or not row[0].strip():
                continue

            if row[0].lstrip().startswith("#"):
                continue

            vals = [float(v.strip()) for v in row if v.strip()]

            if len(vals) != 3:
                raise ValueError(
                    f"Expected three coordinates: {row}"
                )

            rows.append(vals)

    if not rows:
        raise ValueError("No k-points found.")

    return np.asarray(rows, dtype=float)

# ============================================================
# XML HELPERS
# ============================================================
def local_name(element):
    return element.tag.split("}")[-1]

def read_fermi_ry(xml_file):

    root = ET.parse(xml_file).getroot()

    for elem in root.iter():

        if local_name(elem) in {
            "fermi_energy",
            "two_fermi_energies"
        }:

            values = np.fromstring(
                " ".join(elem.itertext()),
                sep=" "
            )

            if values.size:
                return float(values[0])

    raise ValueError(
        f"No Fermi energy found in {xml_file}"
    )

def read_xml_eigenvalues(xml_file):

    root = ET.parse(xml_file).getroot()

    entries = [
        e for e in root.iter()
        if local_name(e) == "ks_energies"
    ]

    values = []

    for entry in entries:

        eig = next(
            (
                e for e in entry.iter()
                if local_name(e) == "eigenvalues"
            ),
            None
        )

        if eig is not None and eig.text:
            values.append(
                np.fromstring(eig.text, sep=" ")
            )

    if not values:
        raise ValueError(
            f"No eigenvalues found in {xml_file}"
        )

    return values

# ============================================================
# READ QE PORTABLE WAVEFUNCTION
# ============================================================
def read_qe_wfc_dat(filename):

    with open(filename, "rb") as f:

        f.seek(4)

        ik = int(
            np.fromfile(f, dtype=np.int32, count=1)[0]
        )

        xk = np.fromfile(
            f, dtype=np.float64, count=3
        )

        ispin = int(
            np.fromfile(f, dtype=np.int32, count=1)[0]
        )

        gamma_only = bool(
            np.fromfile(f, dtype=np.int32, count=1)[0]
        )

        scale = float(
            np.fromfile(f, dtype=np.float64, count=1)[0]
        )

        f.seek(8, 1)

        ngw, igwx, npol, nbnd = map(
            int,
            np.fromfile(f, dtype=np.int32, count=4)
        )

        f.seek(8, 1)

        b = np.array([
            np.fromfile(f, dtype=np.float64, count=3)
            for _ in range(3)
        ])

        f.seek(8, 1)

        mill = np.fromfile(
            f,
            dtype=np.int32,
            count=3 * igwx
        ).reshape(igwx, 3)

        f.seek(8, 1)

        coeff = np.empty(
            (nbnd, npol * igwx),
            dtype=np.complex128
        )

        for band in range(nbnd):

            coeff[band] = np.fromfile(
                f,
                dtype=np.complex128,
                count=npol * igwx
            )

            f.seek(8, 1)

    return {
        "ik": ik,
        "xk": xk,
        "ispin": ispin,
        "gamma_only": gamma_only,
        "scale": scale,
        "ngw": ngw,
        "igwx": igwx,
        "npol": npol,
        "nbnd": nbnd,
        "b": b,
        "mill": mill,
        "coeff": coeff
    }

# ============================================================
# K-PATH DISTANCE
# ============================================================
def cumulative_distance(kpts, reciprocal_rows):

    cart = kpts @ reciprocal_rows

    delta = np.diff(cart, axis=0)

    return np.r_[
        0.0,
        np.cumsum(
            np.linalg.norm(delta, axis=1)
        )
    ]

# ============================================================
# MAIN
# ============================================================
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--prefix", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--path", required=True)

    parser.add_argument(
        "--matrix",
        required=True,
        nargs=9,
        type=int
    )

    parser.add_argument(
        "--spin",
        choices=["up", "down"],
        default="up"
    )

    parser.add_argument(
        "--fermi-ry",
        type=float,
        default=None
    )

    # GPAW default: 0.002 Hartree FWHM
    parser.add_argument(
        "--width-ev",
        type=float,
        default=0.002 * HA_TO_EV
    )

    # GPAW default: 10000 energy points
    parser.add_argument(
        "--nenergy",
        type=int,
        default=10000
    )

    parser.add_argument(
        "--emin",
        type=float,
        default=-5.0
    )

    parser.add_argument(
        "--emax",
        type=float,
        default=4.5
    )

    parser.add_argument(
        "--tick-indices",
        nargs="*",
        type=int,
        default=None
    )

    parser.add_argument(
        "--tick-labels",
        nargs="*",
        default=None
    )

    parser.add_argument(
        "--output",
        default=None
    )

    args = parser.parse_args()

    M = np.asarray(
        args.matrix,
        dtype=int
    ).reshape(3, 3)

    if round(np.linalg.det(M)) == 0:
        parser.error("Matrix must be nonsingular.")

    MinvT = np.linalg.inv(M).T

    k_prim = read_path(args.path)

    save_dir = (
        Path(args.outdir) /
        f"{args.prefix}.save"
    )

    xml_file = save_dir / "data-file-schema.xml"

    if not xml_file.exists():
        parser.error(f"QE XML not found: {xml_file}")

    files = glob.glob(
        str(save_dir / "wfc*.dat")
    )

    files = sorted(
        files,
        key=lambda p: int(Path(p).stem[3:])
        if Path(p).stem[3:].isdigit()
        else 0
    )

    if not files:
        parser.error("No QE portable wavefunctions found.")

    by_index = {
        int(Path(p).stem[3:]): p
        for p in files
        if Path(p).stem[3:].isdigit()
    }

    xml_eigs = read_xml_eigenvalues(xml_file)

    if (
        len(xml_eigs) != len(k_prim)
        or len(by_index) < len(k_prim)
    ):

        parser.error(
            f"Path: {len(k_prim)}, "
            f"XML k-points: {len(xml_eigs)}, "
            f"WFC files: {len(by_index)}. "
            "Check the QE band-path calculation."
        )

    # ========================================================
    # READ SPIN-RESOLVED EIGENVALUES AND UNFOLDING WEIGHTS
    # ========================================================
    k_super = k_prim @ M.T

    weights = []
    energies_ha = []

    b = None

    wanted_spin = (
        1 if args.spin == "up" else 2
    )

    for ik, (kp, Kraw) in enumerate(
        zip(k_prim, k_super)
    ):

        filename = by_index.get(ik + 1)

        if filename is None:
            raise ValueError(
                f"Missing WFC file for k-point {ik + 1}"
            )

        wfc = read_qe_wfc_dat(filename)

        if wfc["ik"] != ik + 1:
            raise ValueError(
                f"Incorrect k-point index in {filename}"
            )

        if wfc["npol"] != 1 or wfc["gamma_only"]:
            raise ValueError(
                "Requires non-spinor, non-Gamma-only WFC files."
            )

        if wfc["ispin"] != wanted_spin:
            raise ValueError(
                f"Wrong spin channel in {filename}"
            )

        if b is None:
            b = wfc["b"]

        elif not np.allclose(b, wfc["b"], atol=1e-8):
            raise ValueError(
                "Reciprocal lattice changed between WFC files."
            )

        wfc_k = wfc["xk"] @ np.linalg.inv(b)

        translation = Kraw - wfc_k

        T = np.rint(translation).astype(int)

        if np.linalg.norm(translation - T) > 2e-5:
            raise ValueError(
                f"Incorrect mapped k-point in {filename}"
            )

        coeff = wfc["coeff"]
        mill = wfc["mill"]

        e_ha = xml_eigs[ik] * RY_TO_HA

        if len(e_ha) != coeff.shape[0]:
            raise ValueError(
                f"Band count mismatch at k-point {ik + 1}"
            )

        # Same plane-wave projection principle as GPAW
        g_prim = (mill - T) @ MinvT

        mask = np.all(
            np.isclose(
                g_prim,
                np.rint(g_prim),
                atol=1e-6
            ),
            axis=1
        )

        norm = np.sum(
            np.abs(coeff) ** 2,
            axis=1
        )

        w = np.sum(
            np.abs(coeff[:, mask]) ** 2,
            axis=1
        ) / np.maximum(norm, 1e-30)

        energies_ha.append(e_ha)
        weights.append(w)

        print(
            f"k {ik + 1:4d}/{len(k_prim)} "
            f"selected G vectors: {mask.sum():5d} "
            f"weight range: {w.min():.5f} - {w.max():.5f}"
        )

    energies_ha = np.asarray(energies_ha)
    weights = np.asarray(weights)

    # ========================================================
    # FERMI ENERGY
    # ========================================================
    if args.fermi_ry is None:
        ef_ry = read_fermi_ry(xml_file)
    else:
        ef_ry = args.fermi_ry

    ef_ha = ef_ry * RY_TO_HA

    # ========================================================
    # GPAW-COMPATIBLE SPECTRAL FUNCTION
    # ========================================================
    width_ha = args.width_ev / HA_TO_EV

    # GPAW uses the complete eigenvalue range +/- 5 widths
    emin_ha = np.min(energies_ha) - 5 * width_ha
    emax_ha = np.max(energies_ha) + 5 * width_ha

    E_ha = np.linspace(
        emin_ha,
        emax_ha,
        args.nenergy
    )

    A = np.zeros(
        (len(k_prim), len(E_ha)),
        dtype=float
    )

    eta = width_ha / 2

    for ik in range(len(k_prim)):

        for e0, wt in zip(
            energies_ha[ik],
            weights[ik]
        ):

            lorentzian = (
                eta / np.pi
            ) / (
                (E_ha - e0) ** 2 + eta ** 2
            )

            A[ik] += wt * lorentzian

    # A is now in 1/Hartree, matching GPAW's stored convention.
    E_ev = (E_ha - ef_ha) * HA_TO_EV

    print("\n========== SPECTRAL FUNCTION ==========")
    print(f"Spin channel       : {args.spin}")
    print(f"Fermi energy       : {ef_ha * HA_TO_EV:.6f} eV")
    print(f"Lorentzian FWHM    : {args.width_ev:.6f} eV")
    print(f"Energy points      : {args.nenergy}")
    print(f"Spectral max       : {A.max():.6f} 1/Ha")
    print(f"Spectral min       : {A.min():.6f} 1/Ha")
    print("=======================================\n")

    # ========================================================
    # K-PATH AXIS
    # ========================================================
    Bp = b @ M.T

    x = cumulative_distance(
        k_prim,
        Bp
    )

    # ========================================================
    # PLOT: GPAW-MATCHED STYLE
    # ========================================================
    fig, ax = plt.subplots(
        figsize=(6, 6)
    )

    extent = [
        x[0],
        x[-1],
        E_ev[0],
        E_ev[-1]
    ]

    im = ax.imshow(
        A.T,
        aspect="auto",
        extent=extent,
        origin="lower",
        cmap="viridis",
        norm=PowerNorm(gamma=0.4),
        interpolation="bicubic"
    )

    ax.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=0.8,
        alpha=0.6
    )

    # Optional high-symmetry labels
    if args.tick_indices is not None:

        if args.tick_labels is None:
            parser.error(
                "--tick-indices requires --tick-labels"
            )

        if len(args.tick_indices) != len(args.tick_labels):
            parser.error(
                "Tick indices and labels must have equal lengths."
            )

        if any(
            i < 0 or i >= len(x)
            for i in args.tick_indices
        ):
            parser.error(
                "A tick index is outside the k-path."
            )

        tick_positions = x[args.tick_indices]

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(
            args.tick_labels,
            fontsize=16
        )

        for xtick in tick_positions:
            ax.axvline(
                xtick,
                color="black",
                linestyle="--",
                linewidth=0.7,
                alpha=0.4
            )

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(args.emin, args.emax)

    ax.tick_params(
        axis="both",
        labelsize=16
    )

    ax.set_ylabel(
        r"Energy - $E_F$ (eV)",
        fontsize=16
    )

    ax.set_title(
        f"CsPbBr$_3$ Spin-{args.spin.capitalize()}",
        fontsize=16
    )

    cbar = plt.colorbar(im, ax=ax)

    cbar.set_label(
        r"Spectral Density $A(k,E)$ (1/Ha)",
        fontsize=16
    )

    cbar.ax.tick_params(labelsize=16)

    plt.tight_layout()

    out = args.output or (
        f"unfolded_QE_{args.spin}_GPAW_comparison.png"
    )

    plt.savefig(
        out,
        dpi=300
    )

    # ========================================================
    # SAVE DATA
    # ========================================================
    np.savez(
        f"{Path(out).with_suffix('').as_posix()}_data.npz",

        k_primitive=k_prim,
        path_distance=x,

        energy_ev=E_ev,
        energy_ha=E_ha,

        spectral_function=A,
        weights=weights,

        eigenvalues_ha=energies_ha,
        fermi_ha=ef_ha,

        width_ha=width_ha,
        width_ev=args.width_ev
    )

    print(f"Saved: {out}")
    print("Saved spectral-function data.")

if __name__ == "__main__":
    main()
