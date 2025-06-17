import os
import numpy as np
import matplotlib.pyplot as plt

def PlotBands(bands_file, dat_file, fermi_energy=0.0, output_file="plots/band_structure.png"):
    if not os.path.exists(dat_file):
        raise FileNotFoundError(f"File not found: {dat_file}")

    # Extract path segment labels from the bands input file
    labels = []
    with open(bands_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("! Path:"):
                labels = line.strip().split("Path:")[1].strip().split(" → ")

    # Read the .dat file, skipping header lines
    kpoints = []
    bands = []
    with open(dat_file, "r") as f:
        lines = f.readlines()

    current_k = None
    current_e = []
    for line in lines:
        try:
            parts = list(map(float, line.split()))
            if len(parts) == 3:
                if current_k is not None and current_e:
                    bands.append(current_e)
                current_k = parts
                kpoints.append(len(bands))
                current_e = []
            else:
                current_e.extend(parts)
        except ValueError:
            continue

    if current_e:
        bands.append(current_e)

    bands = np.array(bands)
    x = np.arange(len(bands))

    # Plot
    plt.figure(figsize=(8, 6))
    for i in range(bands.shape[1]):
        plt.plot(x, bands[:, i] - fermi_energy, color='black', lw=1)

    # Add vertical lines at high-symmetry points
    if labels and len(labels) > 1:
        segment_length = len(x) // (len(labels) - 1)
        xticks = [i * segment_length for i in range(len(labels))]
        for xtick in xticks:
            plt.axvline(x=xtick, color='gray', linestyle='--', linewidth=0.5)
        plt.xticks(xticks, labels, fontsize=12)

    plt.ylabel("Energy (eV)", fontsize=13)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.title("Electronic Band Structure", fontsize=15)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()
