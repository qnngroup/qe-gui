# build_k_path.py

def BuildKPath(kpoints_text, path_string, points_per_segment):
    """
    Constructs the string to write into the K_POINTS crystal_b section
    using explicit point-by-point interpolation along high-symmetry segments.
    """
    # Parse high symmetry points
    label_to_coords = {}
    for line in kpoints_text.strip().splitlines():
        tokens = line.split()
        if len(tokens) != 4:
            continue
        label, x, y, z = tokens
        label_to_coords[label] = (float(x), float(y), float(z))

    # Parse path and build segments
    labels = path_string.strip().split()
    all_kpoints = []
    for i in range(len(labels) - 1):
        start = label_to_coords[labels[i]]
        end = label_to_coords[labels[i + 1]]
        for j in range(int(points_per_segment)):
            frac = j / (int(points_per_segment) - 1)
            interp = [
                start[k] + frac * (end[k] - start[k]) for k in range(3)
            ]
            all_kpoints.append(interp)

    # Format output
    k_lines = [f"{len(all_kpoints)}"]
    for kpt in all_kpoints:
        k_lines.append(f"{kpt[0]:.6f} {kpt[1]:.6f} {kpt[2]:.6f} 1")

    return "K_POINTS crystal\n" + "\n".join(k_lines)
