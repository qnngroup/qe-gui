import os
import json
from tabs.tab_functions.build_k_path import BuildKPath


def write_qe_inputs(input_data):
    prefix = input_data["electronic"]["scf"]["prefix"]
    input_dir = f"./inputs/{prefix}"
    tmp_dir = f"./tmp/{prefix}"
    pseudo_dir = f"./pseudo"

    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    def control_block(kind):
        e = input_data["electronic"][kind]
        lines = [
            " &control",
            f"    calculation = '{kind}'",
            f"    prefix = '{e['prefix']}'",
            f"    outdir = '{tmp_dir}'",
            f"    pseudo_dir = '{pseudo_dir}'",
            f"    verbosity = '{e['verbosity']}'",
            f"    wf_collect = {'.true.' if e['wf_collect'] else '.false.'}",
            f"    disk_io = '{e['disk_io']}'",
            " /"
        ]
        return "\n".join(lines)

    def system_block():
        s = input_data["system"]
        lines = [
            " &system",
            f"    ibrav = {s['ibrav']}",
            f"    nat = {s['nat']}",
            f"    ntyp = {s['ntyp']}",
            f"    ecutwfc = {s['ecutwfc']}",
            f"    ecutrho = {s['ecutrho']}",
            f"    occupations = '{s['occupations']}'",
            f"    force_symmorphic = {'.true.' if s['symmorphic'] else '.false.'}",

            ]
        if int(s["ibrav"]) == 4:
            lines.append(f"    celldm(3) = {s['celldm(3)']}")

        if s["occupations"] == "smearing":
            lines.append(f"    smearing = '{s['smearing']}'")
            lines.append(f"    degauss = {s['degauss']}")
        if s["nbnd"]:
            lines.append(f"    nbnd = {s['nbnd']}")
       
        lines.append(f"    nspin = {s['nspin']}")
        lines.append(f"    input_dft = '{s['input_dft']}'")

        a_angstrom = float(s.get("celldm_angstrom", "0"))
        a_bohr = a_angstrom / 0.529177
        lines.append(f"    celldm(1) = {a_bohr:.6f}")

        if s["ibrav"] == "4" and float(s.get("celldm_ratio", "0")) > 0:
            lines.append(f"    celldm(3) = {s['celldm_ratio']}")

        lines.append(" /")
        return "\n".join(lines)

    def electrons_block():
        g = input_data["general"]
        lines = [
            " &electrons",
            f"    conv_thr = {g['conv_thr']}",
            f"    mixing_beta = {g['mixing_beta']}",
            f"    electron_maxstep = {g['electron_maxstep']}",
            " /"
        ]
        return "\n".join(lines)

    def atomic_species_block(atomic_positions_text):
        psuedo_dir = f"./psuedo"
        lines = atomic_positions_text.strip().splitlines()
        atom_names = [line.split()[0] for line in lines if line.strip()]
        unique_atoms = sorted(set(atom_names), key=atom_names.index)
        with open("./tabs/tab_functions/atomic_masses.json") as f:
            atomic_masses = json.load(f)

        try:
            pseudo_files = os.listdir(pseudo_dir)
        except FileNotFoundError:
            raise RuntimeError(f"⚠ Pseudopotential directory '{pseudo_dir}' not found!")

        block_lines = ["ATOMIC_SPECIES"]
        for atom in unique_atoms:
            matches = [f for f in pseudo_files if f.startswith(f"{atom}_") and f.endswith(".upf")]
            if not matches:
                raise RuntimeError(f"❌ No pseudopotential found for element '{atom}' in '{pseudo_dir}'")
            block_lines.append(f"{atom} {atomic_masses[atom]} {matches[0]}")

        return "\n".join(block_lines)

    def atomic_positions_block():
        pos = input_data["general"]["atomic_positions"].strip()
        return "ATOMIC_POSITIONS crystal\n" + pos if pos else ""

    def cell_parameters_block():
        cell = input_data["general"]["cell_parameters"].strip()
        return "CELL_PARAMETERS angstrom\n" + cell if cell else ""

    def k_points_block(kind):
        kdata = input_data["kpoints"]
        if kind == "scf" or kdata["mode"] == "automatic":
            m = kdata["mesh"]
            return f"K_POINTS automatic\n{m['nk1']} {m['nk2']} {m['nk3']} {m['shift1']} {m['shift2']} {m['shift3']}"
        else:
            return BuildKPath(kdata["high_symmetry_points"], kdata["path"], kdata["points_per_segment"])

    def build_file(kind):
        blocks = [
            control_block(kind),
            system_block(),
            electrons_block(),
            atomic_species_block(
                input_data["general"]["atomic_positions"],
            ),
            atomic_positions_block(),
            cell_parameters_block(),
            k_points_block(kind)
        ]
        return "\n\n".join([b for b in blocks if b.strip()])

    scf_content = build_file("scf")
    nscf_content = build_file("nscf")

    with open(f"{input_dir}/{prefix}.scf.in", "w") as f:
        f.write(scf_content)

    with open(f"{input_dir}/{prefix}.nscf.in", "w") as f:
        f.write(nscf_content)

    print(f"✔ SCF input written to: {input_dir}/{prefix}.scf.in")
    print(f"✔ NSCF input written to: {input_dir}/{prefix}.nscf.in")
    return f"{input_dir}/{prefix}.scf.in", f"{input_dir}/{prefix}.nscf.in"


def write_bands_in_file(prefix, path_data, hs_points_text, points_per_segment):
    hs_map = {}
    for line in hs_points_text.strip().splitlines():
        if line.strip():
            tokens = line.strip().split()
            label, coords = tokens[0], list(map(float, tokens[1:]))
            hs_map[label] = coords

    segments = path_data.strip().split()
    if len(segments) < 2:
        return

    input_dir = f"./inputs/{prefix}"
    output_dir = f"./outputs/{prefix}"
    tmp_dir = f"./tmp/{prefix}"

    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    filepath = f"{input_dir}/{prefix}.bands.in"

    with open(filepath, "w") as f:
        f.write("&bands\n")
        f.write(f"   prefix = '{prefix}',\n")
        f.write(f"   outdir = '{tmp_dir}',\n")
        f.write(f"   filband = '{output_dir}/{prefix}.bands.dat'\n")
        f.write("/\n\n")
        f.write(BuildKPath(hs_points_text, path_data, points_per_segment))
        f.write("\n! Path: " + " → ".join(segments) + "\n")

    print(f"✔ K-Path input written to: {filepath}")


def write_yambo_inputs(values):
    sim_type = values["simulation_type"]
    if sim_type == "RPA Dielectric Function":
        write_rpa_input(values)
    elif sim_type == "GW (G0W0)":
        write_gw_input(values)
    elif sim_type == "BSE (Excitons)":
        write_bse_input(values)
    else:
        raise ValueError(f"Unknown simulation type: {sim_type}")

def write_rpa_input(values, filename="rpa.in"):
    path = values["qe_output_path"]
    if not path:
        raise ValueError("QE output path not provided.")

    er_lo, er_hi = values["energy_range"].split()
    dr_lo, dr_hi = values["damping"].split()

    lines = [
        "chi                             # [R][CHI] Linear response calculation"
    ]

    if values.get("include_optics", False):
        lines += [
            "optics                          # [R OPT] Optical properties",
            "dipoles                         # [R DIP] Compute oscillator strengths",
            f"% QpntsRXd\n {values['q_start']} | {values['q_end']} |\n%",
            f"% EnRngeXd\n {er_lo} | {er_hi} | eV\n%",
            f"% DmRngeXd\n {dr_lo} | {dr_hi} | eV\n%",
            f"% BndsRnXd\n {values['band_min']} | {values['band_max']} |\n%",
            f"ETStpsXd = {values['nstps']}",
            f"% LongDrXd\n {values['polarization']} |\n%"
        ]

    lines += [
        f"% QpntsRX\n {values['q_start']} | {values['q_end']} |\n%",
        f"% BndsRnXp\n {values['band_min']} | {values['band_max']} |\n%",
        f"NGsBlkXp= {values['g_cutoff']}        RL      # [Xp] Response block size",
        f"% LongDrXp\n {values['polarization']} |\n%",
        f"% EnRnge\n {er_lo} | {er_hi} | eV\n%",
        f"% DmRnge\n {dr_lo} | {dr_hi} | eV\n%",
        f"Chimod= \"{values['dielectric_model']}\"     # [Xp] Dielectric model",
        f"ElecTemp= {values['electronic_temp']}     eV  # Electronic temperature",
    ]

    out_path = os.path.join(path, filename)
    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"\u2714 Wrote RPA input file to {out_path}")


def write_gw_input(values, filename="gw.in"):
    path = values["qe_output_path"]
    if not path:
        raise ValueError("QE output path not provided.")

    lines = [
        "gw0                            # [R GW] GW approximation",
        f"% BndsRnXp\n {values['gw_band_min']} | {values['gw_band_max']} |\n%",
        f"EXXRLvcs= {values['exx_cutoff']} Ry    # Exchange cutoff",
        f"PPAPntXp= {values['ppa_energy']}       eV   # Plasmon-pole energy",
        f"QPkrange= {values['qpkrange']}         # k-point and band range for QP",
        f"GTermKind= \"{values['gterm_kind']}\"     # G-term correction"
    ]

    out_path = os.path.join(path, filename)
    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"✔ Wrote GW input file to {out_path}")
def write_bse_input(values, filename="bse.in"):
    path = values["qe_output_path"]
    if not path:
        raise ValueError("QE output path not provided.")

    lines = [
        "bse                            # [R BSE] Bethe–Salpeter equation",
        f"% BSEBands\n {values['bse_band_min']} | {values['bse_band_max']} |\n%",
        f"BSKmod= \"{values['bse_kernel']}\"        # BSE kernel type",
        f"BSENGBlk= {values['bseng_block']} RL     # G-vector cutoff block",
        f"KfnQP_E= \"{values['kfnqp_input']}\"       # Input QP file from GW",
        f"%Krange\n {values['krange_input']} \n%"  # k-point range
    ]

    out_path = os.path.join(path, filename)
    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"✔ Wrote BSE input file to {out_path}")




