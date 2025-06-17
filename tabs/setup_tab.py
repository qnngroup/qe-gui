# setup_tab.py

from PySide6.QtWidgets import QWidget, QTabWidget, QFormLayout, QLineEdit, QComboBox, QCheckBox, QVBoxLayout, QGroupBox, QTextEdit, QRadioButton, QHBoxLayout, QDialog, QLabel
from PySide6.QtCore import Qt

from tabs.tab_functions.qe_writer import write_bands_in_file

class ElectronicTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        pre_box = QGroupBox("General")
        pre_form = QFormLayout()

        self.prefix = QLineEdit("InP")
        self.prefix.setToolTip("Prefix for all output files (used by QE internally).")

        self.wf_collect = QCheckBox("Store full wavefunctions")
        self.wf_collect.setToolTip("Enable this if you want to post-process wavefunctions (e.g., for Yambo).")
        self.wf_collect.setChecked(True)

        self.verbosity = QComboBox()
        self.verbosity.addItems(["low", "medium", "high"])
        self.verbosity.setToolTip("Amount of information QE prints during SCF.")

        self.disk_io = QComboBox()
        self.disk_io.addItems(["default", "low", "high", "none", "hdf5"])
        self.disk_io.setToolTip("Controls how QE stores data on disk.")

        pre_form.addRow("File Prefix:", self.prefix)
        pre_form.addRow(self.wf_collect)
        pre_form.addRow(self.verbosity)
        pre_form.addRow(self.disk_io)


        pre_box.setLayout(pre_form)
        # --- SCF Section ---
        scf_box = QGroupBox("Self-Consistent Field (SCF) Settings")
        scf_form = QFormLayout()


        self.tprnfor = QCheckBox("Print forces")
        self.tprnfor.setToolTip("Check this to output calculated atomic forces.")

        self.tstress = QCheckBox("Print stress tensor")
        self.tstress.setToolTip("Check this to output the full stress tensor.")

        scf_form.addRow(self.tprnfor)
        scf_form.addRow(self.tstress)

        scf_box.setLayout(scf_form)

        # --- NSCF Section ---
        nscf_box = QGroupBox("Non-SCF / Band Structure Settings")
        nscf_form = QFormLayout()

        self.nscf_calculation = QComboBox()
        self.nscf_calculation.addItems(["nscf", "bands"])
        self.nscf_calculation.setToolTip("Choose 'nscf' for density of states or 'bands' for band structure.")

        nscf_form.addRow("Calculation Mode:", self.nscf_calculation)
    
        nscf_box.setLayout(nscf_form)

        # Add to main layout
        layout.addWidget(pre_box)
        layout.addWidget(scf_box)
        layout.addWidget(nscf_box)
        layout.addStretch()
        self.setLayout(layout)

    def get_values(self):
        return {
            "scf": {
                "prefix": self.prefix.text(),
                "verbosity": self.verbosity.currentText(),
                "wf_collect": self.wf_collect.isChecked(),
                "disk_io": self.disk_io.currentText(),
                "tstress": self.tstress.isChecked(),
                "tprnfor": self.tprnfor.isChecked(),
            },
            "nscf": {
                "prefix": self.prefix.text(),
                "calculation": self.nscf_calculation.currentText(),
                "disk_io": self.disk_io.currentText(),
                "verbosity": self.verbosity.currentText(),
                "wf_collect": self.wf_collect.isChecked()
            }
        }


class SystemTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()

        # Bravais lattice dropdown
        self.bravais_lattice = QComboBox()
        self.bravais_lattice.addItems([
            "Free (ibrav=0)",
            "Simple Cubic (ibrav=1)",
            "FCC (ibrav=2)",
            "BCC (ibrav=3)",
            "Simple Tetragonal (ibrav=4)",
            "Body-Centered Tetragonal (ibrav=6)",
            "Hexagonal (ibrav=4 + celldm3)"
        ])
        self.bravais_lattice.setToolTip("Choose the crystal Bravais lattice. This determines the unit cell shape and symmetry.")
        layout.addRow("Bravais Lattice:", self.bravais_lattice)

        # celldm(1) in angstrom with conversion
        self.celldm_angstrom = QLineEdit("5.8688")
        self.celldm_angstrom.setToolTip("Lattice constant a in Ångström. Will be converted to Bohr internally.")
        self.celldm_bohr_label = QLabel("→ 11.090 Bohr")
        self.celldm_bohr_label.setAlignment(Qt.AlignRight)

        def update_bohr_label():
            try:
                a_angstrom = float(self.celldm_angstrom.text())
                a_bohr = a_angstrom / 0.529177
                self.celldm_bohr_label.setText(f"→ {a_bohr:.3f} Bohr")
            except ValueError:
                self.celldm_bohr_label.setText("→ ? Bohr")

        self.celldm_angstrom.textChanged.connect(update_bohr_label)
        update_bohr_label()
        layout.addRow("Lattice Constant a (Å):", self.celldm_angstrom)
        layout.addRow("", self.celldm_bohr_label)

        self.cell_3 = QLineEdit("8")
        self.cell_3.setToolTip("Only used if ibrav=4.")
        layout.addRow("Third Cell Dimension", self.cell_3)

       
        # Atom and type counts
        self.nat = QLineEdit("2")
        self.nat.setToolTip("Total number of atoms in the unit cell. Must match the number of lines in ATOMIC_POSITIONS.")
        self.ntyp = QLineEdit("2")
        self.ntyp.setToolTip("Number of unique atomic species (e.g. In and P = 2). Must match ATOMIC_SPECIES.")
        layout.addRow("# Basis Atoms:", self.nat)
        layout.addRow("# Distinct Atoms:", self.ntyp)

        # Cutoffs
        self.ecutwfc = QLineEdit("50")
        self.ecutwfc.setToolTip("Plane-wave kinetic energy cutoff (Ry). Start with 40–60 Ry.")
        self.ecutrho = QLineEdit("400")
        self.ecutrho.setToolTip("Charge density cutoff (Ry). Typically 4–8× ecutwfc.")
        layout.addRow("ecutwfc (Ry):", self.ecutwfc)
        layout.addRow("ecutrho (Ry):", self.ecutrho)

        # Occupations and smearing
        self.occupations = QComboBox()
        self.occupations.addItems(["fixed", "smearing", "tetrahedra"])
        self.occupations.setToolTip("Choose electron occupation scheme. Use 'fixed' for insulators, 'smearing' for metals.")

        self.smearing = QComboBox()
        self.smearing.addItems(["gaussian", "mp", "mv", "fermi-dirac"])
        self.smearing.setToolTip("Smearing method (only used if occupations = 'smearing').")

        self.degauss = QLineEdit("0.01")
        self.degauss.setToolTip("Smearing width in Ry. Try 0.01–0.05 for metals.")
        layout.addRow("occupations:", self.occupations)
        layout.addRow("smearing:", self.smearing)
        layout.addRow("degauss (Ry):", self.degauss)

        # Functional and spin
        self.input_dft = QComboBox()
        self.input_dft.addItems(["PBE", "LDA", "BLYP", "HSE", "PZ", "revPBE"])
        self.input_dft.setToolTip("Exchange-correlation functional used in DFT.")

        self.nspin = QComboBox()
        self.nspin.addItems(["1 (non-magnetic)", "2 (spin-polarized)"])
        self.nspin.setToolTip("Spin configuration. Use 2 for spin-polarized (e.g. magnetic) materials.")

        self.nbnd = QLineEdit("")
        self.nbnd.setToolTip("Needed for band structure plots. Set 40-60 for Yambo calculations")
        self.nbnd = QLineEdit("20")
        layout.addRow("input_dft:", self.input_dft)
        layout.addRow("nspin:", self.nspin)
        layout.addRow("Number of Bands:", self.nbnd)

        self.symmorphic = QCheckBox("Symmorphic")
        self.symmorphic.setToolTip("forces a symmophic output (e.g., for Yambo).")
        self.symmorphic.setChecked(True)
        layout.addRow("Symmorphic",self.symmorphic)
        


        self.setLayout(layout)

    def get_values(self):
        # Compute the correct ibrav number
        ibrav_map = {
            0: 0,
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 6,  # BCT
            6: 4    # Hexagonal (special case, needs celldm(3) too)
        }

        return {
            "ibrav": str(ibrav_map.get(self.bravais_lattice.currentIndex(), 0)),
            "nat": self.nat.text(),
            "ntyp": self.ntyp.text(),
            "ecutwfc": self.ecutwfc.text(),
            "ecutrho": self.ecutrho.text(),
            "occupations": self.occupations.currentText(),
            "smearing": self.smearing.currentText(),
            "degauss": self.degauss.text(),
            "input_dft": self.input_dft.currentText(),
            "nspin": "1" if "non" in self.nspin.currentText() else "2",
            "nbnd": self.nbnd.text(),
            "celldm_angstrom": self.celldm_angstrom.text(),
            "symmorphic":self.symmorphic.isChecked(),
            "celldm(3)" : self.cell_3.text()
        }

class GeneralTab(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()

        

        # ATOMIC_POSITIONS block — multiline
        self.atomic_positions = QTextEdit()
        self.atomic_positions.setPlaceholderText("Example:\nIn 0.00 0.00 0.00\nP  0.25 0.25 0.25")

        # Optional: CELL_PARAMETERS (if needed)
        self.cell_parameters = QTextEdit()
        self.cell_parameters.setPlaceholderText("Optional — only used if ibrav = 0")

        # Electronic convergence settings
        self.conv_thr = QLineEdit("1.0e-8")
        self.mixing_beta = QLineEdit("0.7")
        self.electron_maxstep = QLineEdit("100")

        layout.addRow("ATOMIC_POSITIONS (crystal):", self.atomic_positions)
        layout.addRow("CELL_PARAMETERS (optional):", self.cell_parameters)
        layout.addRow("Convergence Threshold (conv_thr):", self.conv_thr)
        layout.addRow("Mixing Beta:", self.mixing_beta)
        layout.addRow("Electron Max Steps:", self.electron_maxstep)

        self.setLayout(layout)

    def get_values(self):
            return {
                "atomic_positions": self.atomic_positions.toPlainText(),
                "cell_parameters": self.cell_parameters.toPlainText(),
                "conv_thr": self.conv_thr.text(),
                "mixing_beta": self.mixing_beta.text(),
                "electron_maxstep": self.electron_maxstep.text()
            }   

class KPointsTab(QWidget):
    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()

        # --- SCF Mesh Section ---
        scf_box = QGroupBox("SCF K-Point Mesh")
        scf_layout = QFormLayout()
        self.nk1 = QLineEdit("6")
        self.nk2 = QLineEdit("6")
        self.nk3 = QLineEdit("6")
        self.shift1 = QLineEdit("0")
        self.shift2 = QLineEdit("0")
        self.shift3 = QLineEdit("0")
        scf_layout.addRow("Nk1 × Nk2 × Nk3:", self._row(self.nk1, self.nk2, self.nk3))
        scf_layout.addRow("Shift1 × Shift2 × Shift3:", self._row(self.shift1, self.shift2, self.shift3))
        scf_box.setLayout(scf_layout)

        # --- Mode Selection ---
        mode_box = QGroupBox("K-Point Mode")
        self.auto_mode = QRadioButton("Automatic (mesh)")
        self.path_mode = QRadioButton("Path Mode (crystal)")
        self.path_mode.setChecked(True)
        mode_layout = QVBoxLayout()
        mode_layout.addWidget(self.auto_mode)
        mode_layout.addWidget(self.path_mode)
        mode_box.setLayout(mode_layout)

        # --- High Symmetry Point Definitions ---
        hs_box = QGroupBox("High Symmetry Points")
        self.hs_points = QTextEdit()
        self.hs_points.setPlaceholderText("Example:\nΓ 0.0 0.0 0.0\nX 0.5 0.0 0.0\nW 0.5 0.25 0.75")
        hs_layout = QVBoxLayout()
        hs_layout.addWidget(self.hs_points)
        hs_box.setLayout(hs_layout)

        # --- Path Builder ---
        path_box = QGroupBox("K-Path")
        self.path_line = QLineEdit("Γ X W L Γ")
        self.points_per_segment = QLineEdit("20")
        path_layout = QFormLayout()
        path_layout.addRow("Path (space-separated):", self.path_line)
        path_layout.addRow("Points per segment:", self.points_per_segment)
        path_box.setLayout(path_layout)

        main_layout.addWidget(scf_box)
        main_layout.addWidget(mode_box)
        main_layout.addWidget(hs_box)
        main_layout.addWidget(path_box)
        self.setLayout(main_layout)

        self.load_high_symmetry_points("FCC")

    def _row(self, *widgets):
        hbox = QHBoxLayout()
        for w in widgets:
            hbox.addWidget(w)
        return hbox

    def get_values(self, prefix=None):
        mode = "automatic" if self.auto_mode.isChecked() else "path"
        hs = self.hs_points.toPlainText()
        path = self.path_line.text()
        points = self.points_per_segment.text()

        if mode == "path" and prefix:
            write_bands_in_file(prefix, path, hs, points)

        return {
            "mesh": {
                "nk1": self.nk1.text(),
                "nk2": self.nk2.text(),
                "nk3": self.nk3.text(),
                "shift1": self.shift1.text(),
                "shift2": self.shift2.text(),
                "shift3": self.shift3.text()
            },
            "mode": mode,
            "high_symmetry_points": hs,
            "path": path,
            "points_per_segment": points
        }

    def load_high_symmetry_points(self, lattice_type):
        if lattice_type.upper() == "FCC":
            points = {
                "Γ": (0.0, 0.0, 0.0),
                "X": (0.5, 0.0, 0.5),
                "W": (0.5, 0.25, 0.75),
                "K": (0.375, 0.375, 0.75),
                "L": (0.5, 0.5, 0.5),
                "U": (0.625, 0.25, 0.625)
            }
            lines = [f"{label} {x:.3f} {y:.3f} {z:.3f}" for label, (x, y, z) in points.items()]
            self.hs_points.setPlainText("\n".join(lines))
            self.path_line.setText("Γ X W K Γ L U W L")
        else:
            self.hs_points.setPlainText("")
    
    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()

        # --- SCF Mesh Section ---
        scf_box = QGroupBox("SCF K-Point Mesh")
        scf_layout = QFormLayout()
        self.nk1 = QLineEdit("6")
        self.nk2 = QLineEdit("6")
        self.nk3 = QLineEdit("6")
        self.shift1 = QLineEdit("0")
        self.shift2 = QLineEdit("0")
        self.shift3 = QLineEdit("0")
        scf_layout.addRow("Nk1 × Nk2 × Nk3:", self._row(self.nk1, self.nk2, self.nk3))
        scf_layout.addRow("Shift1 × Shift2 × Shift3:", self._row(self.shift1, self.shift2, self.shift3))
        scf_box.setLayout(scf_layout)

        # --- Mode Selection ---
        mode_box = QGroupBox("K-Point Mode")
        self.auto_mode = QRadioButton("Automatic (mesh)")
        self.path_mode = QRadioButton("Path Mode (crystal_b)")
        self.path_mode.setChecked(True)
        mode_layout = QVBoxLayout()
        mode_layout.addWidget(self.auto_mode)
        mode_layout.addWidget(self.path_mode)
        mode_box.setLayout(mode_layout)

        # --- High Symmetry Point Definitions ---
        hs_box = QGroupBox("High Symmetry Points")
        self.hs_points = QTextEdit()
        self.hs_points.setPlaceholderText("Example:\nΓ 0.0 0.0 0.0\nX 0.5 0.0 0.0\nW 0.5 0.25 0.75")
        hs_layout = QVBoxLayout()
        hs_layout.addWidget(self.hs_points)
        hs_box.setLayout(hs_layout)

        # --- Path Builder ---
        path_box = QGroupBox("K-Path")
        self.path_line = QLineEdit("Γ X W L Γ")
        self.points_per_segment = QLineEdit("20")
        path_layout = QFormLayout()
        path_layout.addRow("Path (space-separated):", self.path_line)
        path_layout.addRow("Points per segment:", self.points_per_segment)
        path_box.setLayout(path_layout)

        # Assemble layout
        main_layout.addWidget(scf_box)
        main_layout.addWidget(mode_box)
        main_layout.addWidget(hs_box)
        main_layout.addWidget(path_box)
        self.setLayout(main_layout)

        # Load default high-symmetry points for FCC lattice
        self.load_high_symmetry_points("FCC")

    def _row(self, *widgets):
        hbox = QHBoxLayout()
        for w in widgets:
            hbox.addWidget(w)
        return hbox

    def get_values(self, prefix=None):
        mode = "automatic" if self.auto_mode.isChecked() else "path"
        hs = self.hs_points.toPlainText()
        path = self.path_line.text()
        points = self.points_per_segment.text()

        if mode == "path" and prefix:
            write_bands_in_file(prefix, path, hs, points)

        return {
            "mesh": {
                "nk1": self.nk1.text(),
                "nk2": self.nk2.text(),
                "nk3": self.nk3.text(),
                "shift1": self.shift1.text(),
                "shift2": self.shift2.text(),
                "shift3": self.shift3.text()
            },
            "mode": mode,
            "high_symmetry_points": hs,
            "path": path,
            "points_per_segment": points
        }

    def load_high_symmetry_points(self, lattice_type):
        if lattice_type.upper() == "FCC":
            points = {
                "Γ": (0.0, 0.0, 0.0),
                "X": (0.5, 0.0, 0.5),
                "W": (0.5, 0.25, 0.75),
                "K": (0.375, 0.375, 0.75),
                "L": (0.5, 0.5, 0.5),
                "U": (0.625, 0.25, 0.625)
            }
            lines = [f"{label} {x:.3f} {y:.3f} {z:.3f}" for label, (x, y, z) in points.items()]
            self.hs_points.setPlainText("\n".join(lines))
            self.path_line.setText("Γ X W K Γ L U W L")
        else:
            self.hs_points.setPlainText("")

class SetupTabs(QTabWidget):
    def __init__(self):
        super().__init__()
        self.electronic_tab = ElectronicTab()
        self.system_tab = SystemTab()
        self.general_tab = GeneralTab()
        self.kpoints_tab = KPointsTab()

        self.addTab(self.electronic_tab, "ELECTRON")
        self.addTab(self.system_tab, "SYSTEM")
        self.addTab(self.general_tab, "LATTICE")
        self.addTab(self.kpoints_tab, "K-PATH")
        # Add other tabs like &ELECTRONS here later

