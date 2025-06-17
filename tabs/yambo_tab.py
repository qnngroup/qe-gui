from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QPushButton, QFileDialog,
    QFormLayout, QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTabWidget, QWidget, QHBoxLayout, QCheckBox
)
import subprocess
import os

class YamboTab(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Yambo Simulation Setup")
        self.setMinimumWidth(600)

        main_layout = QVBoxLayout()

        # --- Directory Selection ---
        self.label = QLabel("Select a Quantum ESPRESSO output directory (prefix.save):")
        main_layout.addWidget(self.label)

        self.select_button = QPushButton("Select Folder")
        self.select_button.clicked.connect(self.select_folder)
        main_layout.addWidget(self.select_button)

        # --- Simulation Type Dropdown ---
        self.sim_type_combo = QComboBox()
        self.sim_type_combo.addItems(["RPA Dielectric Function", "GW (G0W0)", "BSE (Excitons)"])
        self.sim_type_combo.currentTextChanged.connect(self.switch_tab)
        main_layout.addWidget(QLabel("Select simulation type:"))
        main_layout.addWidget(self.sim_type_combo)

        # --- Stacked Tabs for Simulation Settings ---
        self.tabs = QTabWidget()
        self.tabs.setTabBarAutoHide(True)

        self.rpa_tab = self.build_rpa_tab()
        self.gw_tab = self.build_gw_tab()
        self.bse_tab = self.build_bse_tab()

        self.tabs.addTab(self.rpa_tab, "RPA")
        self.tabs.addTab(self.gw_tab, "GW")
        self.tabs.addTab(self.bse_tab, "BSE")

        main_layout.addWidget(self.tabs)

        # --- Run Buttons ---
        self.run_button = QPushButton("Generate Input File")
        self.run_button.clicked.connect(self.run_yambo)
        self.run_button.setEnabled(False)
        main_layout.addWidget(self.run_button)

        self.setLayout(main_layout)

    def select_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select QE Output Directory")
        if folder:
            self.qe_output_path = folder
            self.label.setText(f"Selected: {folder}")
            self.run_button.setEnabled(True)

    def switch_tab(self, text):
        index = {"RPA Dielectric Function": 0, "GW (G0W0)": 1, "BSE (Excitons)": 2}[text]
        self.tabs.setCurrentIndex(index)

    def run_yambo(self):
        self.input_data = self.get_values()
        self.accept()

    def build_rpa_tab(self):
        tab = QWidget()
        layout = QFormLayout()

        self.include_optics = QCheckBox("Include optics module (ε₁, ε₂)")
        layout.addRow("", self.include_optics)

        self.q_start = QSpinBox()
        self.q_end = QSpinBox()
        self.q_start.setRange(1, 9999)
        self.q_end.setRange(1, 9999)
        self.q_start.setValue(1)
        self.q_end.setValue(217)
        q_widget = QWidget()
        q_layout = QHBoxLayout()
        q_layout.addWidget(self.q_start)
        q_layout.addWidget(QLabel("to"))
        q_layout.addWidget(self.q_end)
        q_widget.setLayout(q_layout)
        layout.addRow("Q-point index range:", q_widget)

        self.band_min = QSpinBox()
        self.band_max = QSpinBox()
        self.band_min.setRange(1, 9999)
        self.band_max.setRange(1, 9999)
        self.band_min.setValue(1)
        self.band_max.setValue(20)
        band_widget = QWidget()
        band_layout = QHBoxLayout()
        band_layout.addWidget(self.band_min)
        band_layout.addWidget(QLabel("to"))
        band_layout.addWidget(self.band_max)
        band_widget.setLayout(band_layout)
        layout.addRow("Bands used in polarization function:", band_widget)

        self.g_cutoff = QSpinBox()
        self.g_cutoff.setValue(1)
        layout.addRow("G-vector cutoff (RL):", self.g_cutoff)

        self.polarization = QLineEdit("1.0 0.0 0.0")
        self.polarization.setToolTip("Direction of external electric field (x y z)")
        layout.addRow("Polarization direction:", self.polarization)

        self.energy_range = QLineEdit("0.0 20.0")
        self.energy_steps = QLineEdit("1000")
        self.damping = QLineEdit("0.01 0.01")
        layout.addRow("Energy range (eV):", self.energy_range)
        layout.addRow("Number of E Steps", self.energy_steps)
        layout.addRow("Damping value (eV):", self.damping)

        self.model = QComboBox()
        self.model.addItems(["RPA", "RE", "ALDA"])
        layout.addRow("Dielectric model:", self.model)

        self.temp = QDoubleSpinBox()
        self.temp.setValue(0.0)
        layout.addRow("Electronic temperature (eV):", self.temp)

        tab.setLayout(layout)
        return tab

    def build_gw_tab(self):
        tab = QWidget()
        layout = QFormLayout()

        self.gw_band_min = QSpinBox()
        self.gw_band_max = QSpinBox()
        self.gw_band_min.setRange(1, 9999)
        self.gw_band_max.setRange(1, 9999)
        self.gw_band_min.setValue(1)
        self.gw_band_max.setValue(100)
        band_widget = QWidget()
        band_layout = QHBoxLayout()
        band_layout.addWidget(self.gw_band_min)
        band_layout.addWidget(QLabel("to"))
        band_layout.addWidget(self.gw_band_max)
        band_widget.setLayout(band_layout)
        layout.addRow("Bands for GW correction:", band_widget)

        self.exx_cutoff = QSpinBox()
        self.exx_cutoff.setValue(40)
        layout.addRow("Exchange cutoff (EXXRLvcs, Ry):", self.exx_cutoff)

        self.ppa_energy = QDoubleSpinBox()
        self.ppa_energy.setDecimals(2)
        self.ppa_energy.setValue(12.0)
        layout.addRow("Plasmon-pole energy (eV):", self.ppa_energy)

        self.qpkrange = QLineEdit("1|10|1|10")
        self.qpkrange.setToolTip("Format: kstart|kend|bandstart|bandend")
        layout.addRow("QP correction range:", self.qpkrange)

        self.gterm_kind = QComboBox()
        self.gterm_kind.addItems(["none", "BG"])
        layout.addRow("G-term kind:", self.gterm_kind)

        tab.setLayout(layout)
        return tab

    def build_bse_tab(self):
        tab = QWidget()
        layout = QFormLayout()

        self.bse_band_min = QSpinBox()
        self.bse_band_max = QSpinBox()
        self.bse_band_min.setRange(1, 9999)
        self.bse_band_max.setRange(1, 9999)
        self.bse_band_min.setValue(1)
        self.bse_band_max.setValue(10)
        band_widget = QWidget()
        band_layout = QHBoxLayout()
        band_layout.addWidget(self.bse_band_min)
        band_layout.addWidget(QLabel("to"))
        band_layout.addWidget(self.bse_band_max)
        band_widget.setLayout(band_layout)
        layout.addRow("Bands used in BSE:", band_widget)

        self.bse_kernel = QComboBox()
        self.bse_kernel.addItems(["IP", "ALDA", "SEX"])
        layout.addRow("BSE kernel type:", self.bse_kernel)

        self.bseng_block = QSpinBox()
        self.bseng_block.setValue(1)
        layout.addRow("BSE G-vector block (RL):", self.bseng_block)

        self.kfnqp_input = QLineEdit("gwqps")
        self.kfnqp_input.setToolTip("Input file with GW quasiparticle corrections")
        layout.addRow("GW correction file:", self.kfnqp_input)

        self.krange_input = QLineEdit("1|100")
        self.krange_input.setToolTip("K-point range in format: kstart|kend")
        layout.addRow("K-point range:", self.krange_input)

        tab.setLayout(layout)
        return tab

    def get_values(self):
        sim_type = self.sim_type_combo.currentText()
        values = {
            "simulation_type": sim_type,
            "qe_output_path": getattr(self, "qe_output_path", None)
        }

        if sim_type == "RPA Dielectric Function":
            values.update({
                "include_optics": self.include_optics.isChecked(),
                "q_start": self.q_start.value(),
                "q_end": self.q_end.value(),
                "band_min": self.band_min.value(),
                "band_max": self.band_max.value(),
                "g_cutoff": self.g_cutoff.value(),
                "polarization": self.polarization.text().strip(),
                "energy_range": self.energy_range.text().strip(),
                "damping": self.damping.text().strip(),
                "dielectric_model": self.model.currentText(),
                "electronic_temp": self.temp.value(),
                "nstps": self.energy_steps.text()
            })

        elif sim_type == "GW (G0W0)":
            values.update({
                "gw_band_min": self.gw_band_min.value(),
                "gw_band_max": self.gw_band_max.value(),
                "exx_cutoff": self.exx_cutoff.value(),
                "ppa_energy": self.ppa_energy.value(),
                "qpkrange": self.qpkrange.text().strip(),
                "gterm_kind": self.gterm_kind.currentText()
            })

        elif sim_type == "BSE (Excitons)":
            values.update({
                "bse_band_min": self.bse_band_min.value(),
                "bse_band_max": self.bse_band_max.value(),
                "bse_kernel": self.bse_kernel.currentText(),
                "bseng_block": self.bseng_block.value(),
                "kfnqp_input": self.kfnqp_input.text().strip(),
                "krange_input": self.krange_input.text().strip()
            })

        return values
