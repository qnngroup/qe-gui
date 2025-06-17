# tabs/run_yambo_tab.py

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QPushButton, QFileDialog,
    QComboBox, QHBoxLayout, QMessageBox
)
import subprocess
import os

class RunYamboTab(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Run Yambo Simulation")
        self.setMinimumWidth(500)

        layout = QVBoxLayout()

        self.label = QLabel("Select QE output directory (./tmp/prefix/prefix.save):")
        layout.addWidget(self.label)

        self.select_button = QPushButton("Select Folder")
        self.select_button.clicked.connect(self.select_folder)
        layout.addWidget(self.select_button)

        self.sim_type_combo = QComboBox()
        self.sim_type_combo.addItems(["RPA", "GW (G0W0)", "BSE (Excitons)"])
        layout.addWidget(QLabel("Simulation type:"))
        layout.addWidget(self.sim_type_combo)

        self.run_button = QPushButton("Run p2y + Yambo")
        self.run_button.setEnabled(False)
        self.run_button.clicked.connect(self.run_yambo)
        layout.addWidget(self.run_button)

        self.setLayout(layout)

    def select_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select QE Output Directory")
        if folder:
            self.qe_output_path = folder
            self.label.setText(f"Selected: {folder}")
            self.run_button.setEnabled(True)

    def run_yambo(self):
    

        sim_type = self.sim_type_combo.currentText()
        print(sim_type)
        input_files = {
            "RPA": "rpa.in",
            "GW (G0W0)": "gw.in",
            "BSE (Excitons)": "bse.in"
        }
        input_file = input_files.get(sim_type)
        print(input_file)

        if not input_file:
            QMessageBox.critical(self, "Invalid Simulation Type", "Unrecognized simulation type selected.")
            return

        in_path = os.path.join(self.qe_output_path, input_file)
        job_name = os.path.splitext(input_file)[0]

        if not os.path.isfile(in_path):
            QMessageBox.critical(self, "Input File Not Found", f"Expected input file not found: {in_path}Please generate it first.")
            return

        try:
            subprocess.run(["p2y"], cwd=self.qe_output_path, check=True)
            subprocess.run(["yambo", "-F", input_file, "-J", job_name], cwd=self.qe_output_path, check=True)
            print("✔ Yambo run complete.")
        except subprocess.CalledProcessError as e:
            QMessageBox.critical(self, "Yambo Error", f"Yambo failed to run:{e}")

   
          



