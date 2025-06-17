from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton,
    QFormLayout, QMessageBox, QFileDialog, QTabWidget
)
import subprocess
import os


class RunDFTWidget(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DFT Utilities")
        self.setMinimumSize(400, 200)

        layout = QVBoxLayout()
        self.tabs = QTabWidget()

        self.tabs.addTab(self.create_pw_tab(), "Run SCF + NSCF")
        self.tabs.addTab(self.create_bands_tab(), "Run bands.x")

        layout.addWidget(self.tabs)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        layout.addWidget(close_button)

        self.setLayout(layout)

    def create_pw_tab(self):
        pw_tab = QVBoxLayout()
        form = QFormLayout()

        # SCF input
        self.scf_input = QLineEdit()
        scf_browse = QPushButton("Browse")
        scf_browse.clicked.connect(lambda: self.select_file(self.scf_input, "inputs", "Input Files (*.in);;All Files (*)"))
        scf_row = QHBoxLayout()
        scf_row.addWidget(self.scf_input)
        scf_row.addWidget(scf_browse)
        form.addRow("SCF Input File:", scf_row)

        # NSCF input
        self.nscf_input = QLineEdit()
        nscf_browse = QPushButton("Browse")
        nscf_browse.clicked.connect(lambda: self.select_file(self.nscf_input, "inputs", "Input Files (*.in);;All Files (*)"))
        nscf_row = QHBoxLayout()
        nscf_row.addWidget(self.nscf_input)
        nscf_row.addWidget(nscf_browse)
        form.addRow("NSCF Input File:", nscf_row)

        run_button = QPushButton("▶ Run pw.x on SCF + NSCF")
        run_button.clicked.connect(self.run_dft)

        pw_tab.addLayout(form)
        pw_tab.addWidget(run_button)

        container = QDialog()
        container.setLayout(pw_tab)
        return container

    def create_bands_tab(self):
        bands_tab = QVBoxLayout()
        form = QFormLayout()

        self.bands_in = QLineEdit()
        bands_browse = QPushButton("Browse")
        bands_browse.clicked.connect(lambda: self.select_file(self.bands_in, "inputs", "Input Files (*.in);;All Files (*)"))
        bands_row = QHBoxLayout()
        bands_row.addWidget(self.bands_in)
        bands_row.addWidget(bands_browse)
        form.addRow("bands.in File:", bands_row)

        run_button = QPushButton("▶ Run bands.x")
        run_button.clicked.connect(self.run_bands_x)

        bands_tab.addLayout(form)
        bands_tab.addWidget(run_button)

        container = QDialog()
        container.setLayout(bands_tab)
        return container

    def select_file(self, line_edit, subfolder, file_filter):
        start_dir = os.path.join(os.getcwd(), subfolder)
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", start_dir, file_filter)
        if file_path:
            line_edit.setText(file_path)

    def run_dft(self):
        scf = self.scf_input.text().strip()
        nscf = self.nscf_input.text().strip()

        for file in [scf, nscf]:
            if not os.path.isfile(file):
                QMessageBox.critical(self, "Error", f"Missing input file:\n{file}")
                return

        try:
            subprocess.run(f"pw.x < {scf} > {scf.replace('.in', '.out')}", shell=True, check=True)
            subprocess.run(f"pw.x < {nscf} > {nscf.replace('.in', '.out')}", shell=True, check=True)
            QMessageBox.information(self, "Success", "SCF and NSCF calculations completed.")
        except subprocess.CalledProcessError:
            QMessageBox.critical(self, "Error", "One or both pw.x runs failed.")

    def run_bands_x(self):
        bands_in_path = self.bands_in.text().strip()
        if not os.path.isfile(bands_in_path):
            QMessageBox.critical(self, "Missing File", f"Cannot find bands.in file:\n{bands_in_path}")
            return

        try:
            subprocess.run(f"bands.x < {bands_in_path} > {bands_in_path.replace('.in', '.out')}", shell=True, check=True)
            QMessageBox.information(self, "Success", "bands.x completed successfully.")
        except subprocess.CalledProcessError:
            QMessageBox.critical(self, "Error", "bands.x failed to run.")
