from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QLineEdit, QPushButton,
    QTabWidget, QFormLayout, QMessageBox, QFileDialog, QHBoxLayout
)
import os
from tabs.tab_functions.plot_bands import PlotBands

class PlottingWidget(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Plot Band Structure")
        self.setMinimumSize(400, 200)

        layout = QVBoxLayout()
        tabs = QTabWidget()

        # ---- Only Tab: Plot bands from .dat file ----
        plot_tab = QVBoxLayout()
        plot_form = QFormLayout()

        self.dat_file = QLineEdit()
        dat_browse = QPushButton("Browse")
        dat_browse.clicked.connect(lambda: self.select_file(self.dat_file, "outputs", "Data Files (*.dat);;All Files (*)"))
        dat_row = QHBoxLayout()
        dat_row.addWidget(self.dat_file)
        dat_row.addWidget(dat_browse)
        plot_form.addRow("bands.dat file:", dat_row)

        self.bands_in = QLineEdit()
        bands_browse = QPushButton("Browse")
        bands_browse.clicked.connect(lambda: self.select_file(self.bands_in, "inputs", "Input Files (*.in);;All Files (*)"))
        bands_row = QHBoxLayout()
        bands_row.addWidget(self.bands_in)
        bands_row.addWidget(bands_browse)
        plot_form.addRow("bands.in file:", bands_row)

        plot_button = QPushButton("📊 Plot Band Structure")
        plot_button.clicked.connect(self.run_plot_bands)
        plot_tab.addLayout(plot_form)
        plot_tab.addWidget(plot_button)

        plot_tab_container = QDialog()
        plot_tab_container.setLayout(plot_tab)
        tabs.addTab(plot_tab_container, "Plot bands")

        layout.addWidget(tabs)
        self.setLayout(layout)

    def select_file(self, line_edit, subfolder, file_filter):
        start_dir = os.path.join(os.getcwd(), subfolder)
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", start_dir, file_filter)
        if file_path:
            line_edit.setText(file_path)

    def run_plot_bands(self):
        dat_path = self.dat_file.text().strip()
        bands_path = self.bands_in.text().strip()

        if not os.path.isfile(dat_path):
            QMessageBox.critical(self, "Missing File", f"Cannot find bands.dat file: {dat_path}")
            return
        if not os.path.isfile(bands_path):
            QMessageBox.critical(self, "Missing File", f"Cannot find bands.in file: {bands_path}")
            return

        try:
            PlotBands(bands_file=bands_path, dat_file=dat_path)
            QMessageBox.information(self, "Success", "Band structure plot generated.")
            self.accept()  
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Plotting failed:\n{e}")


