# main.py

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QVBoxLayout,
    QWidget, QDialog, QLabel, QDialogButtonBox, QFormLayout 
)


import sys
import subprocess
from tabs.setup_tab import SetupTabs  
from tabs.run_dft_tab import RunDFTWidget
from tabs.plotting_tab import PlottingWidget
from tabs.tab_functions.qe_writer import write_qe_inputs
from tabs.tab_functions.qe_writer import write_yambo_inputs
from tabs.run_yambo_tab import RunYamboTab
from tabs.yambo_tab import YamboTab

#Widget for building Inputs
class SetupDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("QE Setup")
        self.tabs = SetupTabs()

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)   # Accepts and closes dialog
        buttons.rejected.connect(self.reject)   # Cancels and closes dialog

        layout = QVBoxLayout()
        layout.addWidget(self.tabs)
        layout.addWidget(buttons)
        self.setLayout(layout)

    
    def get_all_inputs(self):
        prefix = self.tabs.electronic_tab.prefix.text()  # or scf_prefix if you prefer
        return {
            "electronic": self.tabs.electronic_tab.get_values(),
            "system": self.tabs.system_tab.get_values(),
            "general": self.tabs.general_tab.get_values(),
            "kpoints": self.tabs.kpoints_tab.get_values(prefix)
        }



#main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Quantum ESPRESSO Launcher")
        self.setMinimumWidth(400)

        layout = QVBoxLayout()

        open_btn = QPushButton("Quantum Espresso Input Builder")
        open_btn.clicked.connect(self.open_setup)
        layout.addWidget(open_btn)

        run_btn = QPushButton("DFT Functions")
        run_btn.clicked.connect(self.launch_run_dialog)
        layout.addWidget(run_btn)

        plot_btn = QPushButton("Band Structure Plotter")
        plot_btn.clicked.connect(self.launch_plot_dialog)
        layout.addWidget(plot_btn)

        self.yambo_button = QPushButton("Yambo Input Builder")
        self.yambo_button.clicked.connect(self.open_yambo_tab)
        layout.addWidget(self.yambo_button)

        self.yambo_run_button = QPushButton("Run Yambo Calculations")
        self.yambo_run_button.clicked.connect(self.open_run_yambo_tab)
        layout.addWidget(self.yambo_run_button)

        central = QWidget()
        central.setLayout(layout)
        self.setCentralWidget(central)
        self.input_data = None

    def open_setup(self):
        dlg = SetupDialog()
        if dlg.exec():
            self.input_data = dlg.get_all_inputs()
            write_qe_inputs(self.input_data)  # Automatically write files after setup

    def launch_run_dialog(self):
        dlg = RunDFTWidget()
        dlg.exec()

    def launch_plot_dialog(self):
        dlg = PlottingWidget()
        dlg.exec()
    def open_yambo_tab(self):
        dlg = YamboTab(self)
        if dlg.exec():
            self.input_data = dlg.input_data
            write_yambo_inputs(self.input_data)  
    def open_run_yambo_tab(self):
        dlg = RunYamboTab(self)
        dlg.exec()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.resize(500, 300)
    window.show()
    sys.exit(app.exec())
