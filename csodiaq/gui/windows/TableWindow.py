from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QTabWidget,
)

from csodiaq.gui.windows.tabs import IdentificationTabWindow
from csodiaq.gui.windows.tabs import ScoringTabWindow
from csodiaq.gui.windows.tabs import TargetedReanalysisTabWindow

class TableWindow(QWidget):
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)

        self.tabs = QTabWidget()
        self.idTab = IdentificationTabWindow()
        self.tabs.addTab(self.idTab, "Peptide Identification")
        self.scoreTab = ScoringTabWindow()
        self.tabs.addTab(self.scoreTab, "Scoring and FDR Filtering")
        self.reanalysisTab = TargetedReanalysisTabWindow()
        self.tabs.addTab(self.reanalysisTab, "Targeted Peptide Reanalysis")

        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)