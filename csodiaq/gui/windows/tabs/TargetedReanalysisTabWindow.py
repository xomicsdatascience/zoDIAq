from csodiaq.gui.windows.tabs.TabWindow import TabWindow
from PyQt5.QtWidgets import (
    QLabel,
    QPushButton,
    QLineEdit,
    QCheckBox,
)
class TargetedReanalysisTabWindow(TabWindow):
    def set_file_layout(self, fileLayout) -> None:
        self.add_identification_output_directory_field(fileLayout)

    def set_setting_layout(self, settingLayout) -> None:
        self.add_max_peptide_per_protein_setting_field(settingLayout)
        self.add_bin_value_proximity_field(settingLayout)
        self.add_heavy_isotope_field(settingLayout)

    def set_args(self) -> list:
        return []

    def add_identification_output_directory_field(self, fileLayout):
        self.scoringOutputDirText = QLabel('CsoDIAq Scoring Output Directory:')
        self.scoringOutputDirBtn = QPushButton('Browse')
        self.scoringOutputDirBtn.clicked.connect(lambda: self.open_browser_to_find_file_or_directory(self.libFile))
        self.scoringOutputDirFile = QLineEdit()
        fileLayout.addRow(self.scoringOutputDirText, self.scoringOutputDirBtn)
        fileLayout.addRow(self.scoringOutputDirFile)

    def add_max_peptide_per_protein_setting_field(self, settingLayout):
        self.proteinText = QLabel(
            'Maximum number of peptides per protein:')
        self.protein = QLineEdit()
        self.protein.setPlaceholderText("no protein analysis")
        settingLayout.addRow(self.proteinText, self.protein)

    def add_bin_value_proximity_field(self, settingLayout):
        self.binValueText = QLabel(
            'Proximity m/z values should be to a bin value:')
        self.binValue = QLineEdit()
        self.binValue.setPlaceholderText("0.75")
        settingLayout.addRow(self.binValueText, self.binValue)

    def add_heavy_isotope_field(self, settingLayout):
        self.heavyCheckBox = QCheckBox()
        settingLayout.addRow('Include Heavy Isotopes for SILAC protocol:', self.heavyCheckBox)
        self.heavyCheckBox.setCheckState(2)




