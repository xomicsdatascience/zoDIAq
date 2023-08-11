from csodiaq.gui.windows.tabs.TabWindow import TabWindow
from PyQt5.QtWidgets import (
    QListWidget,
    QPushButton,
    QLabel,
    QLineEdit,
    QCheckBox,
)
class IdentificationTabWindow(TabWindow):
    def set_file_layout(self, fileLayout):
        self.add_dia_input_file_field(fileLayout)
        self.add_library_file_field(fileLayout)
        self.add_output_directory_field(fileLayout)

    def set_setting_layout(self, settingLayout):
        self.add_match_tolerance_setting_field(settingLayout)
        self.add_correction_setting_fields(settingLayout)

    def set_args(self) -> list:
        return []

    def add_dia_input_file_field(self, fileLayout):
        self.diaFiles = QListWidget()
        self.delBtn = QPushButton('Delete')
        self.delBtn.clicked.connect(self.delete_selected_values)
        self.clrBtn = QPushButton('Clear')
        self.clrBtn.clicked.connect(self.diaFiles.clear)
        self.diaFileText = QLabel('DIA Data File:')
        self.diaBtn = QPushButton('Browse')
        self.diaBtn.clicked.connect(lambda: self.open_browser_to_find_file_or_directory(self.diaFiles))
        self.diaFileText.setToolTip(
            'Data File Requirements:\n-Required, must be populated\n-Must be of type MzXML')
        fileLayout.addRow(self.diaFileText, self.diaBtn)
        fileLayout.addRow(self.diaFiles)
        fileLayout.addRow(self.delBtn, self.clrBtn)


    def add_library_file_field(self, fileLayout):
        self.libFileText = QLabel('Library File:')
        self.libBtn = QPushButton('Browse')
        self.libBtn.clicked.connect(lambda: self.open_browser_to_find_file_or_directory(self.libFile))
        self.libFileText.setToolTip(
            'Library Spectra File Requirements:\n-Required, must be populated\n-Must be in MGF (.mgf) or TraML (.csv or .tsv) format')
        self.libFile = QLineEdit()
        fileLayout.addRow(self.libFileText, self.libBtn)
        fileLayout.addRow(self.libFile)

    def add_output_directory_field(self, fileLayout):
        self.outDirText = QLabel('Outfile Directory:')
        self.dirBtn = QPushButton('Browse')
        self.outDirText.setToolTip(
            'Data File Requirements:\n-Required, must be populated')
        self.outDir = QLineEdit()
        self.dirBtn.clicked.connect(
            lambda: self.open_browser_to_find_file_or_directory(self.outDir, isFile=False))
        fileLayout.addRow(self.outDirText, self.dirBtn)
        fileLayout.addRow(self.outDir)

    def add_match_tolerance_setting_field(self, settingLayout):
        self.matchToleranceText = QLabel(
            'Pre-correction m/z difference tolerance for peak matches (in ppm):')
        self.matchTolerance = QLineEdit()
        self.matchTolerance.setPlaceholderText("30")
        settingLayout.addRow(self.matchToleranceText, self.matchTolerance)

    def add_correction_setting_fields(self, settingLayout):
        self.corr = QLineEdit()
        self.corrCheckBox = QCheckBox()
        self.corrText = QLabel('Corrective Standard Deviations:')
        self.corrText.setToolTip(
            'Corrective Standard Deviations Requirements:\n-Must be blank or a float (decimal value) between or equal to 0.5 and 2.0')
        self.histCheckBox = QCheckBox()
        self.corrCheckBox.stateChanged.connect(
            lambda: self.enable_line_edits_only_if_checkbox_is_checked(self.corrCheckBox, self.corr))
        self.corrCheckBox.stateChanged.connect(
            lambda: self.enable_second_checkbox_only_if_first_checkbox_is_checked(self.corrCheckBox, self.histCheckBox))
        settingLayout.addRow('Correction (recommended):', self.corrCheckBox)
        settingLayout.addRow(self.corrText, self.corr)
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        self.corr.setPlaceholderText('customized')
        self.corrCheckBox.setCheckState(2)
        self.histCheckBox.setCheckState(2)