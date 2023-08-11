from csodiaq.gui.windows.tabs.TabWindow import TabWindow
from PyQt5.QtWidgets import (
    QLabel,
    QPushButton,
    QLineEdit,
)
class ScoringTabWindow(TabWindow):
    def set_file_layout(self, fileLayout):
        self.add_identification_output_directory_field(fileLayout)

    def set_setting_layout(self, settingLayout):
        self.add_no_setting_disclaimer_field(settingLayout)

    def set_args(self) -> list:
        return []

    def add_identification_output_directory_field(self, fileLayout):
        self.idOutputDirText = QLabel('CsoDIAq Identification Output Directory:')
        self.idOutputDirBtn = QPushButton('Browse')
        self.idOutputDirBtn.clicked.connect(lambda: self.open_browser_to_find_file_or_directory(self.libFile))
        self.idOutputDirFile = QLineEdit()
        fileLayout.addRow(self.idOutputDirText, self.idOutputDirBtn)
        fileLayout.addRow(self.idOutputDirFile)

    def add_no_setting_disclaimer_field(self, settingLayout):
        disclaimerText = QLabel('No settings currently implemented for the score step.')
        settingLayout.addRow(disclaimerText)
