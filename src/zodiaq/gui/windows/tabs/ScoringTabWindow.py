from zodiaq.gui.windows.tabs.TabWindow import TabWindow
from zodiaq.zodiaqParser import _IdentificationOutputDirectory
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
        args = ["score"]
        args.extend(self.get_arg_from_text_field_if_present(self.idOutputDir, "-i"))
        return args

    def check_args_for_invalid_input(self, args):
        argsConfirmedChecks = []
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args,
                "-i",
                _IdentificationOutputDirectory(),
                self.idOutputDirText,
                isRequired=True,
            )
        )
        return 0 not in argsConfirmedChecks

    def add_identification_output_directory_field(self, fileLayout):
        self.idOutputDirText = QLabel("zoDIAq Identification Output Directory:")
        self.idOutputDir = QLineEdit()
        self.idOutputDirBtn = QPushButton("Browse")
        self.idOutputDirBtn.clicked.connect(
            lambda: self.open_browser_to_find_file_or_directory(
                self.idOutputDir, isFile=False
            )
        )
        fileLayout.addRow(self.idOutputDirText, self.idOutputDirBtn)
        fileLayout.addRow(self.idOutputDir)

    def add_no_setting_disclaimer_field(self, settingLayout):
        disclaimerText = QLabel("No settings currently implemented for the score step.")
        settingLayout.addRow(disclaimerText)
