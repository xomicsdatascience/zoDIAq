from zodiaq.gui.windows.tabs.TabWindow import TabWindow
from zodiaq.zodiaqParser import (
    _ScoringOutputDirectory,
    _RestrictedInt,
    _RestrictedFloat,
)
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
        args = ["targetedReanalysis"]
        args.extend(
            self.get_arg_from_text_field_if_present(self.scoringOutputDir, "-i")
        )
        args.extend(self.get_arg_from_text_field_if_present(self.protein, "-p"))
        args.extend(
            self.get_flag_from_checkbox_if_checked(self.heavyCheckBox, "-heavy")
        )
        args.extend(self.get_arg_from_text_field_if_present(self.binValue, "-b"))
        return args

    def check_args_for_invalid_input(self, args):
        argsConfirmedChecks = []
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args,
                "-i",
                _ScoringOutputDirectory(),
                self.scoringOutputDirText,
                isRequired=True,
            )
        )
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args, "-p", _RestrictedInt("protein", minValue=1), self.proteinText
            )
        )
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args,
                "-b",
                _RestrictedFloat("binValueProximity", minValue=0.01),
                self.binValueText,
            )
        )
        return 0 not in argsConfirmedChecks

    def add_identification_output_directory_field(self, fileLayout):
        self.scoringOutputDirText = QLabel("zoDIAq Scoring Output Directory:")
        self.scoringOutputDir = QLineEdit()
        self.scoringOutputDirBtn = QPushButton("Browse")
        self.scoringOutputDirBtn.clicked.connect(
            lambda: self.open_browser_to_find_file_or_directory(
                self.scoringOutputDir, isFile=False
            )
        )
        fileLayout.addRow(self.scoringOutputDirText, self.scoringOutputDirBtn)
        fileLayout.addRow(self.scoringOutputDir)

    def add_max_peptide_per_protein_setting_field(self, settingLayout):
        self.proteinText = QLabel("Maximum number of peptides per protein:")
        self.protein = QLineEdit()
        self.protein.setPlaceholderText("no protein analysis")
        settingLayout.addRow(self.proteinText, self.protein)

    def add_bin_value_proximity_field(self, settingLayout):
        self.binValueText = QLabel("Proximity m/z values should be to a bin value:")
        self.binValue = QLineEdit()
        self.binValue.setPlaceholderText("0.75")
        settingLayout.addRow(self.binValueText, self.binValue)

    def add_heavy_isotope_field(self, settingLayout):
        self.heavyCheckBox = QCheckBox()
        settingLayout.addRow(
            "Include Heavy Isotopes for SILAC protocol:", self.heavyCheckBox
        )
        self.heavyCheckBox.setCheckState(2)
