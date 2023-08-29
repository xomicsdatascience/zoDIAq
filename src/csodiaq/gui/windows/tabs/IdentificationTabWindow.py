from csodiaq.gui.windows.tabs.TabWindow import TabWindow
from csodiaq.csodiaqParser import (
    _InputQueryFile,
    _LibraryFile,
    _OutputDirectory,
    _RestrictedFloat,
)
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
        args = ["id"]
        args.extend(self.get_dia_input_file_args())
        args.extend(self.get_arg_from_text_field_if_present(self.libFile, "-l"))
        args.extend(self.get_arg_from_text_field_if_present(self.outDir, "-o"))
        args.extend(self.get_arg_from_text_field_if_present(self.matchTolerance, "-t"))
        args.extend(
            self.get_arg_from_text_field_if_present(self.correctionDegree, "-c")
        )
        args.extend(self.get_flag_from_checkbox_if_checked(self.histCheckBox, "-hist"))
        shouldDoCorrection = self.get_flag_from_checkbox_if_checked(
            self.correctionCheckBox, "-correction"
        )
        if len(shouldDoCorrection) == 0:
            args.append("--noCorrection")
        return args

    def check_args_for_invalid_input(self, args):
        argsConfirmedChecks = []
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args, "-i", _InputQueryFile(), self.diaFileText, isRequired=True
            )
        )
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args, "-l", _LibraryFile(), self.libFileText, isRequired=True
            )
        )
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args, "-o", _OutputDirectory("id"), self.outDirText, isRequired=True
            )
        )
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args,
                "-t",
                _RestrictedFloat("matchTolerance", minValue=1, maxValue=60),
                self.matchToleranceText,
            )
        )
        argsConfirmedChecks.append(
            self.check_if_arg_is_invalid_using_parsing_object(
                args,
                "-c",
                _RestrictedFloat("correctionDegree", minValue=0.5, maxValue=2),
                self.correctionDegreeText,
            )
        )
        return 0 not in argsConfirmedChecks

    def add_dia_input_file_field(self, fileLayout):
        self.diaFiles = QListWidget()
        self.delBtn = QPushButton("Delete")
        self.delBtn.clicked.connect(self.delete_selected_values)
        self.clrBtn = QPushButton("Clear")
        self.clrBtn.clicked.connect(self.diaFiles.clear)
        self.diaFileText = QLabel("DIA Data File:")
        self.diaBtn = QPushButton("Browse")
        self.diaBtn.clicked.connect(
            lambda: self.open_browser_to_find_file_or_directory(self.diaFiles)
        )
        self.diaFileText.setToolTip(
            "Data File Requirements:\n-Required, must be populated\n-Must be of type MzXML"
        )
        fileLayout.addRow(self.diaFileText, self.diaBtn)
        fileLayout.addRow(self.diaFiles)
        fileLayout.addRow(self.delBtn, self.clrBtn)

    def add_library_file_field(self, fileLayout):
        self.libFileText = QLabel("Library File:")
        self.libFile = QLineEdit()
        self.libBtn = QPushButton("Browse")
        self.libBtn.clicked.connect(
            lambda: self.open_browser_to_find_file_or_directory(self.libFile)
        )
        self.libFileText.setToolTip(
            "Library Spectra File Requirements:\n-Required, must be populated\n-Must be in MGF (.mgf) or TraML (.csv or .tsv) format"
        )
        fileLayout.addRow(self.libFileText, self.libBtn)
        fileLayout.addRow(self.libFile)

    def add_output_directory_field(self, fileLayout):
        self.outDirText = QLabel("Outfile Directory:")
        self.dirBtn = QPushButton("Browse")
        self.outDirText.setToolTip(
            "Data File Requirements:\n-Required, must be populated"
        )
        self.outDir = QLineEdit()
        self.dirBtn.clicked.connect(
            lambda: self.open_browser_to_find_file_or_directory(
                self.outDir, isFile=False
            )
        )
        fileLayout.addRow(self.outDirText, self.dirBtn)
        fileLayout.addRow(self.outDir)

    def add_match_tolerance_setting_field(self, settingLayout):
        self.matchToleranceText = QLabel(
            "Pre-correction m/z difference tolerance for peak matches (in ppm):"
        )
        self.matchTolerance = QLineEdit()
        self.matchTolerance.setPlaceholderText("30")
        settingLayout.addRow(self.matchToleranceText, self.matchTolerance)

    def add_correction_setting_fields(self, settingLayout):
        self.correctionDegree = QLineEdit()
        self.correctionCheckBox = QCheckBox()
        self.correctionDegreeText = QLabel("Corrective Standard Deviations:")
        self.correctionDegreeText.setToolTip(
            "Corrective Standard Deviations Requirements:\n-Must be blank or a float (decimal value) between or equal to 0.5 and 2.0"
        )
        self.histCheckBox = QCheckBox()
        self.correctionCheckBox.stateChanged.connect(
            lambda: self.enable_line_edits_only_if_checkbox_is_checked(
                self.correctionCheckBox, self.correctionDegree
            )
        )
        self.correctionCheckBox.stateChanged.connect(
            lambda: self.enable_second_checkbox_only_if_first_checkbox_is_checked(
                self.correctionCheckBox, self.histCheckBox
            )
        )
        settingLayout.addRow("Correction (recommended):", self.correctionCheckBox)
        settingLayout.addRow(self.correctionDegreeText, self.correctionDegree)
        settingLayout.addRow("Create Histogram:", self.histCheckBox)
        self.correctionDegree.setPlaceholderText("customized")
        self.correctionCheckBox.setCheckState(2)
        self.histCheckBox.setCheckState(2)

    def get_dia_input_file_args(self):
        files = [self.diaFiles.item(i).text() for i in range(self.diaFiles.count())]
        return [idArg for file in files for idArg in ("-i", file)]

    def enable_line_edits_only_if_checkbox_is_checked(
        self, checkBox, lineEdit, filledText=""
    ):
        if checkBox.isChecked():
            lineEdit.setText(filledText)
            lineEdit.setPlaceholderText("customized")
            lineEdit.setEnabled(True)
        else:
            lineEdit.clear()
            lineEdit.setPlaceholderText("no correction")
            lineEdit.setEnabled(False)

    def enable_second_checkbox_only_if_first_checkbox_is_checked(
        self, checkBox1, checkBox2
    ):
        if checkBox1.isChecked():
            checkBox2.setEnabled(True)
        else:
            checkBox2.setCheckState(False)
            checkBox2.setEnabled(False)
