from PyQt5.QtWidgets import (
    QLabel,
    QWidget,
    QFormLayout,
    QPushButton,
    QFileDialog,
    QPlainTextEdit,
    QVBoxLayout,
    QListWidget,
    QLineEdit,
)
from PyQt5.QtCore import QProcess
from PyQt5.QtGui import QFont
from abc import ABC, ABCMeta, abstractmethod
import numpy as np


class AbstractTabWindow(ABCMeta, type(QWidget)):
    pass


class TabWindow(QWidget, metaclass=AbstractTabWindow):
    def __init__(self):
        super().__init__()
        self.process = None
        self.VALID_COLOR = "lightgreen"
        self.INVALID_COLOR = "rgb(255,99,71)"
        fullTabLayout = QVBoxLayout()
        self.add_instruction_link_field(fullTabLayout)
        self.add_file_input_fields(fullTabLayout)
        self.add_setting_input_fields(fullTabLayout)
        self.add_process_running_field(fullTabLayout)
        self.setLayout(fullTabLayout)

    def add_instruction_link_field(self, fullTabLayout):
        instructionLabel = QLabel("Link to Instruction Wiki:")
        instructionLabel.setFont(QFont("Times", 17, QFont.Bold))
        self.instructionLink = QLineEdit()
        self.instructionLink.setReadOnly(True)
        self.set_instruction_link(self.instructionLink)
        fullTabLayout.addWidget(instructionLabel)
        fullTabLayout.addWidget(self.instructionLink)

    def add_file_input_fields(self, fullTabLayout):
        fileLayout = QFormLayout()
        fileInputLabel = QLabel("File Input:")
        fileInputLabel.setFont(QFont("Times", 17, QFont.Bold))
        fullTabLayout.addWidget(fileInputLabel)
        self.set_file_layout(fileLayout)
        fullTabLayout.addLayout(fileLayout)

    def add_setting_input_fields(self, fullTabLayout):
        settingsInputLabel = QLabel("Settings:")
        settingsInputLabel.setFont(QFont("Times", 17, QFont.Bold))
        fullTabLayout.addWidget(settingsInputLabel)
        settingLayout = QFormLayout()
        self.set_setting_layout(settingLayout)
        fullTabLayout.addLayout(settingLayout)

    def add_process_running_field(self, fullTabLayout):
        self.runBtn = QPushButton("Execute")
        self.killBtn = QPushButton("Kill Process")
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        fullTabLayout.addWidget(self.runBtn)
        fullTabLayout.addWidget(self.killBtn)
        fullTabLayout.addWidget(self.text)
        self.runBtn.setDefault(True)
        self.killBtn.setEnabled(False)
        self.killBtn.clicked.connect(self.kill_process)
        self.runBtn.clicked.connect(self.start_process)

    def set_instruction_link(self, instructionLinkText: type(QFormLayout)) -> None:
        instructionLinkText.setText("https://github.com/xomicsdatascience/zoDIAq")

    @abstractmethod
    def set_file_layout(self, fileLayout: type(QFormLayout)) -> None:
        pass

    @abstractmethod
    def set_setting_layout(self, settingLayout: QFormLayout) -> None:
        pass

    @abstractmethod
    def set_args(self) -> list:
        pass

    @abstractmethod
    def check_args_for_invalid_input(self, args) -> bool:
        pass

    def message(self, s):
        self.text.appendPlainText(s)

    def kill_process(self):
        self.process.kill()

    def process_finished(self):
        self.message("Process finished.")
        self.killBtn.setEnabled(False)
        self.process = None

    def handle_stderr(self):
        data = self.process.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.process.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

    def start_process(self):
        args = self.set_args()
        if not self.check_if_process_is_running() and self.check_args_for_invalid_input(
            args
        ):
            self.killBtn.setEnabled(True)
            self.message("Executing process")
            self.message("zodiaq " + " ".join(args))
            self.process = QProcess()
            self.process.readyReadStandardOutput.connect(self.handle_stdout)
            self.process.readyReadStandardError.connect(self.handle_stderr)
            self.process.finished.connect(self.process_finished)
            self.process.start("zodiaq", args)

    def check_if_process_is_running(self):
        return self.process

    def delete_selected_values(self):
        for x in self.diaFiles.selectedItems():
            self.diaFiles.takeItem(self.diaFiles.row(x))

    def open_browser_to_find_file_or_directory(self, text, isFile=True):
        dialog = QFileDialog()
        if type(text) is QListWidget and isFile:
            fname = QFileDialog.getOpenFileNames(self, "Open File", "c:\\")
            text.addItems([x for x in fname[0]])
        elif type(text) is QLineEdit and isFile:
            fname = QFileDialog.getOpenFileName(self, "Open File", "c:\\")
            text.setText(fname[0])
        else:
            fname = QFileDialog.getExistingDirectory(self, "Open Directory", "c:\\")
            text.setText(fname)

    def get_arg_from_text_field_if_present(self, textObject, flag):
        if textObject.text():
            return [flag, textObject.text()]
        else:
            return []

    def get_flag_from_checkbox_if_checked(self, checkBox, flag):
        if checkBox.isChecked():
            return [flag]
        else:
            return []

    def check_if_arg_is_invalid_using_parsing_object(
        self, args, targetFlag, parsingObject, textObject, isRequired=False
    ):
        targetIndices = (
            np.array([index for (index, flag) in enumerate(args) if flag == targetFlag])
            + 1
        )
        targets = [args[i] for i in targetIndices]
        try:
            if len(targets) == 0 and isRequired:
                raise ValueError
            for target in targets:
                parsingObject(target)
        except:
            textObject.setStyleSheet("background-color: " + self.INVALID_COLOR)
            return 0
        textObject.setStyleSheet("background-color: " + self.VALID_COLOR)
        return 1
