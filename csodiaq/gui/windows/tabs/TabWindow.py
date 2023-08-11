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


class AbstractTabWindow(ABCMeta, type(QWidget)):
    pass


class TabWindow(QWidget, metaclass=AbstractTabWindow):
    def __init__(self):
        super().__init__()

        self.p = None
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
        instructionLinkText.setText(
            "https://github.com/xomicsdatascience/CsoDIAq"
        )

    @abstractmethod
    def set_file_layout(self, fileLayout: type(QFormLayout)) -> None:
        pass

    @abstractmethod
    def set_setting_layout(self, settingLayout: QFormLayout) -> None:
        pass

    @abstractmethod
    def set_args(self) -> list:
        pass

    def enable_line_edits_only_if_checkbox_is_checked(self, checkBox, lineEdit, filledText=''):
        if checkBox.isChecked():
            lineEdit.setText(filledText)
            lineEdit.setPlaceholderText('customized')
            lineEdit.setEnabled(True)
        else:
            lineEdit.clear()
            lineEdit.setPlaceholderText('no correction')
            lineEdit.setEnabled(False)

    def enable_second_checkbox_only_if_first_checkbox_is_checked(self, checkBox1, checkBox2):
        if checkBox1.isChecked():
            checkBox2.setEnabled(True)
        else:
            checkBox2.setCheckState(False)
            checkBox2.setEnabled(False)



    def message(self, s):
        self.text.appendPlainText(s)

    def kill_process(self):
        self.p.kill()

    def process_finished(self):
        self.message("Process finished.")
        self.killBtn.setEnabled(False)
        self.p = None

    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

    def start_process(self):
        args = []
        args = self.set_args()
        self.killBtn.setEnabled(True)

        if self.p is None:  # No process running.
            self.message("Executing process")
            self.message("csodiaq " + " ".join(args))
            self.p = (
                QProcess()
            )  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.
            self.p.start("csodiaq", args)

    def delete_selected_values(self):
        for x in self.diaFiles.selectedItems():
            self.diaFiles.takeItem(self.diaFiles.row(x))

    def open_browser_to_find_file_or_directory(self, text, isFile=True):
        dialog = QFileDialog()
        if type(text) is QListWidget and isFile:
            fname = QFileDialog.getOpenFileNames(self, 'Open File', 'c:\\')
            text.addItems([x for x in fname[0]])
        elif type(text) is QLineEdit and isFile:
            fname = QFileDialog.getOpenFileName(self, 'Open File', 'c:\\')
            text.setText(fname[0])
        else:
            fname = QFileDialog.getExistingDirectory(
                self, 'Open Directory', 'c:\\')
            text.setText(fname)
