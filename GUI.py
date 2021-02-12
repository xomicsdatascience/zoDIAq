# Filename: dialog.py

"""Dialog-Style application."""

import sys

from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QDialog
from PyQt5.QtWidgets import QDialogButtonBox
from PyQt5.QtWidgets import QFormLayout
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QCheckBox
from PyQt5.QtWidgets import QFileDialog




class IdWindow(QWidget):
    """Dialog."""
    def __init__(self, parent=None):
        """Initializer."""
        super().__init__(parent)
        self.idButton = QPushButton('Peptide/Protein Identification')
        self.quantButton = QPushButton('SILAC Quantification')
        dlgLayout = QVBoxLayout()

        dlgLayout.addWidget(QLabel('Peptide/Protein Identification'))

        modeLayout = QHBoxLayout()
        modeLayout.addWidget(self.idButton)
        modeLayout.addWidget(self.quantButton)
        dlgLayout.addLayout(modeLayout)
        self.idButton.setEnabled(False)

        dlgLayout.addWidget(QLabel('File Input:'))

        self.le = QLineEdit()
        self.btn1 = QPushButton("QFileDialog object")
        self.btn1.clicked.connect(self.getfile)

        fileLayout = QFormLayout()
        fileLayout.addRow('Sample File Choosing Button:', self.btn1)
        fileLayout.addRow('Chosen File:', self.le)
        fileLayout.addRow('DIA Data File:', QLineEdit())
        fileLayout.addRow('Library File:', QLineEdit())
        fileLayout.addRow('Outfile Directory:', QLineEdit())
        dlgLayout.addLayout(fileLayout)

        dlgLayout.addWidget(QLabel('Settings:'))


        self.histCheckBox = QCheckBox()
        self.histOutFile = QLineEdit()
        self.histCheckBox.stateChanged.connect(lambda:self.check_grey(self.histCheckBox, self.histOutFile, 'no histogram'))
        self.protCheckBox = QCheckBox()
        self.protOutFile = QLineEdit()
        self.protCheckBox.stateChanged.connect(lambda:self.check_grey(self.protCheckBox, self.protOutFile, 'peptides, no protein inference','1'))


        settingLayout = QFormLayout()


        settingLayout.addRow('Parent Mass Tolerance (in DA):', QLineEdit())
        settingLayout.addRow('Initial Fragment Mass Tolerance (in ppm):', QLineEdit())
        settingLayout.addRow('Corrective Standard Deviations:', QLineEdit())
        settingLayout.addRow('Number of Scans in One Cycle:', QLineEdit())
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        settingLayout.addRow('Histogram Outfile: ', self.histOutFile)
        settingLayout.addRow('Protein Inference:', self.protCheckBox)
        settingLayout.addRow('Number of Represented Peptides per Protein: ', self.protOutFile)

        dlgLayout.addLayout(settingLayout)
        self.histOutFile.setText('no histogram')
        self.histOutFile.setEnabled(False)
        self.protOutFile.setText('peptides, no protein inference')
        self.protOutFile.setEnabled(False)
        self.le.setEnabled(False)

        btns = QDialogButtonBox()
        btns.setStandardButtons(
            QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        dlgLayout.addWidget(btns)
        self.setLayout(dlgLayout)

    def check_grey(self, checkBox, lineEdit, emptyText, filledText=''):
        if checkBox.isChecked():
            lineEdit.setText(filledText)
            lineEdit.setEnabled(True)
        else:
            lineEdit.setText(emptyText)
            lineEdit.setEnabled(False)

    def getfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file',
            'c:\\')
        self.le.setText(fname[0])

class quantWindow(QWidget):
    """Dialog."""
    def __init__(self, parent=None):
        """Initializer."""
        super().__init__(parent)
        self.idButton = QPushButton('Peptide/Protein Identification')
        self.quantButton = QPushButton('SILAC Quantification')
        dlgLayout = QVBoxLayout()

        dlgLayout.addWidget(QLabel('SILAC Quantification'))

        modeLayout = QHBoxLayout()
        modeLayout.addWidget(self.idButton)
        modeLayout.addWidget(self.quantButton)
        dlgLayout.addLayout(modeLayout)
        self.quantButton.setEnabled(False)

        formLayout = QFormLayout()
        formLayout.addRow('Outfile Directory and Header:', QLineEdit())
        formLayout.addRow('Library Spectrum Path:', QLineEdit())
        formLayout.addRow('FDR Output:', QLineEdit())
        formLayout.addRow('DISPA Targetted Re-Analysis Directory:', QLineEdit())
        dlgLayout.addLayout(formLayout)
        btns = QDialogButtonBox()
        btns.setStandardButtons(
            QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        dlgLayout.addWidget(btns)
        self.setLayout(dlgLayout)

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.startIdWindow()

    def startIdWindow(self):
        self.IdWindow = IdWindow(self)
        self.setWindowTitle('CsoDIAq')
        self.setCentralWidget(self.IdWindow)
        self.IdWindow.quantButton.clicked.connect(self.startQuantWindow)
        self.show()

    def startQuantWindow(self):
        self.quantWindow = quantWindow(self)
        self.setWindowTitle('CsoDIAq')
        self.setCentralWidget(self.quantWindow)
        self.quantWindow.idButton.clicked.connect(self.startIdWindow)
        self.show()
if __name__ == '__main__':
    app = QApplication(sys.argv)
    dlg = MainWindow()
    dlg.show()
    sys.exit(app.exec_())
