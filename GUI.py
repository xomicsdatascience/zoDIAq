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
from PyQt5.QtWidgets import QTabWidget
from PyQt5.QtWidgets import QListWidget
from PyQt5.QtCore import Qt

__version__ = '0.1'
__author__ = 'Caleb Webster Cranney'


class IdWindow(QWidget):
    """Dialog."""
    def __init__(self, parent=None):
        """Initializer."""
        super().__init__(parent)

        dlgLayout = QVBoxLayout()

        dlgLayout.addWidget(QLabel('File Input:'))
        self.btn1 = QPushButton('Browse')
        self.diaFiles = QListWidget()
        self.delBtn = QPushButton('Delete')
        self.delBtn.clicked.connect(self.delete_selected_values)
        self.clrBtn = QPushButton('Clear')
        self.clrBtn.clicked.connect(self.diaFiles.clear)
        self.btn1.clicked.connect(lambda:self.getfile(self.diaFiles))
        self.btn2 = QPushButton('Browse')
        self.libFile = QLineEdit()
        self.btn2.clicked.connect(lambda:self.getfile(self.libFile))
        self.btn3 = QPushButton('Browse')
        self.outDir = QLineEdit()
        self.btn3.clicked.connect(lambda:self.getfile(self.outDir, isFile=False))
        fileLayout = QFormLayout()
        fileLayout.addRow('DIA Data File:', self.btn1)
        fileLayout.addRow(self.diaFiles)
        fileLayout.addRow(self.delBtn, self.clrBtn)
        fileLayout.addRow('Library File:', self.btn2)
        fileLayout.addRow(self.libFile)
        fileLayout.addRow('Outfile Directory:', self.btn3)
        fileLayout.addRow(self.outDir)
        dlgLayout.addLayout(fileLayout)

        dlgLayout.addWidget(QLabel('Settings:'))
        self.fragMassTol = QLineEdit()
        self.corrStDev = QLineEdit()
        self.histCheckBox = QCheckBox()
        self.histFile = QLineEdit()
        self.histCheckBox.stateChanged.connect(lambda:self.check_grey(self.histCheckBox, self.histFile))
        self.protCheckBox = QCheckBox()
        self.protTarg = QLineEdit()
        self.protCheckBox.stateChanged.connect(lambda:self.check_grey(self.protCheckBox, self.protTarg,'1'))
        settingLayout = QFormLayout()
        settingLayout.addRow('Initial Fragment Mass Tolerance (in ppm):', self.fragMassTol)
        settingLayout.addRow('Corrective Standard Deviations:', self.corrStDev)
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        settingLayout.addRow('Histogram Outfile: ', self.histFile)
        settingLayout.addRow('Protein Inference:', self.protCheckBox)
        settingLayout.addRow('Number of Target Peptides per Protein: ', self.protTarg)
        dlgLayout.addLayout(settingLayout)
        self.fragMassTol.setPlaceholderText('20')
        self.corrStDev.setPlaceholderText('customized')
        self.histFile.setPlaceholderText('no histogram')
        self.histFile.setEnabled(False)
        self.protTarg.setPlaceholderText('peptides, no protein inference')
        self.protTarg.setEnabled(False)
        self.libFile.setEnabled(False)
        self.outDir.setEnabled(False)

        self.runBtn = QPushButton('Run')
        dlgLayout.addWidget(self.runBtn)
        self.runBtn.setDefault(True)
#        self.runBtns.button(QDialogButtonBox.Ok).setText('Run')

        self.setLayout(dlgLayout)

    def check_grey(self, checkBox, lineEdit, filledText=''):
        if checkBox.isChecked():
            lineEdit.setText(filledText)
            lineEdit.setEnabled(True)
        else:
            lineEdit.clear()
            lineEdit.setEnabled(False)

    def getfile(self, text, isFile=True):
        dialog = QFileDialog()
        if type(text) is QListWidget and isFile: fname = QFileDialog.getOpenFileNames(self, 'Open File', 'c:\\'); text.addItems([x for x in fname[0]])
        elif type(text) is QLineEdit and isFile: fname = QFileDialog.getOpenFileName(self, 'Open File', 'c:\\'); text.setText(fname[0])
        else: fname = QFileDialog.getExistingDirectory(self, 'Open Directory', 'c:\\'); text.setText(fname)

    def delete_selected_values(self):
        for x in self.diaFiles.selectedItems():
            self.diaFiles.takeItem(self.diaFiles.row(x))

    def return_values(self):
        diaFiles = []
        for i in range(self.diaFiles.count()): diaFiles.append(self.diaFiles.item(i).text())
        libFile = self.libFile.text()
        outDir = self.outDir.text()
        fragMassTol = self.fragMassTol.text()
        corrStDev = self.corrStDev.text()
        histFile = self.histFile.text()
        protTarg = self.protTarg.text()
        return [diaFiles, libFile, outDir, fragMassTol, corrStDev, histFile, protTarg]

class quantWindow(QWidget):
    """Dialog."""
    def __init__(self, parent=None):
        """Initializer."""
        super().__init__(parent)
        dlgLayout = QVBoxLayout()

        formLayout = QFormLayout()
        formLayout.addRow('Outfile Directory and Header:', QLineEdit())
        formLayout.addRow('Library Spectrum Path:', QLineEdit())
        formLayout.addRow('FDR Output:', QLineEdit())
        formLayout.addRow('DISPA Targetted Re-Analysis Directory:', QLineEdit())
        dlgLayout.addLayout(formLayout)
        self.runBtn = QPushButton('Run')
        dlgLayout.addWidget(self.runBtn)
        self.runBtn.setDefault(True)

        self.setLayout(dlgLayout)


    def return_values(self):
        return False

class MyTableWidget(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab1 = IdWindow()
        self.tab2 = quantWindow()
        self.tabs.addTab(self.tab1,'Peptide/Protein Identification')
        self.tabs.addTab(self.tab2,'SILAC Quantification')

        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)

class Controller:

    def __init__(self, view):
        """Controller initializer."""
        self._view = view
        # Connect signals and slots
        self._connectSignals()


    def _runIdentification(self):
        ret = self._view.table_widget.tab1.return_values()
        for x in ret: print(x)

    def _runQuantification(self):
        ret = self._view.table_widget.tab1.return_values()
        print(ret)

    def _connectSignals(self):
        self._view.table_widget.tab1.runBtn.clicked.connect(lambda: self._runIdentification())
        self._view.table_widget.tab2.runBtn.clicked.connect(lambda: self._runQuantification())


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    view = MainWindow()
    view.show()

    Controller(view=view)

    sys.exit(app.exec_())
