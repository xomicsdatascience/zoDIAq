# Filename: dialog.py

"""Dialog-Style application."""

import sys
import os

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

import csodiaq_base_functions as cbf
import csodiaq_menu_functions as menu

#!/usr/bin/env python3

import sys

from PyQt5.QtCore import pyqtSignal, pyqtSlot, QProcess, QTextCodec
from PyQt5.QtGui import QTextCursor
from PyQt5.QtWidgets import QApplication, QPlainTextEdit


__version__ = '0.1'
__author__ = 'Caleb Webster Cranney'


class IdWindow(QWidget):
    """Dialog."""
    def __init__(self, parent=None):
        """Initializer."""
        super().__init__(parent)

        self.p = None


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
        self.protCheckBox = QCheckBox()
        self.protTarg = QLineEdit()
        self.protCheckBox.stateChanged.connect(lambda:self.check_grey(self.protCheckBox, self.protTarg,'1'))
        settingLayout = QFormLayout()
        settingLayout.addRow('Initial Fragment Mass Tolerance (in ppm):', self.fragMassTol)
        settingLayout.addRow('Corrective Standard Deviations:', self.corrStDev)
        settingLayout.addRow('Protein Inference:', self.protCheckBox)
        settingLayout.addRow('Number of Target Peptides per Protein: ', self.protTarg)
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        dlgLayout.addLayout(settingLayout)
        self.fragMassTol.setPlaceholderText('20')
        self.corrStDev.setPlaceholderText('customized')
        self.protTarg.setPlaceholderText('peptides, no protein inference')
        self.protTarg.setEnabled(False)
        self.libFile.setEnabled(False)
        self.outDir.setEnabled(False)

        self.runBtn = QPushButton('Execute')
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        dlgLayout.addWidget(self.runBtn)
        dlgLayout.addWidget(self.text)
        self.runBtn.setDefault(True)
#        self.runBtns.button(QDialogButtonBox.Ok).setText('Run')

        self.setLayout(dlgLayout)

    def check_grey(self, checkBox, lineEdit, filledText):
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
        hist = self.histCheckBox.isChecked()
        protTarg = self.protTarg.text()
        #return False
        return {'diaFiles':diaFiles,
                'libFile':libFile,
                'outDir': outDir,
                'fragMassTol':fragMassTol,
                'corrStDev':corrStDev,
                'hist':hist,
                'protTarg':protTarg}

    def message(self, s):
        self.text.appendPlainText(s)

    def start_process(self):
        if self.p is None:  # No process running.
            self.message("Executing process")
            self.p = QProcess()  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.stateChanged.connect(self.handle_state)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.
            self.p.start("python3", ['dummy_script.py'])

    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

    def handle_state(self, state):
        states = {
            QProcess.NotRunning: 'Not running',
            QProcess.Starting: 'Starting',
            QProcess.Running: 'Running',
        }
        state_name = states[state]
        self.message(f"State changed: {state_name}")

    def process_finished(self):
        self.message("Process finished.")
        self.p = None

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
        self.runBtn = QPushButton('Execute')
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
        #guiValues = self._view.table_widget.tab1.return_values()
        self._view.table_widget.tab1.start_process()
        '''
        #if not guiValues: return
        guiValues = {'diaFiles': ['/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/20190405_MCF7_FAIMS_18_2.mzXML', '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/20190411_DI2A_1to16_n1b.mzXML'], 'libFile': '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv', 'outDir': '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/GUIOutput', 'fragMassTol': '10', 'corrStDev': '1', 'hist': True, 'protTarg': '1'}


        lib = cbf.library_file_to_dict(guiValues['libFile'])
        for i in range(len(guiValues['diaFiles'])):
            expFile = guiValues['diaFiles'][i]
            outFileHeader = expFile.split('/')[-1].split('.')[0]
            outFile = guiValues['outDir']+'/CsoDIAq-file' + str(i) + '_' + outFileHeader + '.csv'
            print(outFile)
            menu.write_csodiaq_output(lib, expFile, outFile)
            menu.write_ppm_offset_tolerance(outFile, hist=guiValues['hist'])

            menu.write_csodiaq_output(lib, expFile, outFile, corrected=True)
            menu.write_csodiaq_fdr_outputs(outFile, corrected=True)
            menu.write_DISPA_targeted_reanalysis_files(outFile, proteins = int(guiValues['protTarg']))
        '''

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
    print('works?')
    view = MainWindow()
    view.show()

    Controller(view=view)

    sys.exit(app.exec_())
