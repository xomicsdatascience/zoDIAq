# Filename: dialog.py

"""Dialog-Style application."""

import sys
import os
import csv

from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QAction
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
        self.dict = None
        self.INVALID = -1
        self.VALID_COLOR = 'lightgreen'
        self.INVALID_COLOR = 'rgb(255,99,71)'


        dlgLayout = QVBoxLayout()

        dlgLayout.addWidget(QLabel('File Input:'))
        self.diaFiles = QListWidget()
        self.delBtn = QPushButton('Delete')
        self.delBtn.clicked.connect(self.delete_selected_values)
        self.clrBtn = QPushButton('Clear')
        self.clrBtn.clicked.connect(self.diaFiles.clear)
        self.diaFileText = QLabel('*DIA Data File:')
        self.btn1 = QPushButton('Browse')
        self.btn1.clicked.connect(lambda:self.getfile(self.diaFiles))
        self.libFileText = QLabel('*Library File:')
        self.btn2 = QPushButton('Browse')
        self.btn2.clicked.connect(lambda:self.getfile(self.libFile))
        self.libFile = QLineEdit()
        self.outDirText = QLabel('*Outfile Directory:')
        self.btn3 = QPushButton('Browse')
        self.outDir = QLineEdit()
        self.btn3.clicked.connect(lambda:self.getfile(self.outDir, isFile=False))
        fileLayout = QFormLayout()
        fileLayout.addRow(self.diaFileText, self.btn1)
        fileLayout.addRow(self.diaFiles)
        fileLayout.addRow(self.delBtn, self.clrBtn)
        fileLayout.addRow(self.libFileText, self.btn2)
        fileLayout.addRow(self.libFile)
        fileLayout.addRow(self.outDirText, self.btn3)
        fileLayout.addRow(self.outDir)
        dlgLayout.addLayout(fileLayout)

        dlgLayout.addWidget(QLabel('Settings:'))
        self.fragMassTolText = QLabel('Initial Fragment Mass Tolerance (in ppm):')
        self.fragMassTol = QLineEdit()
        self.corr = QLineEdit()
        self.histCheckBox = QCheckBox()
        self.protCheckBox = QCheckBox()
        self.corrCheckBox = QCheckBox()
        self.corrText = QLabel('Corrective Standard Deviations:')
        self.protTarg = QLineEdit()
        self.protTargText = QLabel('Number of Target Peptides per Protein: ')
        self.protCheckBox.stateChanged.connect(lambda:self.check_grey(self.protCheckBox, self.protTarg,filledText='1'))
        self.corrCheckBox.stateChanged.connect(lambda:self.check_grey(self.corrCheckBox, self.corr))
        self.corrCheckBox.stateChanged.connect(lambda:self.enable_second_box(self.corrCheckBox, self.histCheckBox))
        settingLayout = QFormLayout()
        settingLayout.addRow(self.fragMassTolText, self.fragMassTol)
        settingLayout.addRow('Correction:', self.corrCheckBox)
        settingLayout.addRow(self.corrText, self.corr)
        settingLayout.addRow('Protein Inference:', self.protCheckBox)
        settingLayout.addRow(self.protTargText, self.protTarg)
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        dlgLayout.addLayout(settingLayout)
        self.fragMassTol.setPlaceholderText('20')
        self.corr.setPlaceholderText('customized')
        self.corr.setEnabled(False)
        self.protTarg.setPlaceholderText('untargeted peptide analysis')
        self.protTarg.setEnabled(False)
        self.libFile.setEnabled(False)
        self.outDir.setEnabled(False)
        self.histCheckBox.setEnabled(False)


        self.runBtn = QPushButton('Execute')
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        dlgLayout.addWidget(self.runBtn)
        dlgLayout.addWidget(self.text)
        self.runBtn.setDefault(True)
        self.runBtn.clicked.connect(self.start_process)
#        self.runBtns.button(QDialogButtonBox.Ok).setText('Run')

        self.setLayout(dlgLayout)

    def check_grey(self, checkBox, lineEdit, filledText=''):
        if checkBox.isChecked():
            lineEdit.setText(filledText)
            lineEdit.setEnabled(True)
        else:
            lineEdit.clear()
            lineEdit.setEnabled(False)

    def enable_second_box(self, checkBox1, checkBox2):
        if checkBox1.isChecked():
            checkBox2.setEnabled(True)
        else:
            checkBox2.setCheckState(False)
            checkBox2.setEnabled(False)

    def getfile(self, text, isFile=True):
        dialog = QFileDialog()
        if type(text) is QListWidget and isFile: fname = QFileDialog.getOpenFileNames(self, 'Open File', 'c:\\'); text.addItems([x for x in fname[0]])
        elif type(text) is QLineEdit and isFile: fname = QFileDialog.getOpenFileName(self, 'Open File', 'c:\\'); text.setText(fname[0])
        else: fname = QFileDialog.getExistingDirectory(self, 'Open Directory', 'c:\\'); text.setText(fname)

    def delete_selected_values(self):
        for x in self.diaFiles.selectedItems():
            self.diaFiles.takeItem(self.diaFiles.row(x))

    def set_dict(self):
        tempDict = {}
        #print(self.libFile.text())
        tempDict['diaFiles'] = self.return_dia_file_values(self.diaFiles, self.diaFileText, permittedTypes=['mzxml'])
        tempDict['libFile'] = self.return_file_path_value(self.libFile.text(), self.libFileText, permittedTypes=['mgf','csv','tsv'])
        tempDict['outDir'] = self.return_file_path_value(self.outDir.text(), self.outDirText)
        tempDict['fragMassTol'] = self.return_integer_above_0(self.fragMassTol.text(), self.fragMassTolText)
        if self.corrCheckBox.isChecked(): tempDict['corr'] = self.return_float_between_P5_2(self.corr.text(), self.corrText)
        else: tempDict['corr'] = self.return_valid(self.corrText, False)
        if self.histCheckBox.isChecked(): tempDict['hist'] = True
        else: tempDict['hist'] = False
        if self.protCheckBox.isChecked(): tempDict['protTarg'] = self.return_integer_above_0(self.protTarg.text(), self.protTargText)
        else: tempDict['protTarg'] = self.return_valid(self.protTargText, False)


        if -1 in list(tempDict.values()): self.dict = None; return
        else: self.dict = tempDict

    def return_dia_file_values(self, filesObject, text, permittedTypes=[]):
        if filesObject.count() == 0:
            return self.return_valid(text)
        files = [filesObject.item(i).text() for i in range(filesObject.count())]
        if len(permittedTypes) != 0:
            for file in files:
                if file.split('.')[-1].lower() not in permittedTypes: return self.return_valid(text)
        return self.return_valid(text, files)

    def return_file_path_value(self, path, text, permittedTypes=[]):
        if len(path)==0:
            return self.return_valid(text)
        if len(permittedTypes)!=0 and path.split('.')[-1].lower() not in permittedTypes: return self.return_valid(text)
        return self.return_valid(text, path)

    def return_integer_above_0(self, x, text):
        if x=='': return self.return_valid(text, 'default')
        try:
            v = int(x)
            if v < 1:
                return self.return_valid(text)
            return self.return_valid(text, x)
        except ValueError:
            return self.return_valid(text)

    def return_float_between_P5_2(self, x, text):
        if x=='': return self.return_valid(text, 'default')
        try:
            v = float(x)
            if v < 0.5 or v > 2.0:
                return self.return_valid(text)
            return self.return_valid(text, x)
        except ValueError:
            return self.return_valid(text)


    def return_valid(self, text, x=-1):
        if x==-1:
            text.setStyleSheet('background-color: ' + self.INVALID_COLOR)
            return -1
        else:
            text.setStyleSheet('background-color: ' + self.VALID_COLOR)
            return x

    def message(self, s):
        self.text.appendPlainText(s)

    def start_process(self):
        self.set_dict()
        if self.dict==None: return
        '''
        self.dict = {'diaFiles': ['/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/20190411_DI2A_1to16_n1b.mzXML'],
        'libFile': '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv',
        'outDir': '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/GUIOutput',
        'fragMassTol': '50',
        'corr': '1',
        'hist': True,
        'protTarg': '1'}
        '''

        if self.p is None:  # No process running.

            args = ['id']
            for file in self.dict['diaFiles']: args += ['-f', file]

            args += ['-l', self.dict['libFile']]
            args += ['-o', self.dict['outDir']+'/']
            if self.dict['fragMassTol']!= 'default': args += ['-m', self.dict['fragMassTol']]
            if self.dict['corr'] == 'default': args += ['-c']
            elif self.dict['corr']: args += ['-c', self.dict['corr']]
            if self.dict['hist']: args += ['-hist']
            if self.dict['protTarg']: args += ['-p', self.dict['protTarg']]

            '''
            #print(args)
            dummyFile = '/Users/calebcranney/Desktop/a4551d2b-9e7b-4f90-99ca-dfd0e711a922.jpg'
            dummyMGF = '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/human.faims.fixed.decoy.mgf'
            dummyCSV = '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/human_31peaks_noloss_400to2000_pt2mz_subset.csv'

            args = [ 'id',
                '-f', '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/20190411_DI2A_1to16_n1b.mzXML',
                '-l', '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/Input/human_31peaks_noloss_400to2000_pt2mz.tsv',
                '-o', '/Users/calebcranney/Desktop/Meyer Lab Project/csoDIAq/Data/GUIOutput/',
                '-m', '10',
                '-c', '1',
                '-hist',
                '-p', '1'
            ]
            '''
            self.message("Executing process")
            self.p = QProcess()  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.

            self.p.start("python3", ['csodiaq.py'] + args)

    def handle_stderr(self):
        data = self.p.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        self.message(stderr)

    def handle_stdout(self):
        data = self.p.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.message(stdout)

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


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.table_widget = MyTableWidget(self)
        self.setCentralWidget(self.table_widget)
        quit = QAction("Quit", self)
        quit.triggered.connect(self.closeEvent)

    def closeEvent(self, event):
        if self.table_widget.tab1.p: self.table_widget.tab1.p.kill()

def main():
    app = QApplication(sys.argv)
    view = MainWindow()
    view.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
