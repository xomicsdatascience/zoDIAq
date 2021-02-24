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
from PyQt5.QtWidgets import QComboBox
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


class csodiaqWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.p = None
        self.dict = None
        self.INVALID = -1
        self.VALID_COLOR = 'lightgreen'
        self.INVALID_COLOR = 'rgb(255,99,71)'
        self.debug_dicts = []


        dlgLayout = QVBoxLayout()
        fileLayout = QFormLayout()
        settingLayout = QFormLayout()

        self.set_files_parent(fileLayout)
        self.set_files_child(fileLayout)
        dlgLayout.addWidget(QLabel('File Input:'))
        dlgLayout.addLayout(fileLayout)

        self.set_settings_parent(settingLayout)
        self.set_settings_child(settingLayout)
        dlgLayout.addWidget(QLabel('Settings:'))
        dlgLayout.addLayout(settingLayout)

        self.runBtn = QPushButton('Execute')
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        dlgLayout.addWidget(self.runBtn)
        dlgLayout.addWidget(self.text)
        self.runBtn.setDefault(True)
#        self.runBtn.clicked.connect(self.debug)
        self.runBtn.clicked.connect(self.start_process)

        self.setLayout(dlgLayout)

    def set_files_parent(self, fileLayout):
        self.diaFiles = QListWidget()
        self.delBtn = QPushButton('Delete')
        self.delBtn.clicked.connect(self.delete_selected_values)
        self.clrBtn = QPushButton('Clear')
        self.clrBtn.clicked.connect(self.diaFiles.clear)
        self.diaFileText = QLabel('DIA Data File:')
        self.diaBtn = QPushButton('Browse')
        self.diaBtn.clicked.connect(lambda:self.getfile(self.diaFiles))
        self.libFileText = QLabel('Library File:')
        self.libBtn = QPushButton('Browse')
        self.libBtn.clicked.connect(lambda:self.getfile(self.libFile))
        self.libFile = QLineEdit()
        self.outDirText = QLabel('Outfile Directory:')
        self.dirBtn = QPushButton('Browse')
        self.outDir = QLineEdit()
        self.dirBtn.clicked.connect(lambda:self.getfile(self.outDir, isFile=False))
        fileLayout.addRow(self.diaFileText, self.diaBtn)
        fileLayout.addRow(self.diaFiles)
        fileLayout.addRow(self.delBtn, self.clrBtn)
        fileLayout.addRow(self.libFileText, self.libBtn)
        fileLayout.addRow(self.libFile)
        fileLayout.addRow(self.outDirText, self.dirBtn)
        fileLayout.addRow(self.outDir)

    def set_files_child(self, fileLayout):
        pass

    def set_settings_parent(self, settingLayout):
        self.fragMassTolText = QLabel('Initial Fragment Mass Tolerance (in ppm):')
        self.fragMassTol = QLineEdit()
        self.corr = QLineEdit()
        self.corrCheckBox = QCheckBox()
        self.corrText = QLabel('Corrective Standard Deviations:')
        self.histCheckBox = QCheckBox()

        self.corrCheckBox.stateChanged.connect(lambda:self.check_grey(self.corrCheckBox, self.corr))
        self.corrCheckBox.stateChanged.connect(lambda:self.enable_second_box(self.corrCheckBox, self.histCheckBox))
        settingLayout.addRow(self.fragMassTolText, self.fragMassTol)
        settingLayout.addRow('Correction:', self.corrCheckBox)
        settingLayout.addRow(self.corrText, self.corr)
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        self.fragMassTol.setPlaceholderText('20')
        self.corr.setPlaceholderText('customized')
        self.corr.setEnabled(False)
        self.libFile.setEnabled(False)
        self.outDir.setEnabled(False)
        self.histCheckBox.setEnabled(False)


    def set_settings_child(self, settingLayout):
        pass

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
        self.set_dict_parent(tempDict)
        self.set_dict_child(tempDict)

        if -1 in list(tempDict.values()): self.dict = None; return
        else: self.dict = tempDict

    def set_dict_parent(self, tempDict):
        tempDict['diaFiles'] = self.return_dia_file_values(self.diaFiles, self.diaFileText, permittedTypes=['mzxml'])
        tempDict['libFile'] = self.return_file_path_value(self.libFile.text(), self.libFileText, permittedTypes=['mgf','csv','tsv'])
        tempDict['outDir'] = self.return_file_path_value(self.outDir.text(), self.outDirText)
        tempDict['fragMassTol'] = self.return_integer_above_0(self.fragMassTol.text(), self.fragMassTolText)
        if self.corrCheckBox.isChecked(): tempDict['corr'] = self.return_float_between_P5_2(self.corr.text(), self.corrText)
        else: tempDict['corr'] = self.return_valid(self.corrText, False)
        if self.histCheckBox.isChecked(): tempDict['hist'] = True
        else: tempDict['hist'] = False


    def set_dict_child(self, tempDict):
        pass

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
        #self.set_dict()
        #if self.dict==None: return
        #print(self.dict)
        self.dict = {'diaFiles': ['/Users/calebcranney/Desktop/0_DataFiles/ID1.mzXML'], 'libFile': '/Users/calebcranney/Desktop/0_DataFiles/lib_tsv.tsv', 'outDir': '/Users/calebcranney/Desktop/0_DataFiles/GUIOutput', 'fragMassTol': 'default', 'corr': 'default', 'hist': True, 'protTarg': '1'}
        #self.dict = {'diaFiles': ['/Users/calebcranney/Desktop/0_DataFiles/quant1.mzXML', '/Users/calebcranney/Desktop/0_DataFiles/quant2.mzXML'], 'libFile': '/Users/calebcranney/Desktop/0_DataFiles/lib_tsv.tsv', 'outDir': '/Users/calebcranney/Desktop/0_DataFiles/GUIOutput', 'fragMassTol': 'default', 'corr': 'default', 'hist': True, 'idFile': '/Users/calebcranney/Desktop/0_DataFiles/CsoDIAq-allCVs.csv', 'libPeaks': '3', 'minMatch': '1', 'ratioType': 'median'}


        if self.p is None:  # No process running.
            args = []
            self.set_args_parent(args)
            self.set_args_child(args)

            self.message("Executing process")
            self.p = QProcess()  # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            self.p.finished.connect(self.process_finished)  # Clean up once complete.

            self.p.start('csodiaq', args)

    def set_args_parent(self, args):
        for file in self.dict['diaFiles']: args += ['-f', file]
        args += ['-l', self.dict['libFile']]
        args += ['-o', self.dict['outDir']+'/']
        if self.dict['fragMassTol']!= 'default': args += ['-t', self.dict['fragMassTol']]
        if self.dict['corr'] == 'default': args += ['-c']
        elif self.dict['corr']: args += ['-c', self.dict['corr']]
        if self.dict['hist']: args += ['-hist']

    def set_args_child(self, args):
        pass

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

    def set_variables_debug(self, tempDict):

        self.__init__()

        if 'diaFiles' in tempDict and tempDict['diaFiles'] != 'default': self.diaFiles.addItems([x for x in tempDict['diaFiles']])
        if 'libFile' in tempDict and tempDict['libFile'] != 'default': self.libFile.setText(tempDict['libFile'])
        if 'outDir' in tempDict and tempDict['outDir'] != 'default': self.outDir.setText(tempDict['outDir'])
        if 'fragMassTol' in tempDict and tempDict['fragMassTol'] != 'default': self.fragMassTol.setText(tempDict['fragMassTol'])
        if 'corr' in tempDict:
            self.corrCheckBox.setCheckState(True)
            if tempDict['corr'] != 'default': self.corr.setText(tempDict['corr'])
        if 'hist' in tempDict and tempDict['hist'] == 'default': self.histCheckBox.setCheckState(True)

    def debug(self):
        for tempDict in self.debug_dicts:
            self.dict = None

            print(tempDict['title'] + ': ' + tempDict['result'])
            gui = False

            if tempDict['result']=='CMD Only (Fail)': print('GUI: N/A')
            else:
                self.set_variables_debug(tempDict)
                self.set_dict()
                if self.dict==None and tempDict['result']=='Fail': print('GUI: PASS')
                elif self.dict!=None and tempDict['result']=='Succeed': print('GUI: PASS'); gui=True
                else: print('GUI: FAIL XXXXXXXXXXXXXXXXXX')

            if gui:
                args = []
                self.set_args_parent(args)
                self.set_args_child(args)
            else:
                args = self.set_args_debug(tempDict)

            noError = [True]
            tempP = QProcess()
            tempP.start('csodiaq', args)
            tempP.waitForFinished(5000)

            if len(tempP.readAllStandardError())==0:
                if tempDict['result']=='Succeed': print('CMD: PASS')
                else: print('CMD: FAIL XXXXXXXXXXXXXXXXXX')
                tempP.kill()
            else:
                if tempDict['result']=='Fail' or tempDict['result']=='CMD Only (Fail)': print('CMD: PASS')
                else: print('CMD: FAIL XXXXXXXXXXXXXXXXXX')
            print('\n')



    def set_args_debug(self, tempDict):
        args = []
        if 'diaFiles' in tempDict:
            if tempDict['diaFiles'] == 'default': args += ['-f']
            else:
                for file in tempDict['diaFiles']: args += ['-f', file]
        if 'libFile' in tempDict:
            if tempDict['libFile'] == 'default': args += ['-l']
            else: args += ['-l', tempDict['libFile']]
        if 'outDir' in tempDict:
            if tempDict['outDir'] == 'default': args += ['-o']
            else: args += ['-o', tempDict['outDir']+'/']
        if 'fragMassTol' in tempDict:
            if tempDict['fragMassTol'] == 'default': args += ['-t']
            else: args += ['-t', tempDict['fragMassTol']]
        if 'corr' in tempDict:
            if tempDict['corr'] == 'default': args += ['-c']
            else: args += ['-c', tempDict['corr']]
        if 'hist' in tempDict:
            if tempDict['hist'] == 'default': args += ['-hist']
            else: args += ['-hist', tempDict['hist']]
        return args

class IdWindow(csodiaqWindow):
    def __init__(self):
        super().__init__()
        self.debug_baseline_dict = {
            'title': 'Baseline',
            'result': 'Succeed',
            'diaFiles': ['/Users/calebcranney/Desktop/0_DataFiles/ID1.mzXML'],
            'libFile': '/Users/calebcranney/Desktop/0_DataFiles/lib_tsv.tsv',
            'outDir': '/Users/calebcranney/Desktop/0_DataFiles/GUIOutput',
            'fragMassTol': '20',
            'corr': '1',
            'hist': 'default',
            'protTarg': '1'
        }
        self.debug_dicts = [
            {
                'title': self.debug_baseline_dict['title'],
                'result': self.debug_baseline_dict['result'],
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-f) no flag passed',
                'result': 'Fail',
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-f) no argument provided',
                'result': 'Fail',
                'diaFiles': 'default',
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-f) incorrect file type',
                'result': 'CMD Only (Fail)',
                'diaFiles': ['/Users/calebcranney/Desktop/test.txt'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-f) file does not exist',
                'result': 'CMD Only (Fail)',
                'diaFiles': ['/Users/calebcranney/Desktop/t.txt'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-l) no flag raised',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-l) no argument provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': 'default',
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-l) file type MGF',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': '/Users/calebcranney/Desktop/0_DataFiles/lib_mgf.mgf',
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-l) file type CSV',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': '/Users/calebcranney/Desktop/0_DataFiles/lib_csv.csv',
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-l) Wrong file type',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': '/Users/calebcranney/Desktop/test.txt',
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-l) File does not exist',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': '/Users/calebcranney/Desktop/t.txt',
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-o) no flag raised',
                'result': 'Fail',
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-o) no argument provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': 'default',
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-o) not a working directory',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': '/Users/calebcranney/Desktop/noDir/',
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-o) a file, not a directory',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': '/Users/calebcranney/Desktop/test.txt',
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-m) no flag raised',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-m) no argument provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': 'default',
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-m) non-integer provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': 'test',
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-m) 0 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': '0',
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-m) Value less than 0 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': '-1',
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },

            {
                'title': '(-c) no flag raised (and no histogram flag raised)',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) no flag raised (and histogram flag raised)',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) no argument provided',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': 'default',
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) non-float provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': 'test',
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) 0.5 provided',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': '0.5',
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) 2 provided',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': '2',
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) 0.4999 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': '0.4999',
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-c) 2.001 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': '2.001',
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-h) no flag raised',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-h) argument provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': 'test',
                'protTarg': self.debug_baseline_dict['protTarg']
            },
            {
                'title': '(-p) no flag raised',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
            },
            {
                'title': '(-p) no argument provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': 'default'
            },
            {
                'title': '(-p) non-integer provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': 'test'
            },
            {
                'title': '(-p) 0 provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': '0'
            },
            {
                'title': '(-p) Value less than 0 provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'protTarg': '-1'
            },
        ]

    def set_settings_child(self, settingLayout):
        self.protTarg = QLineEdit()
        self.protTargText = QLabel('Number of Target Peptides per Protein: ')
        self.protCheckBox = QCheckBox()
        settingLayout.addRow(self.protTargText, self.protTarg)
        settingLayout.addRow('Protein Inference:', self.protCheckBox)
        self.protCheckBox.stateChanged.connect(lambda:self.check_grey(self.protCheckBox, self.protTarg,filledText='1'))
        self.protTarg.setPlaceholderText('untargeted peptide analysis')
        self.protTarg.setEnabled(False)

    def set_dict_child(self, tempDict):
        if self.protCheckBox.isChecked(): tempDict['protTarg'] = self.return_integer_above_0(self.protTarg.text(), self.protTargText)
        else: tempDict['protTarg'] = self.return_valid(self.protTargText, False)

    def set_args_child(self, args):
        args.insert(0,'id')
        if self.dict['protTarg']: args += ['-p', self.dict['protTarg']]

    def set_variables_debug(self, tempDict):
        super().set_variables_debug(tempDict)
        if 'protTarg' in tempDict and tempDict['protTarg'] != 'default': self.protTarg.setText(tempDict['protTarg'])

    def set_args_debug(self, tempDict):
        args = super().set_args_debug(tempDict)
        args.insert(0,'id')
        if 'protTarg' in tempDict:
            if tempDict['protTarg'] == 'default': args += ['-p']
            else: args += ['-p', tempDict['protTarg']]
        return args

class quantWindow(csodiaqWindow):
    def __init__(self):
        super().__init__()
        self.debug_baseline_dict = {
            'title': 'Baseline',
            'result': 'Succeed',
            'diaFiles': ['/Users/calebcranney/Desktop/0_DataFiles/quant1.mzXML'],
            'libFile': '/Users/calebcranney/Desktop/0_DataFiles/lib_tsv.tsv',
            'outDir': '/Users/calebcranney/Desktop/0_DataFiles/GUIOutput',
            'fragMassTol': '20',
            'corr': '1',
            'hist': 'default',
            'idFile': '/Users/calebcranney/Desktop/0_DataFiles/CsoDIAq-allCVs.csv',
            'libPeaks': '3',
            'minMatch': '1',
            'ratioType': 'median'
        }
        self.debug_dicts = [
            {
                'title': self.debug_baseline_dict['title'],
                'result': self.debug_baseline_dict['result'],
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-i) no flag raised',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-i) no argument given',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': 'default',
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-i) wrong file type',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': '/Users/calebcranney/Desktop/test.txt',
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-i) file does not exist',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': '/Users/calebcranney/Desktop/t.txt',
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-p) no flag raised',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-p) no argument provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': 'default',
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-p) non-integer provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': 'test',
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-p) 0 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': '0',
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-p) value less than 0 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': '-1',
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-m) no flag raised',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-m) no argumnt provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': 'default',
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-m) non-integer provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': 'test',
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-m) 0 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': '0',
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-m) value less than 0 provided',
                'result': 'Fail',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': '-1',
                'ratioType': self.debug_baseline_dict['ratioType']
            },
            {
                'title': '(-r) no flag raised',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
            },
            {
                'title': '(-r) no argument provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': 'default'
            },
            {
                'title': '(-r) "median" provided',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': 'median'
            },
            {
                'title': '(-r) "mean" provided',
                'result': 'Succeed',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': 'mean'
            },
            {
                'title': '(-r) Invalid value provided',
                'result': 'CMD Only (Fail)',
                'diaFiles': self.debug_baseline_dict['diaFiles'],
                'libFile': self.debug_baseline_dict['libFile'],
                'outDir': self.debug_baseline_dict['outDir'],
                'fragMassTol': self.debug_baseline_dict['fragMassTol'],
                'corr': self.debug_baseline_dict['corr'],
                'hist': self.debug_baseline_dict['hist'],
                'idFile': self.debug_baseline_dict['idFile'],
                'libPeaks': self.debug_baseline_dict['libPeaks'],
                'minMatch': self.debug_baseline_dict['minMatch'],
                'ratioType': 'test'
            },
    ]

    def set_files_child(self, fileLayout):
        self.idFileText = QLabel('CsoDIAq ID Output File:')
        self.idFile = QLineEdit()
        self.idBtn = QPushButton('Browse')
        self.idBtn.clicked.connect(lambda:self.getfile(self.idFile))
        fileLayout.addRow(self.idFileText, self.idBtn)
        fileLayout.addRow(self.idFile)


    def set_settings_child(self, settingLayout):
        self.libPeaks = QLineEdit()
        self.libPeaks.setPlaceholderText('all spectra peaks')
        self.libPeaksText = QLabel('Number of Max Peaks per Library Spectra: ')
        self.minMatch = QLineEdit()
        self.minMatch.setPlaceholderText('default: 1 of 3 most intense peaks')
        self.minMatchText = QLabel('Number of Min Peak Matches Required: ')
        self.ratioType = QComboBox()
        self.ratioType.addItem('median')
        self.ratioType.addItem('mean')
        settingLayout.addRow(self.libPeaksText, self.libPeaks)
        settingLayout.addRow(self.minMatchText, self.minMatch)
        settingLayout.addRow(QLabel('Ratio Selection Method:'), self.ratioType)

    def set_dict_child(self, tempDict):
        tempDict['idFile'] = self.return_file_path_value(self.idFile.text(), self.idFileText, permittedTypes=['csv'])
        tempDict['libPeaks'] = self.return_integer_above_0(self.libPeaks.text(), self.libPeaksText)
        tempDict['minMatch'] = self.return_integer_above_0(self.minMatch.text(), self.minMatchText)
        tempDict['ratioType'] = self.ratioType.currentText()


    def set_args_child(self, args):
        args.insert(0,'quant')
        args += ['-i', self.dict['idFile']]
        if self.dict['libPeaks']!= 'default': args += ['-p', self.dict['libPeaks']]
        if self.dict['minMatch']!= 'default': args += ['-m', self.dict['minMatch']]
        args += ['-r', self.dict['ratioType']]

    def set_variables_debug(self, tempDict):
        super().set_variables_debug(tempDict)
        if 'idFile' in tempDict: self.idFile.setText(tempDict['idFile'])
        if 'libPeaks' in tempDict: self.libPeaks.setText(tempDict['libPeaks'])
        if 'minMatch' in tempDict: self.minMatch.setText(tempDict['minMatch'])
        if 'ratioType' in tempDict: self.ratioType.setCurrentText(tempDict['ratioType'])

    def set_args_debug(self, tempDict):
        args = super().set_args_debug(tempDict)
        args.insert(0,'quant')
        if 'idFile' in tempDict:
            if tempDict['idFile'] == 'default': args += ['-i']
            else: args += ['-i', tempDict['idFile']]
        args.insert(0,'id')
        if 'libPeaks' in tempDict:
            if tempDict['libPeaks'] == 'default': args += ['-p']
            else: args += ['-p', tempDict['libPeaks']]
        if 'minMatch' in tempDict:
            if tempDict['minMatch'] == 'default': args += ['-m']
            else: args += ['-m', tempDict['minMatch']]
        if 'ratioType' in tempDict:
            if tempDict['ratioType'] == 'default': args += ['-r']
            else: args += ['-r', tempDict['ratioType']]
        return args

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
