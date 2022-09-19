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
import numpy as np
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
        self.killBtn = QPushButton('Kill Process')
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        dlgLayout.addWidget(self.runBtn)
        dlgLayout.addWidget(self.killBtn)
        dlgLayout.addWidget(self.text)
        self.runBtn.setDefault(True)
        self.killBtn.setEnabled(False)
        self.killBtn.clicked.connect(self.kill_process)
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
        self.diaBtn.clicked.connect(lambda: self.getfile(self.diaFiles))
        self.diaFileText.setToolTip(
            'Data File Requirements:\n-Required, must be populated\n-Must be of type MzXML')
        self.libFileText = QLabel('Library File:')
        self.libBtn = QPushButton('Browse')
        self.libBtn.clicked.connect(lambda: self.getfile(self.libFile))
        self.libFileText.setToolTip(
            'Library Spectra File Requirements:\n-Required, must be populated\n-Must be in MGF (.mgf) or TraML (.csv or .tsv) format')
        self.libFile = QLineEdit()
        self.outDirText = QLabel('Outfile Directory:')
        self.dirBtn = QPushButton('Browse')
        self.outDirText.setToolTip(
            'Data File Requirements:\n-Required, must be populated')
        self.outDir = QLineEdit()
        self.dirBtn.clicked.connect(
            lambda: self.getfile(self.outDir, isFile=False))
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
        self.fragMassTolText = QLabel(
            'Initial Fragment Mass Tolerance (in ppm):')
        self.fragMassTolText.setToolTip(
            'Fragment Mass Tolerance Requirements:\n-Must be blank or an integer greater than 0')
        self.fragMassTol = QLineEdit()
        self.corr = QLineEdit()
        self.corrCheckBox = QCheckBox()
        self.corrText = QLabel('Corrective Standard Deviations:')
        self.corrText.setToolTip(
            'Corrective Standard Deviations Requirements:\n-Must be blank or a float (decimal value) between or equal to 0.5 and 2.0')
        self.histCheckBox = QCheckBox()

        self.corrCheckBox.stateChanged.connect(
            lambda: self.check_grey(self.corrCheckBox, self.corr))
        self.corrCheckBox.stateChanged.connect(
            lambda: self.enable_second_box(self.corrCheckBox, self.histCheckBox))
        settingLayout.addRow(self.fragMassTolText, self.fragMassTol)
        settingLayout.addRow('Correction:', self.corrCheckBox)
        settingLayout.addRow(self.corrText, self.corr)
        settingLayout.addRow('Create Histogram:', self.histCheckBox)
        self.fragMassTol.setPlaceholderText('30')
        self.corr.setPlaceholderText('customized')
        self.corrCheckBox.setCheckState(2)
        self.libFile.setEnabled(False)
        self.outDir.setEnabled(False)
        self.histCheckBox.setCheckState(2)

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

    def delete_selected_values(self):
        for x in self.diaFiles.selectedItems():
            self.diaFiles.takeItem(self.diaFiles.row(x))

    def set_dict(self):
        tempDict = {}
        self.set_dict_parent(tempDict)
        self.set_dict_child(tempDict)

        if -1 in list(tempDict.values()):
            self.dict = None
            return
        else:
            self.dict = tempDict

    def set_dict_parent(self, tempDict):
        tempDict['diaFiles'] = self.return_dia_file_values(
            self.diaFiles, self.diaFileText, permittedTypes=['mzxml'])
        tempDict['libFile'] = self.return_file_path_value(
            self.libFile.text(), self.libFileText, permittedTypes=['mgf', 'csv', 'tsv'])
        tempDict['outDir'] = self.return_file_path_value(
            self.outDir.text(), self.outDirText)
        tempDict['fragMassTol'] = self.return_integer_above_0(
            self.fragMassTol.text(), self.fragMassTolText)

        tempDict['num_peaks'] = self.validate_positive_integer(
            self.num_peaks_for_quantification.text(), self.num_peaks_text)
        if(tempDict['num_peaks'] == 0):
            tempDict['num_peaks'] = np.inf

            # Check for common peptide/protein
        if self.do_common_peptide_quantification.isChecked():
            tempDict['do_common_peptide_quantification'] = True
        else:
            tempDict['do_common_peptide_quantification'] = False
        if self.do_common_protein_quantification.isChecked():
            tempDict['do_common_protein_quantification'] = True
        else:
            tempDict['do_common_protein_quantification'] = False

        if self.corrCheckBox.isChecked():
            tempDict['corr'] = self.return_float_between_P5_2(
                self.corr.text(), self.corrText)
        else:
            tempDict['corr'] = self.return_valid(self.corrText, False)
        if self.histCheckBox.isChecked():
            tempDict['hist'] = True
        else:
            tempDict['hist'] = False

    def set_dict_child(self, tempDict):
        pass

    def return_dia_file_values(self, filesObject, text, permittedTypes=[]):
        if filesObject.count() == 0:
            return self.return_valid(text)
        files = [filesObject.item(i).text()
                 for i in range(filesObject.count())]
        if len(permittedTypes) != 0:
            for file in files:
                if file.split('.')[-1].lower() not in permittedTypes:
                    return self.return_valid(text)
        return self.return_valid(text, files)

    def return_file_path_value(self, path, text, permittedTypes=[]):
        if len(path) == 0:
            return self.return_valid(text)
        if len(permittedTypes) != 0 and path.split('.')[-1].lower() not in permittedTypes:
            return self.return_valid(text)
        return self.return_valid(text, path)

    def return_integer_above_0(self, x, text):
        if x == '':
            return self.return_valid(text, 'default')
        try:
            v = int(x)
            if v < 1:
                return self.return_valid(text)
            return self.return_valid(text, x)
        except ValueError:
            return self.return_valid(text)

    def validate_positive_integer(self,
                                  x: str,
                                  label_to_color: QLabel = None):
        '''
        Validates that 'x' is a positive integer. Optionally color the label red/green if value is invalid/valid.
        Parameters
        ----------
        x : str
            String to validate and convert to int.
        label_to_color : QLabel
            Optional.  Label to color.

        Returns
        -------
        int
            Converted integer.
        '''
        if x == '':
            return self.return_valid(label_to_color, x)
        try:
            int_x = int(x)
            if int_x < 0:
                return self.return_valid(label_to_color)
            return self.return_valid(label_to_color, x)
        except ValueError:
            return self.return_valid(label_to_color)


    def return_float_between_P5_2(self, x, text):
        if x == '':
            return self.return_valid(text, 'default')
        try:
            v = float(x)
            if v < 0.5 or v > 2.0:
                return self.return_valid(text)
            return self.return_valid(text, x)
        except ValueError:
            return self.return_valid(text)

    def return_valid(self, text, x=-1):
        if x == -1:
            text.setStyleSheet('background-color: ' + self.INVALID_COLOR)
            return -1
        else:
            text.setStyleSheet('background-color: ' + self.VALID_COLOR)
            return x

    def message(self, s):
        self.text.appendPlainText(s)

    def start_process(self):
        self.set_dict()
        if self.dict == None:
            return
        self.killBtn.setEnabled(True)

        if self.p is None:  # No process running.
            args = []
            self.set_args_parent(args)
            self.set_args_child(args)

            self.message("Executing process")
            # Keep a reference to the QProcess (e.g. on self) while it's running.
            self.p = QProcess()
            self.p.readyReadStandardOutput.connect(self.handle_stdout)
            self.p.readyReadStandardError.connect(self.handle_stderr)
            # Clean up once complete.
            self.p.finished.connect(self.process_finished)
            self.p.start('csodiaq', args)

    def set_args_parent(self, args):
        for file in self.dict['diaFiles']:
            args += ['-f', file]
        args += ['-l', self.dict['libFile']]
        args += ['-o', self.dict['outDir']+'/']
        if self.dict['fragMassTol'] != 'default':
            args += ['-t', self.dict['fragMassTol']]
        if self.dict['corr'] == 'default':
            args += ['-c']
        elif self.dict['corr']:
            args += ['-c', self.dict['corr']]
        if self.dict['hist']:
            args += ['-hist']
        if self.dict['num_peaks']:
            args += ['--peaks', self.dict['num_peaks']]
        if self.dict['do_common_peptide_quantification']:
            args += ['--commonpeptide']
        if self.dict['do_common_protein_quantification']:
            args += ['--commonprotein']

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

    def kill_process(self):
        self.p.kill()

    def process_finished(self):
        self.message("Process finished.")
        self.killBtn.setEnabled(False)
        self.p = None


class IdWindow(csodiaqWindow):
    def __init__(self):
        super().__init__()

    def set_settings_child(self, settingLayout):
        self.protTarg = QLineEdit()
        self.protTargText = QLabel('Number of Target Peptides per Protein: ')
        self.protTargText.setToolTip(
            'Target Peptide Requirements:\n-Must be blank or an integer greater than 0')
        self.protCheckBox = QCheckBox()
        self.query = QLineEdit()
        self.queryText = QLabel('Maximum Number of Query Spectra to Pool: ')
        self.queryText.setToolTip(
            'Query Spectra Pooling Requirements:\n-Must be blank or an integer greater than 0')
        self.heavyCheckBox = QCheckBox()

        self.num_peaks_for_quantification = QLineEdit()
        self.num_peaks_for_quantification.setPlaceholderText('0')
        self.num_peaks_text = QLabel('# peaks for quantification:')

        self.do_common_peptide_quantification = QCheckBox()
        self.do_common_peptide_quantification_text = QLabel("Get shared peptides across input files")
        self.do_common_protein_quantification = QCheckBox()
        self.do_common_protein_quantification_text = QLabel("Get shared protein across input files")

        settingLayout.addRow(self.protTargText, self.protTarg)
        settingLayout.addRow('Protein Inference:', self.protCheckBox)
        settingLayout.addRow(self.queryText, self.query)
        settingLayout.addRow(
            'Permit Heavy Targets in Re-Analysis File:', self.heavyCheckBox)

        # Number of fragments/peaks
        settingLayout.addRow(self.num_peaks_text, self.num_peaks_for_quantification)

        # Checkboxes for post-analysis quantification
        settingLayout.addRow(self.do_common_peptide_quantification_text, self.do_common_peptide_quantification)
        settingLayout.addRow(self.do_common_protein_quantification_text, self.do_common_protein_quantification)
        self.do_common_protein_quantification.setCheckState(2)
        self.do_common_peptide_quantification.setCheckState(2)

        self.protCheckBox.stateChanged.connect(lambda: self.check_grey(
            self.protCheckBox, self.protTarg, filledText='1'))
        self.protTarg.setPlaceholderText('untargeted peptide analysis')
        self.protTarg.setEnabled(False)
        self.query.setPlaceholderText('pool all matching query spectra')
        self.heavyCheckBox.setCheckState(2)

    def set_dict_child(self, tempDict):
        if self.protCheckBox.isChecked():
            tempDict['protTarg'] = self.return_integer_above_0(
                self.protTarg.text(), self.protTargText)
        else:
            tempDict['protTarg'] = self.return_valid(self.protTargText, False)
        tempDict['query'] = self.return_integer_above_0(
            self.query.text(), self.queryText)
        tempDict['heavy'] = self.heavyCheckBox.isChecked()

    def set_args_child(self, args):
        args.insert(0, 'id')
        if self.dict['protTarg']:
            args += ['-p', self.dict['protTarg']]
        if self.dict['query'] != 'default':
            args += ['-q', self.dict['query']]
        if self.dict['heavy']:
            args += ['-heavy']


class quantWindow(csodiaqWindow):
    def __init__(self):
        super().__init__()

    def set_files_child(self, fileLayout):
        self.idFileText = QLabel('CsoDIAq ID Output File:')
        self.idFileText.setToolTip(
            'CsoDIAq ID File Requirements:\n-Required, must be populated')
        self.idFile = QLineEdit()
        self.idBtn = QPushButton('Browse')
        self.idBtn.clicked.connect(lambda: self.getfile(self.idFile))
        fileLayout.addRow(self.idFileText, self.idBtn)
        fileLayout.addRow(self.idFile)

    def set_settings_child(self, settingLayout):
        self.libPeaks = QLineEdit()
        self.libPeaks.setPlaceholderText('all spectra peaks')
        self.libPeaksText = QLabel('Number of Max Peaks per Library Spectra: ')
        self.libPeaksText.setToolTip(
            'Library Peaks Requirements:\n-Must be blank or an integer greater than 0')
        self.minMatch = QLineEdit()
        self.minMatch.setPlaceholderText('default: 1 of 3 most intense peaks')
        self.minMatchText = QLabel('Number of Min Peak Matches Required: ')
        self.minMatchText.setToolTip(
            'Minimum Peak Match Requirements:\n-Must be blank or an integer greater than 0')
        self.ratioType = QComboBox()
        self.ratioType.addItem('median')
        self.ratioType.addItem('mean')
        settingLayout.addRow(self.libPeaksText, self.libPeaks)
        settingLayout.addRow(self.minMatchText, self.minMatch)
        settingLayout.addRow(QLabel('Ratio Selection Method:'), self.ratioType)

    def set_dict_child(self, tempDict):
        tempDict['idFile'] = self.return_file_path_value(
            self.idFile.text(), self.idFileText, permittedTypes=['csv'])
        tempDict['libPeaks'] = self.return_integer_above_0(
            self.libPeaks.text(), self.libPeaksText)
        tempDict['minMatch'] = self.return_integer_above_0(
            self.minMatch.text(), self.minMatchText)
        tempDict['ratioType'] = self.ratioType.currentText()

    def set_args_child(self, args):
        args.insert(0, 'quant')
        args += ['-i', self.dict['idFile']]
        if self.dict['libPeaks'] != 'default':
            args += ['-p', self.dict['libPeaks']]
        if self.dict['minMatch'] != 'default':
            args += ['-m', self.dict['minMatch']]
        args += ['-r', self.dict['ratioType']]


class MyTableWidget(QWidget):

    def __init__(self, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)

        self.tabs = QTabWidget()
        self.tab1 = IdWindow()
        self.tab2 = quantWindow()
        self.tabs.addTab(self.tab1, 'Peptide/Protein Identification')
        self.tabs.addTab(self.tab2, 'SILAC Quantification')

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
        if self.table_widget.tab1.p:
            self.table_widget.tab1.p.kill()


def main():
    app = QApplication(sys.argv)
    view = MainWindow()
    view.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
