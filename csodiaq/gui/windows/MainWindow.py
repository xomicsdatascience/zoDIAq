from PyQt5.QtWidgets import (
    QMainWindow,
)
from PyQt5.QtCore import Qt
from csodiaq.gui.windows.TableWindow import TableWindow


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("CsoDIAq")
        self.tableWindow = TableWindow(self)
        self.setCentralWidget(self.tableWindow)
