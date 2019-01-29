from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit, QInputDialog
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
import sys
from time import time
from reverseSplice import *

class WorkerSignals(QObject):
    """
    Signals class that is used for the GUI when emitting custom signals
    """

    finished = pyqtSignal()

class OutputGenerator(QRunnable):
    """

    """

    def __init__(self, fn, *args, **kwargs):
        super(OutputGenerator, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        self.fn(*self.args)
        self.signals.finished.emit()


class Example(QWidget):

    def __init__(self):
        super().__init__()

        self.initUI()
        # global variables initialised here here
        self.proteinFile = ""
        self.peptideFile = ""
        # initialise threadpool
        self.threadpool = QThreadPool()

    def initUI(self):

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.initialiseWidgets()

        self.move(300, 150)
        self.setWindowTitle('Reverse Splicer')
        self.show()

    def initialiseWidgets(self):
        self.importProtein = QPushButton('Import Protein Fasta')
        self.grid.addWidget(self.importProtein, 1, 1)
        self.importProtein.clicked.connect(self.uploadFile)

        self.importPeptide = QPushButton('Import Peptide Fasta')
        self.grid.addWidget(self.importPeptide, 2, 1)
        self.importPeptide.clicked.connect(self.uploadFile)

        self.generateOutput = QPushButton('Generate Output')
        self.generateOutput.setEnabled(False)
        self.grid.addWidget(self.generateOutput, 6, 1)
        self.generateOutput.clicked.connect(self.outputCheck)

        self.linCheckbox = QCheckBox('Linear')
        self.linCheckbox.setEnabled(False)
        self.linCheckbox.stateChanged.connect(self.enableOutput)
        self.cisCheckbox = QCheckBox('Cis')
        self.cisCheckbox.setEnabled(False)
        self.cisCheckbox.stateChanged.connect(self.enableOutput)
        self.transCheckbox = QCheckBox('Trans')
        self.transCheckbox.setEnabled(False)
        self.transCheckbox.stateChanged.connect(self.enableOutput)
        self.grid.addWidget(self.linCheckbox, 3, 1)
        self.grid.addWidget(self.cisCheckbox, 4, 1)
        self.grid.addWidget(self.transCheckbox, 5, 1)

    def uploadFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home/')
        if fname[0][-5:] == 'fasta':
            if self.sender() == self.importProtein:
                self.proteinFile = fname[0]
                QMessageBox.about(self, 'Message', 'Protein input fasta successfully uploaded!')
            else:
                self.peptideFile = fname[0]
                QMessageBox.about(self, 'Message', 'Peptide input fasta successfully uploaded!')
        if self.proteinFile == "" or self.peptideFile == "":
            self.linCheckbox.setEnabled(False)
            self.cisCheckbox.setEnabled(False)
        else:
            self.linCheckbox.setEnabled(True)
            self.cisCheckbox.setEnabled(True)

    def getOutputPath(self):

        """
        Called after generate output is clicked. Opens a window to select a file location to save the output to.
        Returns False if no path is selected, otherwise returns the selected path.
        """

        outputFile = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        if outputFile == '':
            return False
        else:
            text, ok = QInputDialog.getText(self, 'Input Dialog',
                                            'Enter your file name:')

            if ok:
                outputPath = outputFile + '/' + text
            else:
                return False
        print(outputPath)
        return outputPath

    def outputCheck(self):
        if self.proteinFile == "" or self.peptideFile == "":
            QMessageBox.about(self, 'Message', 'Please Upload a Fasta File before generating output!')
        else:

            reply = QMessageBox.question(self, 'Message', 'Are these the correct files you wish to run?' + "\n" +
                                                            "Please ensure you haven't switched the protein and peptide input." + '\n'+ '\n' +
                                                            'Protein File: ' + self.proteinFile + '\n' + '\n' +
                                                            'Peptide File: ' + self.peptideFile)
            if reply == QMessageBox.Yes:
                start = time()
                outputPath = self.getOutputPath()
                if outputPath is not False:

                    self.outputGen = OutputGenerator(self.createOutput, outputPath, self.proteinFile, self.peptideFile)
                    self.outputGen.signals.finished.connect(self.outputFinished)
                    self.threadpool.start(self.outputGen)
                    self.outputLabel = QLabel("Generating Output. Please Wait!")
                    self.grid.addWidget(self.outputLabel,4,1)
                    #generateOutputNew(outputPath, self.minPeptideLen, self.inputFile)

    def enableOutput(self):
        if self.linCheckbox.isChecked() or self.cisCheckbox.isChecked() or self.transCheckbox.isChecked():
            self.generateOutput.setEnabled(True)
        else:
            self.generateOutput.setEnabled(False)

    def createOutput(self, outputPath, proteinFile, peptideFile):
        generateOutput(outputPath, proteinFile, peptideFile)

    def outputFinished(self):
        QMessageBox.about(self, "Message", "All done!")
        self.grid.removeWidget(self.outputLabel)
        self.outputLabel.deleteLater()
        self.outputLabel = None



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())