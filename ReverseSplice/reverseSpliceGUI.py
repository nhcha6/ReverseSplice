from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QTabWidget, QVBoxLayout, \
    QFileDialog, QGridLayout, QLabel, QComboBox, QCheckBox, QMessageBox, QDesktopWidget, \
    QProgressBar, QLineEdit, QInputDialog, QGroupBox, QFormLayout
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtCore import *
from PyQt5.QtCore import pyqtSlot
import sys
from time import time
from reverseSplice import *
import platform
import os

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

        self.outputPath = None

    def initUI(self):

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.initialiseWidgets()

        self.move(300, 150)
        self.setWindowTitle('Reverse Splicer')
        self.show()

    def closeEvent(self, event):
        print('closed')
        # windows close command
        if platform.system() == 'Windows':
            os.system('taskkill /f /fi "WINDOWTITLE eq Reverse Splicer" /t')
        # mac close command
        else:
            os.system("ps aux |grep reverseSpliceGUI | grep -v 'pattern_of_process_you_dont_want_to_kill' | awk '{print $2}' |xargs kill")

    def initialiseWidgets(self):
        self.importProtein = QPushButton('Import Protein Fasta')
        self.grid.addWidget(self.importProtein, 1, 1)
        self.importProtein.clicked.connect(self.uploadFile)

        self.importPeptide = QPushButton('Import Peptide Fasta')
        self.grid.addWidget(self.importPeptide, 2, 1)
        self.importPeptide.clicked.connect(self.uploadFile)

        self.generateOutput = QPushButton('Generate Output')
        self.generateOutput.setEnabled(False)
        self.grid.addWidget(self.generateOutput, 7, 1)
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
        self.minTransLen = QComboBox()
        for i in range(2,9):
            self.minTransLen.addItem(str(i))
        self.minTransLen.setEnabled(False)
        self.grid.addWidget(self.linCheckbox, 3, 1)
        self.grid.addWidget(self.cisCheckbox, 4, 1)
        self.grid.addWidget(self.transCheckbox, 5, 1)
        self.grid.addWidget(self.minTransLen, 6, 1)

    def uploadFile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open File', '/home/')
        if fname[0][-5:] == 'fasta':
            if self.sender() == self.importProtein:
                # check that the input contains suitable sequence lengths for a protein
                with open(fname[0], "rU") as handle:
                    counter = 0
                    totalLen = 0
                    for record in SeqIO.parse(handle, 'fasta'):
                        if counter == 100:
                            break
                        lngth = len(str(record.seq))
                        totalLen += lngth
                        counter += 1
                # Calculate aveLen of initial sequences and ask if it is the correct input file if aveLen is less than 30.
                aveLen = totalLen/counter
                if aveLen < 30:
                    response = QMessageBox.question(self, 'Message', 'The average length of the initial protein sequences is only: ' + str(round(aveLen, 1)) + '. Are you sure this is the correct protein file?')
                    if response == QMessageBox.Yes:
                        pass
                    else:
                        return
                # If it passes the aveLen check, we reach here and fasta file us uploaded.
                self.proteinFile = fname[0]
                QMessageBox.about(self, 'Message', 'Protein input fasta successfully uploaded!')
            else:
                # check that the input contains suitable sequence lengths for a peptide
                with open(fname[0], "rU") as handle:
                    counter = 0
                    totalLen = 0
                    for record in SeqIO.parse(handle, 'fasta'):
                        if counter == 100:
                            break
                        lngth = len(str(record.seq))
                        totalLen += lngth
                        counter += 1
                # Calculate aveLen of initial sequences and ask if it is the correct input file if aveLen is greater than 20.
                aveLen = totalLen / counter
                if aveLen > 20:
                    response = QMessageBox.question(self, 'Message',
                                         'The average length of the initial peptide sequences is: ' + str(round(aveLen, 1)) + '. Are you sure this is the correct peptide file?')
                    if response == QMessageBox.Yes:
                        pass
                    else:
                        return
                # If it passes the aveLen check, we reach here and fasta file us uploaded.
                self.peptideFile = fname[0]
                QMessageBox.about(self, 'Message', 'Peptide input fasta successfully uploaded!')
        if self.proteinFile == "" or self.peptideFile == "":
            self.linCheckbox.setEnabled(False)
            self.cisCheckbox.setEnabled(False)
            self.transCheckbox.setEnabled(False)
        else:
            self.linCheckbox.setEnabled(True)
            self.cisCheckbox.setEnabled(True)
            self.transCheckbox.setEnabled(True)

    def getOutputPath(self):

        """
        Called after generate output is clicked and the user confirms their input. Opens a window to select a file location
        to save the output to, and if valid opens a window to input the file name.
        """
        # opens a window to select file location.
        self.outputPath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        # if no outout path is returned, simply return to the main GUI and the user can choose to recommence the file location
        # selection process if they desire.
        if self.outputPath == '':
            return
        # else if a valid path is selected, bring up a dialog to input the file name
        else:
            self.filePathDialog()

    def filePathDialog(self):
        """
        This function initialises and shows the filing naming popup.
        """
        self.outputNameBox = QGroupBox('Output Name')
        self.outputNameLayout = QFormLayout()
        self.outputNameLayout.addRow(QLabel("Add a name for the output file."))
        self.outputNameLayout.addRow(QLabel('Banned characters: \ / : * " < > |'))
        self.fileName = QLineEdit()
        self.fileName.textChanged[str].connect(self.nameChecker)
        self.button = QPushButton("Create Output")
        self.valid = QLabel("Valid")
        self.button.clicked.connect(self.returnPath)
        self.outputNameLayout.addRow(self.fileName, self.valid)
        self.outputNameLayout.addRow(self.button)
        self.outputNameBox.setLayout(self.outputNameLayout)
        self.outputNameBox.show()

    def nameChecker(self, input):
        """
        This function is called every time the file name lineEdit is updated. It takes the param input, which is the
        text in the lineEdit, and checks if it is a valid file name.
        :param input:
        """
        # assign bannedCharacters to variables.
        bannedCharacters = set('\/:*"<>|')
        # if the input has no intersection with the banned characters it is valid. If so, update the label validity label
        # and set ensure the generate output button is concurrently enabled/disabled.
        if len(set(input).intersection(bannedCharacters)) == 0:
            self.valid.setText("Valid")
            self.button.setEnabled(True)
        else:
            self.valid.setText("Invalid")
            self.button.setEnabled(False)

    def returnPath(self):
        """
        Called when create output button in the file name dialog is clicked. It takes self.outputPath and adds the
        name input by the user. It then creates the specific names for lin/cis/trans, adds the time and adds these
        output paths to a dictionary. It lastly calls the code to create the output.
        :return:
        """
        # create the base output file name which will be used to create the specific names for lin/cis/tras
        outputFile = self.outputPath + '/' + self.fileName.text()
        print(outputFile)
        self.outputGen = OutputGenerator(self.createOutput, outputFile, self.proteinFile, self.peptideFile, self.linFlag,
                                         self.cisFlag, self.transFlag, int(self.minTransLength))
        self.outputGen.signals.finished.connect(self.outputFinished)
        self.threadpool.start(self.outputGen)
        self.outputLabel = QLabel("Generating Output. Please Wait!")
        self.grid.addWidget(self.outputLabel, 8, 1)
        # close the output name box.
        self.outputNameBox.close()

    def outputCheck(self):
        if self.proteinFile == "" or self.peptideFile == "":
            QMessageBox.about(self, 'Message', 'Please Upload a Fasta File before generating output!')
        else:

            # declare flags globally so they can be reached by returnPath() which now creates the computational process
            # to initiate the code.
            self.linFlag = self.linCheckbox.isChecked()
            self.cisFlag = self.cisCheckbox.isChecked()
            self.transFlag = self.transCheckbox.isChecked()
            self.minTransLength = self.minTransLen.currentText()

            reply = QMessageBox.question(self, 'Message', 'Are these the correct files you wish to run?' + "\n" +
                                                            "Please ensure you haven't switched the protein and peptide input." + '\n'+ '\n' +
                                                            'Protein File: ' + self.proteinFile + '\n' + '\n' +
                                                            'Peptide File: ' + self.peptideFile + '\n' + '\n' +
                                                            'Min Cleavage Length: ' + self.minTransLength + '\n' +
                                                            'Linear: ' + str(self.linFlag) + ', Cis: ' + str(self.cisFlag) + ', Trans: ' + str(self.transFlag))

            if reply == QMessageBox.Yes:
                self.getOutputPath()

    def enableOutput(self):
        if self.transCheckbox.isChecked():
            self.minTransLen.setEnabled(True)
        else:
            self.minTransLen.setEnabled(False)
        if self.linCheckbox.isChecked() or self.cisCheckbox.isChecked() or self.transCheckbox.isChecked():
            self.generateOutput.setEnabled(True)
        else:
            self.generateOutput.setEnabled(False)

    def createOutput(self, outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag, minTransFlag):
        generateOutput(outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag, minTransFlag)

    def outputFinished(self):
        QMessageBox.about(self, "Message", "All done!")
        self.grid.removeWidget(self.outputLabel)
        self.outputLabel.deleteLater()
        self.outputLabel = None



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())