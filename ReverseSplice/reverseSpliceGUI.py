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
    Signals class that is used to call function from within a QRunnable task.
    """

    finished = pyqtSignal()

class OutputGenerator(QRunnable):
    """
    Used in self.returnPath() to run self.createOutput() in a separate thread to that which handles the GUI.
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
        """
        This function is run when the thread is started, and simply runs the input function self.createOutput() before
        emitting the finished signal which is linked to self.outputFinished().
        :return:
        """
        self.fn(*self.args)
        self.signals.finished.emit()


class Example(QWidget):
    """
    Example is the main window widget.
    """

    def __init__(self):
        """
        Called when the Example class is created. This function simply calls self.initUI() before initialising the
        following global variables.
        """
        super().__init__()

        self.initUI()
        # global variables initialised here here
        self.proteinFile = ""
        self.peptideFile = ""
        # initialise threadpool
        self.threadpool = QThreadPool()

        self.outputPath = None

    def initUI(self):
        """
        Called from the initialiser function this function creates the grid layout for the window, calls the
        function to create and add all the widgets to the window and displays the window in the centre of the display.
        :return:
        """

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.initialiseWidgets()

        self.move(300, 150)
        self.setWindowTitle('Reverse Splicer')
        self.show()

    def closeEvent(self, event):
        """
        Called automatically on click of the exit button. This function is configured to test if the system is being
        run on windows or mac, and uses the appropriate command to close the window and kill all processes which were
        generated by the program.

        :param event: details that close event is taking place.
        :return:
        """
        print('closed')
        # windows close command
        if platform.system() == 'Windows':
            os.system('taskkill /f /fi "WINDOWTITLE eq Reverse Splicer" /t')
        # mac close command
        else:
            os.system("ps aux |grep reverseSpliceGUI | grep -v 'pattern_of_process_you_dont_want_to_kill' | awk '{print $2}' |xargs kill")

    def initialiseWidgets(self):
        """
        Called from self.initUI(), this function creates all the required widgets and adds them to the main window's
        grid layout.
        :return:
        """
        self.importProtein = QPushButton('Import Protein Fasta')
        self.grid.addWidget(self.importProtein, 1, 1)
        self.importProtein.clicked.connect(self.uploadFile)

        self.importPeptide = QPushButton('Import Peptide Fasta')
        self.grid.addWidget(self.importPeptide, 2, 1)
        self.importPeptide.clicked.connect(self.uploadFile)

        self.generateOutput = QPushButton('Generate Output')
        self.generateOutput.setEnabled(False)
        self.grid.addWidget(self.generateOutput, 11, 1)
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
        self.overlapCheckbox = QCheckBox('Cis Overlap Off')
        self.overlapCheckbox.setEnabled(False)
        self.minTransLen = QComboBox()
        self.minTransLabel = QLabel("Min Trans Cleavage Length:")
        for i in range(2,9):
            self.minTransLen.addItem(str(i))
        self.minTransLen.setEnabled(False)
        self.minTransLabel.setEnabled(False)
        self.maxDist = QComboBox()
        self.maxDistLabel = QLabel("Max Distance Cis Cleavages: ")
        self.maxDist.addItem('None')
        for i in range(1,26):
            self.maxDist.addItem(str(i))
        self.maxDist.setEnabled(False)
        self.maxDistLabel.setEnabled(False)
        self.grid.addWidget(self.linCheckbox, 3, 1)
        self.grid.addWidget(self.cisCheckbox, 4, 1)
        self.grid.addWidget(self.transCheckbox, 5, 1)
        self.grid.addWidget(self.overlapCheckbox, 6,1)
        self.grid.addWidget(self.maxDistLabel, 7, 1)
        self.grid.addWidget(self.maxDist, 8, 1)
        self.grid.addWidget(self.minTransLen, 10, 1)
        self.grid.addWidget(self.minTransLabel, 9, 1)


    def uploadFile(self):
        """
        Called on click of the self.importProtein or self.importPeptide buttons, this function opens a window for the
        user to select a file to upload to the program. The function checks that a Fasta file has been uploaded, and
        then checks that the sequences in the file resemble that of a protein or peptide file. If the files contain
        sequence of unexpected length, a message box is presented informing the user to check their input.
        Once a successful input is selected, the path is stored in self.peptideFile or self.proteinFile depending on
        the nature of the input. Lastly, the function checks if both required files have been uploaded, and if so
        enables the remaining user inputs.
        :return:
        """
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
        Called after the self.generateOutput button is clicked and the user confirms their input. This function opens a
        window to select a file location to save the output to, and if valid path is selected opens a window to input
        the file name.
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
        Called by self.getOutputPath(), this function initialises and shows the file name pop-up box.
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
        Called by self.filePathDialog, this function is called every time the file name lineEdit is updated. It takes
        the param input, which is the text in the lineEdit, and checks if it is a valid file name. It then updates
        the validity label and enables/disables the confirm name button accordingly.
        :param input: the text in the line-edit.
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
        Called when the create output button, self.button, in the file name dialog is clicked. It takes self.outputPath
        and adds the name input by the user. It lastly create the QRunnable instance required to run self.createOutput()
        in a new thread in the threadpool, starts the new thread and closes the output name dialog box.
        :return:
        """
        # create the base output file name which will be used to create the specific names for lin/cis/tras
        outputFile = self.outputPath + '/' + self.fileName.text()
        print(outputFile)
        self.outputGen = OutputGenerator(self.createOutput, outputFile, self.proteinFile, self.peptideFile, self.linFlag,
                                         self.cisFlag, self.transFlag, self.overlapFlag, int(self.minTransLength),
                                         int(self.maxDistance))
        self.outputGen.signals.finished.connect(self.outputFinished)
        self.threadpool.start(self.outputGen)
        self.outputLabel = QLabel("Generating Output. Please Wait!")
        self.grid.addWidget(self.outputLabel, 12, 1)
        # close the output name box.
        self.outputNameBox.close()

    def outputCheck(self):
        """
        Called on click of the self.generateOutput button, this function checks that valid peptide and protein files
        have been input by the user, and if so extracts all user input and presents them back to the user for
        confirmation. If the user confirms the inputs, self.getOutputPath() is called.
        :return:
        """
        if self.proteinFile == "" or self.peptideFile == "":
            QMessageBox.about(self, 'Message', 'Please Upload a Fasta File before generating output!')
        else:

            # declare flags globally so they can be reached by returnPath() which now creates the computational process
            # to initiate the code.
            self.linFlag = self.linCheckbox.isChecked()
            self.cisFlag = self.cisCheckbox.isChecked()
            self.transFlag = self.transCheckbox.isChecked()
            self.overlapFlag = self.overlapCheckbox.isChecked()
            self.maxDistance = self.maxDist.currentText()
            self.minTransLength = self.minTransLen.currentText()


            reply = QMessageBox.question(self, 'Message', 'Are these the correct files you wish to run?' + "\n" +
                                                            "Please ensure you haven't switched the protein and peptide input." + '\n'+ '\n' +
                                                            'Protein File: ' + self.proteinFile + '\n' + '\n' +
                                                            'Peptide File: ' + self.peptideFile + '\n' + '\n' +
                                                            'Min Cleavage Length: ' + self.minTransLength + '\n' +
                                                            'Max Distance Cis: ' + self.maxDistance + '\n' +
                                                            'Linear: ' + str(self.linFlag) + ', Cis: ' + str(self.cisFlag) + ', Trans: ' + str(self.transFlag) + '\n' +
                                                            'Overlap Off: ' + str(self.overlapFlag))

            if reply == QMessageBox.Yes:
                self.getOutputPath()

    def enableOutput(self):
        """
        Called each time one of the splice type checkboxes changes state. This function enables the min trans length
        widgets if trans is selected, the cis overlap and max distance widgets if cis is selected and the
        self.generateOutput button if one of cis, trans or linear is selected.
        :return:
        """
        if self.transCheckbox.isChecked():
            self.minTransLen.setEnabled(True)
            self.minTransLabel.setEnabled(True)
        else:
            self.minTransLen.setEnabled(False)
            self.minTransLabel.setEnabled(False)
        if self.cisCheckbox.isChecked():
            self.overlapCheckbox.setEnabled(True)
            self.maxDist.setEnabled(True)
            self.maxDistLabel.setEnabled(True)
        else:
            self.overlapCheckbox.setEnabled(False)
            self.maxDist.setEnabled(False)
            self.maxDistLabel.setEnabled(False)
        if self.linCheckbox.isChecked() or self.cisCheckbox.isChecked() or self.transCheckbox.isChecked():
            self.generateOutput.setEnabled(True)
        else:
            self.generateOutput.setEnabled(False)


    def createOutput(self, outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag, overlapFlag, minTrans,
                     maxDistance):
        """
        Called as the worker function to QRunnable class OutputGenerator to initiate computation. This function
        merely passes the user inputs to generateOutput() in reverseSplice.py. This must be done within a thread to
        ensure that the GUI remains responsive while the output is being computed.

        :param outputPath: the path of the output file as selected by the user.
        :param proteinFile: the path of the protein sequence file.
        :param peptideFile: the path of the peptide sequence file.
        :param linFlag: True if the user wishes to return a file with information on where each peptide could have been
        generated from via linear splicing.
        :param cisFlag: True if the user wishes to return a file with information on where each peptide could have been
        generated from via cis splicing.
        :param transFlag: True if the user wishes to return a file with information on where each peptide could have been
        generated from via trans splicing.
        :param overlapFlag: True if the user has selected no overlap when running cis splicing.
        :param minTransLen: the minimum length a cleavage must be for it to be reported in the output file as an origin
        for a trans spliced peptide.
        :param maxDistance: the maximum distance two cleavages can be away from each other for them to combined in cis
        splicing. Will be 'None' if the user wants no maximum.
        :return:
        """
        generateOutput(outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag, overlapFlag, minTrans, maxDistance)

    def outputFinished(self):
        """
        Called when the finished signal is emitted from with the OutputGenerator QRunnable class. This function simply
        removes the label from the GUI the which communicates to the user that the output is being computed, and then
        displays a message detailing that the output has finished.
        :return:
        """
        QMessageBox.about(self, "Message", "All done!")
        self.grid.removeWidget(self.outputLabel)
        self.outputLabel.deleteLater()
        self.outputLabel = None



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())