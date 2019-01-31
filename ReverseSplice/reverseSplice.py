from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from datetime import datetime
import multiprocessing
from multiprocessing import Queue
import logging
import re

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

# Define maxLen
maxLen = 20
STOP = "STOP"

def protFastaToDict(protFile):
    """
    :param protFile: A Fasta file path (containing all the proteins which could serve as an origin location)
    :return protDict: A dictionary which contains a protein sequence as the key, and the protein name as the value.
    """

    protDict = {}
    with open(protFile, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seq = str(record.seq)
            name = str(record.name).split('|')[1]
            if seq in protDict.keys():
                protDict[seq].append(name)
            else:
                protDict[seq] = name
    return protDict


def generateOrigins(protDict, pepFile, outputPath, linFlag, cisFlag, transFlag):
    """
    :param protData: A dictionary containing protein sequences as the key with their origin as the value
    :param pepFile: A file with a list of peptides that you want to find the origin locations for
    :param outputPath: Where you want to store the output file
    :param linFlag: Boolean flag to determine checking for linear peptides
    :param cisFlag: Boolean flag to determine checking for cis peptides
    :param transFlag: Boolean flag to determine checking for trans peptides

    The purpose of this function is to find the aforementioned origins based on the flags passed through and then write
    the output to a Fasta file
    """

    if linFlag:
        # Find linear origins
        findLinOrigins(protDict, pepFile, outputPath)
    if cisFlag:
        # Find cis origins
        findCisOrigins(protDict, pepFile, outputPath)
    if transFlag:
        # Find trans origins > work in progress
        print('Trans')


def findLinOrigins(protDict, pepFile, outputPath):
    """
    :param protDict: A dictionary containing protein sequences as the key with their origin as the value
    :param pepFile: A file containing a list of peptides that you want to find the linear origin locations for
    :return linOriginsDict: Has the peptide as a key and a list of tuples of the form (originProtName, locations).
                            Locations store information on where the corresponding peptide could have been generated
                            from in the relevant origin protein.
    :Data structure summary: linOriginsDict[peptide] = [(proteinName, locations),(proteinName, locations)]
    """
    numWorkers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    outputPath = outputPath + '_' + 'Linear' + '-' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'

    pool = multiprocessing.Pool(processes=numWorkers, initializer=processLinInitArgs,
                                initargs=(toWriteQueue,))

    writerProcess = multiprocessing.Process(target=linearWriter, args=(toWriteQueue, outputPath, 'Linear'))
    writerProcess.start()

    # iterate through each peptide
    with open(pepFile, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            pep = str(record.seq)
            logging.info('Process started for: ' + str(pep))
            pool.apply_async(linearOrigin, args=(pep, protDict))
        pool.close()
        pool.join()
    logging.info("Pool joined")
    toWriteQueue.put(STOP)
    writerProcess.join()


def linearOrigin(pep, protDict):

    linOriginDict = {}
    # initialise the key as an empty list in the outputDict
    linOriginDict[pep] = []
    # iterate through each protSeq in the keys of protDict
    for protSeq in protDict.keys():

        # initialise the locations holder
        locations = []
        # re.finditer(substring, string) returns an iterable with the start and finish indices of all the locations
        # of the substring in the string. The iterator is empty if the susbset does not appear at all.
        # We thus iterate through all instances of the subset and append it to a list of locations.
        for x in re.finditer(pep, protSeq):
            # convert the start index and end index to a list of indexes and then append to locations
            # locations structure is a list of lists: [[1,2,3,4],[5,6,7,8]]
            locations.append([i for i in range(x.start(), x.end())])
        # if nothing is added to locations, it means that the peptide was not found in the protein, and we continue
        # iterating through proteins.
        if locations:
            linOriginDict[pep].append((protDict[protSeq], locations))
        # otherwise if we have added to locations, we append the protName/location tup to linOriginDict

    linearOrigin.toWriteQueue.put(linOriginDict)
    logging.info('Process complete for: ' + str(pep))


def linearWriter(toWriteQueue, outputPath, spliceType):
    finalLinOriginDict = {}

    while True:
        linOriginDict = toWriteQueue.get()
        if linOriginDict == STOP:
            logging.info("STOPPED writer queue")
            break

        for key, value in linOriginDict.items():
            if key in finalLinOriginDict.keys():
                finalLinOriginDict[key].append(linOriginDict[key])
            else:
                finalLinOriginDict[key] = linOriginDict[key]
    writeToFasta(finalLinOriginDict, outputPath, spliceType)


def processLinInitArgs(toWriteQueue):
    """
    Designed to initialise arguments for main process
    """

    linearOrigin.toWriteQueue = toWriteQueue


# as with findLinOrigins(), this takes the protDict and pepList and outputs a dictionary with the peptides as
# keys and a list of tuples containing origin peptide names and the splitsLocation data as the value.
# Data structure summary: cisOriginsDict[peptide] = [(proteinName, locations),(proteinName, locations)]
def findCisOrigins(protDict, pepFile, outputPath):
    """
    :param protDict: A dictionary containing protein sequences as the key with their origin as the value
    :param pepFile: A file containing a list of peptides that you want to find the linear origin locations for
    :return:
    """
    cisOriginDict = {}
    outputPath = outputPath + '_' + 'Cis' + '-' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'

    numWorkers = multiprocessing.cpu_count()
    toWriteQueue = multiprocessing.Queue()
    pool = multiprocessing.Pool(processes=numWorkers, initializer=processCisInitArgs,
                                initargs=(toWriteQueue,))
    writerProcess = multiprocessing.Process(target=linearWriter, args=(toWriteQueue, outputPath, 'Cis'))
    writerProcess.start()

    # iterate through each pep in pepList
    with open(pepFile, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            pep = str(record.seq)
            logging.info('Cis Process started for: ' + str(pep))

            pool.apply_async(cisOrigin, args=(pep, protDict))
        pool.close()
        pool.join()
    logging.info("Pool joined")
    toWriteQueue.put(STOP)
    writerProcess.join()

def cisOrigin(pep, protDict):

    cisOriginDict = {}

    # initialise that key in the dictionary
    cisOriginDict[pep] = []
    # find the splits which could be combined to create the peptide using Cis splicing.
    # cisSplits is a list of tups, where each tuple contains two complimentary splits.
    # cisSplits = [('A', 'BCD'),('AB', 'CD'),('ABC', 'D')]
    cisSplits = findCisSplits(pep)
    # iterate through each protein in protDict.keys()
    for protSeq in protDict.keys():
        # ignore that protSeq if the linear splice comes from it already.
        if pep in protSeq:
            continue
        # find the location data corresponding the current protSeq.
        locations = findCisIndexes(cisSplits, protSeq)
        # if there no pairs of splits are located in protSeq, findCisIndexes() will return an empty list.
        # If so, continue.
        if locations == []:
            continue
        # if it is possible to create the peptide using cis splicing from the current protein, add the protName
        # and locations tuple to cisOriginsDict
        cisOriginDict[pep].append((protDict[protSeq], locations))
    cisOrigin.toWriteQueue.put(cisOriginDict)
    logging.info('Cis Process complete for: ' + str(pep))


def processCisInitArgs(toWriteQueue):
    cisOrigin.toWriteQueue = toWriteQueue

# takes a peptide and returns a list of tuples, where each tuple is a possible pair of subsequences which
# could be combined to make the peptide.
def findCisSplits(pep):
    cisSplits = []
    lngth = len(pep)
    for i in range(1,lngth):
        split1 = pep[0:i]
        split2 = pep[i:lngth+1]
        cisSplits.append((split1,split2))
    return cisSplits

# returns the location data of where different pairs of cisSplits (which combine to form a given peptide)
# exist in a protein sequence. The data structure is a series of embedded lists.
# The outer list contains lists relating to each cisSplit combination which has both splits found in the protSeq.
# The list relating to each cisSplit combination contains two lists; a list of all the places a first split can
# occur and a list of all the places the second split can occur.
# Thus our final data structure looks like:
# outerList[firstSplitCombo[firstSplitList[[1,2,3,4],[5,6,7,8]], secondSplitList[[9,10,11,12],[13,14,15,16]], secondSplitCombo[etc....]]
def findCisIndexes(cisSplits, protSeq):

    totalLocations = []
    # iterate through each splitTup
    for splitTup in cisSplits:
        splitLocations = []
        # assign split1 and split 2
        split1 = splitTup[0]
        split2 = splitTup[1]
        splitLoc1 = []
        splitLoc2 = []
        # append the location of any split1 occurences in protSeq to splitLoc1
        for x in re.finditer(split1, protSeq):
            splitLoc1.append([i for i in range(x.start(), x.end())])
        # append the location of any split2 occurences in protSeq to splitLoc2
        for x in re.finditer(split2, protSeq):
            splitLoc2.append([i for i in range(x.start(), x.end())])
        # if either of the splitLocs are empty, the peptide cannot be formed using this split combo from the given
        # protSeq. Thus continue if either are empty.
        if splitLoc1 == [] or splitLoc2 == []:
            continue
        # if both splits exist, we check if the length of either of the splits is 1. If so, there are likely to be heaos
        # of locations where this split exists. Thus, if the length is 1 we change the location to be simply the amino
        # acid instead of all the positions it is located at within the peptide.
        if len(split1) == 1:
            splitLoc1 = split1
        if len(split2) == 1:
            splitLoc2 = split2
        # we append all the stored up location data of split1 and split2 to the splitLocations list, which
        # stores data relevant only to the specific combnation of splits in question. We then append
        # this list to the totalLocations list which stores all the location data across all potential split combinations.
        splitLocations.append(splitLoc1)
        splitLocations.append(splitLoc2)
        totalLocations.append(splitLocations)
    return totalLocations


# takes the origin dict of the form originDict[peptide] = [(proteinName, locations),(proteinName, locations)...] and writes
# it to the filePath given. Also takes the spliceType argument to know how to format the csv and name it.
def writeToFasta(originDict, outputPath, spliceType):
    with open(outputPath, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        # write a title with the splice type for the entire document, and a blank row underneath.
        title = spliceType + ' Peptide Locations'
        writer.writerow([title])
        writer.writerow([])
        # iterate throught the originDict to write each entry to file.
        for pep, origins in originDict.items():
            # write the peptide to file in a row
            infoRow = ['Peptide', pep]
            writer.writerow(infoRow)
            # write the header for the origin information.
            pepHeader = ['Prot Name', 'Location']
            writer.writerow(pepHeader)
            # if origins = [], the peptide could not be found in the prot.fasta file using the given splice type.
            if origins == []:
                writer.writerow(['None', []])
            # if origins contains data, an origin location was found and we must format the data accordingly.
            else:
                # call the data configuration function relevant to the splice type.
                if spliceType == 'Linear':
                    dataRows = linDataRow(origins)
                elif spliceType == 'Cis':
                    dataRows = cisDataRow(origins)
                # write the formated data to the row.
                for protData in dataRows:
                    writer.writerow(protData)
            # write a blank row and repeat for the next peptide.
            writer.writerow([])

# takes the linear origin data for a given peptide and formats it for writing in a line in the csv.
def linDataRow(origins):
    # iterate through each tuple (there is a tuple for every protein that was found to produce the peptide)
    dataRows = []
    for tuple in origins:
        # initialise the first column of the row to be the name of the protein.
        dataRow = [tuple[0]]
        # the second element of the tuple stores the location data. For each location found in the current protein,
        # add the information to the next column of the row.
        for location in tuple[1]:
            dataRow.append(location)
        # append dataRow to dataRows
        dataRows.append(dataRow)
    # return dataRows to be written to csv.
    return dataRows

# takes the cis origin data for a given peptide and formats it for writing in a line in the csv.
def cisDataRow(origins):
    dataRows = []
    # iterate through each tuple (there is a tuple for every protein that was found to produce a peptide)
    for tuple in origins:
        # initialise the first column of the data row to be the origin protein name.
        dataRow = [tuple[0]]
        # our location data is stored in tuple[1], it is in the form of the imbedded lists created findCisIndexes().
        locations = tuple[1]
        # locations contains a list for each pair of splits which can produce the peptide from the current protein.
        for splitCombo in locations:
            # we have a location string for each splits pair that produced the peptide. The string has the form:
            # split1Location1 or split1Location2.... with split2Location1 or split2Location2......
            # [1,2,3,4] or [5,6,7,8] with [9.10.11.12] or [13,14,15,16]
            strng = " "
            # each splitCombo list contains two lists: one for split1 locations one for split2 locations.
            for splitOptions in splitCombo:
                # if the string is not in its initialised form " ", we have switched from split1 to split2 locations and
                # need to add the "with" to the string
                if strng != " ":
                    strng += " with "
                # iterate through the locations of the specific split and add them to the strings with ors between them.
                for split in splitOptions:
                    if strng[-1] == ']':
                        strng += " or "
                        strng += str(split)
                    # no "or" required if it is the first split to be added to the string, or if it follows the "with"
                    else:
                        strng += str(split)
            # remove the initial space input at the start of the string
            strng = strng[1:]
            # append the string which has been built up for the given split combination to the next column of the dataRow
            dataRow.append(strng)
        # append the build up dataRow to dataRows
        dataRows.append(dataRow)
    # return dataRows
    return dataRows


def generateOutput(outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag):
    protDict = protFastaToDict(proteinFile)
    generateOrigins(protDict, peptideFile, outputPath, linFlag, cisFlag, transFlag)
