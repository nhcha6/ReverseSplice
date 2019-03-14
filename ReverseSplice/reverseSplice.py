from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from datetime import datetime
import multiprocessing
from multiprocessing import Queue
import logging
import re
import io
import traceback

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

# Define maxLen
maxLen = 20
STOP = "STOP"

#LTLWTGNN

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
            protDict[name] = seq
    return protDict


def generateOrigins(protDict, pepFile, outputPath, linFlag, cisFlag, transFlag, minTransLen):
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
        findTransOrigins(protDict, pepFile, outputPath, minTransLen)


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

    writerProcess = multiprocessing.Process(target=linearWriter, args=(toWriteQueue, outputPath, 'Linear', protDict))
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
    try:
        linOriginDict = {}
        # initialise the key as an empty list in the outputDict
        linOriginDict[pep] = []
        # iterate through each protSeq in the keys of protDict
        for key, value in protDict.items():

            # initialise the locations holder
            locations = []
            # change the Is to Ls for both the pep and prot, as they have identical masses and
            # are indeciferable on mass spec.
            alteredPep = pep.replace('I', 'L')
            alteredProt = value.replace('I', 'L')
            # re.finditer(substring, string) returns an iterable with the start and finish indices of all the locations
            # of the substring in the string. The iterator is empty if the susbset does not appear at all.
            # We thus iterate through all instances of the subset and append it to a list of locations.
            for x in re.finditer(alteredPep, alteredProt):
                # convert the start index and end index to a list of indexes and then append to locations
                # locations structure is a list of lists: [[1,2,3,4],[5,6,7,8]]
                locations.append([x.start(), x.end() - 1])
            # if nothing is added to locations, it means that the peptide was not found in the protein, and we continue
            # iterating through proteins.
            if locations:
                # otherwise if we have added to locations, we append the protName/location tup to linOriginDict
                linOriginDict[pep].append((key, locations))
        linearOrigin.toWriteQueue.put(linOriginDict)
        logging.info('Process complete for: ' + str(pep))

    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in worker process: ' + pep + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e


def linearWriter(toWriteQueue, outputPath, spliceType, protDict):

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
    writeToFasta(finalLinOriginDict, outputPath, spliceType, protDict)

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
    writerProcess = multiprocessing.Process(target=linearWriter, args=(toWriteQueue, outputPath, 'Cis', protDict))
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

    try:
        cisOriginDict = {}

        # initialise that key in the dictionary
        cisOriginDict[pep] = []
        # create the altered pep:
        alteredPep = pep.replace('I', 'L')
        # find the splits which could be combined to create the peptide using Cis splicing.
        # cisSplits is a list of tups, where each tuple contains two complimentary splits.
        # cisSplits = [('A', 'BCD'),('AB', 'CD'),('ABC', 'D')]
        cisSplits = findCisSplits(alteredPep)
        # iterate through each protein in protDict.keys()
        for protName, protSeq in protDict.items():
            # replace all Is with Js as they are indeciferable on mass spec.
            alteredProt = protSeq.replace('I', 'L')
            # ignore that protSeq if the linear splice comes from it already.
            if alteredPep in alteredProt:
                continue
            # find the location data corresponding the current protSeq.
            locations = findCisIndexes(cisSplits, alteredProt)
            # if there no pairs of splits are located in protSeq, findCisIndexes() will return an empty list.
            # If so, continue.
            if locations == []:
                continue
            # if it is possible to create the peptide using cis splicing from the current protein, add the protName
            # and locations tuple to cisOriginsDict
            cisOriginDict[pep].append((protName, locations))
        cisOrigin.toWriteQueue.put(cisOriginDict)
        logging.info('Cis Process complete for: ' + str(pep))

    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in worker process: ' + pep + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e

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
            splitLoc1.append([x.start(), x.end() - 1])
        # append the location of any split2 occurences in protSeq to splitLoc2
        for x in re.finditer(split2, protSeq):
            splitLoc2.append([x.start(), x.end() - 1])
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
def writeToFasta(originDict, outputPath, spliceType, protDict):
    with open(outputPath, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        # write a title with the splice type for the entire document, and a blank row underneath.
        title = spliceType + ' Peptide Locations'
        writer.writerow([title])
        writer.writerow([])
        # there is a single listing for every unique location/combination that a peptide can come from.
        # they are listed under the format of the header
        header = ['Prot Name', 'Peptide', 'Pep in Prot', 'Location']
        writer.writerow(header)
        # iterate throught the originDict to write each entry to file.
        for pep, origins in originDict.items():
            # if origins = [], the peptide could not be found in the prot.fasta file using the given splice type.
            if origins == []:
                continue
            # if origins contains data, an origin location was found and we must format the data accordingly.
            else:
                # call the data configuration function relevant to the splice type.
                if spliceType == 'Linear':
                    dataRows = linDataRow(origins, pep, protDict)
                elif spliceType == 'Cis':
                    dataRows = cisDataRowNew(origins, pep, protDict)
                else:
                    dataRows = transDataRow(origins, pep, protDict)
                # write the formated data to the row.
                for protData in dataRows:
                    writer.writerow(protData)
            # write a blank row and repeat for the next peptide.
            writer.writerow([])

# takes the linear origin data for a given peptide and formats it for writing in a line in the csv.
def linDataRow(origins, pep, protDict):
    # iterate through each tuple (there is a tuple for every protein that was found to produce the peptide)
    dataRows = []
    for tuple in origins:
        # initialise the first column of the row to be the name of the protein.
        firstHalf = [tuple[0]]
        # append the peptide as it exist
        firstHalf.append(pep)
        # the second element of the tuple stores the location data. For each location found in the current protein,
        # add the information to the next column of the row.
        for location in tuple[1]:
            protPep = protDict[tuple[0]][location[0]:location[1]+1]
            secondHalf = [protPep, location]
            dataRow = firstHalf + secondHalf
            # append dataRow to dataRows
            dataRows.append(dataRow)
    # return dataRows to be written to csv.
    return dataRows

# takes the cis origin data for a given peptide and formats it for writing in a line in the csv
def cisDataRowNew(origins, pep, protDict):
    dataRows = []
    for tuple in origins:
        dataRow = []
        protName = tuple[0]
        locationData = tuple[1]
        firstHalf = [protName]
        firstHalf.append(pep)
        #print(firstHalf)
        for splitCombo in locationData:
            split1List = splitCombo[0]
            split2List = splitCombo[1]
            for split1 in split1List:
                for split2 in split2List:
                    prot = protDict[protName]
                    # check that they combine in the correct order
                    # check for a split of len1
                    if len(split1) == 1:
                        pepInProt = split1 + prot[split2[0]:split2[1]+1]
                    elif len(split2) == 1:
                        pepInProt = prot[split1[0]:split1[1] + 1] + split2
                    else:
                        pepInProt = prot[split1[0]:split1[1]+1] + prot[split2[0]:split2[1]+1]
                    location = str(split1) + ' and ' + str(split2)
                    secondHalf = [pepInProt]
                    secondHalf.append(location)
                    dataRow = firstHalf + secondHalf
                    dataRows.append(dataRow)
    return dataRows

def findTransOrigins(protDict, pepFile, outputPath, minTransLen):
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
    outputPath = outputPath + '_' + 'Trans' + '-' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'

    pool = multiprocessing.Pool(processes=numWorkers, initializer=processTransInitArgs,
                                initargs=(toWriteQueue,))

    writerProcess = multiprocessing.Process(target=linearWriter, args=(toWriteQueue, outputPath, 'Trans', protDict))
    writerProcess.start()

    # iterate through each peptide
    with open(pepFile, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            pep = str(record.seq)
            logging.info('Process started for: ' + str(pep))
            pool.apply_async(transOrigin, args=(pep, protDict, minTransLen))
        pool.close()
        pool.join()
    logging.info("Pool joined")
    toWriteQueue.put(STOP)
    writerProcess.join()

def transOrigin(pep,protDict, minTransLen):
    try:
        transOriginDict = {}

        # initialise that key in the dictionary
        transOriginDict[pep] = []
        # create the altered pep:
        alteredPep = pep.replace('I', 'L')
        # find the splits which could be combined to create the peptide using Cis splicing.
        # cisSplits is a list of tups, where each tuple contains two complimentary splits.
        # cisSplits = [('A', 'BCD'),('AB', 'CD'),('ABC', 'D')]
        transSplits = findCisSplits(alteredPep)

        # format splits so it iterates in the order we want it to.
        transSplits = editTransSplits(transSplits, minTransLen)

        # pepFound bool allows us to skip all splits with both entries under MIN_TRANS_LEN if
        # it has already been found.
        pepFound = False

        # iterate through transSlits
        for splitCombo in transSplits:
            # declare the two splits in the combo
            split1 = splitCombo[0]
            split2 = splitCombo[1]

            # if the first entry is less than the min lengths and the pep has already been found
            # we can break
            if len(split1) < minTransLen and pepFound:
                break

            # declare holder for split1 location
            if len(split1) >= minTransLen:
                splitLoc1 = []
            else:
                splitLoc1 = False
            # declare holder for split2 location
            if len(split2) >= minTransLen:
                splitLoc2 = []
            else:
                splitLoc2 = False

            # iterate through each protein in protDict.keys()
            for protName, protSeq in protDict.items():
                # replace all Is with Js as they are indeciferable on mass spec.
                alteredProt = protSeq.replace('I', 'L')
                # check for the presence of split1
                if splitLoc1 == True:
                    pass
                else:
                    for x in re.finditer(split1, alteredProt):
                        if splitLoc1 == False:
                            splitLoc1 = True
                            break
                        else:
                            splitLoc1.append([protName, x.start(), x.end() - 1])
                # check for the presence of split2
                if splitLoc2 == True:
                    pass
                else:
                    for x in re.finditer(split2, alteredProt):
                        if splitLoc2 == False:
                            splitLoc2 = True
                            break
                        else:
                            splitLoc2.append([protName, x.start(), x.end() - 1])

            if splitLoc1 == False or splitLoc2 == False:
                continue
            if splitLoc1 == True and splitLoc2 == True:
                pepFound = True
                continue

            if splitLoc1 != [] and splitLoc2 == True:
                toAppend = splitLoc1
            elif splitLoc2 != [] and splitLoc1 == True:
                toAppend = splitLoc2
            elif splitLoc1 != [] and splitLoc2 != []:
                toAppend = splitLoc1 + splitLoc2
            else:
                continue

            transOriginDict[pep] += toAppend

        if transOriginDict[pep] == [] and pepFound:
            transOriginDict[pep] = [True]

        transOrigin.toWriteQueue.put(transOriginDict)

    except Exception as e:

        exc_buffer = io.StringIO()

        traceback.print_exc(file=exc_buffer)

        errorString = 'Uncaught exception in worker process: ' + pep + '\n%s'

        logging.error(

            errorString,

            exc_buffer.getvalue())

        raise e

def editTransSplits(splits, minTransLen):
    splits1 = []
    splits2 = []
    for tuple in splits:
        # sort the tuple so that the longest split appears first so that it is checked first
        tuple = sorted(tuple)
        # we want the tuples which have both splits < MIN_TRANS_LEN to be at the end. We only run
        # them if none of the previous tuples have been found.
        if len(tuple[1]) < minTransLen:
            splits2.append(tuple)
        else:
            splits1.append(tuple)
    splitsNew = splits1 + splits2
    return splitsNew

def transDataRow(origins, pep, protDict):
    dataRows = []
    for location in origins:
        if location == True:
             dataRow = [pep, "Formed only by cleavages under max length."]
             dataRows.append(dataRow)
        else:
            protName = location[0]
            startRef = location[1]
            endRef = location[2] + 1
            pepInProt = protDict[protName][startRef:endRef]
            dataRow = [location[0], pep, pepInProt, [startRef + 1, endRef]]
            dataRows.append(dataRow)
    return dataRows


def processTransInitArgs(toWriteQueue):
    transOrigin.toWriteQueue = toWriteQueue


def generateOutput(outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag, minTransLen):
    protDict = protFastaToDict(proteinFile)
    generateOrigins(protDict, peptideFile, outputPath, linFlag, cisFlag, transFlag, minTransLen)
