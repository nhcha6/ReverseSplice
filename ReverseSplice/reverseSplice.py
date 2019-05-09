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

PROT_THRESH = 2
STOP = "STOP"

def generateOutput(outputPath, proteinFile, peptideFile, linFlag, cisFlag, transFlag, minTransLen):
    """
    Called from Example.createOutput() in reverseSpliceGUI.py, this function recieves the input parameters from the
    GUI code, creates the protein dictionary from the input path and then calls generateOrigins to begin the
    origin computation process.

    :param outputPath: the path of the output file as selected by the user.
    :param proteinFile: the path of the protein sequence file.
    :param peptideFile: the path of the peptide sequence file.
    :param linFlag: True if the user wishes to return a file with information on where each peptide could have been
    generated from via linear splicing.
    :param cisFlag: True if the user wishes to return a file with information on where each peptide could have been
    generated from via cis splicing.
    :param transFlag: True if the user wishes to return a file with information on where each peptide could have been
    generated from via trans splicing.
    :param minTransLen: the minimum length a cleavage must be for it to be reported in the output file as an origin
    for a trans spliced peptide.
    :return:
    """
    protDictList = protFastaToDict(proteinFile)
    generateOrigins(protDictList, peptideFile, outputPath, linFlag, cisFlag, transFlag, minTransLen)

def protFastaToDict(protFile):
    """
    Called by generateOutput(), this function stores all the sequences in a fasta file in a dictionary.

    :param protFile: A Fasta file path (containing all the proteins which could serve as an origin location)
    :return protDictList: A list of dictionaries which contain a protein sequence as the key, and the protein name as the
    value. The list is creates to ensure the dictionary is small enough to be passed into the pool as a global variable.
    """
    protDictList = []
    protDict = {}
    with open(protFile, "rU") as handle:
        counter = 1
        for record in SeqIO.parse(handle, 'fasta'):
            seq = str(record.seq)
            name = 'rec' + str(counter)
            protDict[name] = seq
            if counter % PROT_THRESH == 0:
                protDictList.append(protDict)
                protDict = {}
            counter+=1
        protDictList.append(protDict)
    return protDictList

def generateOrigins(protDictList, pepFile, outputPath, linFlag, cisFlag, transFlag, minTransLen):
    """
    Called by generateOutput(), this function controls the calling of a separate function for linear, cis and trans
    splicing.

    :return protDictList: A list of dictionaries which contain a protein sequence as the key, and the protein name as the
    value.
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
        findLinOrigins(protDictList, pepFile, outputPath)
    if cisFlag:
        # Find cis origins
        findCisOrigins(protDictList, pepFile, outputPath)
    if transFlag:
        # Find trans origins > work in progress
        findTransOrigins(protDictList, pepFile, outputPath, minTransLen)

def findLinOrigins(protDictList, pepFile, outputPath):
    """
    Called by generateOrigins(), this function takes the protDict and pepList and creates processes which compute
    where within the protein list each peptide could have been generated from via linear splicing. Each process it
    creates outputs a dictionary with the peptides as keys and a list of tuples containing origin peptide names and the
    splitsLocation data for linear splicing as the value.

    :return protDictList: A list of dictionaries which contain a protein sequence as the key, and the protein name as the
    value.
    :param pepFile: A file containing a list of peptides that you want to find the linear origin locations for
    :return linOriginsDict: Has the peptide as a key and a list of tuples of the form (originProtName, locations).
                            Locations store information on where the corresponding peptide could have been generated
                            from in the relevant origin protein.
    :Data structure summary: linOriginsDict[peptide] = [(proteinName, locations),(proteinName, locations)]
    """
    outputPath = outputPath + '_' + 'Linear' + '-' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'

    for protDict in protDictList:
        numWorkers = multiprocessing.cpu_count()
        toWriteQueue = multiprocessing.Queue()

        pool = multiprocessing.Pool(processes=numWorkers, initializer=processLinInitArgs,
                                    initargs=(toWriteQueue,protDict))

        writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, 'Linear', protDict))
        writerProcess.start()

        # iterate through each peptide for each split up of the input
        with open(pepFile, "rU") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                pep = str(record.seq)
                logging.info('Process started for: ' + str(pep))
                pool.apply_async(linearOrigin, args=(pep,))
        pool.close()
        pool.join()
        logging.info("Pool joined")
        toWriteQueue.put(STOP)
        writerProcess.join()

def linearOrigin(pep):
    """
    Called as the worker function to the pool in findLinOrigins(), this function takes an individual peptide and a
    dictionary of protein sequences, and returns the proteins and locations within from which the peptide could be
    generated via linear splicing. Once the origin data is compiled, it is added to the toWriteQueue so that it
    can be processed by the linearWriter() function.

    :param pep: the peptide which the user wishes to find potential linear splicing origins of.
    :param protDict: a dictionary containing a portion (possibly all) of the input protein sequences.
    :return:
    """
    try:

        print(pep)
        linOriginDict = {}
        # initialise the key as an empty list in the outputDict
        linOriginDict[pep] = []
        # iterate through each protSeq in the keys of protDict
        for key, value in proteinDict.items():

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

def processLinInitArgs(toWriteQueue, protDict):
    """
    Called from findLinOrigins() when the pool is initialised, this function simply gives linearOrigin() (the worker
    function for each process in the pool) access to the toWriteQueue.
    """

    linearOrigin.toWriteQueue = toWriteQueue
    global proteinDict
    proteinDict = protDict

def findCisOrigins(protDictList, pepFile, outputPath):
    """
    Called by generateOrigins(), this function takes the protDict and pepList and creates processes which compute
    where within the protein list each peptide could have been generated from via cis splicing. Each process it
    creates outputs a dictionary with the peptides as keys and a list of tuples containing origin peptide names and the
    splitsLocation data for cis splicing as the value.

    Data structure summary: cisOriginsDict[peptide] = [(proteinName, locations),(proteinName, locations)]

    :return protDictList: A list of dictionaries which contain a protein sequence as the key, and the protein name as the
    value.
    :param pepFile: A file containing a list of peptides that you want to find the linear origin locations for
    :return:
    """
    outputPath = outputPath + '_' + 'Cis' + '-' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'
    for protDict in protDictList:
        numWorkers = multiprocessing.cpu_count()
        toWriteQueue = multiprocessing.Queue()
        pool = multiprocessing.Pool(processes=numWorkers, initializer=processCisInitArgs,
                                    initargs=(toWriteQueue,))
        writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, 'Cis', protDict))
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
    """
    Called as the worker process to the pool in findCisOrigins(), this function takes an individual peptide and a
    dictionary of protein sequences, and returns the proteins and locations within them from which the peptide could be
    generated via cis splicing. Once the origin data is compiled, it is added to the toWriteQueue so that it
    can be processed by the writer() function.

    :param pep: the peptide which the user wishes to find potential linear splicing origins of.
    :param protDict: a dictionary containing all the input protein sequences.
    :return:
    """
    try:

        linFlag = False
        cisOriginDict = {}

        # initialise that key in the dictionary
        cisOriginDict[pep] = []
        # create the altered pep:
        alteredPep = pep.replace('I', 'L')
        # find the splits which could be combined to create the peptide using Cis splicing.
        # cisSplits is a list of tups, where each tuple contains two complimentary splits.
        # cisSplits = [('A', 'BCD'),('AB', 'CD'),('ABC', 'D')]
        cisSplits = findSplits(alteredPep)
        # iterate through each protein in protDict.keys()
        for protName, protSeq in protDict.items():
            # replace all Is with Js as they are indeciferable on mass spec.
            alteredProt = protSeq.replace('I', 'L')
            # change the linFlag if the pep exists as a linear peptide somewhere in the input
            if alteredPep in alteredProt:
                linFlag = True
            # find the location data corresponding the current protSeq.
            locations = findCisIndexes(cisSplits, alteredProt)
            # if there no pairs of splits are located in protSeq, findCisIndexes() will return an empty list.
            # If so, continue.
            if locations == []:
                continue
            # if it is possible to create the peptide using cis splicing from the current protein, add the protName
            # and locations tuple to cisOriginsDict
            cisOriginDict[pep].append((protName, locations))
        if linFlag:
            cisOrigin.toWriteQueue.put({})
        else:
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
    """
    Called from findLinOrigins() when the pool is initialised, this function simply gives linearOrigin() (the worker
    function for each process in the pool) access to the toWriteQueue.
    """

    cisOrigin.toWriteQueue = toWriteQueue

def findSplits(pep):
    """
    Called by cisOrigins() and transOrigins(), this function takes a peptide and returns a list of tuples, where each
    tuple is a possible pair of subsequences which could be combined to make the peptide.

    :param pep: the input peptide. From this peptide, a list of all the pair of cleavages which could be combined to
    make the peptide is returned.
    :return cisSplits: a list of tuples, where each tuple is a possible pair of subsequences which could be combined
    to make the peptide.
    """
    cisSplits = []
    lngth = len(pep)
    for i in range(1,lngth):
        split1 = pep[0:i]
        split2 = pep[i:lngth+1]
        cisSplits.append((split1,split2))
    return cisSplits

def findCisIndexes(cisSplits, protSeq):
    """
    Called by cisOrigins(), this function returns the location data of where different pairs of cisSplits (which
    combine to form a given peptide) exist in a protein sequence. The output data structure is a series of embedded
    lists. The outer list contains lists relating to each cisSplit combination which has both splits found in the protSeq.
    The list relating to each cisSplit combination contains two lists; a list of all the places a first split can
    occur and a list of all the places the second split can occur.
    Thus our final data structure looks like:
    outerList[firstSplitCombo[firstSplitList[[1,2,3,4],[5,6,7,8]], secondSplitList[[9,10,11,12],[13,14,15,16]], secondSplitCombo[etc....]]

    :param cisSplits: a list of all the pairs of cleavages which could be combined to form a given peptide via cis
    splicing.
    :param protSeq: the proteins sequence within which these pairs are being searched for. If both pairs are found,
    the location data of where they were found within the protein is added to the data structure described above.
    :return:
    """
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

def findTransOrigins(protDictList, pepFile, outputPath, minTransLen):
    """
    Called by generateOrigins(), this function takes the protDict and pepList and creates processes which compute
    where within the protein list each peptide could have been generated from via trans splicing. Each process it
    creates outputs a dictionary with the peptides as keys and a list of tuples containing origin peptide names and the
    splitsLocation data for trans splicing as the value.
    Data structure summary: transOriginsDict[peptide] = [(proteinName, locations),(proteinName, locations)]

    :return protDictList: A list of dictionaries which contain a protein sequence as the key, and the protein name as the
    value.
    :param pepFile: A file containing a list of peptides that you want to find the linear origin locations for
    :return:
    """
    outputPath = outputPath + '_' + 'Trans' + '-' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'
    for protDict in protDictList:
        numWorkers = multiprocessing.cpu_count()
        toWriteQueue = multiprocessing.Queue()

        pool = multiprocessing.Pool(processes=numWorkers, initializer=processTransInitArgs,
                                    initargs=(toWriteQueue,))

        writerProcess = multiprocessing.Process(target=writer, args=(toWriteQueue, outputPath, 'Trans', protDict))
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
    """
    Called as the worker process to the pool in findTransOrigins(), this function takes an individual peptide and a
    dictionary of protein sequences, and returns the proteins and locations within them from which the peptide could be
    generated via trans splicing. Once the origin data is compiled, it is added to the toWriteQueue so that it
    can be processed by the writer() function.

    :param pep: the peptide which the user wishes to find potential linear splicing origins of.
    :param protDict: a dictionary containing all the input protein sequences.
    :param minTransLen: the minimum length that a cleavage must be for its location to be reported in the origin data
    of a trans spliced peptide.
    :return:

    """
    try:
        # initialise the dictionary
        transOriginDict = {}
        # initialise that key in the dictionary
        transOriginDict[pep] = []
        # create the altered pep:
        alteredPep = pep.replace('I', 'L')
        # find the splits which could be combined to create the peptide using trans or cis splicing.
        # transSplits is a list of tups, where each tuple contains two complimentary splits.
        # transSplits = [('A', 'BCD'),('AB', 'CD'),('ABC', 'D')]
        transSplits = findSplits(alteredPep)

        # format splits so it iterates in the order we want it to.
        transSplits = editTransSplits(transSplits, minTransLen)

        # pepFound bool allows us to skip all splits with both entries under minTransLen if
        # it has already been found.
        pepFound = False

        # iterate through each combination in transSplits
        for splitCombo in transSplits:
            # declare the two splits in the combo
            split1 = splitCombo[0]
            split2 = splitCombo[1]

            # if the first entry is less than the min lengths and the pep has already been found
            # we can break
            if len(split2) < minTransLen and (pepFound or transOriginDict[pep]!=[]):
                break

            # declare holder for split1 location. If the corresponding split is greater than minTransLen, we want to
            # initialise it as a list which will store possible locations. If it is smaller than minTransLen, we don't
            # care about the location data just if it can be found or not. Thus we initialise it as False, and later
            # will change it to True if it is founds anywhere.
            if len(split1) >= minTransLen:
                splitLoc1 = []
            else:
                splitLoc1 = False
            # declare holder for split2 location in the same format as split1
            if len(split2) >= minTransLen:
                splitLoc2 = []
            else:
                splitLoc2 = False

            # iterate through each protein in protDict.keys()
            for protName, protSeq in protDict.items():
                # replace all Is with Js as they are indeciferable on mass spec.
                alteredProt = protSeq.replace('I', 'L')

                # check for the presence of split1 and split2 in the same protSeq. If both exist, it is a cis or lin
                # peptide so should be ignored from the trans output.
                if split1 in alteredProt and split2 in alteredProt:
                    transOrigin.toWriteQueue.put(transOriginDict)
                    return

                # check for presence of split1 in the current protein
                # if splitLoc1 == True, we know that this split has been found and we can
                # continue through without checking the current protein for it.
                if splitLoc1 == True:
                    pass
                # if it is not True, we check for its presence in alteredProt. We build up an iterable
                # which stores the locations of the split in protein. If it can't be found it won't enter the
                # for loop. If at least one can be found, splitLoc1 is either set to True, or the location
                # information is appended to it.
                else:
                    for x in re.finditer(split1, alteredProt):
                        if splitLoc1 == False:
                            splitLoc1 = True
                            break
                        else:
                            splitLoc1.append([protName, x.start(), x.end() - 1])

                # check for the presence of split2, exactly the same as we have in splitLoc1
                if splitLoc2 == True:
                    pass
                else:
                    for x in re.finditer(split2, alteredProt):
                        if splitLoc2 == False:
                            splitLoc2 = True
                            break
                        else:
                            splitLoc2.append([protName, x.start(), x.end() - 1])

            # When we get to here we have completed an iteration through every protein for the current combination of
            # splits. We know need to decide what to do with splitLoc1 based on what data it holds.
            # if either splitLoc variables equal False (one of the splits was less than minLenTrans and not found),
            # we want to continue to the next iteration without doing anything.
            if splitLoc1 == False or splitLoc2 == False:
                continue
            # if both splitLoc variables are True (both splits were found but both are smaller than minLenTrans),
            # we want to update pepFound and continue without adding to transOriginsDict.
            if splitLoc1 == True and splitLoc2 == True:
                pepFound = True
                continue

            # we reach here if at least one of the splits was longer than minTransLen and thus the corresponding
            # splitLoc was initialised as a list.
            # if splitLoc1 is not an empty list and splitLoc2 is True, this split combination can be used to create
            # the trans splicing of the pep. Thus, we update toAppend with the location data stored in splitLoc1.
            if splitLoc1 != [] and splitLoc2 == True:
                toAppend = splitLoc1
            # same scenario as previous except splitLoc1 and 2 are switched.
            elif splitLoc2 != [] and splitLoc1 == True:
                toAppend = splitLoc2
            # if we get to here both splitLocs must be lists. If they both aren't empty, both splits were found and
            # thus both splitLoc data must be added to toAppend.
            elif splitLoc1 != [] and splitLoc2 != []:
                toAppend = splitLoc1 + splitLoc2
            else:
                continue

            # add toAppend to transOriginDict. toAppend will be empty if no transSplicing was found involving
            # splits of length greater than minTransLen.
            transOriginDict[pep] += toAppend

        # if after iterating through all splitCombos and all proteins for each combo no location data has been added,
        # we must check if pepFound has been set to True to ensure no splitCombos with split lengths < minTransLen
        # were found. If it is, simplt add True to the dictionary.
        if transOriginDict[pep] == [] and pepFound:
            transOriginDict[pep] = [True]

        # add this transOriginDict (related to the given pep) to the writer function.
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
    """
    Called by transOrigin(), this function takes the splits that could be combined to create the current peptide,
    and sorts them in the required order for the rest of the algorithm to work. We want the longest split to be first
    in each pair, and we want the tuple pairs which have both splits with length less than minTransLen to be last
    in the list.
    :param splits: a list of tuples, where each tuple is a possible pair of subsequences which could be combined to
    make the peptide.
    :param minTransLen: the minimum length that a cleavage must be for its location to be reported in the origin data
    of a trans spliced peptide.
    :return splitsNew: splits sorted so that the longest split in each tuple appears first, and the tuples with both
    splits less than minTransLen at the end of the list.
    """
    #print(splits)
    splits1 = []
    splits2 = []
    for tuple in splits:
        # sort the tuple so that the longest split appears first so that it is checked first
        tuple = sorted(tuple, key=len)
        # we want the tuples which have both splits < MIN_TRANS_LEN to be at the end. We only run
        # them if none of the previous tuples have been found.
        if len(tuple[1]) < minTransLen:
            splits2.append(tuple)
        else:
            splits1.append(tuple)
    splitsNew = splits1 + splits2
    #print(splitsNew)
    return splitsNew

def processTransInitArgs(toWriteQueue):
    """
    Called from findTransOrigins() when the pool is initialised, this function simply gives transOrigin() (the worker
    function for each process in the pool) access to the toWriteQueue.
    """
    transOrigin.toWriteQueue = toWriteQueue

def writer(toWriteQueue, outputPath, spliceType, protDict):
    """
    The writer function for linear spliced origin computation. This function simply pulls origin data to be written
    to file from the toWriteQueue, stores it in the finalLinOriginDict and once all processes are finished, writes
    the final output csv. finalLinOriginsData structure: { inputPeptide: locationDataStructure }
    Note that the locationDataStructure differs for trans, linear and cis splicing.

    :param toWriteQueue: the multiprocessing.Queue which the origin data constructed by linearOrigin() is pushed to
    at the end of each process.
    :param outputPath: the path of the linear output csv file.
    :param spliceType: the type of splicing being run in the current iteration.
    :param protDict: the dictionary containing all the input proteins, which is required when writing to file.
    :return:
    """

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

def writeToFasta(originDict, outputPath, spliceType, protDict):
    """
    This function is called by writer(), and takes the ouptut data dictionary and writes it to the filePath given.
    Also takes the spliceType argument to know how to format the csv and name it.

    :param originDict: the dictionary containing the input peptides as keys and the location data as values. Structure:
    originDict[peptide] = [(proteinName, locations),(proteinName, locations)...]
    :param outputPath: the file location and file name to which the output is to be written.
    :param spliceType: the type of splicing being computed. Each splice type has a different output syntax, and the
    type of splicing also needs to be added to the output file name.
    :param protDict: a dictionary containing the input protein data. This is needed to return slight differences in the
    peptide and origin due to the program not treating I/J differently.
    :return:
    """
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

def linDataRow(origins, pep, protDict):
    """
    Called by writeToFasta, this function takes the linear origin data for a given peptide and formats it for writing
    to the csv file.

    :param origins: the data structure containing information on where a given peptide was found within the input
    proteins. Has the form: [('ProtName1', [[location1]]), ('ProtName2', [[location1], [location2]])]
    :param pep: the peptide which the origin data related to.
    :param protDict: a dictionary containing the input protein data. This is needed to return slight differences in the
    peptide and origin due to the program not treating I/J differently.
    :return dataRows: a list of lists, where each sublist is a row of data which is to be written to file. Each sublist
    has the format: [protName, peptide, pepInProt, location]:
    """
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

def cisDataRowNew(origins, pep, protDict):
    """
    Called by writeToFasta(), this function takes the cis origin data for a given peptide and formats it for writing
    to the csv file.

    :param origins: the data structure containing information on where a given peptide was found within the input
    proteins. Has the form:
    [('ProtName1', [[[b1 location1, b1 location2][y1 location1, y1 location2]],[[b3 location1],[y3 location1, y3 location2]]]), ('ProtName2', ....)]
    :param pep: the peptide which the origin data related to.
    :param protDict: a dictionary containing the input protein data. This is needed to return slight differences in the
    peptide and origin due to the program not treating I/J differently.
    :return dataRows: a list of lists, where each sublist is a row of data which is to be written to file. Each sublist
    has the format: [protName, peptide, pepInProt, location]
    """
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

def transDataRow(origins, pep, protDict):
    """
    Called by writeToFasta(), this function takes the trans origin data for a given peptide and formats it for writing
    to the csv file.

    :param origins: the data structure containing information on where a given peptide was found within the input
    proteins. Has the form: [[protName, startIndex, endIndex]..] where each sublist refers to the location of an
    individual cleavage which can be combined with another cleavage somewhere in the protein file to create the peptide.
    :param pep: the peptide which the origin data related to.
    :param protDict: a dictionary containing the input protein data. This is needed to return slight differences in the
    peptide and origin due to the program not treating I/J differently.
    :return dataRows: a list of lists, where each sublist is a row of data which is to be written to file. Each sublist
    has the format: [protName, peptide, pepInProt, location]
    """
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