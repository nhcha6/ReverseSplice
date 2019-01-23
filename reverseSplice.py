from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from datetime import datetime
import re

#added comment

# cis:  generate all split pairs possible to create a peptide, and
#       then search through all proteins in second file to see if
#       any proteins contain both splits.

# lin: search if peptide is in any proteins.

# trans: generate all splits of 5+ possible, and search for them in the file.

def protFastaToDict(protFile):
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

def pepFastaToList(pepFile):
    pepList = []
    with open(pepFile, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seq = str(record.seq)
            pepList.append(seq)
    return pepList

def generateOrigins(protDict, pepList, linFlag, cisFlag, transFlag):
    if linFlag:
        #lin
        linOriginDict = findLinOrigins(protDict, pepList)
        writeToFasta(linOriginDict, 'linOrigins', 'Linear')
    if cisFlag:
        #cis
        cisOriginDict = findCisOrigins(protDict, pepList)
        writeToFasta(cisOriginDict, 'cisOrigins', 'Cis')
    if transFlag:
        #trans
        print('Trans')
    return cisOriginDict

def findLinOrigins(protDict, pepList):
    linOriginDict = {}
    for pep in pepList:
        linOriginDict[pep] = []
        for protSeq in protDict.keys():
            locations = []
            for x in re.finditer(pep, protSeq):
                locations.append([i for i in range(x.start(), x.end())])
            if locations == []:
                continue
            linOriginDict[pep].append((protDict[protSeq], locations))
    return linOriginDict

def findCisOrigins(protDict, pepList):
    cisOriginDict = {}
    for pep in pepList:
        cisOriginDict[pep] = []
        cisSplits = findCisSplits(pep)
        for protSeq in protDict.keys():
            # ignore that protSeq if the linear splice comes from it already.
            if pep in protSeq:
                continue
            locations = findCisIndexes(cisSplits, protSeq)
            if locations == []:
                continue
            cisOriginDict[pep].append((protDict[protSeq], locations))
    return cisOriginDict

def findCisSplits(pep):
    cisSplits = []
    lngth = len(pep)
    for i in range(1,lngth):
        split1 = pep[0:i]
        split2 = pep[i:lngth+1]
        cisSplits.append((split1,split2))
    return cisSplits

def findCisIndexes(cisSplits, protSeq):
    totalLocations = []
    for splitTup in cisSplits:
        splitLocations = []
        split1 = splitTup[0]
        split2 = splitTup[1]
        splitLoc1 = []
        splitLoc2 = []
        for x in re.finditer(split1, protSeq):
            splitLoc1.append([i for i in range(x.start(), x.end())])
        for x in re.finditer(split2, protSeq):
            splitLoc2.append([i for i in range(x.start(), x.end())])
        if splitLoc1 == [] or splitLoc2 == []:
            continue
        if len(split1) == 1:
            splitLoc1 = split1
        if len(split2) == 1:
            splitLoc2 = split2
        splitLocations.append(splitLoc1)
        splitLocations.append(splitLoc2)
        totalLocations.append(splitLocations)
    return totalLocations

def writeToFasta(originDict, outputPath, spliceType):
    finalPath = outputPath + '_' + datetime.now().strftime("%d%m%y_%H%M") + '.csv'
    with open(finalPath, 'a', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        title = spliceType + ' Peptide Locations'
        writer.writerow([title])
        writer.writerow([])
        for pep, origins in originDict.items():
            infoRow = ['Peptide', pep]
            writer.writerow(infoRow)
            pepHeader = ['Prot Name', 'Location']
            writer.writerow(pepHeader)
            if origins == []:
                writer.writerow(['None', []])
            else:
                if spliceType == 'Linear':
                    dataRow = linDataRow(origins)
                elif spliceType == 'Cis':
                    dataRow = cisDataRow(origins)
                writer.writerow(dataRow)
            writer.writerow([])

def linDataRow(origins):
    for tuple in origins:
        dataRow = [tuple[0]]
        for location in tuple[1]:
            dataRow.append(location)
        return dataRow

def cisDataRow(origins):
    dataRow = []
    for tuple in origins:
        dataRow = [tuple[0]]
        locations = tuple[1]
        for splitCombo in locations:
            strng = " "
            for splitOptions in splitCombo:
                if strng != " ":
                    strng += " with "
                for split in splitOptions:
                    if strng[-1] == ']':
                        strng += " or "
                        strng += str(split)
                    else:
                        strng += str(split)
            dataRow.append(strng)
    return dataRow



protDict = protFastaToDict('Example.fasta')

pepList = pepFastaToList('samplePeps.fasta')

print(generateOrigins(protDict, pepList, True, True, False))
