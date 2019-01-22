from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv

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
        print('Cis')
    if transFlag:
        #trans
        print('Trans')
    return linOriginDict

def findLinOrigins(protDict, pepList):
    linOriginDict = {}
    for pep in pepList:
        linOriginDict[pep] = []
        for protSeq in protDict.keys():
            if pep in protSeq:
                linOriginDict[pep].append((protDict[protSeq], protSeq))
    return linOriginDict

def writeToFasta(originDict, outputPath, spliceType):
    with open(outputPath + '.csv', 'a', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        title = spliceType + ' Peptide Locations'
        writer.writerow([title])
        writer.writerow([])
        for pep, origins in originDict.items():
            infoRow = ['Peptide', pep]
            writer.writerow(infoRow)
            pepHeader = ['Prot Name', 'Sequence']
            writer.writerow(pepHeader)
            if originDict[pep] == []:
                writer.writerow(['None', 'None'])
            else:
                for tuple in origins:
                    writer.writerow([tuple[0],tuple[1]])
            writer.writerow([])


protDict = protFastaToDict('Example.fasta')

pepList = pepFastaToList('samplePeps.fasta')

linOriginDict = generateOrigins(protDict, pepList, True, False, False)
print(linOriginDict)