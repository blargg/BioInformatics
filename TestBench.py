import os
import time
from collections import namedtuple
import Motif as M
from FileFormat import (readSequences, readSites, readMotifLength,
                        writePositionWeightMatrix, writeSites)

PerformanceResults = namedtuple("PerformanceResults",
                                ['time', 'algorithmScore', 'motifScore'])


def runOnAllDatasets(algorithm, datasetfolder):
    for folder in os.listdir(datasetfolder):
        subfolder = os.path.join(datasetfolder, folder)
        print(folder)
        runOnGroup(algorithm, subfolder)


def runOnGroup(algorithm, folder):
    def folderName(i):
        return os.path.join(folder, "dataset" + str(i))

    results = [test(algorithm, folderName(i)) for i in range(10)]

    for i, p in enumerate(results):
        print("run " + str(i))
        print(p)
        print()

    print('average time = ' +
          str(float(sum(r.time for r in results)) / len(results)))
    print('average algorithm score = ' +
	      str(float(sum(r.algorithmScore for r in results)) / len(results)))
    print('average motif score = ' +
	      str(float(sum(r.motifScore for r in results)) / len(results)))


def test(algorithm, datasetfolder):
    """algorithm is a function that takes a motif length and a set of sequences
    and returns a list of sequences,positions of the best motif

    datasetfolder is the name of the dataset to test on
    returns PerformanceResults of the run"""

    seq = readSequences(os.path.join(datasetfolder, "sequences.fa"))
    motifLength = readMotifLength(os.path.join(datasetfolder,
                                               "motiflength.txt"))
    actualLocations = readSites(os.path.join(datasetfolder, "sites.txt"))
    assert len(seq) == len(actualLocations), "there should be a site location"
    "for every sequence"

    startTime = time.clock()
    result = algorithm(motifLength, seq)
    endTime = time.clock()

    elapsedTime = endTime - startTime
    runScore = M.informationInSelection(motifLength, result)
    actualScore = M.informationInSelection(motifLength,
                                           zip(seq, actualLocations))

    pwfilename = os.path.join(datasetfolder, "predictedmotif.txt")
    writePositionWeightMatrix(pwfilename, motifLength, result)

    sitesfilename = os.path.join(datasetfolder, "predictedsites.txt")
    writeSites(sitesfilename, [site for (_, site) in result])

    return PerformanceResults(time=elapsedTime,
                              algorithmScore=runScore,
                              motifScore=actualScore)
