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


def overlappingPos(ml, actual, predicted):
    return max(0, ml - abs(actual - predicted))


def totalOverlap(ml, actuals, predicteds):
    assert len(actuals) == len(predicteds), "must be same length"
    return sum(overlappingPos(ml, a, p) for (a, p) in zip(actuals, predicteds))


def results(datasetfolder):
    seq = readSequences(os.path.join(datasetfolder, "sequences.fa"))
    motifLength = readMotifLength(os.path.join(datasetfolder,
                                               "motiflength.txt"))
    actualLocations = readSites(os.path.join(datasetfolder, "sites.txt"))
    predictedLocations = readSites(os.path.join(datasetfolder,
                                                "predictedsites.txt"))
    assert len(seq) == len(actualLocations), "one location for every seq"
    assert len(actualLocations) == len(predictedLocations), "same number of"
    " predictions"

    possibleOverlap = len(seq) * motifLength
    return (totalOverlap(motifLength, actualLocations, predictedLocations),
            possibleOverlap)


def resultsGroup(folder):
    totalResults =\
        [results(os.path.join(folder, sub)) for sub in os.listdir(folder)]
    positionResults = [overlap for (overlap, _) in totalResults]
    possibleOverlap = totalResults[0][1]
    assert all([possibleOverlap == po for (_, po) in totalResults])
    print(folder)
    print("possible overlap = " + str(possibleOverlap))
    print("average overlap = " + str(sum(positionResults) /
          len(positionResults)))


def resultsAll(folder):
    for sub in os.listdir(folder):
        path = os.path.join(folder, sub)
        print(folder)
        resultsGroup(path)
