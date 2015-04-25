import collections
from math import log
from GenerateSequence import bases


def informationInMotif(selectedSubsequences):
    """selectedSubsequences is a list of subsequences that are being considered
    for a motif
        returns the information gain of this selection (higher is better)"""

    # check if all the subseqences are the same length
    lengths = [len(s) for s in selectedSubsequences]
    assert(all(l == lengths[0] for l in lengths))

    length = lengths[0]

    counters = [collections.Counter() for x in range(length)]
    for seq in selectedSubsequences:
        for (i, c) in enumerate(seq):
            counters[i].update([c])

    total = 0.0
    for i in range(length):
        for base in bases:
            p = counters[i][base] / sum(counters[i].values())
            if p == 0.0:
                total += 0.0
            else:
                total += p * log(p / 0.25)

    return total


def extractSelections(motifLength, selections):
    subsequences = [s[i:i+motifLength] for (s, i) in selections]
    assert(all(len(s) == motifLength for s in subsequences))
    return subsequences


def informationInSelection(motifLength, selections):
    """ motif length is the lenth of the motif
    selections is a list of (sequence, position) where
    sequence is the whole sequence
    position is the sarting position for the subsequence to select
    returns the information in that selection"""
    subsequences = [s[i:i+motifLength] for (s, i) in selections]
    assert(all(len(s) == motifLength for s in subsequences))
    return informationInMotif(subsequences)
