import Motif as M
import random


def uniformRandomSampling(iterations, motifLength, sequences):
    """Finds motifs by selecting a random placement, and checking the
    fitness"""
    # chose an arbitrary best to start
    bestSelection = [(s, 0) for s in sequences]

    for _ in range(iterations):
        selections =\
            [(s, random.randrange(len(s) - motifLength + 1))
             for s in sequences]
        if M.informationInSelection(motifLength, selections) >\
                M.informationInSelection(motifLength, bestSelection):
            bestSelection = selections

    return bestSelection


def SGDpossibleMutations(maxIndex, currentSelection):
    mutations = []
    for i, v in enumerate(currentSelection):
        assert v < maxIndex, "currentSelection is invalid"
        for j in range(maxIndex):
            if j is not v:
                m = currentSelection[:]
                m[i] = j
                mutations.append(m)

    return mutations


def SGDimprove(evalPositions, motifLength, selections, sequences):
    seqLen = len(sequences[0])
    maxIndex = seqLen - motifLength
    possiblePos = SGDpossibleMutations(maxIndex, selections)

    return max(possiblePos, key=evalPositions)


def SGD(iterations, motifLength, sequences):
    """Selects a random assignment of positions, then repeadedly makes small
    changes until it reaches a local maximum, then repeats"""
    bestSelection = [0 for s in sequences]

    def evalPositions(pos):
        return M.infoInPositions(motifLength, sequences, pos)

    for _ in range(iterations):
        lastSelection = [random.randrange(len(s) - motifLength)
                         for s in sequences]
        currentSelection = SGDimprove(evalPositions, motifLength,
                                      lastSelection, sequences)
        while abs(evalPositions(currentSelection) -
                  evalPositions(lastSelection)) < 0.000001:
            assert len(lastSelection) == len(sequences),\
                    "there should be exactly one selection per seq"
            lastSelection = currentSelection
            currentSelection = SGDimprove(evalPositions,
                                          motifLength,
                                          lastSelection,
                                          sequences)
            print(currentSelection)
            print(evalPositions(currentSelection))
        bestSelection = max(bestSelection, currentSelection, key=evalPositions)
    return zip(sequences, bestSelection)


def greedy(motifLength, sequences):

    def evalWithSeq(seq):
        return lambda pos: M.infoInPositions(motifLength, seq, pos)

    maxIndex = len(sequences[0]) - motifLength
    initialoptions = [[i, j] for i in range(maxIndex) for j in range(maxIndex)]
    cumulative = max(initialoptions, key=evalWithSeq(sequences[:2]))

    for i in range(2, len(sequences)):
        options = []
        for j in range(maxIndex):
            copyList = cumulative[:]
            copyList.append(j)
            options.append(copyList)
        cumulative = max(options, key=evalWithSeq(sequences[:i+1]))

    return list(zip(sequences, cumulative))
