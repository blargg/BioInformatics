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
