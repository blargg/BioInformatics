import random


def randomBase():
    return random.choice(['a', 'c', 'g', 't'])


def randomSequence(length):
    seq = ""
    for i in range(length):
        seq += randomBase()
    return seq


def randomMotif(length, numVaried):
    # 0 <= numVaried <= length
    numVaried = min(length, numVaried)
    numVaried = max(0, numVaried)

    motifList = list(randomSequence(length))
    varyPositions = random.sample(range(length), numVaried)
    for pos in varyPositions:
        motifList[pos] = '*'
    return "".join(motifList)


def sequenceFromMotif(motif):
    def changeStarToBase(s):
        if s == "*":
            return randomBase()
        else:
            return s
    return "".join([changeStarToBase(x) for x in list(motif)])


def overWriteAt(s1, s2, location):
    return s1[:location] + s2 + s1[location + len(s2):]


def plantAt(seq, motif, location):
    plantSeq = sequenceFromMotif(motif)
    return overWriteAt(seq, plantSeq, location)


def plantRandom(seq, motif):
    plantLocation = random.randrange(len(seq) - len(motif) + 1)
    return (plantAt(seq, motif, plantLocation), plantLocation)
