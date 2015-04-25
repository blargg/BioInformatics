import GenerateSequence as Gen
import Motif as M
import os


def motifFile(name, motif):
    return name + "\t" + str(len(motif)) + "\t" + motif.upper()


def writeMotif(filename, motifname, motif):
    with open(filename, 'w') as f:
        f.write(motifFile(motifname, motif))


def readMotifFile(filename):
    with open(filename, 'r') as f:
        return f.readline().split('\t')[2].lower()


def writeMotifLength(filename, motif):
    with open(filename, 'w') as f:
        f.write(str(len(motif)))


def readMotifLength(filename):
    with open(filename, 'r') as f:
        return int(f.readline())


def writeSequences(filename, sequences):
    with open(filename, 'w') as f:
        for (i, s) in enumerate(sequences):
            f.write("Sequence" + str(i) + ">\n")
            f.write(s.upper() + "\n")


def readSequences(filename):
    sequences = []
    with open(filename, 'r') as f:
        while True:
            f.readline()
            seqline = f.readline().lower().strip()
            if not seqline:
                break

            if Gen.isSeq(seqline):
                sequences.append(seqline)
            else:
                print("read a line that was not a sequence:\n" + seqline)
    return sequences


def writeSites(filename, sites):
    with open(filename, 'w') as f:
        for site in sites:
            f.write(str(site) + "\n")


def readSites(filename):
    sites = []
    with open(filename, 'r') as f:
        for line in f:
            sites.append(int(line.strip()))
    return sites


def writePositionWeightMatrix(filename, motifLength, selections):
    subseq = M.extractSelections(motifLength, selections)

    counts = [{base: 0 for base in Gen.bases} for i in range(motifLength)]

    for seq in subseq:
        for i, char in enumerate(seq):
            counts[i][char] += 1

    with open(filename, 'w') as f:
        f.write(">PMOTIF\t" + str(motifLength) + "\n")
        for pos in counts:
            line = "\t".join(str(pos[c]) for c in Gen.bases)
            f.write(line + "\n")
        f.write("<")


def generateDataset(name, ml, nm, sl, sc):
    folder = name

    motif = Gen.randomMotif(ml, nm)
    seqs = [Gen.randomSequence(sl) for i in range(sc)]
    plants = [Gen.plantRandom(seq, motif) for seq in seqs]
    plantedseqs = [x for (x, _) in plants]
    plantLocation = [x for (_, x) in plants]

    if not os.path.exists(folder):
        os.makedirs(folder)
    writeSequences(folder + "/sequences.fa", plantedseqs)
    writeMotifLength(folder + "/motiflength.txt", motif)
    writeMotif(folder + "/motif.txt", "MOTIF1", motif)
    writeSites(folder + "/sites.txt", plantLocation)


def generatePerameterSettings():
    allsettings = []
    default_args = [8, 1, 500, 10]
    allsettings.append(default_args)

    generateDataset(*default_args)

    for nm in [0, 2]:
        args = default_args[:]
        args[1] = nm
        allsettings.append(args)

    for ml in [6, 7]:
        args = default_args[:]
        args[0] = ml
        allsettings.append(args)

    for sc in [5, 20]:
        args = default_args[:]
        args[3] = sc
        allsettings.append(args)


def generateAllDatasets():
    for args in generatePerameterSettings():
        for i in range(10):
            name = "datasets/dataset_" + "_".join(str(x) for x in args)\
                   + "/dataset" + str(i)
            generateDataset(name, *args)


def generateWithParameters(name):
    default_args = [name, 8, 1, 500, 10]

    generateDataset(*default_args)

    for nm in [0, 2]:
        args = default_args[:]
        args[2] = nm
        generateDataset(*args)

    for ml in [6, 7]:
        args = default_args[:]
        args[1] = ml
        generateDataset(*args)

    for sc in [5, 20]:
        args = default_args[:]
        args[4] = sc
        generateDataset(*args)
