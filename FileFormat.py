import GenerateSequence as Gen
import os


def motifFile(name, motif):
    return name + "\t" + str(len(motif)) + "\t" + motif.upper()


def writeMotif(filename, motifname, motif):
    with open(filename, 'w') as f:
        f.write(motifFile(motifname, motif))


def writeMotifLength(filename, motif):
    with open(filename, 'w') as f:
        f.write(str(len(motif)))


def writeSequences(filename, sequences):
    with open(filename, 'w') as f:
        for (i, s) in enumerate(sequences):
            f.write("Sequence" + str(i) + ">\n")
            f.write(s.upper() + "\n")


def writeSites(filename, sites):
    with open(filename, 'w') as f:
        for site in sites:
            f.write(str(site) + "\n")


def generateDataset(name, ml, nm, sl, sc):
    folder = name + "_" + "_".join(str(x) for x in [ml, nm, sl, sc])

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


def generateAllDatasets():
    for i in range(10):
        name = "dataset" + str(i) + "/set"
        generateWithParameters(name)


def generateWithParameters(name):
    default_args = [name, 8, 1, 500, 10]
    for nm in [0, 1, 2]:
        args = default_args[:]
        args[2] = nm
        generateDataset(*args)

    for ml in [6, 7, 8]:
        args = default_args[:]
        args[1] = ml
        generateDataset(*args)

    for sc in [5, 10, 20]:
        args = default_args[:]
        args[4] = sc
        generateDataset(*args)
