#!/usr/bin/env python3

import sys
import random
import numpy as np
import numpy.random as npr
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def initData(seqLength):
    listePos = list()
    for i in seqLength:
        listePos.append({"seqSize":i,"freq":0.0,"nb":0}) #size,freq muta,nb muta
    return(listePos)

def migrate(freq, param):
    if param["txmigra"] != 0:
        print(freq)
        print("Migra L")
        freq *= (1-(npr.binomial(2*param["popSize"],param["txmigra"])/(2*param["popSize"])))
        print(freq)
    elif random.random() < param["probMigra"]:
        print("Migra F")
        freq *= (1-param["indMigra"])

    return(freq)

def evolve(pos, param):
    probMuta = param["txMut"] * pos["seqSize"] * param["freqMutaStop"] * 2*param["popSize"] * (1-pos["freq"])
    if(random.random() < probMuta):
        pos["freq"] += 1/(2*param["popSize"])
        pos["nb"] += 1

    if(pos["freq"] > 0):
        pos["freq"] = migrate(pos["freq"], param)

    pos["freq"] = npr.binomial(2*param["popSize"],pos["freq"]) / (2*param["popSize"])
    return(pos)

def mutationCounter(listeFreq):
    return listeFreq.count(1)

parser = argparse.ArgumentParser(description="Simulate gene loss by stop mutations")
parser.add_argument('-g', '--generations', type=int, required=True, help="Number of generations")
parser.add_argument('-r', '--repeats', type=int, required=True, help="Number of repeats")
parser.add_argument('-m', '--mu', type=float, required=True, help="Mutation rate")
parser.add_argument('-N', '--Ne', type=int, required=True, help="Effective population size")
parser.add_argument('-s', '--stop', type=float, required=True, help="Frequency of stop mutations")
parser.add_argument('-f', '--freqmigra', type=float, required=False, help="Frequency of migration event", default=0)
parser.add_argument('-F', '--migrant', type=float, required=False, help="Frequency of migrant fish", default=0)
parser.add_argument('-l', '--tauxmigra', type=float, required=False,help="Combination of frequency of migration event and of migrant fish", default=0)

args = parser.parse_args()

if args.tauxmigra == 0 and ((args.freqmigra == 0 and args.migrant != 0) or (args.freqmigra != 0 and args.migrant == 0)):
    args.freqmigra = 0
    args.migrant = 0
    print("Assuming no migration")

if args.tauxmigra != 0 and (args.freqmigra != 0 or args.migrant != 0):
    print("ERROR : option -l can't be used together with -F and -f. Quitting")
    sys.exit()


param = {"nbGeneration":args.generations,
"nbRepeats":args.repeats,
"txMut":args.mu,
"freqMutaStop":args.stop,
"popSize":args.Ne,
"probMigra":args.freqmigra,
"indMigra":args.migrant,
"txmigra":args.tauxmigra
}


seqLength = (1050,1011,1056,1158,1059,1218,1056,966,861,945,1101,1221,1071,2802,1653,1227,2520,1062,1011,1146,1164,1227,2958,1197,1047,1881,891,882,1743)

listePos = initData(seqLength)



dimension = (((param["nbGeneration"]//10)),(len(seqLength)+1))
countMutaMatrix = np.zeros(dimension)
#print(countMutaMatrix)
tps1 = time.clock()
for j in range(param["nbRepeats"]):
    sys.stderr.write("\rDoing repeat {} on {}".format(j,param["nbRepeats"]))
    sys.stderr.flush()
    listePos = initData(seqLength)
    for i in range(param["nbGeneration"]):
        for pos in listePos:
            pos = evolve(pos, param)
        if(i % 10 == 0):
            #print(i)
            b = [el["freq"] for el in listePos]
            #print(mutationCounter(b))
            countMutaMatrix[(i//10),mutationCounter(b)] += 1


with PdfPages("simu_stop_{}_generations_{}_repeats_{}_mu_{}_ind_{}_migra_{}_migrants.pdf".format(param["nbGeneration"], param["nbRepeats"], param["txMut"], param["popSize"],param["probMigra"], param["indMigra"])) as pdf:
    for i in range((len(seqLength)+1)):
        countMutaMatrix /= param["nbRepeats"]
        plt.plot(countMutaMatrix[:,i], 'k-')
        plt.ylim([-0.05,1.05])
        plt.title("Probability of {} pseudogenes".format(i))
        pdf.savefig()
        plt.close()
