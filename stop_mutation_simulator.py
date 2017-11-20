#!/usr/bin/env python3

import sys
import random
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def initData(seqLength):
    listePos = list()
    for i in seqLength:
        listePos.append({"seqSize":i,"freq":0.0,"nb":0}) #size,freq muta,nb muta
    return(listePos)

def evolve(pos, param):
    probMuta = param["txMut"] * pos["seqSize"] * param["freqMutaStop"] * 2*param["popSize"] * (1-pos["freq"])
    if(random.random() < probMuta):
        pos["freq"] += 1/(2*param["popSize"])
        pos["nb"] += 1
    if(random.random() < param["probMigra"]):
        pos["freq"] = pos["freq"]*(1-param["indMigra"])
    #print(pos["freq"])
    pos["freq"] = sum(np.random.binomial(1,pos["freq"], 2*param["popSize"])) / (2*param["popSize"])
    return(pos)

def mutationCounter(listeFreq):
    return listeFreq.count(1)

parser = argparse.ArgumentParser(description="Simulate gene loss by stop mutations")
parser.add_argument('-g', '--generations', type=int, required=True, help="Number of generations")
parser.add_argument('-r', '--repeats', type=int, required=True, help="Number of repeats")
parser.add_argument('-m', '--mu', type=float, required=True, help="Mutation rate")
parser.add_argument('-N', '--Ne', type=int, required=True, help="Effective population size")
parser.add_argument('-s', '--stop', type=float, required=True, help="Frequency of stop mutations")
parser.add_argument('-f', '--freqmigra', type=float, required=True, help="Frequency of migration event")
parser.add_argument('-F', '--migrant', type=float, required=True, help="Frequency of migrant fish")

args = parser.parse_args()

param = {"nbGeneration":args.generations,
"nbRepeats":args.repeats,
"txMut":args.mu,
"freqMutaStop":args.stop,
"popSize":args.Ne,
"probMigra":args.freqmigra,
"indMigra":args.migrant}


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
