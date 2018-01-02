#!/usr/bin/env python3

import sys
import random
import numpy as np
import numpy.random as npr
import time
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def initData(seqLength):
    listePos = list()
    for i in seqLength:
        listePos.append({"seqSize":i,"freq":0.0,"nb":0}) #size,freq muta,nb muta
    return(listePos)

"""def migrate(listePos, param):
    if param["txmigra"] != 0:
        for pos in listePos:
            if pos["freq"] > 0:
                pos["freq"] *= (1-())
    elif random.random() < param["probMigra"]:
        for pos in listePos:
            if pos["freq"] > 0:
                pos["freq"] *= (1-param["indMigra"])
    return(listePos)"""

def drift(freq,popSize):
    return npr.binomial(2*popSize,freq)/(2*popSize)

def selection(freq,param,freqMigra):
    averageSelection = ((1-freq)**2)*(1-freqMigra) + freqMigra + 2*freq*(1-freq)*(1-param["h"]*param["fitness"])*(1-freqMigra) + (freq**2)*(1-param["fitness"])*(1-freqMigra)
    newFreq = ((freq**2)*(1-param["fitness"])*(1-freqMigra) + freq*(1-freq)*(1-param["h"]*param["fitness"])*(1-freqMigra))/averageSelection
    return(newFreq)

def evolve(pos, param, freqMigra):
    pos["freq"] = selection(pos["freq"], param, freqMigra)

    probMuta = param["txMut"] * pos["seqSize"] * param["freqMutaStop"] * 2*param["popSize"] * (1-pos["freq"])
    if(random.random() < probMuta):
        pos["freq"] += 1/(2*param["popSize"])
        pos["nb"] += 1


    pos["freq"] = drift(pos["freq"],param["popSize"])
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
parser.add_argument('-S', '--fitness', type=float, required=False, default=0, help="Fitness of mutated alleles")
parser.add_argument('-H', '--dominance', type=float, required=False, default=0, help="Dominance of mutated alleles")

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
"txmigra":args.tauxmigra,
"fitness":args.fitness,
"h":args.dominance
}


seqLength = (1050,1011,1056,1158,1059,1218,1056,966,861,945,1101,1221,1071,2802,1653,1227,2520,1062,1011,1146,1164,1227,2958,1197,1047,1881,891,882,1743)

#listePos = initData(seqLength)

dimension = (((param["nbGeneration"]//10)),(len(seqLength)+1))
countMutaMatrix = np.zeros(dimension)
tps1 = time.clock()
for j in range(param["nbRepeats"]):
    sys.stderr.write("\rDoing repeat {} on {}".format(j,param["nbRepeats"]))
    sys.stderr.flush()
    listePos = initData(seqLength)
    for i in range(param["nbGeneration"]):
        freqMigra = 0
        if param["txmigra"] != 0:
            freqMigra = npr.binomial(2*param["popSize"],param["txmigra"])/(2*param["popSize"])
        elif random.random() < param["probMigra"]:
            freqMigra = param["indMigra"]

        for pos in listePos:
            pos = evolve(pos, param, freqMigra)

        if(i % 10 == 0):
            #print(i)
            b = [el["freq"] for el in listePos]
            countMutaMatrix[(i//10),mutationCounter(b)] += 1


print(countMutaMatrix)

with PdfPages("graphe_simu_stop_{}_generations_{}_repeats_{}_mu_{}_ind_{}_migra_{}_migrants_{}_mig_{}_h_{}_s.pdf".format(param["nbGeneration"], param["nbRepeats"], param["txMut"], param["popSize"],param["probMigra"], param["indMigra"],param["txmigra"], param["h"], param["fitness"])) as pdf:
    for i in range((len(seqLength)+1)):
        plt.plot(countMutaMatrix[:,i]/param["nbRepeats"], 'k-')
        plt.plot()
        plt.ylim([-0.05,1.05])
        plt.title("Probability of {} pseudogenes".format(i))
        pdf.savefig()
        plt.close()

np.savetxt("out_simu_stop_{}_generations_{}_repeats_{}_mu_{}_ind_{}_migra_{}_migrants_{}_mig_{}_h_{}_s.csv".format(param["nbGeneration"], param["nbRepeats"], param["txMut"], param["popSize"],param["probMigra"], param["indMigra"],param["txmigra"], param["h"], param["fitness"]), countMutaMatrix, delimiter=",")
