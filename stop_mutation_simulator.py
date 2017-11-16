#!/usr/bin/env python3

import sys
import random
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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

seqLength = (1050,1011,1056,1158,1059,1218,1056,966,861,945,1101,1221,1071,2802,1653,1227,2520,1062,1011,1146,1164,1227,2958,1197,1047,1881,891,882,1743)

listePos = initData(seqLength)

param = {"nbGeneration":1000,"txMut":1e-1,"freqMutaStop":0.036,"popSize":313,"probMigra":0,"indMigra":0.01}

dimension = (((param["nbGeneration"]//10)+1),(len(seqLength)+1))
countMutaMatrix = np.zeros(dimension, np.int8)
#print(countMutaMatrix)

for j in range(100):
    sys.stdout.write("\rDoing repeat %i on 1000" % j)
    sys.stdout.flush()

    listePos = initData(seqLength)
    for i in range(param["nbGeneration"]):
        for pos in listePos:
            pos = evolve(pos, param)

        if(i % 10 == 0):
            #print(i)
            b = [el["freq"] for el in listePos]
            #print(mutationCounter(b))
            countMutaMatrix[(i//10),mutationCounter(b)] += 1



print(countMutaMatrix)
