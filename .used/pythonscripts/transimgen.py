import matplotlib.pyplot as plt
from numpy import exp, linspace
from math import log10
import pandas as pd

def plotvalues(line):
    line = line.split()
    return tuple(map(float, (line[0], line[1])))

def process(data):
    return list(map(plotvalues, data.split('\n')[:-1]))

data = open('../datafiles/transmission.dat', 'r').read().split('end\n')[:-1]
data = list(map(process, data))

zeros = int(log10(len(data))) + 1
for i in range(len(data)):
    t, F = zip(*data[i])
    plt.ylim(0, 1.05)
    plt.xlabel('t')
    plt.ylabel('F')
    plt.plot(t, F, linewidth = 1)
    plt.savefig('images/frame' + str(i + 1).zfill(zeros))
    plt.clf()
