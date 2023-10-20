import matplotlib.pyplot as plt
from numpy import exp, linspace
from math import log10
import pandas as pd

def plotvalues(line):
    line = line.split()
    return tuple(map(float, (line[0], line[1])))

def process(data):
    return list(map(plotvalues, data.split('\n')[:-1]))

data = open('../datafiles/wave.dat', 'r').read().split('end\n')[:-1]
data = list(map(process, data))

zeros = int(log10(len(data))) + 1
for i in range(len(data)):
    x, psi = zip(*data[i])
    plt.ylim(0, 0.05)
    plt.text(300, 0.045, 'time step = ' + str(i))
    plt.xlabel('x')
    plt.ylabel('norm(\u03A8)')
    plt.plot(x, psi, linewidth = 1)
    plt.savefig('images/frame' + str(i + 1).zfill(zeros))
    plt.clf()
